/*  UPDATE:  Update system for 1 cycle of num_step */
/***********************************************************************
ASLB : Lattice-Boltzmann simulation code for deformable particle+polymer+
lattice Boltzmann fluid
Copyright (C) 2019 Yeng-Long Chen

based on 
Susp3D: Lattice-Boltzmann simulation code for particle-fluid suspensions
Copyright (C) 2003 Tony Ladd

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
***********************************************************************/

/* calls lbe_update to calc velcs_df  */
/* calls n_list to create neighbor lists */
/* calls velcs_update to update velocity */
/* calls hs3d to check for collisions */
/* accumulates particle properties */
/* writes checkpoint file */

#include  "header.h"
#define UPDATE_DEBUG 0
int max_cl_size;

int update (struct object *objects, Float ***velcs_df, int **node_map, struct vector f_ext, double y_prop[MAX_Y][Num_Prop],
	    int num_cycle, int n_cycle, int num_step, int num_modes, int backflow_flag, char *modes_file, char *chk_p_file, 
	    char *chk_f_file, char *chk_c_file, struct DP *DPs, struct DP *cylinders, struct monomer *monomers, struct face *faces, struct DP_param *sphere_pm, struct DP_param *polym_pm, struct DP_param *cyl_pm, VSLStreamStatePtr rngstream, char *work_dir)
{
  extern int    seed;
  extern int    max_x, max_y, max_z;
  extern int    num_x, x_min, x_max;
  extern int    addsteps;
  extern int    num_sph, num_obj;
  extern int    n_proc;
  extern int    wall_flag;
  static struct vector f_mom, f_tot, backflow;
  static struct vector p_lbe, p_sph, p_tot, p_old;
  static double vol, vol_frac;
  static int    init_flag=1, num_list=0;
  double sum_y[MAX_Y][16];
  double den_fac;
  double dx, dy, dz, dt;
  double yn, y1, y2, rad, phi;
  int    n_step, n_sph, n_obj, n_mode, list_flag;
  int    n_dir, n, num_stop, num_coll, max_coll;
  int    Nbead;
  int    x, y, y_front;
  int    i, j, d, temp;
  int    type, start, start0;
  int    errcode;
  void   *io_ptr;
  char   filename[MAX_NL];
  FILE   *file_ptr=0;

  if (init_flag) 
    {
      for (n_sph = 0; n_sph < num_sph; n_sph++)                   /* Zero collision counters */
	objects[n_sph].num_coll = 0;
      
      vol = max_x*max_y*max_z;                                    /* Total volume */
      f_tot.x = f_ext.x;  f_tot.y = f_ext.y;  f_tot.z = f_ext.z;  /* Set external force density */
      
      /*  Read checkpoint file  */
      file_ptr = fopen (chk_p_file, "rb");
      if (file_ptr != 0)
	{
	  io_ptr = (void *) &addsteps;
	  
	  fread (io_ptr, sizeof(int), 1, file_ptr);
	  for (n_sph = 0; n_sph < num_sph; n_sph++)
	    {
	      io_ptr = (void *)    &objects[n_sph].r;
	      fread (io_ptr, sizeof(objects[n_sph].r), 1, file_ptr);
	      io_ptr = (void *)    &objects[n_sph].e;
	      fread (io_ptr, sizeof(objects[n_sph].e), 1, file_ptr);
	      io_ptr = (void *)    &objects[n_sph].u;
	      fread (io_ptr, sizeof(objects[n_sph].u), 1, file_ptr);
	      io_ptr = (void *)    &objects[n_sph].w;
	      fread (io_ptr, sizeof(objects[n_sph].w), 1, file_ptr);
	      io_ptr = (void *)    &objects[n_sph].r_map;
	      fread (io_ptr, sizeof(objects[n_sph].r_map), 1, file_ptr);
	      io_ptr = (void *)    &objects[n_sph].stop_flag;
	      fread (io_ptr, sizeof(objects[n_sph].stop_flag), 1, file_ptr);
	      if (objects[n_sph].stop_flag)
		{
		  objects[n_sph].mass_flag = 0;
		  objects[n_sph].move_flag = 0;
		}
	    }

	  io_ptr = (void *)    &p_old;

	  fread (io_ptr, sizeof(p_old), 1, file_ptr);
	  if (backflow_flag > 1)
	    {
	      io_ptr = (void *)    &f_tot;
	      fread (io_ptr, sizeof(f_tot), 1, file_ptr);
	    }
	  fclose (file_ptr);  file_ptr = 0;
	}

      file_ptr = fopen (chk_f_file, "rb");
      if (file_ptr != 0)
	{
	  io_ptr = (void *) &x_min;
	  fread (io_ptr, sizeof(int), 1, file_ptr);
	  io_ptr = (void *) &x_max;
	  fread (io_ptr, sizeof(int), 1, file_ptr);
	  for (x = 1; x <= num_x; x++)
	    for (y = 1; y <=  max_y; y++)
	      for (n_dir = 0; n_dir < Num_Dir; n_dir++)
		{
		  n = x*(max_y+2) + y;
		  io_ptr = (void *) velcs_df[n][n_dir];
		  fread (io_ptr, sizeof(Float), max_z+2, file_ptr);
		}
	  fclose (file_ptr);  file_ptr = 0;
	}

      file_ptr = fopen (chk_c_file, "rb");
      if (file_ptr != 0)
	{
	  fscanf(file_ptr, "%d %d %d %d %d\n", &sphere_pm->Ntype[0], &sphere_pm->Ntype[1], &sphere_pm->nlevel[0], &sphere_pm->nlevel[1], &sphere_pm->num_beads);
	  fscanf(file_ptr, "%d %d %d %d %d\n", &polym_pm->Ntype[0], &polym_pm->Ntype[1], &polym_pm->N_per_DP[0], &polym_pm->N_per_DP[1], &polym_pm->num_beads);
	  for(x=0; x< sphere_pm->NDP; x++)
	    fscanf(file_ptr, "%le %le %le %le\n", &DPs[x].volume, &DPs[x].area, &DPs[x].disp2, &DPs[x].dr2);
	  for(x=sphere_pm->NDP; x< sphere_pm->NDP+polym_pm->NDP; x++)
	    fscanf(file_ptr, "%le %le %le\n", &DPs[x].Rg2, &DPs[x].disp2, &DPs[x].dr2);
	  /*
	  if(cyl_pm->NDP > 0)
	    for(x=0; x< cyl_pm->NDP; x++)
	      fscanf(file_ptr, "%le %le\n", &cylinders[x].disp2, &cylinders[x].dr2);
	  */
	  for (x = 0; x < sphere_pm->num_beads+polym_pm->num_beads+cyl_pm->num_beads; x++)
	    fscanf(file_ptr, "%le %le %le %le %le %le %le %le %le %le\n", &monomers[x].pos[0], &monomers[x].pos[1], &monomers[x].pos[2], &monomers[x].pos_pbc[0], &monomers[x].pos_pbc[1], &monomers[x].pos_pbc[2], &monomers[x].vel[0], &monomers[x].vel[1], &monomers[x].vel[2], &monomers[x].dr2);
	  fclose (file_ptr);  file_ptr = 0;
	  
	  for(x=0; x<sphere_pm->num_beads+polym_pm->num_beads+cyl_pm->num_beads; x++)
	    printf("#monomer %d at (%le %le %le) \n", x, monomers[x].pos_pbc[0], monomers[x].pos_pbc[1], monomers[x].pos_pbc[2]);
	  
	  for (j=0 ; j<NTYPES ; j++) {
	    if(sphere_pm->nlevel[j] == -1) {
	      sphere_pm->N_per_DP[j] = 162;
	      sphere_pm->face_per_DP[j] = 320;
	    }
	    else if(sphere_pm->nlevel[j] == -2) {
	      sphere_pm->N_per_DP[j] = 642;
	      sphere_pm->face_per_DP[j] = 1280;
	    }
	    else {
	      sphere_pm->N_per_DP[j] = 12;
	      sphere_pm->face_per_DP[j] = 20;
	      temp = 30;
	      for (i=0 ; i<sphere_pm->nlevel[j] ; i++) {
		sphere_pm->N_per_DP[j] += temp;
		sphere_pm->face_per_DP[j] *= 4;
		temp *= 4;
	      }
	    }
	  }
	
	  /* calculate sphere com */
	  if(sphere_pm->NDP > 0) {
	    for(i=0; i<sphere_pm->NDP; i++) {
	      type = (i<sphere_pm->Ntype[0] ? 0 :1);
	      Nbead = sphere_pm->N_per_DP[type];
	      start = (type==0 ? i*Nbead : sphere_pm->Ntype[0]*sphere_pm->N_per_DP[0] + (i-sphere_pm->Ntype[0])*Nbead);

	      for(d=0; d<DIMS; d++) 
		DPs[i].com[d] = 0.0;
	      
	      for(x=0; x<Nbead; x++) {
		for(d=0; d<DIMS; d++)
		  DPs[i].com[d] += monomers[start+x].pos_pbc[d];
	      }
	      for(d=0; d<DIMS; d++) {
		DPs[i].com[d] /= Nbead;
		DPs[i].com0[d] = DPs[i].com[d];
		DPs[i].comold[d] = DPs[i].com[d];
	      }
	    }
	  }
	  if(polym_pm->NDP > 0) {
	    for(i=0; i<polym_pm->NDP; i++) {
	      type = (i<polym_pm->Ntype[0] ? 0 :1);
	      Nbead = polym_pm->N_per_DP[type];
	      start = (type==0 ? i*Nbead+sphere_pm->num_beads : polym_pm->Ntype[0]*polym_pm->N_per_DP[0] + (i-polym_pm->Ntype[0])*Nbead + sphere_pm->num_beads);

	      for(d=0; d<DIMS; d++) 
		DPs[i].com[d] = 0.0;
	      
	      for(x=0; x<Nbead; x++) {
		for(d=0; d<DIMS; d++)
		  DPs[i].com[d] += monomers[start+x].pos_pbc[d];
	      }
	      for(d=0; d<DIMS; d++) {
		DPs[i].com[d] /= Nbead;
		DPs[i].com0[d] = DPs[i].com[d];
		DPs[i].comold[d] = DPs[i].com[d];
	      }
	    }
	  }
	}
    
  /*     /\*read in bond list*\/ */
/*       if(sphere_pm->NDP > 0) { */
/* 	sprintf(filename, "%s/init/bond.dat", work_dir); */
/* 	file_ptr = fopen(filename, "r"); */
/* 	if(file_ptr != 0) { */
/* 	  if(sphere_pm->Ntype[0] > 0) { */
/* 	    for(i=0; i<sphere_pm->N_per_DP[0]; i++) { */
/* 	      fscanf(file_ptr, "%*s %d\n",  &monomers[i].blist[0][0]); */
/* 	      for(j=1; j<=monomers[i].blist[0][0]; j++) */
/* 		fscanf(file_ptr, "%d %d %d\n", &monomers[i].blist[j][0], &monomers[i].blist[j][1], &monomers[i].blist[j][2]); */
/* 	    } */
/* 	  } */
/* 	  if(sphere_pm->Ntype[1] > 0) { */
/* 	    for(i=0; i<sphere_pm->N_per_DP[1]; i++) { */
/* 	      n=sphere_pm->Ntype[0]*sphere_pm->N_per_DP[0]+i; */
/* 	      fscanf(file_ptr, "%*s %d\n",  &monomers[n].blist[0][0]); */
/* 	      for(j=1; j<=monomers[n].blist[0][0]; j++) */
/* 		fscanf(file_ptr, "%d %d %d\n", &monomers[n].blist[j][0], &monomers[n].blist[j][1], &monomers[n].blist[j][2]); */
/* 	    } */
/* 	  } */
/* 	  fclose(file_ptr); file_ptr = 0; */
	  
/* 	  for(i=1; i<sphere_pm->NDP; i++) { */
/* 	    type = (i<sphere_pm->Ntype[0] ? 0 : 1); */
/* 	    start0 = (type == 0 ? 0 : sphere_pm->Ntype[0]*sphere_pm->N_per_DP[0]); */
/* 	    start = (type == 0 ? i*sphere_pm->N_per_DP[0] : sphere_pm->Ntype[0]*sphere_pm->N_per_DP[0] + (i-sphere_pm->Ntype[0])*sphere_pm->N_per_DP[1]); */
/* 	    for(j=0; j<sphere_pm->N_per_DP[type]; j++) { */
/* 	      monomers[start+j].blist[0][0] = monomers[start0+j].blist[0][0]; */
/* 	      for(n=1; n<=monomers[j].blist[0][0]; n++) { */
/* 		monomers[start+j].blist[n][0] = start+(monomers[start0+j].blist[n][0]-start0); */
/* 		monomers[start+j].blist[n][1] = start+(monomers[start0+j].blist[n][1]-start0); */
/* 		monomers[start+j].blist[n][2] = start+(monomers[start0+j].blist[n][2]-start0); */
/* 	      } */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */

      init_flag = 0;
      
      if(cyl_pm->NDP > 0 || sphere_pm->num_beads > 200) {
	sphere_pm->nlist_flag = 1;
	n_list_mon(sphere_pm, polym_pm, cyl_pm, monomers, DPs);                        /* Set neighbor lists */
	num_list++;
      }

      n_cycle = addsteps / num_step;
      printf("addsteps %d, ncycle %d, numcycle %d\n", addsteps, n_cycle, num_cycle);

      if(addsteps == 0) n_cycle = 0;
      else sphere_pm->numsteps = addsteps;

      error_chk ();
    }
	   	
  /*  Loop over 1 cycle of num_step  */
  for (n_step = 0; n_step < num_step; n_step++)
    { 
      Write_Output(n_step, num_step, n_cycle, num_sph, objects, DPs, cylinders, monomers, faces, sphere_pm, polym_pm, cyl_pm, velcs_df, node_map, work_dir);

      if(UPDATE_DEBUG == 1)
	printf("Output done\n");

      if(sphere_pm->nlist_flag == 1 ) {
	for (n_sph = 0, list_flag = 0; n_sph < sphere_pm->num_beads; n_sph++)
	  {
	    dx = n_image(monomers[n_sph].pos[0] - monomers[n_sph].pos_lst[0], max_x);
	    dy = n_image(monomers[n_sph].pos[1] - monomers[n_sph].pos_lst[1], max_y);
	    dz = n_image(monomers[n_sph].pos[2] - monomers[n_sph].pos_lst[2], max_z);
	    
	    if (dx*dx+dy*dy+dz*dz > 0.05)  {
	      list_flag = 1;
	      break;
	    }
	  }
	if (list_flag)  n_list_mon(sphere_pm, polym_pm, cyl_pm, monomers, DPs);                                              /* Update neighbor lists */
	num_list += list_flag;
      }

      if(UPDATE_DEBUG == 1)
	printf("nlist done\n");

      if (num_modes > 0)                                                            /* Write velocity field */
	{
	  if (((n_cycle+1)%num_modes == 0) && (n_step == num_step-1))  file_ptr = fopen(modes_file, "a");
	}
				
      for (y = 1; y <=  max_y; y++)                                                  /* Zero sum_y array */
	for (n_mode = 0; n_mode < 16; n_mode++)
	  sum_y[y][n_mode] = 0.0;
      
      /* Main LBE update routine for particle and fluid moves */
      lbe_update (objects, sphere_pm, polym_pm, cyl_pm, DPs, cylinders, monomers, velcs_df, node_map, f_tot, &p_lbe, sum_y, num_step, n_step, file_ptr, rngstream); 

      if(UPDATE_DEBUG == 1)
	printf("lbeupdate done\n");
      
      /* Update polymer beads with MD */
      if(sphere_pm->NDP+polym_pm->NDP+cyl_pm->NDP > 0) {
	for(n=0; n< sphere_pm->MD_steps; n++) {

	  verlet_update(monomers, faces, DPs, sphere_pm, polym_pm, cyl_pm, velcs_df, node_map, n_step, rngstream);
	  sync_procs();
	}
      }

      if(UPDATE_DEBUG == 1)
	printf("verletupdate done\n");

      /* update velocity distributions and ensure no hard sphere overlaps */
      for (n = 0; n <  Num_Step; n++)
	{
	  dt  = 1.0/(double) Num_Step;
	  velcs_update (objects, dt);                                               /* Update velocities */
	  if(num_obj > 0)
	    hs3d (objects, dt);                                                       /* Update coordinates */
	}

      for (y = 1; y <=  max_y; y++)
	for (n_mode = 0; n_mode < 16; n_mode++)
	  y_prop[y][n_mode] += sum_y[y+1][n_mode];
		
      for (n_sph = 0; n_sph < num_sph; n_sph++)                                     /* Particle props vs. y */
	{
	  if (range(box(objects[n_sph].r.x, max_x), x_min, x_max))                  /* Totals counted once */
	    {
	      y = (int) box(objects[n_sph].r.y, max_y);
	      y_prop[y][16] += objects[n_sph].mass;                                 /* Mass */
	      y_prop[y][17] += objects[n_sph].mass*objects[n_sph].u.x;              /* Momentum */
	      y_prop[y][18] += objects[n_sph].mass*objects[n_sph].u.y;
	      y_prop[y][19] += objects[n_sph].mass*objects[n_sph].u.z;
	      y_prop[y][20] += objects[n_sph].inertia*objects[n_sph].w.x;
	      y_prop[y][21] += objects[n_sph].inertia*objects[n_sph].w.y;
	      y_prop[y][22] += objects[n_sph].inertia*objects[n_sph].w.z;
	      y_prop[y][23] -= objects[n_sph].p_str.xx;                             /* Particle stress */
	      y_prop[y][24] -= objects[n_sph].p_str.yy;
	      y_prop[y][25] -= objects[n_sph].p_str.zz;
	      y_prop[y][26] -= objects[n_sph].p_str.yz;
	      y_prop[y][27] -= objects[n_sph].p_str.zx;
	      y_prop[y][28] -= objects[n_sph].p_str.xy;
	      y_prop[y][29] -= objects[n_sph].pc_str.xx;                            /* Collisional stress */
	      y_prop[y][30] -= objects[n_sph].pc_str.yy;
	      y_prop[y][31] -= objects[n_sph].pc_str.zz;
	      y_prop[y][32] -= objects[n_sph].pc_str.yz;
	      y_prop[y][33] -= objects[n_sph].pc_str.zx;
	      y_prop[y][34] -= objects[n_sph].pc_str.xy;
	      y_prop[y][35] -= objects[n_sph].pf_lub.xx;                            /* Lubrication stress */
	      y_prop[y][36] -= objects[n_sph].pf_lub.yy;
	      y_prop[y][37] -= objects[n_sph].pf_lub.zz;
	      y_prop[y][38] -= objects[n_sph].pf_lub.yz;
	      y_prop[y][39] -= objects[n_sph].pf_lub.zx;
	      y_prop[y][40] -= objects[n_sph].pf_lub.xy;
	      y_prop[y][41] += objects[n_sph].e_diss.t;                             /* Energy dissipation */
	      y_prop[y][42] += objects[n_sph].e_diss.r;
	      y_prop[y][43] += objects[n_sph].e_diss.c;
	    }
	}
      for (n_sph = 0; n_sph < num_sph; n_sph++)                                     /* Particle props vs. y */
	{
	  y = (int) box(objects[n_sph].r.y, max_y);                                 /* Sum particle fluid stresses */
	  y_prop[y][35] -= objects[n_sph].pf_str.xx;
	  y_prop[y][36] -= objects[n_sph].pf_str.yy;
	  y_prop[y][37] -= objects[n_sph].pf_str.zz;
	  y_prop[y][38] -= objects[n_sph].pf_str.yz;
	  y_prop[y][39] -= objects[n_sph].pf_str.zx;
	  y_prop[y][40] -= objects[n_sph].pf_str.xy;

	  objects[n_sph].p_str.xx  = 0.0;                                           /* Zero particle stress tensors */
	  objects[n_sph].p_str.yy  = 0.0;
	  objects[n_sph].p_str.zz  = 0.0;
	  objects[n_sph].p_str.yz  = 0.0;
	  objects[n_sph].p_str.zx  = 0.0;
	  objects[n_sph].p_str.xy  = 0.0;
	  objects[n_sph].pc_str.xx = 0.0;
	  objects[n_sph].pc_str.yy = 0.0;
	  objects[n_sph].pc_str.zz = 0.0;
	  objects[n_sph].pc_str.yz = 0.0;
	  objects[n_sph].pc_str.zx = 0.0;
	  objects[n_sph].pc_str.xy = 0.0;
	  objects[n_sph].pf_str.xx = 0.0;
	  objects[n_sph].pf_str.yy = 0.0;
	  objects[n_sph].pf_str.zz = 0.0;
	  objects[n_sph].pf_str.yz = 0.0;
	  objects[n_sph].pf_str.zx = 0.0;
	  objects[n_sph].pf_str.xy = 0.0;
	  objects[n_sph].pf_lub.xx = 0.0;
	  objects[n_sph].pf_lub.yy = 0.0;
	  objects[n_sph].pf_lub.zz = 0.0;
	  objects[n_sph].pf_lub.yz = 0.0;
	  objects[n_sph].pf_lub.zx = 0.0;
	  objects[n_sph].pf_lub.xy = 0.0;
	  objects[n_sph].e_diss.t  = 0.0;                                           /* Zero energy dissipation */
	  objects[n_sph].e_diss.r  = 0.0;
	  objects[n_sph].e_diss.c  = 0.0;
	}
    }

  f_mom.x = f_mom.y = f_mom.z = 0.0;
  for (y = 1; y <=  max_y; y++)                                                      /* Global sums of y_prop */
    {
      n_mode = Num_Prop;
      global_sum (y_prop[y], n_mode);
      f_mom.x += y_prop[y][1]/(num_step*Rho_Fl);                                    /* Calculate total fluid velocity */
      f_mom.y += y_prop[y][2]/(num_step*Rho_Fl);
      f_mom.z += y_prop[y][3]/(num_step*Rho_Fl);
    }
  for (n_sph = 0, vol_frac = 0; n_sph < num_sph; n_sph++)                           /* Add particle phase velocity */
    {
      vol_frac+= objects[n_sph].vol/vol;
      f_mom.x += objects[n_sph].u.x*objects[n_sph].vol;
      f_mom.y += objects[n_sph].u.y*objects[n_sph].vol;
      f_mom.z += objects[n_sph].u.z*objects[n_sph].vol;
    }
  f_mom.x *= Rho_Fl/vol;  f_mom.y *= Rho_Fl/vol;  f_mom.z *= Rho_Fl/vol;

  p_sph.x = p_sph.y = p_sph.z = 0.0;
  for (n_obj = 0; n_obj < num_obj; n_obj++)                                         /* Drag force */
    {
      p_sph.x += objects[n_obj].p.x;
      p_sph.y += objects[n_obj].p.y;
      p_sph.z += objects[n_obj].p.z;
    }

  for (n_sph = 0; n_sph < num_sph; n_sph++)                                         /* External particle force */
    {
      p_sph.x -= objects[n_sph].f_ext.x*num_step;
      p_sph.y -= objects[n_sph].f_ext.y*num_step;
      p_sph.z -= objects[n_sph].f_ext.z*num_step;
    }
  p_lbe.x += 0.5*f_tot.x*(1.0 - vol_frac)*vol;                                      /* Start up correction */
  p_lbe.y += 0.5*f_tot.y*(1.0 - vol_frac)*vol;
  p_lbe.z += 0.5*f_tot.z*(1.0 - vol_frac)*vol;
  p_tot.x  = p_lbe.x - p_old.x - f_tot.x*num_step*vol;
  p_tot.y  = p_lbe.y - p_old.y - f_tot.y*num_step*vol;
  p_tot.z  = p_lbe.z - p_old.z - f_tot.z*num_step*vol;
  p_old.x  = p_lbe.x;
  p_old.y  = p_lbe.y;
  p_old.z  = p_lbe.z;
  p_lbe.x  = p_tot.x;
  p_lbe.y  = p_tot.y;
  p_lbe.z  = p_tot.z;
  p_tot.x  = p_lbe.x + p_sph.x;
  p_tot.y  = p_lbe.y + p_sph.y;
  p_tot.z  = p_lbe.z + p_sph.z;

  if (backflow_flag > 1)                                                            /* Calculate extra force to set fluid at rest */
    {		
      switch (backflow_flag)
	{
	case 2:
	  backflow.x = 1;
	  backflow.y = 1;
	  backflow.z = 1;
	  break;
	case 3:
	  backflow.x = 1;
	  backflow.y = 0;
	  backflow.z = 0;
	  break;
	case 4:
	  backflow.x = 0;
	  backflow.y = 1;
	  backflow.z = 0;
	  break;
	case 5:
	  backflow.x = 0;
	  backflow.y = 0;
	  backflow.z = 1;
	  break;
	}
		
      f_tot.x = (f_ext.x - f_mom.x/num_step)*backflow.x;                            /* Calculate forcing for zero fluid momentum */
      f_tot.y = (f_ext.y - f_mom.y/num_step)*backflow.y;
      f_tot.z = (f_ext.z - f_mom.z/num_step)*backflow.z;
      for (n_obj = 0; n_obj < num_obj; n_obj++)                                     /* Correct for momentum lost to fixed objects */
	{
	  if (objects[n_obj].mass_flag == 0) 
	    {
	      f_tot.x += (objects[n_obj].p.x/(vol*num_step))*backflow.x;
	      f_tot.y += (objects[n_obj].p.y/(vol*num_step))*backflow.y;
	      f_tot.z += (objects[n_obj].p.z/(vol*num_step))*backflow.z;
	    }
	}
      if (n_proc == 0)  fprintf (stdout, "    Force density     % .5e % .5e % .5e\n", f_tot.x,  f_tot.y,  f_tot.z);
    }
	
  for (n_sph = 0, den_fac = 0.0; n_sph < num_sph; n_sph++)
    den_fac += objects[n_sph].mass/objects[n_sph].vol;                                           /* Mean particle mass density */
  den_fac /= num_sph;

  for (y = 1, y_front = 0; y <= max_y; y++)
    {
      for (n_sph = 0, phi = 0.0; n_sph < num_sph; n_sph++)
	{
	  yn  = box(objects[n_sph].r.y,max_y);
	  rad = objects[n_sph].r_a;
	  if (range(yn,y-rad,y+rad+1))
	    {
	      y1 = y - yn;
	      y2 = y - yn + 1;
	      if (y1 < -rad)  y1 = -rad;
	      if (y2 >  rad)  y2 =  rad;
	      phi += Pi* (double) (rad*rad*(y2-y1)-((y2*y2*y2)-(y1*y1*y1))/3);
	    }
	}
      phi /= (double)(max_x*max_z);
      if (phi > Phi_Sed)  y_front = y;                                                             /* Locate sedimentation front */
    }
  if (y_front > 0 && n_proc == 0)  fprintf (stdout, "    Front location%6d\n", y_front);
  for (n_sph = 0, num_stop = 0, num_coll = 0, max_coll = 0; n_sph < num_sph; n_sph++)
    {
      if (y_front > 0 && objects[n_sph].r.y < (y_front+1) && objects[n_sph].mass_flag)             /* Check for stopped particles */
	{
	  objects[n_sph].mass_flag = 0;
	  objects[n_sph].move_flag = 0;
	  objects[n_sph].stop_flag = 1;
	  objects[n_sph].u.x = 0.0;
	  objects[n_sph].u.y = 0.0;
	  objects[n_sph].u.z = 0.0;
	  objects[n_sph].w.x = 0.0;
	  objects[n_sph].w.y = 0.0;
	  objects[n_sph].w.z = 0.0;
	  if (n_proc == 0)  fprintf (stdout, "    Particle %6d stopped at % .3e % .3e % .3e\n", n_sph, objects[n_sph].r.x, objects[n_sph].r.y, objects[n_sph].r.z);
	}
      num_stop += objects[n_sph].stop_flag;
      num_coll += objects[n_sph].num_coll;
      max_coll  = max(objects[n_sph].num_coll, max_coll);
      objects[n_sph].num_coll  = 0;
    }

  if (n_proc == 0)
    {
      if (num_stop > 0)
	fprintf (stdout, "    Total stopped %6d\n", num_stop);
      fprintf (stdout, "    Max cluster   %6d\n", max_cl_size);
      fprintf (stdout, "    List updates  %6d\n", num_list);
      fprintf (stdout, "    Collision rate    % .3e\n", num_coll/(double) (max(num_step*num_sph,1)));
      fprintf (stdout, "    Max collisions    % .3e\n", max_coll/(double) num_step);
      fprintf (stdout, "    Total momenta     % .5e % .5e % .5e\n", p_tot.x, p_tot.y, p_tot.z);
      fflush  (stdout);

    }
  num_list = 0;
  max_cl_size = 0;

  /*  Write checkpoint file  */
  
  error_chk();                                                                      /* Exit if fatal error */   
	
  n_cycle++;                                                                        /* Increment cycle counter */
  addsteps += num_step;
  
  if (n_proc == 0)  fprintf (stdout, "Begin checkpoint  %6d\n", n_cycle);
  file_ptr = fopen (chk_p_file, "wb");
  if (file_ptr != 0)
    {
      io_ptr = (void *) &addsteps;
      fwrite (io_ptr, sizeof(int), 1, file_ptr);
      for (n_sph = 0; n_sph < num_sph; n_sph++)
	{
	  io_ptr = (void *)     &objects[n_sph].r;
	  fwrite (io_ptr, sizeof(objects[n_sph].r), 1, file_ptr);
	  io_ptr = (void *)     &objects[n_sph].e;
	  fwrite (io_ptr, sizeof(objects[n_sph].e), 1, file_ptr);
	  io_ptr = (void *)     &objects[n_sph].u;
	  fwrite (io_ptr, sizeof(objects[n_sph].u), 1, file_ptr);
	  io_ptr = (void *)     &objects[n_sph].w;
	  fwrite (io_ptr, sizeof(objects[n_sph].w), 1, file_ptr);
	  io_ptr = (void *)     &objects[n_sph].r_map;
	  fwrite (io_ptr, sizeof(objects[n_sph].r_map), 1, file_ptr);
	  io_ptr = (void *)     &objects[n_sph].stop_flag;
	  fwrite (io_ptr, sizeof(objects[n_sph].stop_flag), 1, file_ptr);
	}
      io_ptr = (void *)    &p_old;
      fwrite (io_ptr, sizeof(p_lbe), 1, file_ptr);
      io_ptr = (void *)    &f_tot;
      fwrite (io_ptr, sizeof(f_tot), 1, file_ptr);
      fclose (file_ptr);  file_ptr = 0;
    }
			
  file_ptr = fopen (chk_f_file, "wb");
  if (file_ptr != 0)
    {
      io_ptr = (void *) &x_min;
      fwrite (io_ptr, sizeof(int), 1, file_ptr);
      io_ptr = (void *) &x_max;
      fwrite (io_ptr, sizeof(int), 1, file_ptr);
      for (x = 1; x <= num_x; x++)
	for (y = 1; y <=  max_y; y++)
	  for (n_dir = 0; n_dir < Num_Dir; n_dir++)
	    {
	      n = x*(max_y+2) + y;
	      io_ptr = (void *) velcs_df[n][n_dir];
	      fwrite (io_ptr, sizeof(Float), max_z+2, file_ptr);
	    }
      fclose (file_ptr);  file_ptr = 0;
    }

  file_ptr = fopen (chk_c_file, "wb");
  if (file_ptr != 0)
  {
    if(sphere_pm->NDP+polym_pm->NDP > 0) 
      sphere_props(sphere_pm, polym_pm, DPs, monomers, faces);
    fprintf(file_ptr, "%d %d %d %d %d\n", sphere_pm->Ntype[0], sphere_pm->Ntype[1], sphere_pm->nlevel[0], sphere_pm->nlevel[1], sphere_pm->num_beads);
    fprintf(file_ptr, "%d %d %d %d %d\n", polym_pm->Ntype[0], polym_pm->Ntype[1], polym_pm->N_per_DP[0], polym_pm->N_per_DP[1], polym_pm->num_beads);

    for(x=0; x< sphere_pm->NDP; x++)
      fprintf(file_ptr, "%le %le %le %le\n", DPs[x].volume, DPs[x].area, DPs[x].disp2, DPs[x].dr2);
    for(x=sphere_pm->NDP; x< sphere_pm->NDP+polym_pm->NDP; x++)
      fprintf(file_ptr, "%le %le %le\n", DPs[x].Rg2, DPs[x].disp2, DPs[x].dr2);
    /*    for(x=0; x< cyl_pm->NDP; x++)
      fprintf(file_ptr, "%le %le\n", cylinders[x].disp2, cylinders[x].dr2);
    */
    for (x = 0; x < sphere_pm->num_beads+polym_pm->num_beads+cyl_pm->num_beads; x++)
      fprintf(file_ptr, "%le %le %le %le %le %le %le %le %le %le\n", monomers[x].pos[0], monomers[x].pos[1], monomers[x].pos[2], monomers[x].pos_pbc[0], monomers[x].pos_pbc[1], monomers[x].pos_pbc[2], monomers[x].vel[0], monomers[x].vel[1], monomers[x].vel[2], monomers[x].dr2);
    fclose (file_ptr);  file_ptr = 0;
  }
  
  sync_procs ();                                                                    /* Synchronize processors */
  if (n_proc == 0)
    {
      fprintf (stdout, "End checkpoint    %6d\n", n_cycle);
      fflush  (stdout);
    }

  return (n_cycle);
}
