/*  LBE_UPDATE:  Update lattice-Boltzmann for 1 step  */
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

#include  "header.h"
#define LBE_DEBUG 0

void lbe_update (struct object *objects, struct DP_param *sphere_pm, struct DP_param *polym_pm, struct DP_param *cyl_pm, struct DP *spheres, struct DP *cylinders, struct monomer *monomers, Float ***velcs_df, int ** node_map, struct vector f_ext, struct vector *p_lbe, double sum_y[MAX_Y][16], int num_step, int n_step, FILE *file_ptr, VSLStreamStatePtr rngstream) 
{
  extern int    add_noise;
  extern int    max_x, max_y, max_z;
  extern int    num_x;
  extern int    num_sph, num_obj;
  extern int    n_proc;
  static int    num_map = 0;
  double        sum_zcol[16], f_mom[4], vol;
  double        dx, dy, dz;
  double        rannum[max_y*max_z*Num_Dir_NConserve];
  int           i, j, k, q, nxy;
  int           x, y, m, n, n_sph, map_flag;
  int           errcode;
  FILE          *stream;

  map_flag = 0;
  if (n_step == 0)  map_flag = 1;

  /* if the particles move a distance > Del_Sq, update the node_map */
  for (n_sph = 0, map_flag = 0; n_sph < num_sph; n_sph++)
    {
      dx = objects[n_sph].r.x - objects[n_sph].r_map.x;
      dy = objects[n_sph].r.y - objects[n_sph].r_map.y;
      dz = objects[n_sph].r.z - objects[n_sph].r_map.z;
      if (dx*dx+dy*dy+dz*dz > Del_Sq)  map_flag = 1;
    }
  
  if (map_flag)                                                     /* Update boundary maps and apply fluid-solid boundary conditions */
    {
      bnodes_add  (objects, velcs_df);                              /* Create fluid in old boundary node regions */
      bnodes_init (objects, node_map);                              /* Create boundary node maps */
      bnodes_del  (objects, velcs_df);                              /* Delete fluid in new boundary node regions */
      bnodes_mom  (objects, node_map);                              /* Calculate moments of boundary nodes */
    }
  num_map += map_flag;

  if(LBE_DEBUG == 1) printf("Bnodes done\n");

  //if(cyl_pm->NDP > 0 && (cyl_pm->nlist_flag == 1 || map_flag==1))

  if(cyl_pm->NDP >0) {
    /* cyl_pm->nlist_flag needs to be set */
    /* No updating node_map if the tube is not moving */
    //    if(cyl_pm->nlist_flag ==1 && map_flag == 1)   
    if(cyl_pm->nlist_flag ==1 && map_flag == 1)   
      bnodes_tube(velcs_df, node_map, cyl_pm, cylinders, monomers, sphere_pm->num_beads+polym_pm->num_beads);         /*Tube boundary init */
  }

  if(LBE_DEBUG == 1) printf("bnodes_tube done\n");

  if(sphere_pm->NDP > 0)
    bnodes_dp(node_map, sphere_pm, spheres, monomers, n_step);                /* Deformable particle boundary init */

  if(LBE_DEBUG == 1) printf("bnodes_dp done\n");
                                                                      /* The box bc must go before bounce back or it will overwrite */
  lbe_bconds(velcs_df);	                                            /* Box boundary conditions across multi-processors*/
  if(LBE_DEBUG == 1) printf("bnodes_bconds done\n");

  if(cyl_pm->NDP > 0)
    bnodes_tube_bb(velcs_df, node_map, cyl_pm);

  if(LBE_DEBUG == 1) printf("bnodes_tube bounce back done\n");

  if(num_obj > 0) {
    bnodes_wall  (objects, velcs_df, node_map);                       /* Wall-fluid boundary conditions */
    bnodes_sph_1 (objects, velcs_df, node_map);                       /* Particle-fluid forces */
    force_sum (objects);                                              /* Sum forces over processors */
    implicit_force (objects, f_ext);                                  /* Implicit force calculation and particle velocity update */
    bnodes_sph_2 (objects, velcs_df, node_map);                       /* Particle-fluid boundary conditions */
  }

  if(LBE_DEBUG == 1) printf("bnodes_wall done\n");

  lbe_move (velcs_df);                                                /* XY propagation */

  if(LBE_DEBUG == 1) printf("lbe_move done\n");

  /* check_vel_conservation(velcs_df, n_step); */
  /*  Generate random number array */
  /* Collisions */
  for(i=0; i<max_y*max_z*Num_Dir_NConserve; i++)
    rannum[i]=0.0;
  for (x = 1; x <= num_x; x++) {                  /* Update LBE by columns */
    if(add_noise > 1) {
      errcode = vdRngGaussian( METHOD, rngstream, max_y*max_z*(Num_Dir_NConserve), rannum, 0.0e0, 1.0e0);
      CheckVslError(errcode);
    }
    
    for (y = 1; y <= max_y; y++) {
      n = x*(max_y+2)+y;
  
      lbe_zcol (velcs_df[n], node_map[n], f_ext, y, rannum);
    }
  }

  f_mom[0] = f_mom[1] = f_mom[2] = f_mom[3] = 0.0;

  if (file_ptr)	                                                  /* Write fluid modes */
    {
      for (x = 1; x <= num_x; x++)
	for (y = 1; y <= max_y; y++)
	  {
	    n = x*(max_y+2) + y;
	    modes_write(velcs_df[n], node_map[n], f_ext, file_ptr);
	  }
      fclose (file_ptr);
      file_ptr = 0;
    }
  n_step++;

  if (n_step%num_step == 0)                                         /* End of cycle */
    {
      global_sum (f_mom, 4);
      p_lbe->x = f_mom[1];
      p_lbe->y = f_mom[2];
      p_lbe->z = f_mom[3];
      for (n_sph = 0, vol = 0; n_sph < num_sph; n_sph++)
	vol += objects[n_sph].vol;
      f_mom[0] += (vol - max_x*max_y*max_z)*Rho_Fl;
      if(num_sph == 0) f_mom[0] = 0.0;
		
      if (n_proc == 0)
	{
	  fprintf (stdout, "    Total mass        % .5e\n", (Float)(max_x*max_y*max_z)*Rho_Fl);
	  fprintf (stdout, "    Fluid mass        % .5e\n", f_mom[0]);
	  fprintf (stdout, "    Particle mass   % .5e\n", vol*Rho_Fl);
	  fprintf (stdout, "    Map updates   %6d\n", num_map);
	  fflush (stdout);
	}
      num_map  = 0;
      if (fabs(f_mom[0]) > 10*Tol)                                  /* Mass conservation error */
	fatal_err("mass conservation error during step", n_step);                                             
      error_chk();                                                  /* Exit if Fatal Error */   
    }
}
