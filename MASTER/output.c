/***********************************************************************
 * ASLB : Lattice-Boltzmann simulation code for deformable particle+polymer+
 * lattice Boltzmann fluid
 * Copyright (C) 2019 Yeng-Long Chen
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * ***********************************************************************/


#include "header.h"

void sym33eigen_values(double mx[3][3], double lambda[3]);
void eigen_vector(double mx[3][3], double lambda, double v[3]);

void Write_Output(int n_step, int num_step, int n_cycle, int num_sph, struct object *objects, struct DP *DPs, struct DP *cylinders, struct monomer *monomers, struct face *faces, struct DP_param *sphere_pm, struct DP_param *polym_pm, struct DP_param *cyl_pm, Float ***velcs_df, int **node_map, char *work_dir) 
{
  extern int addsteps;
  int n = (addsteps+n_step)/sphere_pm->write_config;
  
  if(n_step % 1000 == 0)
    Write_time(n_step+addsteps, work_dir);
  if(num_sph != 0)
    Write_Particle(objects, sphere_pm, cyl_pm, n_step, work_dir);
  if(sphere_pm->num_beads+polym_pm->num_beads+cyl_pm->num_beads !=0) {
    Write_Monomer(DPs, cylinders, monomers, sphere_pm, polym_pm, cyl_pm, n_step, work_dir);
    Write_DP(DPs, monomers, faces, sphere_pm, polym_pm, n_step, work_dir);
    if((addsteps+n_step)%sphere_pm->write_config == 0) {
      if(sphere_pm->NDP > 0)
	Write_Sphereconfig(DPs, monomers, faces, sphere_pm, n, work_dir);
      if(polym_pm->NDP > 0)
	Write_Polymconfig(DPs, monomers, polym_pm, n, work_dir, sphere_pm->num_beads);
      if(cyl_pm->NDP > 0 && n_step == 0)
	Write_Cylconfig(cylinders, monomers, cyl_pm,  n, work_dir, sphere_pm->num_beads+polym_pm->num_beads);
    }
    //    Write_RDF(DPs, monomers, sphere_pm, n_step, n_cycle*num_step, num_step, work_dir, 0);
    //    Write_RDF(DPs, monomers, polym_pm, n_step, n_cycle*num_step, num_step, work_dir, sphere_pm->num_beads);
  }

  Write_Stress(DPs, monomers, sphere_pm, n_step, work_dir);
  /* Write_stat(velcs_df, n_step, sphere_pm->kT, work_dir); */ 

  if(n_step % sphere_pm->write_fluid == 0 && sphere_pm->write_fluid < num_step) {
    Write_Fluid(velcs_df, node_map, addsteps + n_step, work_dir); 
    Write_Nodemap(node_map, addsteps + n_step, work_dir);
    Write_Velslice(monomers, velcs_df, node_map, addsteps+n_step, work_dir); 
  }  
}

void Write_Particle(struct object *DP, struct DP_param *DP_pm, struct DP_param *cyl_pm, int n_step, char *work_dir)
{
  extern int max_x, max_y, max_z;
  extern int num_sph;

  static double x0[MAX_COLL], y0[MAX_COLL], z0[MAX_COLL];
  static int buff_nstep[WRITE_BUFFER];
  static double buff_avgdisp2[WRITE_BUFFER], buff_avgdr2[WRITE_BUFFER], buff_avgvel[WRITE_BUFFER];
  static double xp[MAX_COLL][WRITE_BUFFER], yp[MAX_COLL][WRITE_BUFFER], zp[MAX_COLL][WRITE_BUFFER];
  static double ux[MAX_COLL][WRITE_BUFFER], uy[MAX_COLL][WRITE_BUFFER], uz[MAX_COLL][WRITE_BUFFER];
  static double wx[MAX_COLL][WRITE_BUFFER], wy[MAX_COLL][WRITE_BUFFER], wz[MAX_COLL][WRITE_BUFFER];

  int i,j;
  int num_step;
  double dx, dy, dz;
  char filename[MAX_NL];
  double avg_disp2;
  double avg_dr2;
  double avg_vel;
  FILE *stream;

  if(n_step == 0) {
    for(i=0; i<num_sph; i++) {
      x0[i]=DP[i].r.x;
      y0[i]=DP[i].r.y;
      z0[i]=DP[i].r.z;
      DP[i].rold.x=DP[i].r.x;
      DP[i].rold.y=DP[i].r.y;
      DP[i].rold.z=DP[i].r.z;
      DP[i].dr2 = 0.0;
    }
    
    for(j=0; j<WRITE_BUFFER; j++) {
      buff_nstep[j]=0;
      buff_avgdisp2[j]=0.0;
      buff_avgvel[j]=0.0;

      for(i=0; i<MAX_COLL; i++) {
	xp[i][j]=0.0; yp[i][j]=0.0; zp[i][j]=0.0;
	ux[i][j]=0.0; uy[i][j]=0.0; uz[i][j]=0.0;
	wx[i][j]=0.0; wy[i][j]=0.0; wz[i][j]=0.0;
      }
    }
  }

  num_step = n_step / DP_pm->write_time;
  if(n_step % DP_pm->write_time == 0) {
    /* average mean squared displacement */
    for(i=0; i<num_sph; i++) {
      dx = n_image(DP[i].r.x - x0[i], max_x);
      dy = n_image(DP[i].r.y - y0[i], max_y);
      dz = n_image(DP[i].r.z - z0[i], max_z);
      DP[i].disp2 = dx*dx+dy*dy+dz*dz;
      
      dx = n_image(DP[i].r.x - DP[i].rold.x, max_x);
      dy = n_image(DP[i].r.y - DP[i].rold.y, max_y);
      dz = n_image(DP[i].r.z - DP[i].rold.z, max_z);
      DP[i].dr2 += dx*dx+dy*dy+dz*dz;
    }

    avg_disp2 = 0.0;
    avg_dr2 = 0.0;
    avg_vel = 0.0;
    for(i=0; i<num_sph; i++) {
      avg_disp2 += DP[i].disp2;
      avg_dr2 += DP[i].dr2;
      avg_vel += DP[i].u.x * DP[i].u.x + DP[i].u.y*DP[i].u.y + DP[i].u.z*DP[i].u.z;
    }
    avg_disp2 /= num_sph;
    avg_dr2 /= num_sph;
    avg_vel /= num_sph;

    /* Store the data into the write buffer */
    buff_nstep[n_step%WRITE_BUFFER]=n_step;
    buff_avgdisp2[n_step%WRITE_BUFFER]=avg_disp2;
    buff_avgdr2[n_step%WRITE_BUFFER]=avg_dr2;
    buff_avgvel[n_step%WRITE_BUFFER]=avg_vel;
    for(i=0; i<num_sph; i++) {
      xp[i][n_step%WRITE_BUFFER]=DP[i].r.x;
      yp[i][n_step%WRITE_BUFFER]=DP[i].r.y;
      zp[i][n_step%WRITE_BUFFER]=DP[i].r.z;
      ux[i][n_step%WRITE_BUFFER]=DP[i].u.x;
      uy[i][n_step%WRITE_BUFFER]=DP[i].u.y;
      uz[i][n_step%WRITE_BUFFER]=DP[i].u.z;
      wx[i][n_step%WRITE_BUFFER]=DP[i].w.x;
      wy[i][n_step%WRITE_BUFFER]=DP[i].w.y;
      wz[i][n_step%WRITE_BUFFER]=DP[i].w.z;
    }

    for(i=0; i<num_sph; i++) {
      DP[i].rold.x=DP[i].r.x;
      DP[i].rold.y=DP[i].r.y;
      DP[i].rold.z=DP[i].r.z;
    }

    if((num_step+1)%WRITE_BUFFER == 0) {
      sprintf(filename, "%s/data/avg_disp2.dat", work_dir);
      stream = fopen(filename, "a");
      for(j=0; j<WRITE_BUFFER; j++)
	fprintf(stream, "%d %le %le %le\n", buff_nstep[j], buff_avgdisp2[j], buff_avgdr2[j], buff_avgvel[j]);
      fclose(stream);
      
      for(i=0; i<10; i++) {
	sprintf(filename, "%s/data/coll%d.dat", work_dir, i);
	stream = fopen(filename, "a");
	for(j=0; j<WRITE_BUFFER; j++)
	  fprintf(stream, "%d %le %le %le %le %le %le\n", buff_nstep[j], xp[i][j], yp[i][j], zp[i][j], ux[i][j], uy[i][j], uz[i][j]);
	fclose(stream);
      }
    }
  }
}


void Write_Fluid(Float ***velcs_df, int **node_map, int nstep, char *work_dir)
{
  extern int max_x, max_y, max_z;
  extern int num_x;

  char filename[MAX_NL];
  double prop[DIMS+1+2*DIMS];  /* properties for fluid density, fluid velocity, and fluid stress */
  double lvx, lvy, lvz;
  FILE *stream;

  long i,j,k, p,q, q1, nxy;

  /* Write out velocity field in vtk format */
  write_velocity_field(nstep, velcs_df, node_map, work_dir);

}

void Write_Stress(struct DP *DPs, struct monomer *monomers, struct DP_param *DP_pm, int n_step, char *work_dir)
{
  extern int addsteps;
  extern int max_x, max_y, max_z;
  extern Float FluidStress[2*DIMS], ParStress[DIMS][DIMS];
  static double buff_stress[WRITE_BUFFER][6*DIMS];
  static int buff_nstep[WRITE_BUFFER];

  int i, j, n1, d, q, start, type;
  int NDP = DP_pm -> NDP;
  int Ntype[NTYPES];
  int Nbead[NTYPES];
  int num_beads[NTYPES];
  int num_step, tot_time;
  char filename[MAX_NL];
  FILE *stream;

  tot_time = (n_step+addsteps)*DP_pm->MD_steps*DP_pm->dt;

  if(addsteps == 0) {
    sprintf(filename, "%s/data/stress.dat", work_dir);
    stream = fopen(filename, "w");
    fprintf(stream, " #nstep fluidstress ParticleStress Totalstress \n");
    fclose(stream);	
  }

  num_step = n_step / DP_pm->write_time;
  if(n_step % DP_pm->write_time == 0) {    
    for (i=0 ; i<NTYPES ; i++) {
      Ntype[i] = DP_pm->Ntype[i];
      Nbead[i] = DP_pm->N_per_DP[i];
      num_beads[i] = Ntype[i]*Nbead[i];
    }
    
    if(Ntype[0]+Ntype[1] > 0) {
      /* calculate particle stress */
      
      for(i=0; i<DIMS; i++)
	for(j=0; j<DIMS; j++)
	  ParStress[i][j] = 0.0;
      
      for(i=0; i<NDP; i++) {
	type = (i<Ntype[0] ? 0 : 1);
	start = (type==0 ? i*Nbead[0] : Ntype[0]*Nbead[0]+(i-Ntype[0])*Nbead[1]);
	
	
	for (j=0 ; j<Nbead[type] ; j++) {
	  n1 = start + j;
	  
	  for(q=0; q<DIMS; q++)
	    for(d=0; d<DIMS; d++) {
	      monomers[n1].stress[q][d]+=(monomers[n1].pos_pbc[q] - DPs[i].com[q])*(monomers[n1].force[d] - monomers[n1].evforce[d] - monomers[n1].fricforce[d]*Rho_Fl);
	      
	      ParStress[q][d] += monomers[n1].stress[q][d];
	    }
	}
      }
      
      for(i=0; i<DIMS; i++)
	for(j=0; j<DIMS; j++) {
	  ParStress[i][j] /= -max_x*max_y*max_z;
	}
    }

    /* calculate the xx, yy, and zz stresses */
    buff_nstep[num_step%WRITE_BUFFER]=tot_time;
    buff_stress[num_step%WRITE_BUFFER][0] = (FluidStress[0]+FluidStress[1])/3.0;
    buff_stress[num_step%WRITE_BUFFER][1] = (FluidStress[0]+FluidStress[2])/2.0-buff_stress[num_step%WRITE_BUFFER][1];
    buff_stress[num_step%WRITE_BUFFER][2] = FluidStress[0]-buff_stress[num_step%WRITE_BUFFER][0]-buff_stress[num_step%WRITE_BUFFER][1];

    for(q=DIMS; q<2*DIMS; q++)
      buff_stress[num_step%WRITE_BUFFER][q]=FluidStress[q];
    
    buff_stress[num_step%WRITE_BUFFER][2*DIMS  ] = ParStress[0][0];
    buff_stress[num_step%WRITE_BUFFER][2*DIMS+1] = ParStress[1][1];
    buff_stress[num_step%WRITE_BUFFER][2*DIMS+2] = ParStress[2][2];
    buff_stress[num_step%WRITE_BUFFER][2*DIMS+3] = ParStress[0][1];
    buff_stress[num_step%WRITE_BUFFER][2*DIMS+4] = ParStress[1][2];
    buff_stress[num_step%WRITE_BUFFER][2*DIMS+5] = ParStress[0][2];

    buff_stress[num_step%WRITE_BUFFER][2*DIMS+6] = ParStress[0][0] + buff_stress[num_step%WRITE_BUFFER][0];
    buff_stress[num_step%WRITE_BUFFER][2*DIMS+7] = ParStress[1][1] + buff_stress[num_step%WRITE_BUFFER][1];
    buff_stress[num_step%WRITE_BUFFER][2*DIMS+8] = ParStress[2][2] + buff_stress[num_step%WRITE_BUFFER][2];
    buff_stress[num_step%WRITE_BUFFER][2*DIMS+9] = ParStress[0][1] + FluidStress[3];
    buff_stress[num_step%WRITE_BUFFER][2*DIMS+10] = ParStress[1][2] + FluidStress[4];
    buff_stress[num_step%WRITE_BUFFER][2*DIMS+11] = ParStress[0][2] + FluidStress[5];
  
    if((num_step+1)%WRITE_BUFFER == 0) {
      sprintf(filename, "%s/data/stress.dat", work_dir);
      stream = fopen(filename, "a");
      for(j=0; j<WRITE_BUFFER; j++) {
	fprintf(stream, "%d ", buff_nstep[j]);
	for(i=0; i<6*DIMS; i++)
	  fprintf(stream, "%le ", buff_stress[j][i]);
	fprintf(stream, "\n");
      }
      fclose(stream);
    }
  }
}


void Write_Monomer(struct DP *DPs, struct DP *cylinders, struct monomer *monomers, struct DP_param *sphere_pm, struct DP_param *polym_pm, struct DP_param *cyl_pm, int n_step, char *work_dir)
{
  extern int addsteps;
  extern int max_x, max_y, max_z;
  static int buff_nstep[WRITE_BUFFER];
  static double buff_avgdisp2[WRITE_BUFFER], buff_avgdr2[WRITE_BUFFER], buff_avgvel[WRITE_BUFFER];
  static double buff_monvel[WRITE_BUFFER][500], buff_tempscale[WRITE_BUFFER];

  char filename[MAX_NL];
  int j;
  int n, n2, d;
  int maxsize[DIMS];
  int num_write, num_step;
  int tot_time, add_time;
  double avg_disp2, avg_vel, avg_E, avg_dr2;
  int num_beads = sphere_pm->num_beads + polym_pm->num_beads;
  FILE *stream;

  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;

  tot_time = (n_step+addsteps)*sphere_pm->MD_steps*sphere_pm->dt;
  add_time = addsteps*sphere_pm->MD_steps*sphere_pm->dt;

  if(n_step == 0) {
    for(n=0; n<num_beads; n++)
      for(d=0; d<DIMS; d++) {
	monomers[n].pos0[d]=monomers[n].pos_pbc[d];
	monomers[n].pos_old[d]=monomers[n].pos_pbc[d];
      }

    for(j=0; j<WRITE_BUFFER; j++) {
      buff_nstep[j]=0;
      buff_avgdisp2[j]=0.0;
      buff_avgdr2[j]=0.0;
      buff_avgvel[j]=0.0;
      buff_tempscale[j]=0.0;
    }
  }

  if(sphere_pm->numsteps > 1001) num_write = sphere_pm->relax_time/(sphere_pm->MD_steps*sphere_pm->dt);
  else if(sphere_pm->num_beads == 1 || polym_pm->num_beads == 1) num_write = 10;
  else  num_write = sphere_pm->write_config;

  if((n_step+addsteps) % num_write == 0) {
    /* Write out monomer position and velocity */
    sprintf(filename, "%s/data/monpos.dat", work_dir);
    stream=fopen(filename, "a");
    if(num_beads < 500) 
      for(n=0; n<num_beads; n++) 
	fprintf(stream, "%d %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le\t %12.10le %12.10le %12.10le\n", tot_time, monomers[n].pos[0], monomers[n].pos[1], monomers[n].pos[2], monomers[n].vel[0], monomers[n].vel[1], monomers[n].vel[2], monomers[n].fluid_vel[0], monomers[n].fluid_vel[1], monomers[n].fluid_vel[2]);
    else
      for(n=0; n<num_beads; n++) 
	fprintf(stream, "%d %12.10le %12.10le %12.10le\n", tot_time, monomers[n].pos[0], monomers[n].pos[1], monomers[n].pos[2]);
    fclose(stream);

    /* write out bonding config */
    if (WRITE_BOND) {
      sprintf(filename, "%s/data/bond.dat", work_dir);
      stream=fopen(filename, "a");
      for (n=0 ; n<num_beads ; n++)
        for(j=1 ; j<=monomers[n].blist[0][0] ; j++) {
          n2 = monomers[n].blist[j][0];
          fprintf(stream, "%d %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le\n", tot_time, monomers[n].pos[0], monomers[n].pos[1], monomers[n].pos[2], monomers[n2].pos[0], monomers[n2].pos[1], monomers[n2].pos[2]);
        }
      fclose(stream);
    }

  }
  
  /* write out monomer velocities */
  if(num_beads < 5)
    for(n=0; n<num_beads; n++) {
      sprintf(filename, "%s/data/monvel_%d.dat", work_dir, n);
      stream=fopen(filename, "a");

      fprintf(stream, "%d %12.10le %12.10le %12.10le %12.10le\n", tot_time, monomers[n].vel[0], monomers[n].vel[1], monomers[n].vel[2], monomers[n].fricforce[0]);
      
      fclose(stream);
    }

  num_step = n_step / sphere_pm->write_time;
  if(n_step % sphere_pm->write_time == 0) {
    avg_disp2 = 0.0;
    avg_vel = 0.0;
    avg_dr2 = 0.0;
    
    for(n=0; n<num_beads; n++) {
      for(d=0; d<DIMS; d++) {
	avg_disp2 += (monomers[n].pos_pbc[d]-monomers[n].pos0[d])*(monomers[n].pos_pbc[d]-monomers[n].pos0[d]);
	avg_vel += monomers[n].vel[d] * monomers[n].vel[d];
	
	/* Calculate monomer mean square displacement */
	monomers[n].dr2 += (monomers[n].pos_pbc[d]-monomers[n].pos_old[d])*(monomers[n].pos_pbc[d]-monomers[n].pos_old[d]);
	monomers[n].pos_old[d]=monomers[n].pos_pbc[d];
      }
      avg_dr2 += monomers[n].dr2;
    }
    avg_disp2 /= num_beads;
    avg_vel /= num_beads;
    avg_dr2 /= num_beads;

    buff_nstep[num_step%WRITE_BUFFER]=tot_time;
    buff_avgdisp2[num_step%WRITE_BUFFER]=avg_disp2;
    buff_avgdr2[num_step%WRITE_BUFFER]=avg_dr2;
    buff_avgvel[num_step%WRITE_BUFFER]=avg_vel*Rho_Fl*sphere_pm->monmass;  
    buff_tempscale[num_step%WRITE_BUFFER]=sphere_pm->tempscale;  
    
    if((num_step+1)%WRITE_BUFFER == 0) {
      sprintf(filename, "%s/data/avg_disp2.dat", work_dir);
      stream = fopen(filename, "a");
      for(j=0; j<WRITE_BUFFER; j++)
	fprintf(stream, "%d %le %le %le\n", buff_nstep[j], buff_avgdisp2[j], buff_avgdr2[j], buff_avgvel[j]);
      fclose(stream);
    }
  } 
}

void Write_RDF(struct DP *DPs, struct monomer *mon, struct DP_param *DP_pm, struct DP_param *cyl_pm, int n_step, int num_step, char *work_dir, int n0)
{
  extern int addsteps;
  extern int max_x, max_y, max_z;
  static int num_relax=0, num_write=0;
  static double gofy[BINS], gcofy[BINS];
  static int buff_gofy[WRITE_BUFFER][BINS][NTYPES], buff_gcofy[WRITE_BUFFER][BINS][NTYPES];
  static double buff_stretchofy[WRITE_BUFFER][BINS][NTYPES];
  static double buff_Rxofy[WRITE_BUFFER][BINS][NTYPES], buff_Ryofy[WRITE_BUFFER][BINS][NTYPES], buff_Rzofy[WRITE_BUFFER][BINS][NTYPES];
  static int buff_nstep[WRITE_BUFFER];

  char filename[MAX_NL];
  int i,j, k;
  int n, d;
  int maxsize[DIMS];
  int Ntype[NTYPES];
  int Nbead[NTYPES];
  int num_beads[NTYPES];
  int tot_time, add_time, relax_step;
  int type, start;
  double *yc;
  //double range=(double)min(max_x, min(max_y, max_z));
  double range=max_y-1.;
  double binsize=range/BINS;
  FILE *stream;

  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;
  for (i=0 ; i<NTYPES ; i++) {
    Ntype[i] = DP_pm->Ntype[i];
    Nbead[i] = DP_pm->N_per_DP[i];
    num_beads[i] = Ntype[i]*Nbead[i];
  }
  tot_time = (n_step+addsteps)*DP_pm->MD_steps*DP_pm->dt;
  add_time = addsteps*DP_pm->MD_steps*DP_pm->dt;
  relax_step = DP_pm->relax_time/(DP_pm->MD_steps*DP_pm->dt);

  if((yc = (double *)malloc(DP_pm->NDP*sizeof(double)))==NULL) exit(1);
      
  /* Calculate the time dependent g(r,t) */
  if(n_step == 0) {
    for(j=0; j<WRITE_BUFFER; j++) {
      buff_nstep[j]=0;
      for(i=0; i<BINS; i++)
        for (k=0 ; k<NTYPES ; k++){
          buff_gofy[j][i][k]=0;
          buff_gcofy[j][i][k]=0;
          buff_stretchofy[j][i][k]=0.0;
          buff_Rxofy[j][i][k]=0.0;
          buff_Ryofy[j][i][k]=0.0;
          buff_Rzofy[j][i][k]=0.0;
        }
    }
  }

  if(relax_step > 10) {
  if((n_step+addsteps) % (relax_step/10) == 0) {
    for(i=0; i<DP_pm->NDP; i++) {
      type = (i<Ntype[0] ? 0 : 1);
      start = (type==0 ? n0+i*Nbead[0] : n0+Ntype[0]*Nbead[0]+(i-Ntype[0])*Nbead[1]);
      yc[i] = 0.0;
      for(j=0; j<Nbead[type]; j++) {
        yc[i]+=mon[start+j].pos[1];
        if (mon[start+j].pos[1] < range)
          buff_gofy[num_write%WRITE_BUFFER][(int)(mon[start+j].pos[1]/binsize)][type]++;
        else
          fprintf(stderr, "monomer %d y position out of range %le!\n", start+j, mon[start+j].pos[1]);
      }
      yc[i] /= Nbead[type];

      if(yc[i] < range) {
        buff_gcofy[num_write%WRITE_BUFFER][(int)(yc[i]/binsize)][type]++;
        buff_stretchofy[num_write%WRITE_BUFFER][(int)(yc[i]/binsize)][type] += sqrt(DPs[i].stretch[0]*DPs[i].stretch[0]+DPs[i].stretch[1]*DPs[i].stretch[1]);
        buff_Rxofy[num_write%WRITE_BUFFER][(int)(yc[i]/binsize)][type] += sqrt(DPs[i].Rx2);
        buff_Ryofy[num_write%WRITE_BUFFER][(int)(yc[i]/binsize)][type] += sqrt(DPs[i].Ry2);
        buff_Rzofy[num_write%WRITE_BUFFER][(int)(yc[i]/binsize)][type] += sqrt(DPs[i].Rz2);
      }
      else
        fprintf(stderr, "DP %d y position out of range %le!\n", i, yc[i]);
    }
    
    buff_nstep[num_write%WRITE_BUFFER]=tot_time;
    
    if((num_write+1)%WRITE_BUFFER == 0)
      for (k=0 ; k<NTYPES ; k++) {
        sprintf(filename, "%s/data/gofyoft-%d.dat", work_dir, k);
        stream=fopen(filename, "a");
        for(j=0; j<WRITE_BUFFER; j++)
          for(i=0; i<BINS; i++)
            fprintf(stream, "%d %le %le\n", buff_nstep[j], (double)i*binsize, (double)buff_gofy[j][i][k]/num_beads[k]);
        fclose(stream);

        sprintf(filename, "%s/data/gcofyoft-%d.dat", work_dir, k);
        stream=fopen(filename, "a");
        for(j=0; j<WRITE_BUFFER; j++)
          for(i=0; i<BINS; i++)
            fprintf(stream, "%d %le %le\n", buff_nstep[j], (double)i*binsize, (double)buff_gcofy[j][i][k]/Ntype[k]);
        fclose(stream);

        sprintf(filename, "%s/data/stretchofyoft-%d.dat", work_dir, k);
        stream=fopen(filename, "a");
        for(j=0; j<WRITE_BUFFER; j++)
          for(i=0; i<BINS; i++) {
            if (buff_gcofy[j][i][k] != 0) {
              buff_stretchofy[j][i][k] /= (double)buff_gcofy[j][i][k];
              buff_Rxofy[j][i][k] /= (double)buff_gcofy[j][i][k];
              buff_Ryofy[j][i][k] /= (double)buff_gcofy[j][i][k];
              buff_Rzofy[j][i][k] /= (double)buff_gcofy[j][i][k];
            }
            fprintf(stream, "%d %le %le %le %le %le\n", buff_nstep[j], (double)i*binsize, buff_stretchofy[j][i][k], buff_Rxofy[j][i][k], buff_Ryofy[j][i][k], buff_Rzofy[j][i][k]);
          }
        fclose(stream);

        for(j=0; j<WRITE_BUFFER; j++)
          for(i=0; i<BINS; i++) {
            buff_gofy[j][i][k]=0;
            buff_gcofy[j][i][k]=0;
            buff_stretchofy[j][i][k]=0.0;
            buff_Rxofy[j][i][k]=0.0;
            buff_Ryofy[j][i][k]=0.0;
            buff_Rzofy[j][i][k]=0.0;
          }
      }
    num_write++;
  }
  }  
  
  /* Calculate the overall g(r) and gc(r) */
  if(n_step ==0) {
    num_relax = addsteps/relax_step;
    printf("num_relax = %d\n", num_relax);
    
    for(i=0; i<BINS; i++) {
      gofy[i]=0.0;
      gcofy[i]=0.0;
    }
    
    if(num_relax > 0){
      if(n0 == 0) {
	sprintf(filename, "%s/data/gofysph.dat", work_dir);
	if((stream = fopen(filename, "r")) != NULL) {
	  fscanf(stream, "%*s %*s %*s %*s %*s %*s\n");
	  for(i=0; i<BINS; i++)
	    fscanf(stream, "%*s %lf\n", &gofy[i]);
	  gofy[i] *= (double)num_relax;
	  fclose(stream);
	}
	sprintf(filename, "%s/data/gcofysph.dat", work_dir);
	if((stream = fopen(filename, "r")) != NULL) {
	  for(i=0; i<BINS; i++)
	    fscanf(stream, "%*s %lf\n", &gcofy[i]);
	  gcofy[i] *= (double)num_relax;
	  fclose(stream);
	}
      }
      else {
	sprintf(filename, "%s/data/gofypoly.dat", work_dir);
	if((stream = fopen(filename, "r")) != NULL) {
	  fscanf(stream, "%*s %*s %*s %*s %*s %*s\n");
	  for(i=0; i<BINS; i++)
	    fscanf(stream, "%*s %lf\n", &gofy[i]);
	  gofy[i] *= (double)num_relax;
	  fclose(stream);
	}
	sprintf(filename, "%s/data/gcofypoly.dat", work_dir);
	if((stream = fopen(filename, "r")) != NULL) {
	  for(i=0; i<BINS; i++)
	    fscanf(stream, "%*s %lf\n", &gcofy[i]);
	  gcofy[i] *= (double)num_relax;
	  fclose(stream);
	}
      }
    }
  }
      
  if((n_step+addsteps) % relax_step == 0) {
    for(i=0; i<DP_pm->NDP; i++) {
      type = (i<Ntype[0] ? 0 : 1);
      start = (type==0 ? i*Nbead[0] : Ntype[0]*Nbead[0]+(i-Ntype[0])*Nbead[1]);
      yc[i] = 0.0;
      
      for(j=0; j<Nbead[type]; j++) {
        yc[i]+=mon[start+j].pos[1];

        gofy[(int)(mon[start+j].pos[1]/binsize)]++;
      }
      
      gcofy[(int)(yc[i]/(binsize*Nbead[type]))]++;
    }
    
    num_relax++;
  }
  
  if(n_step == num_step-1 && n_step > relax_step) {
    if(n0 == 0) {
      sprintf(filename, "%s/data/gofysph.dat", work_dir);
      stream=fopen(filename, "w");
      fprintf(stream, "# %d beads %d Relaxation times \n", DP_pm->num_beads, num_relax-1);
      for(i=0; i<BINS; i++)
	fprintf(stream, "%le %le\n", i*binsize, gofy[i]/(DP_pm->num_beads*num_relax));
      fclose(stream);
      sprintf(filename, "%s/data/gcofysph.dat", work_dir);
      stream=fopen(filename, "w");
      for(i=0; i<BINS; i++)
	fprintf(stream, "%le %le\n", i*binsize, gcofy[i]/(DP_pm->NDP*num_relax));
      fclose(stream);
    }
    else {
      sprintf(filename, "%s/data/gofypoly.dat", work_dir);
      stream=fopen(filename, "w");
      fprintf(stream, "# %d beads %d Relaxation times \n", DP_pm->num_beads, num_relax-1);
      for(i=0; i<BINS; i++)
	fprintf(stream, "%le %le\n", i*binsize, gofy[i]/(DP_pm->num_beads*num_relax));
      fclose(stream);
      sprintf(filename, "%s/data/gcofypoly.dat", work_dir);
      stream=fopen(filename, "w");
      for(i=0; i<BINS; i++)
	fprintf(stream, "%le %le\n", i*binsize, gcofy[i]/(DP_pm->NDP*num_relax));
      fclose(stream);
    }
  }
  free(yc);
}

void Write_DP(struct DP *spheres, struct monomer *monomers, struct face *faces, struct DP_param *sphere_pm, struct DP_param *polym_pm, int n_step, char *work_dir)
{
  extern int addsteps;
  int i, j, k, n, d;
  char filename[MAX_NL];
  int num_spheres = sphere_pm->NDP;
  int num_poly = polym_pm->NDP;
  int num_DP = num_spheres+num_poly;
  int num_step;
  double avg_disp2;
  double  avg_Rg2;
  double avg_stretchx, avg_stretchy, avg_stretchz;
  double avg_dr2;
  double avg_sp_len;
  double avg_asphr;
  double avg_acldr;
  static int   buff_nstep[WRITE_BUFFER];
  static double buff_avgdisp2[WRITE_BUFFER];
  static double buff_avgdr2[WRITE_BUFFER];
  static double buff_Rg2[WRITE_BUFFER];
  static double buff_stretch[WRITE_BUFFER][DIMS];
  static double buff_sp_len[WRITE_BUFFER];
  static double buff_asphr[WRITE_BUFFER];
  static double buff_acldr[WRITE_BUFFER];
  static double buff_sphereprops[WRITE_BUFFER][14][MAX_DP];
  static double buff_sphereshape[WRITE_BUFFER][6][MAX_DP];
  static double buff_spherevalues[WRITE_BUFFER][2][MAX_DP];
  static double buff_spherevolume[WRITE_BUFFER][4][MAX_DP];
  FILE *stream;

  if(n_step == 0) {
    if(addsteps == 0) {
      sprintf(filename, "%s/data/avg_props.dat", work_dir);
      stream = fopen(filename, "w");
      fprintf(stream, " #nstep    disp2   dr2    Rg2     stretch(x,y,z)   sprng_len   asphericity  acylindricity\n");
      fclose(stream);

      for(n=0; n<num_DP; n++) {
        sprintf(filename, "%s/data/DP_props.%d.dat", work_dir, n);
        stream=fopen(filename, "w");
        fprintf(stream, " #nstep   com   Rg2   avg_sprng_len   disp2   dr2    stretch(x,y,z)    volume     area    theta    omega[z]\n");
        fclose(stream);
      }

      for(n=0; n<num_spheres; n++){
        sprintf(filename, "%s/data/sphere_shape.%d.dat", work_dir, n);
        stream=fopen(filename, "w");
        fprintf(stream, "#nstep    Rx2          Ry2          Rz2      Ia     Ib     Ic\n");
        fclose(stream);
      }
    }

    for(i=0; i<WRITE_BUFFER; i++) {
      buff_avgdisp2[i]=0.0;
      buff_avgdr2[i]=0.0;
      buff_Rg2[i]=0.0;
      buff_sp_len[i]=0.0;
      buff_asphr[i]=0.0;
      buff_acldr[i]=0.0;
      for(d=0; d<DIMS; d++)
	buff_stretch[i][d]=0.0;

      for(k=0; k<MAX_DP; k++) {
        for(j=0; j<14; j++)
          buff_sphereprops[i][j][k]=0.0;
        for(j=0; j<6; j++)
          buff_sphereshape[i][j][k]=0.0;
        for(j=0; j<2; j++)
          buff_spherevalues[i][j][k]=0.0;
        for(j=0; j<4; j++)
          buff_spherevolume[i][j][k]=0.0;
      }  
    }
  }
 
  num_step = n_step / sphere_pm->write_time;
  
  if(n_step % sphere_pm->write_time == 0) {
    sphere_props(sphere_pm, polym_pm, spheres, monomers, faces);
    avg_Rg2 = 0.0;
    avg_stretchx = 0.0;
    avg_stretchy = 0.0;
    avg_stretchz = 0.0;
    avg_disp2 = 0.0;
    avg_dr2 = 0.0;
    avg_sp_len = 0.0;
    avg_asphr = 0.0;
    avg_acldr = 0.0;
  
    for(n=0; n<num_DP; n++) {
      avg_disp2 += spheres[n].disp2;
      avg_Rg2   += spheres[n].Rg2;
      avg_stretchx += spheres[n].stretch[0];
      avg_stretchy += spheres[n].stretch[1];
      avg_stretchz += spheres[n].stretch[2];
      avg_dr2 += spheres[n].dr2;
      avg_sp_len += spheres[n].avg_sprng_len;
      avg_asphr += spheres[n].asphericity;
      avg_acldr += spheres[n].acylindricity;

      buff_sphereprops[num_step%WRITE_BUFFER][0][n] = spheres[n].com[0];
      buff_sphereprops[num_step%WRITE_BUFFER][1][n] = spheres[n].com[1];
      buff_sphereprops[num_step%WRITE_BUFFER][2][n] = spheres[n].com[2];
      buff_sphereprops[num_step%WRITE_BUFFER][3][n] = spheres[n].Rg2;
      buff_sphereprops[num_step%WRITE_BUFFER][4][n] = spheres[n].avg_sprng_len;
      buff_sphereprops[num_step%WRITE_BUFFER][5][n] = spheres[n].disp2;
      buff_sphereprops[num_step%WRITE_BUFFER][6][n] = spheres[n].dr2;
      buff_sphereprops[num_step%WRITE_BUFFER][7][n] = spheres[n].stretch[0];
      buff_sphereprops[num_step%WRITE_BUFFER][8][n] = spheres[n].stretch[1];
      buff_sphereprops[num_step%WRITE_BUFFER][9][n] = spheres[n].stretch[2];
      buff_sphereprops[num_step%WRITE_BUFFER][10][n] = spheres[n].volume;
      buff_sphereprops[num_step%WRITE_BUFFER][11][n] = spheres[n].area;
      buff_sphereprops[num_step%WRITE_BUFFER][12][n] = spheres[n].theta;
      buff_sphereprops[num_step%WRITE_BUFFER][13][n] = spheres[n].omega[2];

      buff_sphereshape[num_step%WRITE_BUFFER][0][n] = spheres[n].Rx2;
      buff_sphereshape[num_step%WRITE_BUFFER][1][n] = spheres[n].Ry2;
      buff_sphereshape[num_step%WRITE_BUFFER][2][n] = spheres[n].Rz2;
      buff_sphereshape[num_step%WRITE_BUFFER][3][n] = spheres[n].Ia;
      buff_sphereshape[num_step%WRITE_BUFFER][4][n] = spheres[n].Ib;
      buff_sphereshape[num_step%WRITE_BUFFER][5][n] = spheres[n].Ic;
    }
    
    avg_disp2 /= num_DP;
    avg_Rg2   /= num_DP;
    avg_stretchx /= num_DP;
    avg_stretchy /= num_DP;
    avg_stretchz /= num_DP;
    avg_dr2 /= num_DP;
    avg_sp_len /= num_DP;
    avg_asphr /= num_DP;
    avg_acldr /= num_DP;

    buff_nstep[num_step%WRITE_BUFFER]=(int)((n_step+addsteps)*sphere_pm->MD_steps*sphere_pm->dt);
    buff_avgdisp2[num_step%WRITE_BUFFER]=avg_disp2;
    buff_avgdr2[num_step%WRITE_BUFFER]=avg_dr2;
    buff_Rg2[num_step%WRITE_BUFFER]=avg_Rg2;
    buff_stretch[num_step%WRITE_BUFFER][0]=avg_stretchx;
    buff_stretch[num_step%WRITE_BUFFER][1]=avg_stretchy;
    buff_stretch[num_step%WRITE_BUFFER][2]=avg_stretchz;
    buff_sp_len[num_step%WRITE_BUFFER]=avg_sp_len;
    buff_asphr[num_step%WRITE_BUFFER]=avg_asphr;
    buff_acldr[num_step%WRITE_BUFFER]=avg_acldr;

    if((num_step+1)%WRITE_BUFFER == 0) {
      sprintf(filename, "%s/data/avg_props.dat", work_dir);
      stream = fopen(filename, "a");
      for(i=0; i<WRITE_BUFFER; i++)
        fprintf(stream, "%d %le %le %le %le %le %le %le %le %le\n", buff_nstep[i], buff_avgdisp2[i], buff_avgdr2[i], buff_Rg2[i], buff_stretch[i][0], buff_stretch[i][1], buff_stretch[i][2], buff_sp_len[i], buff_asphr[i], buff_acldr[i]);
      fclose(stream);

      for(n=0; n<num_spheres+num_poly; n++) {
	for(j=0; j<WRITE_BUFFER; j++) {
	  sprintf(filename, "%s/data/DP_props.%d.dat", work_dir, n);
	  stream=fopen(filename, "a");
	  fprintf(stream, "%d ", buff_nstep[j]);
	  for(i=0; i<14; i++)
	    fprintf(stream, "%le ", buff_sphereprops[j][i][n]);
	  fprintf(stream, "\n");
	  fclose(stream);
	}
      }

      for(n=0; n<num_spheres; n++) {
	for(j=0; j<WRITE_BUFFER; j++) {
	  sprintf(filename, "%s/data/DP_shape.%d.dat", work_dir, n);
	  stream=fopen(filename, "a");
	  fprintf(stream, "%d ", buff_nstep[j]);
	  for(i=0; i<6; i++)
	    fprintf(stream, "%le ", buff_sphereshape[j][i][n]);
	  fprintf(stream, "\n");
	  fclose(stream);
	}
      }
    }
  }
}

void Write_Sphereconfig(struct DP *spheres, struct monomer *monomers, struct face *faces, struct DP_param *sphere_pm,  int n, char *work_dir)
{
  extern int max_x, max_y, max_z;
  int i,j,k,type;
  int nbead[NTYPES], nface[NTYPES];
  int bead_0, bead_f, face_0, face_f;
  int nface_temp;
  int spherenum;
  double maxsize[DIMS];
  char filename[MAX_NL];

  FILE *outfile;

  maxsize[0] = max_x;
  maxsize[1] = max_y;
  maxsize[2] = max_z;

  for (type = 0; type < NTYPES; type++) {
    nbead[type] = sphere_pm->N_per_DP[type]*sphere_pm->Ntype[type];
    nface[type] = sphere_pm->face_per_DP[type]*sphere_pm->Ntype[type];

    if(nbead[type] > 0) {
   
      if(type ==0) {
	bead_0 = 0;
	face_0 = 0;
      
	bead_f = nbead[0];
	face_f = nface[0];
      }
      else {
	bead_0 = nbead[type-1];
	face_0 = nface[type-1];

	bead_f = bead_0+nbead[type];
	face_f = face_0+nface[type];
      }

      if(n < 10)
	sprintf(filename, "%s/data/bond%d_t000%d.vtk", work_dir, type, n);
      else if(n < 100)
	sprintf(filename, "%s/data/bond%d_t00%d.vtk", work_dir, type, n);
      else if(n < 1000)
	sprintf(filename, "%s/data/bond%d_t0%d.vtk", work_dir, type, n);
      else if(n >= 1000)
	sprintf(filename, "%s/data/bond%d_t%d.vtk", work_dir, type, n);

      outfile = fopen(filename, "w");
      fprintf(outfile, "# vtk DataFile Version 2.3   \n");
      fprintf(outfile, "title red blood cell data %d \n", n);
      fprintf(outfile, "ASCII                        \n\n");
      fprintf(outfile, "DATASET UNSTRUCTURED_GRID \n");
      fprintf(outfile, "POINTS %d float\n", nbead[type]);

      for(i=bead_0; i<bead_f; i++)
	fprintf(outfile, "%8.6lf %8.6lf %8.6lf\n", monomers[i].pos[0], monomers[i].pos[1], monomers[i].pos[2]);

      /* check if the face spans across periodic boundaries */
      nface_temp = 0;
      for(i=face_0; i<face_f; i++) {
	if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[1]].pos[0]) > 0.0) 
	  if(maxsize[0]/2-fabs(monomers[faces[i].vertices[1]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
	    if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
	      if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[1]].pos[1]) > 0.0) 
		if(maxsize[1]/2-fabs(monomers[faces[i].vertices[1]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0) 
		  if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0)       
		    if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[1]].pos[2]) > 0.0) 
		      if(maxsize[2]/2-fabs(monomers[faces[i].vertices[1]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
			if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
			  nface_temp++;
      }

      fprintf(outfile, "CELLS %d %d\n", nface_temp, nface_temp*4);
    
      for(i=face_0; i<face_f; i++) 
	if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[1]].pos[0]) > 0.0) 
	  if(maxsize[0]/2-fabs(monomers[faces[i].vertices[1]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
	    if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
	      if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[1]].pos[1]) > 0.0) 
		if(maxsize[1]/2-fabs(monomers[faces[i].vertices[1]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0) 
		  if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0)       
		    if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[1]].pos[2]) > 0.0) 
		      if(maxsize[2]/2-fabs(monomers[faces[i].vertices[1]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
			if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
			  fprintf(outfile, "3 %d %d %d\n", faces[i].vertices[0]-bead_0, faces[i].vertices[1]-bead_0, faces[i].vertices[2]-bead_0);
      
      /* Cell type 5 is triangle */
      fprintf(outfile, "CELL_TYPES %d\n", nface_temp);
      for(i=0; i<nface_temp/2; i++)
	fprintf(outfile, "5\n");
      for(i=nface_temp/2; i<nface_temp; i++)
	fprintf(outfile, "5\n");

      /* Add color for stress */
      fprintf(outfile, "POINT_DATA %d\n", bead_f-bead_0); 
      fprintf(outfile, "SCALARS color float 1\n");
      fprintf(outfile, "LOOKUP_TABLE default\n");
      for(i=bead_0; i<bead_f; i+=4) 
	fprintf(outfile, "%ld %ld %ld %ld\n", monomers[i].DP_id, monomers[i+1].DP_id, monomers[i+2].DP_id, monomers[i+3].DP_id);
	//	fprintf(outfile, "%le %le %le %le\n", monomers[i].stress[0][1], monomers[i+1].stress[0][1], monomers[i+2].stress[0][1], monomers[i+3].stress[0][1]);
      fclose(outfile);
    
    
      /* write out in Tecplot format */
      sprintf(filename, "%s/data/bond%d_tec.dat", work_dir, type);
      
      nface_temp = 0;
      for(i=face_0; i<face_f; i++) {
	if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[1]].pos[0]) > 0.0) 
	  if(maxsize[0]/2-fabs(monomers[faces[i].vertices[1]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
	    if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
	      if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[1]].pos[1]) > 0.0) 
		if(maxsize[1]/2-fabs(monomers[faces[i].vertices[1]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0) 
		  if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0)       
		    if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[1]].pos[2]) > 0.0) 
		      if(maxsize[2]/2-fabs(monomers[faces[i].vertices[1]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
			if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
			  nface_temp++;
      }
      
      outfile = fopen(filename, "a");
      fprintf(outfile, "VARIABLES = \"x\", \"y\", \"z\"  \n");
      fprintf(outfile, "ZONE T = \"TRIANGLES\", N=%d, E=%d, F=FEPOINT, ET=TRIANGLE \n", nbead[type], nface_temp);
      
      for(i=bead_0; i<bead_f; i++)
	fprintf(outfile, "%8.6lf %8.6lf %8.6lf\n", monomers[i].pos[0], monomers[i].pos[1], monomers[i].pos[2]);
      
      for(i=face_0; i<face_f; i++) 
	if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[1]].pos[0]) > 0.0) 
	  if(maxsize[0]/2-fabs(monomers[faces[i].vertices[1]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
	    if(maxsize[0]/2-fabs(monomers[faces[i].vertices[0]].pos[0] - monomers[faces[i].vertices[2]].pos[0]) > 0.0) 
	      if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[1]].pos[1]) > 0.0) 
		if(maxsize[1]/2-fabs(monomers[faces[i].vertices[1]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0) 
		  if(maxsize[1]/2-fabs(monomers[faces[i].vertices[0]].pos[1] - monomers[faces[i].vertices[2]].pos[1]) > 0.0)       
		    if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[1]].pos[2]) > 0.0) 
		      if(maxsize[2]/2-fabs(monomers[faces[i].vertices[1]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
			if(maxsize[2]/2-fabs(monomers[faces[i].vertices[0]].pos[2] - monomers[faces[i].vertices[2]].pos[2]) > 0.0) 
			  fprintf(outfile, "%d %d %d\n", faces[i].vertices[0]-bead_0+1, faces[i].vertices[1]-bead_0+1, faces[i].vertices[2]-bead_0+1);
	    
      fclose(outfile);
    }
  }
}

void Write_Polymconfig(struct DP *DPs, struct monomer *monomers, struct DP_param *polym_pm, int n, char *work_dir, int n0)
{
  extern int max_x, max_y, max_z;
  int i,j,d;
  int type;
  int nbead[2], nchain[2], numbeads[2];
  int bead0, beadf;
  int outln_ct;
  int *ln_flag;
  double maxsize[DIMS];
  double dx[DIMS];
  char filename[MAX_NL];
  FILE *outfile;

  maxsize[0] = max_x;
  maxsize[1] = max_y;
  maxsize[2] = max_z;

  outln_ct = 0;

  for (type = 0; type < NTYPES; type++) {
    nchain[type] = polym_pm->Ntype[type];
    nbead[type] = polym_pm->N_per_DP[type];
    numbeads[type] = nchain[type]*nbead[type];

    if(nchain[type] > 0) {
      
      if(type ==0) {
	bead0 = n0;
	beadf = n0+numbeads[0];
      }
      else {
	bead0 += numbeads[type-1];
	beadf = bead0+numbeads[type];
      }
    
      ln_flag = (int *)malloc(numbeads[type]*sizeof(int));
      if(ln_flag == 0) fatal_err("cannot allocate ln_flag in output", -1);
      for(j=0; j<numbeads[type]; j++)
	ln_flag[j] = 0;
      
      if(n < 10)
	sprintf(filename, "%s/data/bondp%d_t000%d.vtk", work_dir, type, n);
      else if(n < 100)
	sprintf(filename, "%s/data/bondp%d_t00%d.vtk", work_dir, type, n);
      else if(n < 1000)
	sprintf(filename, "%s/data/bondp%d_t0%d.vtk", work_dir, type, n);
      else if(n >= 1000)
	sprintf(filename, "%s/data/bondp%d_t%d.vtk", work_dir, type, n);
     
      outfile = fopen(filename, "w");
      fprintf(outfile, "# vtk DataFile Version 2.3   \n");
      fprintf(outfile, "title polymer chain data %d \n", n);
      fprintf(outfile, "ASCII                        \n\n");
      fprintf(outfile, "DATASET POLYDATA\n");
      fprintf(outfile, "POINTS %d float\n", numbeads[type]);

      outln_ct = 0;
      for(j=bead0; j<beadf; j++) {
	fprintf(outfile, "%8.6lf %8.6lf %8.6lf\n", monomers[j].pos[0], monomers[j].pos[1], monomers[j].pos[2]);
	 
	/* if the bond is out of boundary, do not write */
	if(((j-bead0)%nbead[type])<nbead[type]-1) {
	  for(d=0; d<DIMS; d++) {
	    dx[d] = monomers[j+1].pos[d]-monomers[j].pos[d];
	    if(dx[d] > (maxsize[d]/2.0) || dx[d] < -(maxsize[d]/2.0)) {
	      outln_ct++;
	      ln_flag[j-bead0] = 1;
	      break;
	    }
	  }
	}
      }
    
      fprintf(outfile, "LINES %d %d\n", nchain[type]*(nbead[type]-1)-outln_ct, (nchain[type]*(nbead[type]-1)-outln_ct)*3);
     
      for(i=0; i<nchain[type]; i++)
	for(j=0; j<nbead[type]-1; j++) {
	  beadf = i*nbead[type]+j;
	  if(ln_flag[beadf] == 0)
	    fprintf(outfile, "2 %d %d\n", beadf, beadf+1);
	}
     
      fclose(outfile);
      free(ln_flag);
    }
  }
}

void Write_Cylconfig(struct DP *cylinders, struct monomer *monomers, struct DP_param *cyl_pm,  int n, char *work_dir, int n0)
{
  extern int max_x, max_y, max_z;
  int i,j,k,type;
  int nbead, height;
  int bead_0, bead_f, face_0, face_f, Nperring;
  int nface_temp;
  int cylnum;
  double maxsize[DIMS];
  char filename[MAX_NL];
  FILE *outfile;

  maxsize[0] = max_x;
  maxsize[1] = max_y;
  maxsize[2] = max_z;

  if(cyl_pm->initconfig == 1) {
    Nperring = (int)(2.0*M_PI*(cyl_pm->radius+cyl_pm->ramp));
    height = max_x;
  }
  else if(cyl_pm->initconfig == 2) {
    Nperring = cyl_pm->N_per_ring;
    height = cyl_pm->height;
  }
  else if(cyl_pm->initconfig == 3) {
    Nperring = cyl_pm->N_per_ring;
    height = cyl_pm->height;
  }
  else if(cyl_pm->initconfig == 5) {
    Nperring = cyl_pm->N_per_ring;
    height = cyl_pm->height;
  }
  else if(cyl_pm->initconfig == 6) {
    Nperring = cyl_pm->N_per_ring;
    height = cyl_pm->height;
  }
  bead_0 = n0;
  bead_f = bead_0+cyl_pm->num_beads;
  face_0 = 0;
  face_f = (height-1) * Nperring;

  if(n < 10)
    sprintf(filename, "%s/data/bondcyl_t000%d.vtk", work_dir, n);
  else if(n < 100)
    sprintf(filename, "%s/data/bondcyl_t00%d.vtk", work_dir, n);
  else if(n < 1000)
    sprintf(filename, "%s/data/bondcyl_t0%d.vtk", work_dir, n);
  else if(n >= 1000)
    sprintf(filename, "%s/data/bondcyl_t%d.vtk", work_dir, n);

  outfile = fopen(filename, "w");
  fprintf(outfile, "# vtk DataFile Version 2.3   \n");
  fprintf(outfile, "title tube data %d \n", n);
  fprintf(outfile, "ASCII                        \n");
  fprintf(outfile, "DATASET UNSTRUCTURED_GRID \n");
  fprintf(outfile, "POINTS %d float\n", cyl_pm->num_beads);

  for(i=bead_0; i<bead_f; i++)
    fprintf(outfile, "%8.6lf %8.6lf %8.6lf\n", monomers[i].pos[0], monomers[i].pos[1], monomers[i].pos[2]);

  /* check if the face spans across periodic boundaries */
  nface_temp = face_f-face_0;
  fprintf(outfile, "CELLS %d %d\n", nface_temp, nface_temp*5);
    
  for(i=face_0; i<face_f; i++) {
    if((i+1) % Nperring == 0)
      fprintf(outfile, "4 %d %d %d %d\n", i, i+1-Nperring, i+Nperring, i+1);
    else
      fprintf(outfile, "4 %d %d %d %d\n", i, i+1, i+Nperring, i+1+Nperring);
  }      
  fprintf(outfile, "CELL_TYPES %d\n", nface_temp);
  for(i=0; i<nface_temp; i++)
    fprintf(outfile, "5\n");
      
  fclose(outfile);
}

void Write_Velslice(struct monomer *monomers, Float ***velcs_df, int **node_map, int nstep, char *work_dir)
{
  /* Write out the velocity slicing through the monomer position */
  extern int max_x, max_y, max_z;
  int i,j,k,p,q, nxy;
  double lvx, lvy, lvz;
  double prop[DIMS+1];
  char filename[MAX_NL];
  FILE *stream;

  if(nstep < 10)
    sprintf(filename, "%s/data/velslice.000%d.dat", work_dir, nstep);
  else if(nstep < 100)
    sprintf(filename, "%s/data/velslice.00%d.dat", work_dir, nstep);
  else if(nstep < 1000)
    sprintf(filename, "%s/data/velslice.0%d.dat", work_dir, nstep);
  else if(nstep >= 1000)
    sprintf(filename, "%s/data/velslice.%d.dat", work_dir, nstep);

  stream = fopen(filename, "w");

    
  //i = floor(monomers[0].pos[0]); 
  //k = floor(monomers[0].pos[2]); 
  

  i = (max_x+1)/2;
  k = (max_z+1)/2;
  // j = (max_y+1)/2;

  fprintf(stream, "#monpos (%d, j, %d)\n", i, k);
  for(j=1; j <= max_y; j++) {
    
    nxy = i*(max_y+2) + j;
	
    for(p=0; p<DIMS+1; p++)
      prop[p]=0.0;
    
    if (node_map[nxy][k] != 1) {         /* Fluid node */
      lvx =0.0;
      lvy =0.0;
      lvz =0.0;
      for(q=0; q < Num_Dir; q++) {
	prop[0] += evector[0][q]*velcs_df[nxy][q][k];
	lvx += c_x[q]*velcs_df[nxy][q][k];
	lvy += c_y[q]*velcs_df[nxy][q][k];
	lvz += c_z[q]*velcs_df[nxy][q][k];
      }
      
      /*      prop[0] += lvx*lvx+lvy*lvy+lvz*lvz; */
      prop[1] = lvx /Rho_Fl;
      prop[2] = lvy /Rho_Fl;
      prop[3] = lvz /Rho_Fl;
      
      fprintf(stream, "%d %le %le %le %le \n", j, prop[0], prop[1], prop[2], prop[3]);
    }
  }
  
  fclose(stream);
     
}

void Write_stat(Float ***velcs_df, int nstep, double kT, char *work_dir)
{
  extern int max_x, max_y, max_z;
  extern int num_x;
  extern Float std_xi[Num_Dir];

  int i,j,k,q,p;
  int stat_obs[101][101][101];
  long nxy;
  char filename[MAX_NL];
  double prop[Num_Dir];  
  static double stat_mom[101][101][101][Num_Dir];
  static double stat_momsq[101][101][101][Num_Dir];
  double ave_M, ave_sq_M, var_M;

  FILE *stream;

  // Statistics of the moments 

  /*  
      if(nstep == 100) {
      for(i=1; i<=num_x; i++)
      for(j=0; j<max_y; j++)
      for(k=0; k<max_z; k++) {
      stat_obs[i][j][k] = 0;
      
      for(q=0; q<10; q++) {
      stat_mom[i][j][k][q]=0.0;
      stat_momsq[i][j][k][q]=0.0;
      }
      }
      }

      if(nstep >=100) {
      for(i=1; i<=num_x; i++) 
      for(j=0; j<max_y; j++) {
      nxy = i*max_y + j;
      
      for(k=0; k<max_z; k++) {
	for(q=0; q < 10; q++)
	  prop[q]=0.0;
	
	for(q=0; q < Num_Dir; q++) {
	  prop[0] += evector[0][q]*velcs_df[nxy][q][k];
	  prop[1] += evector[1][q]*velcs_df[nxy][q][k];
	  prop[2] += evector[2][q]*velcs_df[nxy][q][k];
	  prop[3] += evector[3][q]*velcs_df[nxy][q][k];
	  prop[4] += evector[4][q]*velcs_df[nxy][q][k];
	  prop[5] += evector[5][q]*velcs_df[nxy][q][k];
	  prop[6] += evector[6][q]*velcs_df[nxy][q][k];
	  prop[7] += evector[7][q]*velcs_df[nxy][q][k];
	  prop[8] += evector[8][q]*velcs_df[nxy][q][k];
	  prop[9] += evector[9][q]*velcs_df[nxy][q][k];
	  
	}
	
	for(q=0; q < 10; q++) {
	  stat_mom[i][j][k][q] +=prop[q];
	  stat_momsq[i][j][k][q] +=prop[q]*prop[q];
	}
	
	stat_obs[i][j][k]++;
	}
	}
	}

  if(nstep == 1500) {
    for(q=0; q < 10; q++) {
      sprintf(filename, "%s/data/moment_%ld.dat", work_dir, q);

      stream = fopen(filename, "w");
      for(i=1; i<=num_x; i++) 
	for(j=0; j<max_y; j++) 
	  for(k=0; k<max_z; k++) {
	    ave_M = stat_mom[i][j][k][q]/(stat_obs[i][j][k]);
	    ave_sq_M = stat_momsq[i][j][k][q]/(stat_obs[i][j][k]);
	    var_M = ave_sq_M - ave_M*ave_M;
	    
	    if(q > 3)
	      fprintf(stream, "%d %d %d\t%le\t%le\t%le\t%le\n", i,j,k, ave_M, ave_sq_M, var_M, var_M/(std_xi[q]*std_xi[q]));
	    if(q <= 3)
	      fprintf(stream, "%d %d %d\t%le\t%le\t%le\n", i,j,k, ave_M, ave_sq_M, var_M);
	      
	  }
      
      fclose(stream);
    }
  }

  */

  for(i=1; i<=num_x; i++)
    for(j=1; j<=max_y; j++)
      for(k=1; k<=max_z; k++) 
	for(q=0; q<Num_Dir; q++) {
	  stat_mom[i][j][k][q]=0.0;
	  stat_momsq[i][j][k][q]=0.0;
	}
  
  for(i=1; i<=num_x; i++) 
    for(j=1; j<=max_y; j++) {
      nxy = i*(max_y+2) + j;
      
      for(k=1; k<=max_z; k++) {
	for(q=0; q < Num_Dir; q++)
	  prop[q]=0.0;

	for(p=0; p < Num_Dir; p++)
	  for(q=0; q < Num_Dir; q++) 
	    prop[p] += evector[p][q]*velcs_df[nxy][q][k];
	
	for(q=0; q < Num_Dir; q++) {
	  stat_mom[i][j][k][q] =prop[q];
	  stat_momsq[i][j][k][q] =prop[q]*prop[q];
	}
      }
    }

  for(q=0; q < Num_Dir; q++) {
    ave_M = 0.0;
    ave_sq_M = 0.0;
    for(i=1; i<=num_x; i++) 
      for(j=1; j<=max_y; j++) 
	for(k=1; k<=max_z; k++) {
	  ave_M += stat_mom[i][j][k][q];
	  ave_sq_M += stat_momsq[i][j][k][q];
	}
    
    ave_M /= (num_x*max_y*max_z);     	  
    ave_sq_M /= (num_x*max_y*max_z); 
    var_M = ave_sq_M - ave_M*ave_M;
    
    sprintf(filename, "%s/data/moment_%ld.dat", work_dir, q);
    stream = fopen(filename, "a");
    
    if(q > 3)
      fprintf(stream, "%d\t%le\t%le\t%le\n", nstep, ave_M, ave_sq_M, var_M/(std_xi[q]*std_xi[q]));
    if(q <= 3)
      fprintf(stream, "%d\t%le\t%le\t%le\n", nstep, ave_M, ave_sq_M/Rho_Fl, var_M/Rho_Fl);
    
    fclose(stream);
  }
}

void check_vel_conservation(Float ***velcs_df, int nstep)
{
  extern int max_x, max_y, max_z;
  int i, j, k, d, q;
  int nxy;
  Float total_deltaf[DIMS];

  i=(1+max_x)/2;
  j=(1+max_y)/2;
  k=(1+max_z)/2;

  for(d=0; d<DIMS; d++) 
    total_deltaf[d]=0.0;
  
  nxy = i*max_y+j;

  for(d=0; d<DIMS; d++) {
    for(q=0; q<Num_Dir; q++) {
      total_deltaf[0] +=velcs_df[nxy][q][k]*c_x[q];
      total_deltaf[1] +=velcs_df[nxy][q][k]*c_y[q];
      total_deltaf[2] +=velcs_df[nxy][q][k]*c_z[q];
    }
  }

  for(d=0; d<DIMS; d++) {
    if(fabs(total_deltaf[d]) > 1e-12) {
      printf(" step %d check_vel  momentum not conserved at (%d %d %d) (%le %le %le)!\n", nstep, i, j, k, total_deltaf[0], total_deltaf[1], total_deltaf[2]);
      break;    
    }
  }
  
}

void sphere_props(struct DP_param *sphere_pm, struct DP_param *polym_pm, struct DP *DPs,  struct monomer *mono, struct face *faces)
{
  extern int max_x, max_y, max_z;
  int i, j, k, d, n1, n2, n3, d1, d2;
  int numbeads = sphere_pm->num_beads;
  int Nsphere = sphere_pm->NDP;
  int Npoly = polym_pm->NDP;
  int NDP = Nsphere + Npoly;
  int Ntype[NTYPES];
  int Nbead[NTYPES];
  int maxsize[3];
  int nspring, nface;
  int start, type;
  int min;
  double sprng_len, avg_sprng_len;
  double temp, tempstretch;
  double g_tensor[DIMS][DIMS];  /* gyration tensor */
  double I_tensor[DIMS][DIMS];  /* inertia tensor */
  double dr[DIMS], lambda[DIMS], Im[DIMS];
  double q12[DIMS], q13[DIMS], normal[DIMS];
  double y_com;
  double Rx2, Ry2, Rz2;
  double costheta, sintheta;
  double omega[DIMS];
  double y[162];
  double v[DIMS];  
  FILE *stream;

  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;

  /* calculate the sphere center-of-mass */
  for(i=0; i<NDP; i++) {
    if(i < Nsphere) {
      for (j=0 ; j<NTYPES ; j++) {
	Ntype[j] = sphere_pm->Ntype[j];
	Nbead[j] = sphere_pm->N_per_DP[j];
      }
      type = (i<Ntype[0] ? 0 : 1);
      start = (type==0 ? i*Nbead[0] : Ntype[0]*Nbead[0] + (i-Ntype[0])*Nbead[1]);
    }
    else {
      for (j=0 ; j<NTYPES ; j++) {
	Ntype[j] = polym_pm->Ntype[j];
	Nbead[j] = polym_pm->N_per_DP[j];
      }
      type = (i<Ntype[0] ? 0 : 1);
      start = (type==0 ? numbeads+(i-Nsphere)*Nbead[0] : numbeads+ Ntype[0]*Nbead[0] + (i-Nsphere-Ntype[0])*Nbead[1]);
    }

    for(d=0; d<DIMS; d++) 
      DPs[i].com[d] = 0.0;
    
    for(j=0; j<Nbead[type]; j++) {
      for(d=0; d<DIMS; d++)
        DPs[i].com[d] += mono[start+j].pos_pbc[d];
    }
  
    DPs[i].disp2 = 0.0;  
    for(d=0; d<DIMS; d++) {
      DPs[i].com[d] = DPs[i].com[d]/Nbead[type];
      DPs[i].disp2 += (DPs[i].com[d]-DPs[i].com0[d])*(DPs[i].com[d]-DPs[i].com0[d]);
      DPs[i].dr2 += (DPs[i].com[d]-DPs[i].comold[d])*(DPs[i].com[d]-DPs[i].comold[d]);
      DPs[i].comold[d]=DPs[i].com[d];
    }
  }
  
  /* calculate the sphere radius of gyration and stretch */
  for(i=0; i<NDP; i++) {
    if(i < Nsphere) {
      for (j=0 ; j<NTYPES ; j++) {
	Ntype[j] = sphere_pm->Ntype[j];
	Nbead[j] = sphere_pm->N_per_DP[j];
      }
      type = (i<Ntype[0] ? 0 : 1);
      start = (type==0 ? i*Nbead[0] : Ntype[0]*Nbead[0] + (i-Ntype[0])*Nbead[1]);
    }
    else {
      for (j=0 ; j<NTYPES ; j++) {
	Ntype[j] = polym_pm->Ntype[j];
	Nbead[j] = polym_pm->N_per_DP[j];
      }
      type = (i<Ntype[0]+Nsphere ? 0 : 1);
      start = (type==0 ? numbeads+(i-Nsphere)*Nbead[0] : numbeads+ Ntype[0]*Nbead[0] + (i-Nsphere-Ntype[0])*Nbead[1]);
    }

    for(d=0; d<DIMS; d++)
    DPs[i].stretch[d] = 0.0;

    //    printf("Test polymer %d com (%le %le %le)\n", i, DPs[i].com[0], DPs[i].com[1], DPs[i].com[2]);

    for (d1=0 ; d1<DIMS ; d1++)
      for (d2=0 ; d2<DIMS ; d2++)
        g_tensor[d1][d2] = 0.0;

    for(j=0; j<Nbead[type]; j++) {

      for(d=0; d<DIMS; d++)
        dr[d] = mono[start+j].pos_pbc[d]-DPs[i].com[d];

      for (d1=0 ; d1<DIMS ; d1++)
        for (d2=0 ; d2<DIMS ; d2++)
          g_tensor[d1][d2] += dr[d1]*dr[d2];

      for(k=0; k<j; k++) {
        tempstretch = 0.0;
        for(d=0; d<DIMS; d++) {
          temp = mono[start+k].pos_pbc[d]-mono[start+j].pos_pbc[d];
          tempstretch = sqrt(temp*temp);
	  if(tempstretch > DPs[i].stretch[d]) 
	    DPs[i].stretch[d] = tempstretch;
        }
      }

      //      printf("Test monomer %d %le %le %le\n", start+j, mono[start+j].pos_pbc[0], mono[start+j].pos_pbc[1], mono[start+j].pos_pbc[2]);
    }

    for (d1=0 ; d1<DIMS ; d1++)
      for (d2=0 ; d2<DIMS ; d2++)
        g_tensor[d1][d2] /= Nbead[type];

    /* convert the gyration tensor to the inertia tensor */
    I_tensor[0][0] = g_tensor[1][1] + g_tensor[2][2];
    I_tensor[1][1] = g_tensor[0][0] + g_tensor[2][2];
    I_tensor[2][2] = g_tensor[0][0] + g_tensor[1][1];
    I_tensor[0][1] = I_tensor[1][0] = -g_tensor[0][1];
    I_tensor[0][2] = I_tensor[2][0] = -g_tensor[0][2];
    I_tensor[1][2] = I_tensor[2][1] = -g_tensor[1][2];
  
    /* the moments of gyration tensor are the eigenvalues of gyration tensor */
    sym33eigen_values(g_tensor, lambda);
    
    /* principal moments of inertia */
    sym33eigen_values(I_tensor, Im);

    /* sort the eigenvalues small to big */
    for (j=0 ; j<DIMS-1 ; j++) {
      min = j;
      for (k=j+1 ; k<DIMS ; k++)
        if (lambda[k] < lambda[min])
          min = k;
      temp = lambda[j];
      lambda[j] = lambda[min];
      lambda[min] = temp;
    }

    for (j=0 ; j<DIMS-1 ; j++) {
      min = j;
      for (k=j+1 ; k<DIMS ; k++)
        if (Im[k] < Im[min])
          min = k;
      temp = Im[j];
      Im[j] = Im[min];
      Im[min] = temp;
    }

    /* take thet2 to be angle between the eigenvector of I and flow */
    eigen_vector(I_tensor, Im[0], v);
    costheta = fabs(v[0]);
    if (costheta > 1.) costheta = 1.;
    if (costheta < -1.) costheta = -1.;
    DPs[i].theta = acos(costheta)*180./M_PI;

    Rx2 = lambda[0];
    Ry2 = lambda[1];
    Rz2 = lambda[2];

    DPs[i].Rg2 = Rx2 + Ry2 + Rz2;
    DPs[i].Rx2 = Rx2;
    DPs[i].Ry2 = Ry2;
    DPs[i].Rz2 = Rz2;
    // DPs[i].asphericity = Rz2 - (Rx2 + Ry2)/2.;
    DPs[i].acylindricity = Ry2 - Rx2;

    DPs[i].asphericity = ((Rx2-Ry2)*(Rx2-Ry2) + (Rx2-Rz2)*(Rx2-Rz2) + (Ry2-Rz2)*(Ry2-Rz2) )
                              / (2.*(Rx2+Ry2+Rz2)*(Rx2+Ry2+Rz2));
    DPs[i].Ia = Im[0];
    DPs[i].Ib = Im[1];
    DPs[i].Ic = Im[2];

    //    printf("Test Rx2 (%le %le %le) Rg %le\n", Rx2, Ry2, Rz2, DPs[i].Rg2);

  /* calculate the average spring length */
    nspring = 0;
    avg_sprng_len = 0.0;
    
    for (j=0 ; j<Nbead[type] ; j++)
      for (k=1 ; k<=mono[start+j].blist[0][0] ; k++) {
	n1 = start+j;
	n2 = mono[n1].blist[k][0];
	
	sprng_len = 0.0;
	for(d=0; d<DIMS; d++) {
	  temp = mono[n1].pos_pbc[d]-mono[n2].pos_pbc[d];
	  sprng_len += temp*temp;
	}
	sprng_len = sqrt(sprng_len);
	
	avg_sprng_len += sprng_len;
	nspring ++;
      }
    avg_sprng_len /= nspring;
    DPs[i].avg_sprng_len = avg_sprng_len;
  }
  
  /* calculate average angular velocity */
  for (i=0 ; i<NDP ; i++) {
    if(i < Nsphere) {
      for (j=0 ; j<NTYPES ; j++) {
	Ntype[j] = sphere_pm->Ntype[j];
	Nbead[j] = sphere_pm->N_per_DP[j];
      }
      type = (i<Ntype[0] ? 0 : 1);
      start = (type==0 ? i*Nbead[0] : Ntype[0]*Nbead[0] + (i-Ntype[0])*Nbead[1]);
    }
    else {
      for (j=0 ; j<NTYPES ; j++) {
	Ntype[j] = polym_pm->Ntype[j];
	Nbead[j] = polym_pm->N_per_DP[j];
      }
      type = (i<Ntype[0]+Nsphere ? 0 : 1);
      start = (type==0 ? numbeads+(i-Nsphere)*Nbead[0] : numbeads+ Ntype[0]*Nbead[0] + (i-Nsphere-Ntype[0])*Nbead[1]);
    }

    for (d=0 ; d<DIMS ; d++)
      DPs[i].omega[d] = 0.;
    for (j=0 ; j<Nbead[type] ; j++) {
      temp = 0.;
      for (d=0 ; d<DIMS ; d++) {
        dr[d] = mono[start+j].pos_pbc[d] - DPs[i].com[d];
        temp += dr[d]*dr[d];
      }
      product(dr, mono[start+j].vel, omega);
      for (d=0 ; d<DIMS ; d++)
        DPs[i].omega[d] += omega[d]/temp;
    }
    for (d=0 ; d<DIMS ; d++)
      DPs[i].omega[d] /= Nbead[type];
    
    DPs[i].area = 0.;
    DPs[i].volume =0.;
    DPs[i].volume_dupin =0.;
  }
  
  /* calculate surface area and volume */
  for(i=0; i<Nsphere; i++) {
    type = (i<Ntype[0] ? 0 : 1);
    start = (type==0 ? i*sphere_pm->face_per_DP[0] : Ntype[0]*sphere_pm->face_per_DP[0] + (i-Ntype[0])*sphere_pm->face_per_DP[1]);

    DPs[i].area = 0.;
    DPs[i].volume =0.;
    DPs[i].volume_dupin =0.;

    /* loop over all faces */
    /* n1, n2, n3 are the three points of the triangle; q12, q13 are two sides; dr is vector from com to n1 */
    
    for(j=0; j<sphere_pm->face_per_DP[type]; j++) {
      nface = start + j;
      n1 = faces[nface].vertices[0];
      for (d=0 ; d<DIMS ; d++)
        dr[d] = mono[n1].pos_pbc[d] - DPs[i].com[d];
      n2 = faces[nface].vertices[1];
      for (d=0 ; d<DIMS ; d++)
	q12[d] = mono[n2].pos_pbc[d] - mono[n1].pos_pbc[d];
      n3 = faces[nface].vertices[2];
      for (d=0 ; d<DIMS ; d++)
	q13[d] = mono[n3].pos_pbc[d] - mono[n1].pos_pbc[d];

      product(q12, q13, normal);
      faces[nface].area = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2])/2.;
      DPs[i].area += faces[nface].area;
      DPs[i].volume += fabs(dr[0]*normal[0]+dr[1]*normal[1]+dr[2]*normal[2])/6.;
	    
      y_com = (mono[n1].pos_pbc[1]+mono[n2].pos_pbc[1]+mono[n3].pos_pbc[1])/3.;
      /* make the y-component of the normal vector point outwards */
      if(y_com - DPs[i].com[1] >= 0.0) 
	DPs[i].volume_dupin += fabs(normal[1]*y_com/2.);  
      else
	DPs[i].volume_dupin -= fabs(normal[1]*y_com/2.);  
    }
  }
}

/* give eigenvales of a 3x3 symmetric matrix. algorithm taken from Smith (1961) */
void sym33eigen_values(double mx[3][3], double lambda[3])
{
  double m, p, q, phi, temp;
  int i, j;

  /* trace / 3 */
  m = (mx[0][0] + mx[1][1] + mx[2][2])/3.;

  mx[0][0] -= m;
  mx[1][1] -= m;
  mx[2][2] -= m;

  p = 0.;
  /* sum the square of elements */
  for (i=0 ; i<3 ; i++)
    for (j=0 ; j<3 ; j++)
      p+= mx[i][j]*mx[i][j];
  p /= 6.;

  /* determinant */
  q = mx[0][0]*mx[1][1]*mx[2][2] + mx[0][1]*mx[1][2]*mx[2][0] + mx[0][2]*mx[1][0]*mx[2][1]
    - (mx[0][0]*mx[1][2]*mx[2][1] + mx[0][1]*mx[1][0]*mx[2][2] + mx[0][2]*mx[1][1]*mx[2][0]);
  q /= 2.;

  temp = p*p*p - q*q;
  if (temp < 0.) {
    printf("rounding %le to 0!\n", temp);
    temp = 0.;
  }

  phi = atan(sqrt(temp)/q)/3.;

  lambda[0] = m + 2.*sqrt(p)*cos(phi);
  lambda[1] = m - sqrt(p)*(cos(phi)+sqrt(3.)*sin(phi));
  lambda[2] = m - sqrt(p)*(cos(phi)-sqrt(3.)*sin(phi));
}

void eigen_vector(double mx[3][3], double lambda, double v[3])
{
  int i;
  double a, b, c, d, delta, p, q, temp;

  a = mx[0][0] - lambda;
  b = mx[0][1];
  c = mx[1][0];
  d = mx[1][1] - lambda;
  delta = a*d - b*c;
  p = mx[0][2];
  q = mx[1][2];

  v[0] = (b*q-d*p)/delta;
  v[1] = (c*p-a*q)/delta;
  v[2] = 1.;

  temp = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  for (i=0;i<3;i++) v[i] /= temp;
}

/* Output velocity field in the xy-plane in vtk format */
void write_velocity_field(int istep, Float ***velcs_df, int **node_map, char *work_dir)
{
  char PREFIX[20]="u";
  char filename[MAX_NL];
  FILE *outfile;

  int i,j,k,l,d;
  int q, nxy;
  double v[3];
  
  extern int max_x;
  extern int max_y;
  extern int max_z;

  if(istep < 10)
    sprintf(filename, "%s/data/%s.000%d.vtk", work_dir, PREFIX, istep);
  else if(istep < 100)
    sprintf(filename, "%s/data/%s.00%d.vtk", work_dir, PREFIX, istep);
  else if(istep < 1000)
    sprintf(filename, "%s/data/%s.0%d.vtk", work_dir, PREFIX, istep);
  else
    sprintf(filename, "%s/data/%s.%d.vtk", work_dir, PREFIX, istep);

  outfile = fopen(filename, "w");
  
  fprintf(outfile, "# vtk DataFile Version 2.3   \n");
  fprintf(outfile, "title u field(istep=%d)      \n", istep);
  fprintf(outfile, "ASCII                        \n");
  fprintf(outfile, "DATASET STRUCTURED_POINTS    \n");
  fprintf(outfile, "DIMENSIONS %d %d %d \n", max_x+2, max_y+2, max_z+2);
  fprintf(outfile, "ORIGIN  %le %le %le\n", 0.0, 0.0, 0.0);
  fprintf(outfile, "SPACING %le %le %le\n", 1.0, 1.0, 1.0);
  fprintf(outfile, "POINT_DATA %d\n", (max_x+2)*(max_y+2)*(max_z+2));
  fprintf(outfile, "VECTORS velocity double      \n");
  fclose(outfile);

  outfile = fopen(filename, "a");
  //  for(k=1; k<=max_z; k++) {
  //    for(j=1; j<=max_y; j++)
  //      for(i=1; i<=max_x; i++) {

  for(k=0; k<=max_z+1; k++) {
    for(j=0; j<=max_y+1; j++)
      for(i=0; i<=max_x+1; i++) {
		
	nxy = i*(max_y+2)+j;
	
	v[0] = 0.0;
	v[1] = 0.0;
	v[2] = 0.0;
	
	if (node_map[nxy][k]!=1)          /* Fluid node */   
	  for(q=0; q<Num_Dir; q++) {
	    v[0] += c_x[q]*velcs_df[nxy][q][k];
	    v[1] += c_y[q]*velcs_df[nxy][q][k];
	    v[2] += c_z[q]*velcs_df[nxy][q][k];
	  }
	
	fprintf(outfile, "%- 12.10lf\t%- 12.10lf\t%- 12.10lf\n", v[0], v[1], v[2]);
      }
  }
  fclose(outfile);


/*   if(istep < 10) */
/*     sprintf(filename, "%s/%s.000%d.dat", work_dir, PREFIX, istep); */
/*   else if(istep < 100) */
/*     sprintf(filename, "%s/%s.00%d.dat", work_dir, PREFIX, istep); */
/*   else if(istep < 1000) */
/*     sprintf(filename, "%s/%s.0%d.dat", work_dir, PREFIX, istep); */
/*   else */
/*     sprintf(filename, "%s/%s.%d.dat", work_dir, PREFIX, istep); */

/*   outfile = fopen(filename, "w"); */
    
/*   fprintf(outfile," VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"W\"\n"); */
/*   fprintf(outfile, "ZONE I=%ld, J=%ld\n", max_x, max_y); */
/*   k=max_z/2; */
/*   for(j=1; j <= max_y; j++) { */
/*     for(i=1; i<= max_x; i++) { */
/* 	nxy = i*(max_y+2) + j; */
	
/* 	if (node_map[nxy][k]==0) {         /\* Fluid node *\/ */
/* 	  v[0]=0.0; */
/* 	  v[1]=0.0; */
/* 	  v[2]=0.0; */
/* 	  for(q=0; q < Num_Dir; q++) { */
/* 	    v[0]+= c_x[q]*velcs_df[nxy][q][k]; */
/* 	    v[1]+= c_y[q]*velcs_df[nxy][q][k]; */
/* 	    v[2]+= c_z[q]*velcs_df[nxy][q][k]; */
/* 	  } */
	  
/* 	  for(d=0; d<DIMS; d++) */
/* 	    v[d]/=Rho_Fl; */

/* 	  fprintf(outfile, "%d %d %le %le %le \n", i,j, v[0], v[1], v[2]); */
/* 	} */
/*     } */
/*   } */
/*   fclose(outfile); */

  /* write out system fluid density and momentum */
  
/*   if(nstep % WRITE_BUFFER == 0){ */
/*     sprintf(filename, "%s/data/stress.dat", work_dir); */
/*     stream = fopen(filename, "a"); */
    
/*     i = (max_x+1)/2; */
/*     j = (max_y+1)/2; */
/*     k = (max_z+1)/2; */
    
/*     for(p=0; p<DIMS+1+2*DIMS; p++) */
/*       prop[p]=0.0; */
    
/*     for(i=1; i<=num_x; i++) { */
/*       for(j=0; j<max_y; j++) { */
/* 	nxy = i*max_y + j; */
	
/* 	for(k=0; k<max_z; k++) { */

/* 	  lvx =0.0; */
/* 	  lvy =0.0; */
/* 	  lvz =0.0; */
/* 	  for(q=0; q < Num_Dir; q++) { */
/* 	  /\*	    prop[0] += evector[0][q]*velcs_df[nxy][q][k]; *\/ */
/* 	    lvx += evector[1][q]*velcs_df[nxy][q][k]; */
/* 	    lvy += evector[2][q]*velcs_df[nxy][q][k]; */
/* 	    lvz += evector[3][q]*velcs_df[nxy][q][k]; */

/* 	    /\* stresses *\/ */
/* 	    prop[4] += evector[4][q]*velcs_df[nxy][q][k]; */
/* 	    prop[5] += evector[5][q]*velcs_df[nxy][q][k]; */
/* 	    prop[6] += evector[6][q]*velcs_df[nxy][q][k]; */
/* 	    prop[7] += evector[7][q]*velcs_df[nxy][q][k]; */
/* 	    prop[8] += evector[8][q]*velcs_df[nxy][q][k]; */
/* 	    prop[9] += evector[9][q]*velcs_df[nxy][q][k]; */

/* 	  } */
	  
/* 	  prop[0] += lvx*lvx+lvy*lvy+lvz*lvz; */
/* 	  prop[1] += lvx; */
/* 	  prop[2] += lvy; */
/* 	  prop[3] += lvz; */
	  
/* 	} */
/*       } */
/*     } */

    
/*     for(p=0; p<=DIMS; p++) */
/*       prop[p] /= (num_x*max_y*max_z); */
    
    
/*     fprintf(stream, "%d %le %le %le %le ", nstep, prop[0], prop[1], prop[2], prop[3]); */
/*     fprintf(stream, "%le %le %le %le %le %le\n", prop[4], prop[5], prop[6], prop[7], prop[8], prop[9]); */
    
/*     fclose(stream); */
/*   } */
  
}


/* Output velocity field in the xy-plane in vtk format */
void Write_Nodemap(int **node_map, int istep, char *work_dir)
{
  char PREFIX[20]="node_map";
  char filename[MAX_NL];
  FILE *outfile;

  int i,j,k,l,d;
  int q, nxy;
  double v[3];
  
  extern int max_x;
  extern int max_y;
  extern int max_z;

  if(istep < 10)
    sprintf(filename, "%s/data/%s.000%d.vtk", work_dir, PREFIX, istep);
  else if(istep < 100)
    sprintf(filename, "%s/data/%s.00%d.vtk", work_dir, PREFIX, istep);
  else if(istep < 1000)
    sprintf(filename, "%s/data/%s.0%d.vtk", work_dir, PREFIX, istep);
  else
    sprintf(filename, "%s/data/%s.%d.vtk", work_dir, PREFIX, istep);

  outfile = fopen(filename, "w");
  
  fprintf(outfile, "# vtk DataFile Version 2.3   \n");
  fprintf(outfile, "title u field(istep=%d)      \n", istep);
  fprintf(outfile, "ASCII                        \n");
  fprintf(outfile, "DATASET STRUCTURED_POINTS    \n");
  fprintf(outfile, "DIMENSIONS %d %d %d \n", max_x+2, max_y+2, max_z+2);
  fprintf(outfile, "ORIGIN  %le %le %le\n", 0.0, 0.0, 0.0);
  fprintf(outfile, "SPACING %le %le %le\n", 1.0, 1.0, 1.0);
  fprintf(outfile, "POINT_DATA %d\n", (max_x+2)*(max_y+2)*(max_z+2));
  fprintf(outfile, "VECTORS velocity double      \n");
  fclose(outfile);

  outfile = fopen(filename, "a");
  for(k=0; k<=max_z+1; k++) {
    for(j=0; j<=max_y+1; j++)
      for(i=0; i<=max_x+1; i++) {
		
	nxy = i*(max_y+2)+j;
	for(d=0; d<DIMS; d++)
	  v[d] = node_map[nxy][k];

	/*	fprintf(outfile, "(%- d\t %- d\t %- d)\t%- 12.10lf\t%- 12.10lf\t%- 12.10lf\n", i,j,k,v[0], v[1], v[2]); */
	fprintf(outfile, "%- 12.10lf\t%- 12.10lf\t%- 12.10lf\n", v[0], v[1], v[2]);
	/*	fprintf(outfile, "%- d\t%- d\t%- d\n", i, j, k); */
      }
  }
  fclose(outfile);
}

void Write_time(int nstep, char *work_dir)
{
  extern int n_proc;
  char filename[MAX_NL];
  FILE *stream;
  
  sprintf(filename, "%s/data/run_time.dat", work_dir); 
  
  stream = fopen(filename, "a");
  fprintf(stream, "%d %d %le\n", n_proc, nstep, wclock());
  fclose(stream);

}
