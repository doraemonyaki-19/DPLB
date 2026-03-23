/*  LBE: Lattice-Boltzmann Solver; Contains
    LBE_BCONDS (boundary conditions)
    LBE_MOVEXY (move velocity pointers for propagation in xy plane)
    LBE_ZCOL   (propagate in z direction and do collisions)  */
/**********************************************************************
ASLB : Lattice-Boltzmann simulation code for deformable particle+polymer+
lattice Boltzmann fluid
Copyright (C) 2019 Yeng-Long Chen

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

#include "header.h"
extern int    max_x, max_y, max_z;
extern int    num_x, x_min, x_max;
extern int    num_proc, n_proc;
extern double tau[2], tau_v[2], tau_g[2];

/* transfrom distribution to modes */
void ftom(double *, double *, int);

/* calculate the local equilibrium distribution */
void equil(double *, double *, double *);

void lbe_bconds (Float ***velcs_df)
{
  extern int wall_flag;
  double send_buf[MAX_B], recv_buf[MAX_B];
  int x, y, z, xy, n, n_dir;
  int send_id, recv_id;
  int num_msg, num_buf, n_buf, buf_size, ptr;

  num_msg = MAX_B/(Num_Dir_X*(max_z+2));           /* # msgs per buffer */
  num_buf = ((max_y+2)-1)/num_msg + 1;             /* # buffers */

  send_id = (n_proc+1)%num_proc;               /* +x-direction (+1) */
  recv_id = (n_proc-1+num_proc)%num_proc;
  for (n_buf = 0; n_buf < num_buf; n_buf++)
    {
      buf_size = min((max_y+2)-n_buf*num_msg, num_msg);
      ptr  = 0;
      for (y = 0; y < buf_size; y++)
	{
	  xy = num_x*(max_y+2) + y + n_buf*num_msg;
	  for (n = 0; n < Num_Dir_X; n++)
	    {
	      n_dir = x_p[n];
	      for (z = 0; z < (max_z+2); z++)
		{
		  send_buf[ptr] = velcs_df[xy][n_dir][z];
		  ptr++;
		}
	    }
	}

      vector_xchg (send_buf, recv_buf, ptr, send_id, recv_id);           
      ptr  = 0;
      for (y = 0; y < buf_size; y++)
	{
	  xy = y + n_buf*num_msg;
	  for (n = 0; n < Num_Dir_X; n++)
	    {
	      n_dir = x_p[n];
	      for (z = 0; z < max_z+2; z++)
		{
		  velcs_df[xy][n_dir][z] = recv_buf[ptr];
		  ptr++;
		}
	    }
	}
    }

  send_id = (n_proc-1+num_proc)%num_proc;      /* -x-direction (-1) */
  recv_id = (n_proc+1)%num_proc;
  for (n_buf = 0; n_buf < num_buf; n_buf++)
    {
      buf_size = min(max_y+2-n_buf*num_msg, num_msg);
      ptr  = 0;
      for (y = 0; y < buf_size; y++)
	{
	  xy = (max_y+2) + y + n_buf*num_msg;
	  for (n = 0; n < Num_Dir_X; n++)
	    {
	      n_dir = x_m[n];
	      for (z = 0; z < max_z+2; z++)
		{
		  send_buf[ptr] = velcs_df[xy][n_dir][z];
		  ptr++;
		}
	    }
	}

      vector_xchg (send_buf, recv_buf, ptr, send_id, recv_id);           
      ptr  = 0;
      for (y = 0; y < buf_size; y++)
	{
	  xy = (num_x+1)*(max_y+2) + y + n_buf*num_msg;
	  for (n = 0; n < Num_Dir_X; n++)
	    {
	      n_dir = x_m[n];
	      for (z = 0; z < max_z+2; z++)
		{
		  velcs_df[xy][n_dir][z] = recv_buf[ptr];
		  ptr++;
		}
	    }
	}
    }

  /* Periodic boundary conditions for y and z*/
  /* +/- y direction */
  if(wall_flag < 1) 
    for(x=0; x<=num_x+1; x++) {
      xy = x*(max_y+2) + max_y; 
      for(n=0; n<Num_Dir_X; n++) {
	n_dir = y_p[n];
	for(z=0; z<=max_z+1; z++)
	  velcs_df[x*(max_y+2)][n_dir][z] = velcs_df[xy][n_dir][z];
      }
      
      xy = x*(max_y+2) + 1;
      for(n=0; n<Num_Dir_X; n++) {
	n_dir = y_m[n];
	for(z=0; z<=max_z+1; z++)
	  velcs_df[x*(max_y+2)+max_y+1][n_dir][z] = velcs_df[xy][n_dir][z];
      }
    }


  /* +/- z direction */
  if(wall_flag < 2)
    for(xy=0; xy<=(num_x+1)*(max_y+2)+max_y+1; xy++) {
      for(n=0; n<Num_Dir_X; n++) {
	n_dir = z_m[n];
	velcs_df[xy][n_dir][max_z+1]=velcs_df[xy][n_dir][1];
	
	n_dir = z_p[n];
	velcs_df[xy][n_dir][0]=velcs_df[xy][n_dir][max_z];
      }
    }
}


/*_MOVEXY: Updates velocity pointers for propagation in xy-plane
    Uses cyclic permutation of pointers */

void lbe_move (Float ***velcs_df)
{
  Float *tmp_p;
  Float ztmp, vztmp;
  int   x, y, z, xo, yo, q; 
  int   xy, xy_new, xy_old, z_new, z_old;
  int   n, n_dir;

  for (n = 0; n < Num_Dir_X; n++)              /* +x-direction */
#pragma omp parallel for schedule(static) private(n_dir, tmp_p, x, xy_new, xy_old) num_threads(NTHREAD)
    for (y = 0; y <= max_y+1; y++)
      {
	n_dir  = x_p[n];
	tmp_p  = velcs_df[num_x*(max_y+2)+y][n_dir];
	for (x = num_x; x >= 1; x--)
	  {
	    xy_new = x*(max_y+2) + y;
	    xy_old = xy_new  - (max_y+2);
	    velcs_df[xy_new][n_dir] = velcs_df[xy_old][n_dir];
	  }
	velcs_df[y][n_dir] = tmp_p;
      }

  for (n = 0; n < Num_Dir_X; n++)              /* -x-direction */
#pragma omp parallel for schedule(static) private(n_dir, tmp_p, x, xy_new, xy_old) num_threads(NTHREAD)
    for (y = 0; y <= max_y+1; y++)
      {
	n_dir  = x_m[n];
	tmp_p  = velcs_df[(max_y+2)+y][n_dir];
	for (x = 1; x <= num_x; x++)
	  {
	    xy_new = x*(max_y+2) + y;
	    xy_old = xy_new  + (max_y+2);
	    velcs_df[xy_new][n_dir] = velcs_df[xy_old][n_dir];
	  }
	velcs_df[(num_x+1)*(max_y+2)+y][n_dir] = tmp_p;
      }

  /* For DP-TUBE, x should go from 1 to num_x, otherwise the solid node velocities will leak into fluid in the buffer */
  for (n = 0; n < Num_Dir_X; n++)              /* +y-direction */
#pragma omp parallel for schedule(static) private(n_dir, tmp_p, y, xy_new, xy_old) num_threads(NTHREAD)
    for (x = 1; x <= num_x; x++)
      {
	n_dir  = y_p[n];
	tmp_p  = velcs_df[x*(max_y+2)+max_y][n_dir];
	for (y = max_y; y >= 1; y--)
	  {
	    xy_new = x*(max_y+2) + y;
	    xy_old = xy_new - 1;
	    velcs_df[xy_new][n_dir] = velcs_df[xy_old][n_dir];
	  }
	velcs_df[x*(max_y+2)][n_dir] = tmp_p;
      }
				
  for (n = 0; n < Num_Dir_X; n++)              /* -y-direction */
#pragma omp parallel for schedule(static) private(n_dir, tmp_p, y, xy_new, xy_old) num_threads(NTHREAD)
    for (x = 1; x <= num_x; x++)
      {
	n_dir  = y_m[n];
  	tmp_p  = velcs_df[x*(max_y+2)+1][n_dir];
	for (y = 1; y <= max_y; y++)
	  {
	    xy_new = x*(max_y+2) + y;
	    xy_old = xy_new + 1;
	    velcs_df[xy_new][n_dir] = velcs_df[xy_old][n_dir];
	  }
	velcs_df[x*(max_y+2)+(max_y+1)][n_dir] = tmp_p;
      }
  
  for (n = 0; n < Num_Dir_X; n++)               /* +z-direction */
#pragma omp parallel for schedule(static) private(y, z, xy, ztmp, n_dir) num_threads(NTHREAD)
    for (x=1; x<=num_x; x++) 
      for(y=1; y<=max_y; y++) {
	xy = x*(max_y+2)+y;
	n_dir  = z_p[n];
	ztmp  = velcs_df[xy][n_dir][max_z];
	for (z = max_z; z >= 1; z--) 
	  velcs_df[xy][n_dir][z] = velcs_df[xy][n_dir][z-1];
	velcs_df[xy][n_dir][0] = ztmp;
      }
  
  for (n = 0; n < Num_Dir_X; n++) {              /* -z-direction */
#pragma omp parallel for schedule(static) private(y, z, xy, ztmp, n_dir) num_threads(NTHREAD)
    for (x=1; x<=num_x; x++) 
      for(y=1; y<=max_y; y++) {
	xy = x*(max_y+2)+y;
	n_dir  = z_m[n];
	ztmp  = velcs_df[xy][n_dir][1];
	for (z = 1; z <= max_z; z++) 
	  velcs_df[xy][n_dir][z] = velcs_df[xy][n_dir][z+1];

	velcs_df[xy][n_dir][max_z+1] = ztmp;
      }
  }
}

/*  LBE_ZCOL: Collision operator in z for 1 column  */
/*  Two relaxation time method : Chun & Ladd (2009) */
void lbe_zcol (Float **velcz_df, int *node_map, struct vector f_ext, int y, double rannum[])
{
  /*extern mt_state twister; */
  extern Float std_xi[Num_Dir], std_xi_1[Num_Dir];
  extern int add_noise;
  extern Float FluidStress[2*DIMS];

  Float  ThreadStress[NTHREAD*2*DIMS];
  double lambda[Num_Dir], xi[Num_Dir], m2;
  double lam[2], lam_v[2], lam_g[2];
  double dmom[4];
  double f[Num_Dir], modes[Num_Dir], meq[Num_Dir], momfluc[Num_Dir_NConserve];
  int    nproc;
  int    posfluc_flag=0, errcode;
  int    z, z_p, z_m, z10, q;
  int    i,j,k,nyz;
  FILE   *stream;

  lam[0]   = - 1.0/tau[0];                      /* Eigenvalues */
  lam_v[0] = - 1.0/tau_v[0];
  lam_g[0] = - 1.0/tau_g[0];
  lam[1]   = - 1.0/tau[1];                      
  lam_v[1] = - 1.0/tau_v[1];
  lam_g[1] = - 1.0/tau_g[1];

  for(q=0; q<Num_Dir_NConserve ; q++) 
    momfluc[q] = 0.0;

  for(q=0; q<2*DIMS; q++)
    FluidStress[q] = 0.0;

  for(i=0; i<NTHREAD; i++)
    for(q=0; q<2*DIMS; q++)
      ThreadStress[i*2*DIMS+q] = 0.0;

  /* normal modes: m0 = mass; m1-m3 = momenta; m4-m9 = stress  */
#pragma omp parallel for schedule(static) private(q, lambda, xi, f, m2, meq, momfluc, modes) num_threads(NTHREAD)
  for (z = 1; z <= max_z; z++) {
    if (node_map[z] != 1) {            /* Skip solid nodes */
      /* Eigenvalues for the conserved modes */
      for(q=0; q<4; q++)
	lambda[q] = 0.0;

      if(node_map[z] == 0) {
	for(q=4; q<10; q++) 
	  lambda[q] = lam[0];                      /* Eigenvalues */
	
	for(q=10; q<Num_Dir; q++)
	  lambda[q] = lam_g[0];

	for(q=0; q<Num_Dir; q++)
	  xi[q] = std_xi[q];
      }
      else {
	for(q=4; q<10; q++) 
	  lambda[q] = lam[node_map[z]-1];                      /* Eigenvalues */

	for(q=10; q<Num_Dir; q++)
	  lambda[q] = lam_g[node_map[z]-1];

	for(q=0; q<Num_Dir; q++)
	  xi[q] = std_xi_1[q];
      } 

      for(q=0; q<Num_Dir; q++)
	f[q] = velcz_df[q][z] ;

      /* add fluctuations */
      if(add_noise > 1) 
	for(q=0; q<Num_Dir_NConserve ; q++) 
	  momfluc[q] = rannum[((y-1)*max_z+(z-1))*(Num_Dir_NConserve)+q]*xi[q+4]; 

      /* Transform from the distribution to the moments */
      ftom(f, modes, 1); 

      m2 = modes[1]*modes[1]+modes[2]*modes[2]+modes[3]*modes[3];	
      meq[0] = modes[0];
      meq[1] = modes[1];
      meq[2] = modes[2];
      meq[3] = modes[3];
      meq[4] = m2/meq[0];
      meq[5] = (3*modes[1]*modes[1]-m2)/meq[0];
      meq[6] = (modes[2]*modes[2]-modes[3]*modes[3])/meq[0];
      meq[7] = modes[1]*modes[2]/meq[0];
      meq[8] = modes[3]*modes[2]/meq[0];
      meq[9] = modes[3]*modes[1]/meq[0];
      for(q=10; q<Num_Dir; q++)
	meq[q]=0.0;

      nproc = omp_get_thread_num();
      for(q=0; q<2*DIMS; q++)
	ThreadStress[nproc*2*DIMS+q] += modes[q+4]/2;

      /* Collision */
      for(q=4; q<Num_Dir; q++) 
	modes[q] += lambda[q]*(modes[q]-meq[q])+momfluc[q-4]; 
      
      /* add external momentum and stress */
      modes[1] += f_ext.x;
      modes[2] += f_ext.y;
      modes[3] += f_ext.z;
      modes[4] += 2.0*(modes[1]*f_ext.x+modes[2]*f_ext.y+modes[3]*f_ext.z)/meq[0];
      modes[5] += (4.0*modes[1]*f_ext.x-2.0*(modes[2]*f_ext.y+modes[3]*f_ext.z))/meq[0];
      modes[6] += 2.0*(modes[2]*f_ext.y-modes[3]*f_ext.z)/meq[0];
      modes[7] += (modes[1]*f_ext.y+modes[2]*f_ext.x)/meq[0];
      modes[8] += (modes[3]*f_ext.y+modes[2]*f_ext.z)/meq[0];
      modes[9] += (modes[3]*f_ext.x+modes[1]*f_ext.z)/meq[0];

      for(q=0; q<2*DIMS; q++)
	ThreadStress[nproc*2*DIMS+q] += modes[q+4]/2;

      /*  Compute new population densities from modes  */
      ftom(f, modes, -1); 

      for(q=0; q<Num_Dir; q++) 
	velcz_df[q][z] = f[q];
    }
    else {
     for(q=0; q<Num_Dir; q++)
	velcz_df[q][z] = 0.0 ;
    }
  }

  for(i=0; i<NTHREAD; i++) 
    for(q=0; q<2*DIMS; q++)
      FluidStress[q] += ThreadStress[i*2*DIMS+q];
}


void ftom(double *f, double *m, int flag)
{
  int i,j;
  double m8;  
  double m2_3p, m2_3m, m1_3p, m1_3m, m1_2p, m1_2m, m3_4p, m3_4m;
  double m5_6p, m5_6m, m7_8p, m7_8m, m9_10p, m9_10m;
  double m0_16p, m0_16m, m1_10m, m2_11m, m3_12m;
  double m4_10m, m4_10p, m5_17p, m5_17m, m11_7p, m11_7m;
  double m6_18p, m6_18m, m11_12p, m11_12m, m12_9p, m12_9m;    
  double m13_14p, m13_14m, m13_15p, m13_15m;
  double m14_15p, m14_15m, m15_16p, m15_16m, m17_18p, m17_18m;
  double m4_5_17;


  if(flag == 1) {
    /*   
    for(i=0; i<Num_Dir; i++) {
      m[i] = 0.0;
      for(j=0; j<Num_Dir; j++) 
	m[i] += evector[i][j]*f[j];
    }
    */    
    
    m1_2p = f[1]+f[2];  m1_2m = f[1]-f[2];
    m3_4p = f[3]+f[4];  m3_4m = f[3]-f[4];
    m5_6p = f[5]+f[6];  m5_6m = f[5]-f[6];
    m7_8p = f[7]+f[8];  m7_8m = f[7]-f[8];
    m9_10p = f[9]+f[10];  m9_10m = f[9]-f[10];
    m11_12p = f[11]+f[12];  m11_12m = f[11]-f[12];
    m13_14p = f[13]+f[14];  m13_14m = f[13]-f[14];
    m15_16p = f[15]+f[16];  m15_16m = f[15]-f[16];
    m17_18p = f[17]+f[18];  m17_18m = f[17]-f[18];

    m[0] = 0;       for(i=0; i<Num_Dir; i++)      m[0] += f[i];
    m[1] = m1_2m+m11_12m-m13_14m+m15_16m+m17_18m;
    m[2] = m3_4m+m7_8m+m9_10m+m15_16m-m17_18m;
    m[3] = m5_6m+m7_8m-m9_10m+m11_12m+m13_14m;
    m[4] =-f[0];    for(i=7; i<Num_Dir; i++)      m[4] += f[i];
    m[5] = 2.0*(m1_2p-m7_8p-m9_10p)-m3_4p-m5_6p;  for(i=11; i<Num_Dir; i++) m[5]+=f[i];
    m[6] = m3_4p-m5_6p-m11_12p-m13_14p+m15_16p+m17_18p;
    m[7] = m15_16p-m17_18p;
    m[8] = m7_8p-m9_10p;;
    m[9] = m11_12p-m13_14p;
    m[10]= 2.0*(-m1_2m)+m11_12m-m13_14m+m15_16m+m17_18m;
    m[11]= 2.0*(-m3_4m)+m7_8m+m9_10m+m15_16m-m17_18m;
    m[12]= 2.0*(-m5_6m)+m7_8m-m9_10m+m11_12m+m13_14m;
    m[13]= -m11_12m+m13_14m+m15_16m+m17_18m;
    m[14]= m7_8m+m9_10m-m15_16m+m17_18m;
    m[15]= -m7_8m+m9_10m+m11_12m+m13_14m;
    m[16]= f[0]; for(i=1; i<=6; i++) m[16]-=2*f[i]; for(i=7; i<Num_Dir; i++) m[16]+=f[i];
    m[17]= -2.0*(m1_2p+m7_8p+m9_10p)+m3_4p+m5_6p; for(i=11; i<Num_Dir; i++) m[17]+=f[i];
    m[18]= -m3_4p+m5_6p-m11_12p-m13_14p+m15_16p+m17_18p;
  }

  else if(flag == -1) {

    
    m8=36.0*m[8];  
    m1_2p = 12.0*(m[1]+m[2]);     m1_2m = 12.0*(m[1]-m[2]);
    m2_3p = 12.0*(m[2]+m[3]);     m2_3m = 12.0*(m[2]-m[3]);
    m1_3p = 12.0*(m[1]+m[3]);     m1_3m = 12.0*(m[1]-m[3]);
    m1_10m = 24.0*(m[1]-m[10]);   m2_11m = 24.0*(m[2]-m[11]);    m3_12m = 24.0*(m[3]-m[12]);
    m0_16p = 4.0*m[0]+2.0*m[16];  m0_16m = 8.0*(m[0]-m[16]);
    m4_10p = 6.0*(m[4]+m[10]);    m4_10m = 6.0*(m[4]-m[10]);
    m6_18p = 9.0*(m[6]+m[18]);    m6_18m = 18.0*(m[6]-m[18]);
    m5_17p = 3.0*(m[5]+m[17]);    m5_17m = 6.0*(m[5]-m[17]);
    m11_7p=  6.0*m[11]+36.0*m[7]; m11_7m=  6.0*m[11]-36.0*m[7];
    m11_12p= 6.0*(m[11]+m[12]);   m11_12m= 6.0*(m[11]-m[12]);
    m12_9p = 6.0*m[12]+36.0*m[9]; m12_9m = 6.0*m[12]-36.0*m[9];
    m13_15p= 18.0*(m[13]+m[15]);  m13_15m=18.0*(m[13]-m[15]);
    m13_14p= 18.0*(m[13]+m[14]);  m13_14m=18.0*(m[13]-m[14]);
    m14_15p= 18.0*(m[14]+m[15]);  m14_15m=18.0*(m[14]-m[15]);
    m4_5_17= 6.0*(m[4]-m[5]-m[17]);

    f[0] = 48.0*m[0]-72.0*m[4]+24.0*m[16];
    f[1] =  8.0*m[0]+m1_10m+2.0*m5_17m-8.0*m[16];
    f[2] =  8.0*m[0]-m1_10m+2.0*m5_17m-8.0*m[16];
    f[3] =  m0_16m+m2_11m-m5_17m+m6_18m;
    f[4] =  m0_16m-m2_11m-m5_17m+m6_18m;
    f[5] =  m0_16m+m3_12m-m5_17m-m6_18m;
    f[6] =  m0_16m-m3_12m-m5_17m-m6_18m;
    f[7] =  m0_16p+m2_3p+m4_5_17+m11_12p+m8+m14_15m;
    f[8] =  m0_16p-m2_3p+m4_5_17-m11_12p+m8-m14_15m;
    f[9] =  m0_16p+m2_3m+m4_5_17+m11_12m-m8+m14_15p;
    f[10]=  m0_16p-m2_3m+m4_5_17-m11_12m-m8-m14_15p;
    f[11]=  m0_16p+m1_3p+m4_10p+m12_9p-m13_15m+m5_17p-m6_18p;
    f[12]=  m0_16p-m1_3p+m4_10m-m12_9m+m13_15m+m5_17p-m6_18p;
    f[13]=  m0_16p-m1_3m+m4_10m+m12_9m+m13_15p+m5_17p-m6_18p;
    f[14]=  m0_16p+m1_3m+m4_10p-m12_9p-m13_15p+m5_17p-m6_18p;
    f[15]=  m0_16p+m1_2p+m4_10p+m11_7p+m13_14m+m5_17p+m6_18p;
    f[16]=  m0_16p-m1_2p+m4_10m-m11_7m-m13_14m+m5_17p+m6_18p;
    f[17]=  m0_16p+m1_2m+m4_10p-m11_7p+m13_14p+m5_17p+m6_18p;
    f[18]=  m0_16p-m1_2m+m4_10m+m11_7m-m13_14p+m5_17p+m6_18p;

    for(i=0; i<Num_Dir; i++)
      f[i] /= (Rho_Fl*4.0);
  
  }
}

void equil(double *m, double *feq, double *meq)
{
  int q;
  double udotc;
  double uucc;
  double m1 = m[1]/Rho_Fl;
  double m2 = m[2]/Rho_Fl;
  double m3 = m[3]/Rho_Fl;
  double m1m1 = m1*m1;
  double m2m2 = m2*m2;
  double m3m3 = m3*m3;
  double m1m2 = 2*m1*m2;
  double m1m3 = 2*m1*m3;
  double m2m3 = 2*m2*m3;

  for(q=0; q< Num_Dir; q++) {
    udotc = m1*c_x[q]+m2*c_y[q]+m3*c_z[q];
    uucc = m1m1*(c_x[q]*c_x[q]-CS2) +
      m2m2*(c_y[q]*c_y[q]-CS2) +
      m3m3*(c_z[q]*c_z[q]-CS2) +
      m1m2*c_x[q]*c_y[q] +
      m1m3*c_x[q]*c_z[q] +
      m2m3*c_y[q]*c_z[q] ;
    feq[q] = fac[q]*m[0]*(1.0+ udotc/CS2 + uucc/(2*CS2*CS2))/Rho_Fl;
  }

  ftom(feq, meq, 1);
}


/*  LBE_ZCOL: Collision operator and propagation in z for 1 column  */
/*  Ladd stress moments method */
void lbe_zcol_old (Float **velcz_df, int *node_map, struct vector f_ext, int y, double rannum[])
{
  /*extern mt_state twister; */
  extern Float std_xi[Num_Dir], std_xi_1[Num_Dir];
  extern int add_noise;

  double modes[MAX_Z*10];
  double nfluc[Num_Dir], xi[Num_Dir];
  double m, m0, m1, m2, m3, m4, m5, m6, m7, m8, m9;
  double m1p, m2p, m3p, m4p, m5p, m6p, m7p, m8p, m9p;
  double m1m, m2m, m3m, m4m, m5m, m6m, m7m, m8m, m9m;
  double m0e, m4e, m5e, m6e, m7e, m8e, m9e;
  double lambda, lambda_v, lam[2], lam_v[2], lam_g[2];
  int    z, z_p, z_m, z10, q;
  int    nyz;
  FILE   *stream;

  lam[0]   = 1.0 -1.0/tau[0];                      /* Eigenvalues */
  lam_v[0] = 1.0 -1.0/tau_v[0];
  lam[1]   = 1.0 -1.0/tau[1];                      
  lam_v[1] = 1.0 -1.0/tau_v[1];

  /*  Compute normal modes: m0 = mass; m1-m3 = momenta; m4-m9 = stress  */
  for (z = 1, z10 = 0; z <= max_z; z++, z10 += 10)
    {
      if (node_map[z] != 1)  {          /* Skip solid nodes */
	if(node_map[z] == 0) {
	  lambda   = lam[0];                      /* Eigenvalues */
	  lambda_v = lam_v[0];
	  for(q=0; q<Num_Dir; q++)
	    xi[q] = std_xi[q];
	}
	else {
	  lambda   = lam[node_map[z]-1];                      /* Eigenvalues */
	  lambda_v = lam_v[node_map[z]-1];
	  for(q=0; q<Num_Dir; q++)
	    xi[q] = std_xi_1[q];
	}

	nyz=((y-1)*max_z+(z-1))*6;

	m0  = velcz_df[ 0][z];
	m1p = velcz_df[ 1][z]   + velcz_df[ 2][z];
	m1m = velcz_df[ 1][z]   - velcz_df[ 2][z];
	m2p = velcz_df[ 3][z]   + velcz_df[ 4][z];
	m2m = velcz_df[ 3][z]   - velcz_df[ 4][z];
	m3p = velcz_df[ 5][z]   + velcz_df[ 6][z];
	m3m = velcz_df[ 5][z]   - velcz_df[ 6][z];
	m4p = velcz_df[ 7][z]   + velcz_df[ 8][z];
	m4m = velcz_df[ 7][z]   - velcz_df[ 8][z];
	m5p = velcz_df[ 9][z]   + velcz_df[10][z];
	m5m = velcz_df[ 9][z]   - velcz_df[10][z];
	m6p = velcz_df[11][z]   + velcz_df[12][z];
	m6m = velcz_df[11][z]   - velcz_df[12][z];
	m7p = velcz_df[13][z]   + velcz_df[14][z];
	m7m = velcz_df[13][z]   - velcz_df[14][z];
	m8p = velcz_df[15][z]   + velcz_df[16][z];
	m8m = velcz_df[15][z]   - velcz_df[16][z];
	m9p = velcz_df[17][z]   + velcz_df[18][z];
	m9m = velcz_df[17][z]   - velcz_df[18][z];

	m0 += m1p + m2p + m3p + m4p + m5p + m6p + m7p + m8p + m9p;
	m1  = m1m + m8m + m9m + m6m - m7m;
	m2  = m2m + m4m + m5m + m8m - m9m;
	m3  = m3m + m6m + m7m + m4m - m5m;
	m4  = m1p + m8p + m9p + m6p + m7p;
	m5  = m2p + m4p + m5p + m8p + m9p;
	m6  = m3p + m6p + m7p + m4p + m5p;
	m7  = m4p - m5p;
	m8  = m6p - m7p;
	m9  = m8p - m9p;

	m0e = m0 + Rho_Fl;                     /* Equilibrium distribution */
	m4e = m1*m1/m0e;
	m5e = m2*m2/m0e;
	m6e = m3*m3/m0e;
	m7e = m2*m3/m0e;
	m8e = m3*m1/m0e;
	m9e = m1*m2/m0e;
	m4 -= m4e + m0*CS2;                    /* Non-equilibrium stress */
	m5 -= m5e + m0*CS2;
	m6 -= m6e + m0*CS2;
	m7 -= m7e;
	m8 -= m8e;
	m9 -= m9e;

	m   = (m4 + m5 + m6)/3.0;              /* Collisions */
	m4 -= m;
	m5 -= m;
	m6 -= m;
	m4 *= lambda;
	m5 *= lambda;
	m6 *= lambda;
	m7 *= lambda;
	m8 *= lambda;
	m9 *= lambda;
	m  *= lambda_v;

	m4 += m;
	m5 += m;
	m6 += m;

	m4 += m4e;                             /* Non-equilibrium stress */
	m5 += m5e;
	m6 += m6e;
	m7 += m7e;
	m8 += m8e;
	m9 += m9e;

	/* Add extern stress */
	m4 +=  m1*f_ext.x*2.0/m0e;
	m5 +=  m2*f_ext.y*2.0/m0e;
	m6 +=  m3*f_ext.z*2.0/m0e;
	m7 += (m2*f_ext.z + m3*f_ext.y)/m0e;
	m8 += (m3*f_ext.x + m1*f_ext.z)/m0e;
	m9 += (m1*f_ext.y + m2*f_ext.x)/m0e;
			
	/* Add fluid fluctuations*/
	  	  
	if(add_noise > 1) {
	  m4 += xi[4]*rannum[nyz];
	  m5 += xi[5]*rannum[nyz+1];
	  m6 += xi[6]*rannum[nyz+2];
	  m7 += xi[7]*rannum[nyz+3];
	  m8 += xi[8]*rannum[nyz+4];
	  m9 += xi[9]*rannum[nyz+5];
	}
		
	m1 += f_ext.x;                     /* Add external force */
	m2 += f_ext.y;
	m3 += f_ext.z;
			
	modes[z10  ] = m0;                     /* Save modes for column */
	modes[z10+1] = m1;
	modes[z10+2] = m2;
	modes[z10+3] = m3;
	modes[z10+4] = m4;
	modes[z10+5] = m5;
	modes[z10+6] = m6;
	modes[z10+7] = m7;
	modes[z10+8] = m8;
	modes[z10+9] = m9;
      }
    }
		
  /*  Compute new population densities from modes  */

  for (z = 1, z10 = 0; z <= max_z; z++, z10 += 10)
    {
      if (node_map[z] != 1)            /* Fluid nodes */
	{
	  m0  = modes[z10  ]/36.0;
	  m1  = modes[z10+1]/12.0;
	  m2  = modes[z10+2]/12.0;
	  m3  = modes[z10+3]/12.0;
	  m4  = modes[z10+4]/8.0;
	  m5  = modes[z10+5]/8.0;
	  m6  = modes[z10+6]/8.0;
	  m7  = modes[z10+7]/4.0;
	  m8  = modes[z10+8]/4.0;
	  m9  = modes[z10+9]/4.0;

	  m0 -=(m4 + m5 + m6)/3.0;
	  m1p =(m0 + m4)*2.0;
	  m2p =(m0 + m5)*2.0;
	  m3p =(m0 + m6)*2.0;
	  m4p = m0 + m5 + m6 + m7;
	  m5p = m0 + m5 + m6 - m7;
	  m6p = m0 + m6 + m4 + m8;
	  m7p = m0 + m6 + m4 - m8;
	  m8p = m0 + m4 + m5 + m9;
	  m9p = m0 + m4 + m5 - m9;
	  m1m = m1*2.0;
	  m2m = m2*2.0;
	  m3m = m3*2.0;
	  m4m = m2 + m3;
	  m5m = m2 - m3;
	  m6m = m3 + m1;
	  m7m = m3 - m1;
	  m8m = m1 + m2;
	  m9m = m1 - m2;

	  velcz_df[ 0][z] = m0*12.0;
	  velcz_df[ 1][z] = m1p + m1m;
	  velcz_df[ 2][z] = m1p - m1m;
	  velcz_df[ 3][z] = m2p + m2m;
	  velcz_df[ 4][z] = m2p - m2m;
	  velcz_df[ 5][z] = m3p + m3m;
	  velcz_df[ 6][z] = m3p - m3m;
	  velcz_df[ 7][z] = m4p + m4m;
	  velcz_df[ 8][z] = m4p - m4m;
	  velcz_df[ 9][z] = m5p + m5m;
	  velcz_df[10][z] = m5p - m5m;
	  velcz_df[11][z] = m6p + m6m;
	  velcz_df[12][z] = m6p - m6m;
	  velcz_df[13][z] = m7p + m7m;
	  velcz_df[14][z] = m7p - m7m;
	  velcz_df[15][z] = m8p + m8m;
	  velcz_df[16][z] = m8p - m8m;
	  velcz_df[17][z] = m9p + m9m;
	  velcz_df[18][z] = m9p - m9m;
	}
      else                                       /* Zero solid nodes */
	{
	  velcz_df[ 1][z] = 0.0;
	  velcz_df[ 2][z] = 0.0;
	  velcz_df[ 3][z] = 0.0;
	  velcz_df[ 4][z] = 0.0;
	  velcz_df[ 5][z] = 0.0;
	  velcz_df[ 6][z] = 0.0;
	  velcz_df[ 7][z] = 0.0;
	  velcz_df[ 8][z] = 0.0;
	  velcz_df[ 9][z] = 0.0;
	  velcz_df[10][z] = 0.0;
	  velcz_df[11][z] = 0.0;
	  velcz_df[12][z] = 0.0;
	  velcz_df[13][z] = 0.0;
	  velcz_df[14][z] = 0.0;
	  velcz_df[15][z] = 0.0;
	  velcz_df[16][z] = 0.0;
	  velcz_df[17][z] = 0.0;
	  velcz_df[18][z] = 0.0;
	}
    }
}
