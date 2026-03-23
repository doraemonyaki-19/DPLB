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

void ludcmp(double **, int, int *, double *);
void lubksb(double **, int, int *, double b[]);
void matinv(double **, double **, int);

void initlbe(double ***velcs_df)
{
  extern int num_x, max_y, max_z;
  int i,j,k;
  int xy, z, q;
  int c2[Num_Dir];  
  double **a, **y;

  a = (double **) calloc(Num_Dir, sizeof(double *));
  if(a == 0) fatal_err("cannot allocate eigenvector matrix a", -1);
  for(i=0; i<Num_Dir; i++) {
    a[i] = (double *) calloc(Num_Dir, sizeof(double));
    if(a[i] == 0) fatal_err("cannot allocate eigenvector matrix a[i]", -1);
  }
  y = (double **) calloc(Num_Dir, sizeof(double *));
  if(y == 0) fatal_err("cannot allocate eigenvector matrix y", -1);
  for(i=0; i<Num_Dir; i++) {
    y[i] = (double *) calloc(Num_Dir, sizeof(double));
    if(y[i] == 0) fatal_err("cannot allocate eigenvector matrix y[i]", -1);
  }

  /* The transformation matrix from Adhikari et al., see Ladd and Duenweg, PRE, 2007 */

  for(i=0; i< Num_Dir; i++) {
    c2[i] = c_x[i]*c_x[i] + c_y[i]*c_y[i] + c_z[i]*c_z[i];
    evector[0][i] = 1;
    evector[1][i] = c_x[i];
    evector[2][i] = c_y[i];
    evector[3][i] = c_z[i];
    evector[4][i] = c2[i] - 1;
    evector[5][i] = 3*c_x[i]*c_x[i] - c2[i];
    evector[6][i] = c_y[i]*c_y[i] - c_z[i]*c_z[i];
    evector[7][i] = c_x[i]*c_y[i];
    evector[8][i] = c_z[i]*c_y[i];
    evector[9][i] = c_x[i]*c_z[i];
    evector[10][i] = (3*c2[i]-5)*c_x[i];
    evector[11][i] = (3*c2[i]-5)*c_y[i];
    evector[12][i] = (3*c2[i]-5)*c_z[i];
    evector[13][i] = (c_y[i]*c_y[i]-c_z[i]*c_z[i])*c_x[i];
    evector[14][i] = (c_z[i]*c_z[i]-c_x[i]*c_x[i])*c_y[i];
    evector[15][i] = (c_x[i]*c_x[i]-c_y[i]*c_y[i])*c_z[i];
    evector[16][i] = 3*c2[i]*c2[i]-6*c2[i]+1;
    evector[17][i] = (2*c2[i]-3)*(3*c_x[i]*c_x[i]-c2[i]);
    evector[18][i] = (2*c2[i]-3)*(c_y[i]*c_y[i]-c_z[i]*c_z[i]);
  }
  /*
  printf("evector =\n");
  for(i=0; i<Num_Dir; i++) {
    printf("[ ");
    for(j=0; j<Num_Dir; j++) 
      printf("%3.0lf ", evector[i][j]);
    printf(" ]\n");
  }
  */
  for(i=0; i<Num_Dir; i++) {
    wk[i] = 0;
    for(k=0; k<Num_Dir; k++)
      wk[i] += evector[i][k]*evector[i][k]*fac[k]; 
  }

/*
  printf("wk= [ ");
  for(i=0; i<Num_Dir; i++)
	printf("%d ", wk[i]);
  printf(" ]\n"); 
*/

  for(i=0; i<Num_Dir; i++) 
    for(j=0; j<Num_Dir; j++) 
      nmevector[i][j]= evector[i][j]*sqrt((double)fac[j]/(double)wk[i]);

  for(i=0; i<Num_Dir; i++) 
    for(j=0; j<Num_Dir; j++) 
      inv_nmevector[i][j] = nmevector[j][i];

  for(i=0; i<Num_Dir; i++) 
    for(j=0; j<Num_Dir; j++) 
      a[i][j] = evector[i][j];

  matinv(a, y, Num_Dir);

  for(i=0; i<Num_Dir; i++) 
    for(j=0; j<Num_Dir; j++) 
      inv_evector[i][j] = y[i][j];
  /*
  printf("inv_evector =\n");
  for(i=0; i<Num_Dir; i++) {
    printf("[ ");
    for(j=0; j<Num_Dir; j++) 
      printf("%3.0lf ", inv_evector[i][j]*Rho_Fl*4);
    printf(" ]\n");
  }
  */ 
  /* Rewrite the eigenvector matrices in vector format */
  /*
  for(i=0; i<Num_Dir; i++) 
    for(j=0; j<Num_Dir; j++) 
      lin_evector[i*Num_Dir+j] = evector[i][j];
 
  for(i=0; i<Num_Dir; i++) 
    for(j=0; j<Num_Dir; j++) 
      lin_invevector[i*Num_Dir+j] = inv_evector[i][j];
  */

  /*
  printf("\nProduct = \n");
  for(i=0; i<Num_Dir; i++) {
    printf("[ ");
    for(j=0; j<Num_Dir; j++) {
      y[i][j]=0.0;
      for(k=0; k<Num_Dir; k++)
	y[i][j] += a[i][k]*inv_evector[k][j];
      printf("%5.3lf ", y[i][j]);
    }
    printf(" ]\n");
  }
  */  

  for(xy=0; xy <=(num_x+1)*(max_y+2)+(max_y+1); xy++) 
    for(z=0; z<=max_z+1; z++) 
      for(q=0; q<Num_Dir; q++) 
	velcs_df[xy][q][z] = (double)fac[q];
   
  for(i=0; i<Num_Dir; i++) {
    free(a[i]);
    free(y[i]);
  }
}

/* inverting square matrix A using numerical recipes */
void matinv(double **a, double **y, int n)
{
  int *indx, i, j;
  double d, *col;
  double **A;

  A=(double **)calloc((n+1),sizeof(double *));
  if(A == 0) fatal_err("cannot allocate matrix A in matinv", -1);
  for(i=0; i<=n; i++) {
    A[i]=(double *)calloc((n+1),sizeof(double));
    if(A[i] == 0) fatal_err("cannot allocate matrix A[i] in matinv", -1);
  }

  col = (double *)calloc((n+1),sizeof(double));
  if(col == 0) fatal_err("cannot allocate col in matinv", -1);
  indx = (int *)calloc((n+1), sizeof(int));
  if(indx == 0) fatal_err("cannot allocate indx in matinv", -1);

  for(j=0; j<n; j++) {
    col[j]=0.0;
    indx[j]=0;
    for(i=0; i<n; i++)
      A[i+1][j+1] = a[i][j];
  }

  ludcmp(A,n,indx,&d);

  for(j=1; j<=n; j++) {
    for(i=1; i<=n; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb(A,n,indx,col);
    for(i=1;i<=n;i++) y[i-1][j-1]=col[i];
  }

  for(i=0; i<=n; i++)
    free(A[i]);
  free(col);
  free(indx);
}


void ludcmp(double **a, int n, int *indx, double *d)
{
  int i, imax, j, k;
  double big, dum, sum, temp;
  double *vv;

  vv = (double *) calloc((n+1),sizeof(double));
  *d=1.0;
  for(i=1; i<=n; i++) {
    big=0.0;
    for(j=1; j<=n; j++)
      if((temp=fabs(a[i][j])) > big) big=temp;
    if(big == 0.0) exit(13);
    vv[i] = 1.0/big;
  }

  for(j=1; j<=n; j++) {
    /* upper matrix */
    for(i=1; i<j; i++) {
      sum = a[i][j];
      for(k=1; k<i; k++) sum-=a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    /* lower matrix */
    for(i=j;i<=n;i++) {
      sum=a[i][j];
      for(k=1;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if((dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }

    if(j!=imax) {
      for(k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if(a[j][j]==0.0) a[j][j]=1e-6;
    if(j!=n) {
      dum=1.0/a[j][j];
      for(i=j+1;i<=n;i++) a[i][j] *=dum;
    }
  }

  free(vv);
}

void lubksb(double **a, int n, int *indx, double b[])
{
  int i, ii=0, ip, j;
  double sum;

  for(i=1; i<=n; i++) {
    ip = indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if(ii)
      for(j=ii; j<=i-1; j++) sum -= a[i][j]*b[j];
    else if (sum) ii=1;
    b[i]=sum;
  }
  for(i=n;i>=1;i--){
    sum=b[i];
    for(j=i+1; j<=n; j++) sum-=a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

