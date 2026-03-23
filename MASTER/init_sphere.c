/* Model sphere as an icosahedron w/ vertices joined by springs */
/* Model RBC as compressed level 2 spheres from template */
/* spheres always inserted before RBC */
/* Coupled with LBE code to add hydrodynamic forces */ 

/***********************************************************************
 * ASLB : Lattice-Boltzmann simulation code for deformable particle+polymer+
 * lattice Boltzmann fluid
 * Copyright (C) 2010 Yeng-Long Chen
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
#include <time.h>

#define NVERTICES 12   /* number of vertices for an icosahedron */
#define NBONDS 5       /* number of bonds for a vertex on an icosahedron */

void setup(double ***v, int blist[NVERTICES][NBONDS+1][3], double radius[NTYPES]);
void setup_template(double **, int ***, int, double, char *);
int check_bond(int k, int m, struct monomer *monomers);
int check_face(int, int, int, struct face *, int);
int get_midpoint(int a, int b, struct monomer *monomers, int **midpoints);

int sphere_init(struct DP_param *sphere_pm, struct DP_param *cyl_pm, struct DP *spheres, struct monomer *monomers, struct face *faces, char *work_dir)
{
  extern int n_proc;
  extern int max_x, max_y, max_z;
  extern int wall_flag;

  int i,j,k,d,m,f, pnum, m_temp;
  int n, level;
  int n1, n2, n3, n4;
  int n12, n13, n14, n23, n24;
  int nx, ny, nz;
  int trial_count=0;
  int overlap, rot_flag;
  int temp;
  int type, type1, start, start0;
  int NDP = sphere_pm->NDP;
  int num_beads = sphere_pm->num_beads;
  int maxNbeads;
  int MRBC;
  int Ntype[NTYPES];
  int nlevel[NTYPES];
  int Nbeads[NTYPES];
  int Nfaces[NTYPES];
  int blist_ico[NVERTICES][NBONDS+1][3];   /* blist for the starting icosahedron */
  int bcount;
  int ***blist_temp, ***blist_rbc;
  int **midpoints;
  int addface;

  double prefactor;
  double maxsize[DIMS], mindim;
  double bondlength, radius[NTYPES];
  double Vfactor[NTYPES], Sfactor[NTYPES];
  double var;
  double dx, dy, dz, r2, rshift;
  double **centers;            /* center position of spheres */
  double theta, phi, psi;
  double axis[DIMS][DIMS];     /* x,y,z axis of sphere */
  double ***v;   /* the bond list */
  double p_mid[DIMS], dr[DIMS];

  char filename[MAX_NL];
  FILE *stream=0;

  int test=TRUE;
  char test_file[MAX_NL] = {"data/init.dat"};
  FILE *file_ptr=0;

  srand((unsigned)time(0));

  maxsize[0]=(double)max_x;
  maxsize[1]=(double)max_y;
  maxsize[2]=(double)max_z;

  for (i=0 ; i<NTYPES ; i++) {
    Ntype[i] = sphere_pm->Ntype[i];
    nlevel[i] = sphere_pm->nlevel[i];
    Nbeads[i] = sphere_pm->N_per_DP[i];
    Nfaces[i] = sphere_pm->face_per_DP[i];
  }

  if(nlevel[0] == -3 || nlevel[1] == -3 || nlevel[0] ==4 || nlevel[1] ==4)
    MRBC = 2562;
  else if(nlevel[0] == -2 || nlevel[1] == -2 || nlevel[0] ==3 || nlevel[1] ==3)
    MRBC = 642;
  else if(nlevel[0] == -1 || nlevel[1] == -1 || nlevel[0] ==2 || nlevel[1] ==2)
    MRBC = 162;
  else {
    printf("set MRBC in init_sphere\n");
    exit(22);
  }

  /* memory allocation */

  centers = (double**) calloc(NDP, sizeof(double*));
  if (centers == 0)  fatal_err ("cannot allocate centers", -1);

  midpoints = (int **) calloc(num_beads, sizeof(int *));
  if (midpoints == 0) fatal_err("cannot allocate midpoints", -1);

  v = (double ***) calloc(NTYPES, sizeof(double **));
  if (v == 0) fatal_err("cannot allocate bond vectors", -1);

  blist_temp = (int ***) calloc(num_beads, sizeof(int **));
  if (blist_temp == 0) fatal_err("cannot allocate blist", -1);

  blist_rbc = (int ***) calloc(MRBC, sizeof(int **));
  if (blist_rbc == 0) fatal_err("cannot allocate blist", -1);

  for (n=0 ; n<NDP ; n++) {
    centers[n] = (double *) calloc(DIMS, sizeof(double));
    if (centers[n] == 0)  fatal_err ("cannot allocate centers", -1);
  }
  for (i=0 ; i<num_beads ; i++) {
    midpoints[i] = (int *) calloc(MAX_BOND+1, sizeof(int));
    if (midpoints[i] == 0) fatal_err("cannot allocate midpoints", -1);

    blist_temp[i] = (int **) calloc(MAX_BOND+1, sizeof(int *));
    if (blist_temp[i] == 0) fatal_err("cannot allocate blist", -1);

    blist_rbc[i] = (int **) calloc(MAX_BOND+1, sizeof(int *));
    if (blist_rbc[i] == 0) fatal_err("cannot allocate blist", -1);

    for (j=0 ; j<=MAX_BOND ; j++) {
      blist_temp[i][j] = (int *) calloc(DIMS, sizeof(int));
      if (blist_temp[i][j] == 0) fatal_err("cannot allocate blist", -1);

      blist_rbc[i][j] = (int *) calloc(DIMS, sizeof(int));
      if (blist_rbc[i][j] == 0) fatal_err("cannot allocate blist", -1);
    }
  }

  maxNbeads = max(Nbeads[0], Nbeads[1]);

  for(type = 0; type <NTYPES; type++) {
    v[type]= (double **) calloc(maxNbeads, sizeof(double *));
    if (v[type] == 0) fatal_err("cannot allocate bond vectors", -1);
    
    for(i=0; i<maxNbeads; i++) {
      v[type][i] = (double *) calloc(DIMS, sizeof(double));
      if (v[type][i] == 0) fatal_err("cannot allocate bond vectors", -1);
    }
  }

  if (test)
    printf("Adding %d spheres, %d particles in total\n", NDP, num_beads);

  /* Initialize monomer properties */
  for(i=0; i<num_beads; i++) {
    if(sphere_pm->ev_type == 0)   /* HS */
      monomers[i].radius = 0.5;
    else if(sphere_pm->ev_type == 1 || sphere_pm->ev_type == 3) {   /* WCA */
      monomers[i].radius = 0.5;
      if(sphere_pm->nlevel[0] == 3 || sphere_pm->nlevel[0] == -2)
	//	monomers[i].radius=0.25;
	monomers[i].radius = 0.5;
    }
    else if(sphere_pm->ev_type == 2)   /* gaussian */
      monomers[i].radius = sphere_pm->Ss;
    else {
      printf("wrong EV type!!\n");
      exit(1);
    }

    monomers[i].rho = 1.0;
    monomers[i].dr2 = 0.0;
    monomers[i].blist[0][0] = 0;
    monomers[i].blist[0][1] = 0;
    monomers[i].blist[0][2] = 0;

    for(j=0; j<DIMS; j++) {
      monomers[i].vel[j]=0.0;
      monomers[i].force[j]=0.0;
      monomers[i].force0[j]=0.0;
      monomers[i].fluid_vel[j]=0.0;
    }
  }

  /* determine bond length (stretched length) */
  if(sphere_pm->spring == 0)  /* FENE */
    // bondlength = sphere_pm->Q_fene*(monomers[0].radius*2.0)*MAX_EXT;
    bondlength = 1.0 * monomers[0].radius *2.0; /* no stretch */
  else if(sphere_pm->spring == 1) /* WLC */
    bondlength = sphere_pm->nks*sphere_pm->sigma_k*MAX_EXT;
  else if(sphere_pm->spring == 2)  /* harmonic */
    bondlength = 1.0 * monomers[0].radius * 2.0;
  else {
    printf("wrong spring type!!\n");
    exit(1);
  }

  for (i=0 ; i<NTYPES ; i++) {
    if (nlevel[i] == 0)
      radius[i] = bondlength/4.*sqrt(10.+2.*sqrt(5.)); /* radius of icosahedron (around 0.95*bondlength) */
    else if (nlevel[i] == 1) 
      radius[i] = 1.71*bondlength;  /* measured relaxation radius */
    else if (nlevel[i] == 2 || nlevel[i] == -1)
      radius[i] = 3.30*bondlength;  /* measured relaxation radius */
    else if (nlevel[i] == 3 || nlevel[i] == -2)
      radius[i] = 5.70*bondlength;  /* speculated radius */
    else
      radius[i] = 5.70*pow(1.98,nlevel[i]-3)*bondlength;
  }
  
  if (test)
    printf("bondlength=%lf, radius[0]=%lf, radius[1]=%lf\n", bondlength, radius[0], radius[1]);

  /* define rest volume and surface area */
  for(i=0; i<NTYPES; i++) {
    if (nlevel[i] <= -1 || nlevel[i] >=2) {
      if(nlevel[i] <= -1)  {
	Vfactor[i] = 0.6 * bondlength*bondlength*bondlength;
	Sfactor[i] = bondlength*bondlength;
      }
      else {
	Vfactor[i] = 1.0 * bondlength*bondlength*bondlength;
	Sfactor[i] = bondlength*bondlength;
      }

      if(sphere_pm->spring == 0) {
	if(sphere_pm->N_per_DP[i] == 162) {
	  var = log(sphere_pm->H_fene[i])/3.0;
	  sphere_pm->V0[i] = (59.6+231.0*exp(-0.29*exp(var))) * Vfactor[i];
	  sphere_pm->A0[i] = (75.4+138.328*exp(-0.25*exp(var))) * Sfactor[i];
	}
	else if(sphere_pm->N_per_DP[i] == 642) {
	  sphere_pm->V0[i] = (1811.5-187.77*log(sphere_pm->H_fene[i])) * Vfactor[i];
	  sphere_pm->A0[i] = (724.14-56.592*log(sphere_pm->H_fene[i])) * Sfactor[i];
	}
      }
      else if(sphere_pm->spring == 2) {
	if(sphere_pm->N_per_DP[i] == 162) {
	  var = log(sphere_pm->H_fene[i])/3.0;
	  sphere_pm->V0[i] = (151.0+155.55*exp(-0.384*exp(var))) * Vfactor[i];
	  sphere_pm->A0[i] = (139.6+83.85*exp(-0.388*exp(var))) * Sfactor[i];
	}
	else if(sphere_pm->N_per_DP[i] == 642) {
	  sphere_pm->V0[i] = (1885.3 - 102.73*log(sphere_pm->H_fene[i])) * Vfactor[i];
	  sphere_pm->A0[i] = (738.43 - 27.447*log(sphere_pm->H_fene[i])) * Sfactor[i];  
	}
      }
      else if(sphere_pm->spring == 3) {
	if(sphere_pm->N_per_DP[i] == 162) {
	  if(sphere_pm->H_fene[i] < 200) {
	    sphere_pm->V0[i] = (266.45-16.611*log(sphere_pm->H_fene[i]))*Vfactor[i];
	    sphere_pm->A0[i] = (201.31-8.8242*log(sphere_pm->H_fene[i]))*Sfactor[i];
	  }
	  else {
	    sphere_pm->V0[i] = (229.84 - 10.078*log(sphere_pm->H_fene[i])) * Vfactor[i];
	    sphere_pm->A0[i] = (179.89 - 5.2169*log(sphere_pm->H_fene[i])) * Sfactor[i];
	  }
	}
	else if(sphere_pm->N_per_DP[i] == 642) {
	  sphere_pm->V0[i] = 2100.0 *Vfactor[i];
	  sphere_pm->A0[i] = 793.862 *Sfactor[i];
	}
      }
    }
    else if (nlevel[i] == 0) {
      sphere_pm->V0[i] = 2.07121;
      sphere_pm->A0[i] = 8.365;
    }
    else if (nlevel[i] == 1) {
      sphere_pm->V0[i] = 17.1369;
      sphere_pm->A0[i] = 33.5969;
    }
    else {
      fprintf(stderr, "rest volume for nlevel=%d is unknown!", nlevel[i]);
      exit(1);
    }

    //    printf("V0[%d] = %le, A0[%d] = %le\n", i, sphere_pm->V0[i], i, sphere_pm->A0[i]);
  }

  /* generate a new config */
  if (sphere_pm->initconfig != 4) {
    /* monomer index */
    m = 0;
    /* start to grow spheres */
    for(n=0; n<NDP; n++) {
      type = (n<Ntype[0] ? 0 : 1);
      start = (type==0 ? n*Nbeads[0] : Ntype[0]*Nbeads[0] + (n-Ntype[0])*Nbeads[1]);

      if (test)
        printf("Adding sphere %d\n", n);

//      printf("%d %d radius=%lf %lf trial=%d\n", n, m, radius[0], radius[1], trial_count);

      /* set up vertices positions for RBC and sphere */
      if(nlevel[type] <= -1)
	setup_template(v[type], blist_rbc, Nbeads[type], radius[type], work_dir);
      else
	/* determine vertices position and bonding list for the icosahedron */
	setup(v, blist_ico, radius);

      if(trial_count > 1000000) {
	n = 0;
	m = 0;
	radius[0] = radius[0]*0.98;
	radius[1] = radius[1]*0.98;
	
	start = 0;
	type = 0;
	
	/* determine vertices position from template file for RBC*/
	if(nlevel[i] <= -1)
	  setup_template(v[type], blist_rbc, MRBC, radius[type], work_dir);
	/* determine vertices position and bonding list for the icosahedron */
	else 
	  setup(v, blist_ico, radius);
      }
      
      trial_count = 0;
      /* find a center with enough space around */
      do {
        if (trial_count > 1000000) {
          if (n_proc == 0) {
            printf("Tried 1000000 times but cannnot insert DP %d!\n", n);
	    
	    break;
	  }
        }

	/* choose center of the sphere */
	/* randomly between 1 and maxsize[j]-2 */
	for(j=0; j<DIMS; j++)
	  centers[n][j]=(double)rand()/((double)(RAND_MAX)+1.0)*(maxsize[j]-2.0)+1.0;
	
	if(sphere_pm->initconfig == 2 && NDP == 1) {
	  for(j=0; j<DIMS; j++)
	    centers[n][j]=(maxsize[j]-1.0-2.0)/2.0+1.0;
	}
	else if(sphere_pm->initconfig == 3) {
	  if(NDP == 1) {
	    for(j=0; j<DIMS; j++)
	    centers[n][j]=(maxsize[j]-1.0-2.0)/2.0+1.0;
	    centers[n][2] = (maxsize[2]-radius[type]-1.5);
	  }
	  else if(NDP == 2) {
	    for(j=0; j<DIMS; j++)
	      centers[n][j]=(maxsize[j]-1.0-2.0)/2.0+1.0;
	    centers[0][1] += 2.0;
	    centers[0][0] -= radius[type]*1.5;
	    centers[1][1] -= 2.0;
	    centers[1][0] += radius[type]*1.5;
	  }
	}
	/* Use regular formed cluster configuration */
	else if(sphere_pm->initconfig == 5) {
	  rshift=0.55;
	  if(NDP == 1) {
	    for(j=1; j<DIMS; j++)
	      centers[n][j]=(maxsize[j]-2.0)/2.0+1.0;
	    centers[n][0]=radius[type]*3.0;
	  }

	  else if(NDP == 2) {
	    for(j=0; j<DIMS; j++)
	      centers[n][j]=(maxsize[j]-2.0)/2.0+1.0;
	    if(n == 0)
	      centers[0][0]=radius[type]*1.75;
	    else if(n == 1)
	      centers[1][0]=centers[0][0]+2.0*(radius[type]+rshift);
	  }
	  /* triangle */
	  else if(NDP == 3) {
	    for(j=0; j<DIMS; j++)
	      centers[n][j]=(maxsize[j]-2.0)/2.0+1.0;
	    centers[n][0]=radius[type]*3.5;
	    if(n == 2) {
	      centers[2][0]+=(radius[type]+rshift);
	      centers[2][1]-=sqrt(3.0)*(radius[type]+rshift);
	    }
	    else if(n==1) {
	      centers[1][0]-=(radius[type]+rshift);
	      centers[1][1]-=sqrt(3.0)*(radius[type]+rshift);
	    }
	  }
	  /* tetrahedral */
	  else if(NDP == 4) {
	    for(j=0; j<DIMS; j++)
	      centers[n][j]=(maxsize[j]-2.0)/2.0+1.0;
	    centers[n][0]=radius[type]*3.0;
	    if(n == 0) {
	      centers[0][1]+= (radius[type]+rshift);
	      centers[0][2]+= ((radius[type]+rshift)/sqrt(2.0));
	    }
	    else if(n==1) {
	      centers[1][0]-= (radius[type]+rshift);
	      centers[1][2]-= ((radius[type]+rshift)/sqrt(2.0));
	    }	      
	    else if(n == 2) {
	      centers[2][0]+= (radius[type]+rshift);
	      centers[2][2]-= ((radius[type]+rshift)/sqrt(2.0));
	    }
	    else if(n == 3) {
	      centers[3][1]-= (radius[type]+rshift);
	      centers[3][2]+= ((radius[type]+rshift)/sqrt(2.0));
	    }
	  }
	  /* cube */
	  else if(NDP == 8) {
	    for(j=1; j<DIMS; j++)
	      centers[n][j]=(maxsize[j]-2.0)/2.0+1.0-(radius[type]+rshift);
	    centers[n][0]=radius[type]+rshift;
	    if(n==1) 
	      centers[1][0]+= 2.0*(radius[type]+rshift);
	    else if(n == 2) 
	      centers[2][1]+= 2.0*(radius[type]+rshift);
	    else if(n == 3) {
	      centers[3][0]+= 2.0*(radius[type]+rshift);
	      centers[3][1]+= 2.0*(radius[type]+rshift);
	    }
	    else if(n == 4) {
	      centers[4][1]+= 2.0*(radius[type]+rshift);
	      centers[4][2]+= 2.0*(radius[type]+rshift);
	    }
	    else if(n == 5) {
	      centers[5][0]+= 2.0*(radius[type]+rshift);
	      centers[5][1]+= 2.0*(radius[type]+rshift);
	      centers[5][2]+= 2.0*(radius[type]+rshift);
	    }
	    else if(n == 6) 
	      centers[6][2]+= 2.0*(radius[type]+rshift);
	    else if(n == 7) {
	      centers[7][0]+= 2.0*(radius[type]+rshift);
	      centers[7][2]+= 2.0*(radius[type]+rshift);
	    }
	  }
	}
	/* Use flattened cluster configuration */
	else if(sphere_pm->initconfig == 6) {
	  prefactor = 2.8;
	  rshift=0.55; 

	  if(NDP == 1) {
	    for(j=1; j<DIMS; j++)
	      centers[n][j]=(maxsize[j]-2.0)/2.0+1.0;
	    centers[n][0]=radius[type]*3.0;
	  }

	  else if(NDP == 2) {
	    for(j=0; j<DIMS; j++)
	      centers[n][j]=(maxsize[j]-2.0)/2.0+1.0;
	    if(n == 0)
	      centers[0][0]=radius[type]*1.75;
	    else if(n == 1)
	      centers[1][0]=centers[0][0]+prefactor*(radius[type]+rshift);
	  }
	  /* triangle */
	  else if(NDP == 3) {
	    for(j=0; j<DIMS; j++)
	      centers[n][j]=(maxsize[j]-2.0)/2.0+1.0;
	    centers[n][0]=radius[type]*1.75;
	    if(n == 2) {
	      centers[2][0]+=(radius[type]+rshift)* prefactor*sqrt(3.0)/2.0;
	      centers[2][1]+=1.0/2.0*(radius[type]+rshift)*prefactor;
	    }
	    else if(n==1) {
	      centers[1][0]+=(radius[type]+rshift)* prefactor*sqrt(3.0)/2.0;
	      centers[1][1]-=1.0/2.0*(radius[type]+rshift)*prefactor;
	    }
	  }
	  /* tetrahedral */
	  else if(NDP == 4) {
	    for(j=0; j<DIMS; j++)
	      centers[n][j]=(maxsize[j]-1.0)/2.0;
	    centers[n][0]=radius[type]*1.5;
	    if(n==1) 
	      centers[1][0]+= prefactor*sqrt(3.0)*(radius[type]+rshift);
	    else if(n == 2) {
	      centers[2][0]+=prefactor*sqrt(3.0)/2.0*(radius[type]+rshift);
	      centers[2][1]+=prefactor/2.0*(radius[type]+rshift);
	    }
	    else if(n == 3) {
	      centers[3][0]+= prefactor*sqrt(3.0)/2.0*(radius[type]+rshift);
	      centers[3][1]-= prefactor/2.0*(radius[type]+rshift);
	    }
	  }
	  /* cube */
	  else if(NDP == 8) {
	    for(j=1; j<DIMS; j++)
	      centers[n][j]=(maxsize[j]-2.0)/2.0+1.0-(radius[type]+rshift);
	    centers[n][0]=(radius[type]+rshift)*4.0;
	    if(n==1) 
	      centers[1][0]+= 3.0*(radius[type]+rshift);
	    else if(n == 2) 
	      centers[2][1]+= 3.0*(radius[type]+rshift);
	    else if(n == 3) {
	      centers[3][0]+= 3.0*(radius[type]+rshift);
	      centers[3][1]+= 3.0*(radius[type]+rshift);
	    }
	    else if(n == 4) 
	      centers[4][0]-= 3.0*(radius[type]+rshift);
	    else if(n == 5) 
	      centers[5][1]-= 3.0*(radius[type]+rshift);
	    else if(n == 6) {
	      centers[6][0]+= 3.0*(radius[type]+rshift);
	      centers[6][1]-= 3.0*(radius[type]+rshift);
	    }
	    else if(n == 7) {
	      centers[7][0]-= 3.0*(radius[type]+rshift);
	      centers[7][1]+= 3.0*(radius[type]+rshift);
	    }
	  }
	}

        overlap = FALSE;
        trial_count++;

        /* check sphere-sphere overlap */
        for (n1=0 ; n1<n ; n1++) {

          type1 = (n1<Ntype[0] ? 0 : 1);

          dx = centers[n][0] - centers[n1][0];
          dy = centers[n][1] - centers[n1][1];
          dz = centers[n][2] - centers[n1][2];

          if (wall_flag < 3)  /* no x-wall */
            dx = n_image(dx,max_x);
          if (wall_flag < 2)  /* no z-wall */
            dz = n_image(dz,max_z);
          if (wall_flag < 1)  /* no y-wall */
            dy = n_image(dy,max_y);

          r2 = dx*dx + dy*dy + dz*dz;

          if (sqrt(r2) <= (radius[type]+radius[type1])*1.0+1.06)
            overlap = TRUE;
        }
	
        /* check wall-overlap */
	if(sphere_pm->initconfig < 5) {
	  if (wall_flag >= 3)  /* x-wall */
	    if (centers[n][0] + radius[type]*1.0 + 1.05 > max_x - 1 || centers[n][0] - radius[type]*1.0 - 1.05 < 0)
	      overlap = TRUE;
	  if (wall_flag >= 2)  /* z-wall */
	    if (centers[n][2] + radius[type]*1.0 + 1.05 > max_z - 1 || centers[n][2] - radius[type]*1.0 - 1.05 < 0)
	      overlap = TRUE;
	  if (wall_flag >= 1)  /* y-wall */
	    if (centers[n][1] + radius[type]*1.0 + 1.05 > max_y - 1 || centers[n][1] - radius[type]*1.0 - 1.05 < 0)
	      overlap = TRUE;
	}

	if(cyl_pm->initconfig == 6) 
	  if ((centers[n][1] + radius[type]*1.0 + 1.05) > (cyl_pm->offset+cyl_pm->ramp*centers[n][2]))
	    overlap = TRUE;
	
      } while (overlap == TRUE);

      if(trial_count > 10000000) 
	continue;

      if (test)
        printf("sphere %d centers at (%le,%le,%le) radius %le\n", n, centers[n][0], centers[n][1], centers[n][2], radius[type]);

      /* randomly rotate the particle, make sure it stays inside the box */
      rot_flag = 1;
      while(rot_flag == 1) {
	/* choose a direction for the z-axis of the sphere */
	phi = (double)rand()/(double)(RAND_MAX)*2.*M_PI;
	theta = (double)rand()/(double)(RAND_MAX)*M_PI;
	psi = (double)rand()/(double)(RAND_MAX)*M_PI;

	phi = M_PI/2.0;
	theta = 0.0;
	psi = 0.0;

	/* /\* z axis *\/ */
	/* axis[2][0] = cos(phi)*sin(theta); */
	/* axis[2][1] = sin(phi)*sin(theta); */
	/* axis[2][2] = cos(theta); */

	/* /\* y axis is chosen arbitrarily by rotating z-axis down with pi/2 *\/ */
	/* axis[1][0] = cos(phi)*sin(theta+M_PI/2.); */
	/* axis[1][1] = sin(phi)*sin(theta+M_PI/2.); */
	/* axis[1][2] = cos(theta+M_PI/2.); */

	/* /\* x axis is dicided by a cross product *\/ */
	/* product(axis[1], axis[2], axis[0]); */

	/* z axis */
	axis[2][0] = -sin(theta);
	axis[2][1] = sin(phi)*cos(theta);
	axis[2][2] = cos(phi)*cos(theta);

	/* y axis */
	axis[1][0] = cos(theta)*sin(psi);
	axis[1][1] = cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi);
	axis[1][2] = -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi);
	
	/* x axis */
	axis[0][0] = cos(theta)*cos(psi);
	axis[0][1] = -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
	axis[0][2] = sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi);
	
	/* normalize again in case there is numerical error */
	for (j=0 ; j<DIMS ; j++) {
	  r2 = 0.;
	  for (k=0 ; k<DIMS ; k++)
	    r2 += axis[j][k]*axis[j][k];
	  r2 = sqrt(r2);
	  for (k=0 ; k<DIMS ; k++)
	    axis[j][k] /= r2;
	}

	/* insert a template */
	if(nlevel[type] <= -1) {
	  for (i=0 ; i<Nbeads[type] ; i++) {
	    m = start+i;
	    for (j=0 ; j<DIMS ; j++) {
	      monomers[m].pos_pbc[j] = centers[n][j];
	      for (k=0 ; k<DIMS ; k++) 
		monomers[m].pos_pbc[j] += v[type][i][k]*axis[k][j];
	      monomers[m].pos[j] = monomers[m].pos_pbc[j];
	    }

	    //	    printf("m=%d/%d pos=(%le %le) (%le %le) rotflag = %d\n", i, Nbeads[type],  monomers[m].pos_pbc[1], monomers[m].pos_pbc[2], phi, theta, rot_flag);

	    monomers[m].pos[0] = box(monomers[m].pos_pbc[0],maxsize[0]);
	    if(wall_flag < 2) 
	      monomers[m].pos[2] = box(monomers[m].pos_pbc[2],maxsize[2]);
	    if(wall_flag < 1)
	      monomers[m].pos[1] = box(monomers[m].pos_pbc[1],maxsize[1]);
	    
	    /* check if the rotated vertex position go out of walls */
	    if(wall_flag ==0)
	      rot_flag = 0;
	    else {
	      for(j=1; j<=wall_flag; j++) {
		if (monomers[m].pos_pbc[j] > maxsize[j] - 1.0 || monomers[m].pos_pbc[j] < 0.0) {
		  rot_flag = 1;
		  break;
		}
		else 
		  rot_flag = 0;
	      }
	    }

	    if(rot_flag == 1) break;

	    monomers[m].blist[0][0] = blist_rbc[i][0][0];
	    for (j=1 ; j<=blist_rbc[i][0][0] ; j++) {
	      monomers[m].blist[j][0] = start + blist_rbc[i][j][0];
	      monomers[m].blist[j][1] = start + blist_rbc[i][j][1];
	      monomers[m].blist[j][2] = start + blist_rbc[i][j][2];
	    }
	    monomers[m].DP_id = n;
	  }
	}
	/* insert an icosahedron */
	else {
	  m = start;
	  rot_flag = 0;

	  for (i=0 ; i<NVERTICES ; i++) {
	    for (j=0 ; j<DIMS ; j++) {
	      monomers[m].pos_pbc[j] = centers[n][j];
	      for (k=0 ; k<DIMS ; k++)
		monomers[m].pos_pbc[j] += v[type][i][k]*axis[k][j];
	      monomers[m].pos[j] = box(monomers[m].pos_pbc[j],maxsize[j]);
	    }
	    bcount = 0;
	    for (j=1 ; j<=NBONDS ; j++) {
	      k = start + blist_ico[i][j][0];  /* m is bonded to k */
	      /* add k to the blist of m only if m is not already in the blist of k */
	      if (!check_bond(k,m,monomers)) {
		bcount++;
		monomers[m].blist[bcount][0] = k;
		monomers[m].blist[bcount][1] = start + blist_ico[i][j][1];
		monomers[m].blist[bcount][2] = start + blist_ico[i][j][2];
	      }
	    }
	    monomers[m].blist[0][0] = bcount;
	    monomers[m].DP_id = n;
	    monomers[m].type = type;
	    //	    printf("n=%d start=%d m=%d monpos=%le %le %le\n", n, start, m, monomers[m].pos_pbc[0], monomers[m].pos_pbc[1], monomers[m].pos_pbc[2]);
	    m++;
	  }

	  /* sphere triangulation */
	  for (level=0 ; level<nlevel[type] ; level++) {
	    m_temp = m;

	    /* refine the mesh by creating additional points */
	    for (n1=start ; n1<m_temp ; n1++)
	      for (j=1 ; j<=monomers[n1].blist[0][0] ; j++) {  /* loop over all bonds */
		n2 = monomers[n1].blist[j][0];
		/* locate the midpoint */
		r2 = 0.;
		for (k=0 ; k<DIMS ; k++) {
		  p_mid[k] = (monomers[n1].pos_pbc[k]+monomers[n2].pos_pbc[k])/2.;
		  dr[k] = p_mid[k] - centers[n][k];
		  r2 += dr[k]*dr[k];
		}
		r2 = sqrt(r2);
		/* push out the midpoint to the sphere */
		for (k=0 ; k<DIMS ; k++) {
		  monomers[m].pos_pbc[k] = centers[n][k] + dr[k]*radius[type]/r2;
		  monomers[m].pos[k] = box(monomers[m].pos_pbc[k], maxsize[k]);
		}
		monomers[m].DP_id = n;
		monomers[m].type = type;
		midpoints[n1][j] = m;
		//		printf("n=%d start=%d m=%d monpos=%le %le %le\n", n, start, m, monomers[m].pos_pbc[0], monomers[m].pos_pbc[1], monomers[m].pos_pbc[2]);
		m++;
	      }

	    /* set up blist for the new points and update blist for the existing points */
	    for (n1=start ; n1<m_temp ; n1++)
	      for (j=1 ; j<=monomers[n1].blist[0][0] ; j++) {  /* loop over all bonds */
		/* This is the relative position of the points (facing outside) */
		n2 = monomers[n1].blist[j][0];                   /*              4               */
		n3 = monomers[n1].blist[j][1];                   /*              /\              */
		n4 = monomers[n1].blist[j][2];                   /*          14 /__\ 24          */
		                                                 /*            /\  /\            */
		n12 = midpoints[n1][j];                          /*        1  /__12__\  2        */
		n13 = get_midpoint(n1, n3, monomers, midpoints); /*           \  /\  /           */
		n14 = get_midpoint(n1, n4, monomers, midpoints); /*            \/__\/ 23         */
		n23 = get_midpoint(n2, n3, monomers, midpoints); /*          13 \  /             */
		n24 = get_midpoint(n2, n4, monomers, midpoints); /*              \/              */
		                                                 /*              3               */            
                                                       		 /*                              */
		/* update blist: the n1-n2 bond now becomes n1-n12 bond */
		/* we still need the old blist for getting midpoints. So we store the update in blist_temp and recover it later */
		blist_temp[n1][j][0] = n12;
		blist_temp[n1][j][1] = n13;
		blist_temp[n1][j][2] = n14;

		/* set up blist for n12 */
		bcount = 1;
		monomers[n12].blist[bcount][0] = n2;
		monomers[n12].blist[bcount][1] = n23;
		monomers[n12].blist[bcount][2] = n24;
		if (!check_bond(n13,n12,monomers)) {
		  bcount++;
		  monomers[n12].blist[bcount][0] = n13;
		  monomers[n12].blist[bcount][1] = n1;
		  monomers[n12].blist[bcount][2] = n23;
		}
		if (!check_bond(n23,n12,monomers)) {
		  bcount++;
		  monomers[n12].blist[bcount][0] = n23;
		  monomers[n12].blist[bcount][1] = n13;
		  monomers[n12].blist[bcount][2] = n2;
		}
		if (!check_bond(n14,n12,monomers)) {
		  bcount++;
		  monomers[n12].blist[bcount][0] = n14;
		  monomers[n12].blist[bcount][1] = n24;
		  monomers[n12].blist[bcount][2] = n1;
		}
		if (!check_bond(n24,n12,monomers)) {
		  bcount++;
		  monomers[n12].blist[bcount][0] = n24;
		  monomers[n12].blist[bcount][1] = n2;
		  monomers[n12].blist[bcount][2] = n14;
		}
		monomers[n12].blist[0][0] = bcount;
	      }

	    /* recover blist for the old points */
	    for (n1=start ; n1<m_temp ; n1++)
	      for (j=1 ; j<=monomers[n1].blist[0][0] ; j++)
		for (k=0 ; k<3 ; k++)
		  monomers[n1].blist[j][k] = blist_temp[n1][j][k];
	  }
	}
      }
    }

    /* write config to a file */
    if (test) {
      file_name (test_file, work_dir, -1);
      file_ptr = fopen (test_file, "w");
      if (file_ptr == 0)  fatal_err("Could not open test.dat", -1);
      for (m=0 ; m<num_beads ; m++)
        fprintf(file_ptr, "%le %le %le\n", monomers[m].pos[0], monomers[m].pos[1], monomers[m].pos[2]);
      fclose(file_ptr);
    }
  }
  /* read in the initial sphere configuration from file */
  else if(sphere_pm->initconfig == 4) {
    sprintf(filename, "%s/init/init.config", work_dir);
    printf("init config file %s\n", filename);
    stream = fopen(filename, "r");      
    fscanf(stream, "%d %d %d %d %d", &sphere_pm->Ntype[0], &sphere_pm->Ntype[1], &sphere_pm->nlevel[0], &sphere_pm->nlevel[1], &sphere_pm->num_beads);

    sphere_pm->NDP = sphere_pm->Ntype[0] + sphere_pm->Ntype[1];
    NDP = sphere_pm->NDP;
    num_beads = sphere_pm->num_beads;
    
    for (j=0 ; j<NTYPES ; j++) {
     if(sphere_pm->nlevel[j] >= 0) {
	sphere_pm->N_per_DP[j] = 12;
	sphere_pm->face_per_DP[j] = 20;
	temp = 30;
	for (i=0 ; i<sphere_pm->nlevel[j] ; i++) {
	  sphere_pm->N_per_DP[j] += temp;
	  sphere_pm->face_per_DP[j] *= 4;
	  temp *= 4;
	}
      }
      else {
	if(sphere_pm->nlevel[j] == -1) {
	  sphere_pm->N_per_DP[j] = 162;
	  sphere_pm->face_per_DP[j] = 320;
	}
	if(sphere_pm->nlevel[j] == -2) {
	  sphere_pm->N_per_DP[j] = 642;
	  sphere_pm->face_per_DP[j] = 1280;
	}
      }
    }
    
    printf("Type %d %d beads/DP %d faces/DP, Type %d %d beads/DP %d faces/DP\n", 0, sphere_pm->N_per_DP[0], sphere_pm->face_per_DP[0], 1, sphere_pm->N_per_DP[1], sphere_pm->face_per_DP[1]);

    for (i=0 ; i<NTYPES ; i++) {
      Ntype[i] = sphere_pm->Ntype[i];
      nlevel[i] = sphere_pm->nlevel[i];
      Nfaces[i] = sphere_pm->face_per_DP[i];
      Nbeads[i] = sphere_pm->N_per_DP[i];
    }

    for(i=0; i<NDP; i++) {
      type = (i<Ntype[0] ? 0 : 1);
      fscanf(stream, "%le %le %le %le", &sphere_pm->V0[type], &sphere_pm->A0[type], &spheres[i].disp2, &spheres[i].dr2);
    }

    for(i=0; i<NDP; i++) {
      type = (i<Ntype[0] ? 0 : 1);
      start = (type == 0 ? i*Nbeads[0] : Ntype[0]*Nbeads[0] + (i-Ntype[0])*Nbeads[1]);
      if(n_proc == 0)
        printf("#sphere %d \n", i);
      for(n=0; n<Nbeads[type]; n++) {
        pnum = start+n;
        fscanf(stream, "%le %le %le %le %le %le %le %le %le %le", &monomers[pnum].pos[0], &monomers[pnum].pos[1], &monomers[pnum].pos[2], &monomers[pnum].pos_pbc[0], &monomers[pnum].pos_pbc[1], &monomers[pnum].pos_pbc[2], &monomers[pnum].vel[0], &monomers[pnum].vel[1], &monomers[pnum].vel[2], &monomers[pnum].dr2);

        monomers[pnum].DP_id=i;
        monomers[pnum].type = type;

        if(n_proc == 0)
          printf("#monomer %d at (%le %le %le) (%le %le %le)\n", pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], monomers[pnum].pos[2], monomers[pnum].pos_pbc[0], monomers[pnum].pos_pbc[1], monomers[pnum].pos_pbc[2]);
      }
    }
    fclose(stream);
  
    /*read in bond list*/
    sprintf(filename, "%s/init/bond.dat", work_dir);
    stream = fopen(filename, "r");
    if(Ntype[0] > 0) {
      for(i=0; i<Nbeads[0]; i++) {
	fscanf(stream, "%*s %d\n",  &monomers[i].blist[0][0]);
	for(j=1; j<=monomers[i].blist[0][0]; j++)
	  fscanf(stream, "%d %d %d\n", &monomers[i].blist[j][0], &monomers[i].blist[j][1], &monomers[i].blist[j][2]);
      }
    }
    if(Ntype[1] > 0) {
      for(i=0; i<Nbeads[1]; i++) {
	n=Ntype[0]*Nbeads[0]+i;
	fscanf(stream, "%*s %d\n",  &monomers[n].blist[0][0]);
	for(j=1; j<=monomers[n].blist[0][0]; j++)
	  fscanf(stream, "%d %d %d\n", &monomers[n].blist[j][0], &monomers[n].blist[j][1], &monomers[n].blist[j][2]);
      }
    }
    fclose(stream);

    for(i=1; i<NDP; i++) {
      type = (i<Ntype[0] ? 0 : 1);
      start0 = (type == 0 ? 0 : Ntype[0]*Nbeads[0]);
      start = (type == 0 ? i*Nbeads[0] : Ntype[0]*Nbeads[0] + (i-Ntype[0])*Nbeads[1]);
      for(j=0; j<Nbeads[type]; j++) {
	monomers[start+j].blist[0][0] = monomers[start0+j].blist[0][0];
	for(n=1; n<=monomers[j].blist[0][0]; n++) {
	  monomers[start+j].blist[n][0] = start+(monomers[start0+j].blist[n][0]-start0);
	  monomers[start+j].blist[n][1] = start+(monomers[start0+j].blist[n][1]-start0);
	  monomers[start+j].blist[n][2] = start+(monomers[start0+j].blist[n][2]-start0);
	}
      }
    }
  }

  sprintf(filename, "%s/init/bond.dat", work_dir);
  stream = fopen(filename, "w");
  if(Ntype[0] > 0) {
    for(i=0; i<Nbeads[0]; i++) {
      fprintf(stream, "%d %d\n",  i, monomers[i].blist[0][0]);
      for(j=1; j<=monomers[i].blist[0][0]; j++)
	fprintf(stream, "%d %d %d\n", monomers[i].blist[j][0], monomers[i].blist[j][1], monomers[i].blist[j][2]);
    }
  }

  if(Ntype[1] > 0) {  
    for(i=0; i<Nbeads[1]; i++) {
      n=Ntype[0]*Nbeads[0]+i;
      fprintf(stream, "%d %d\n",  n, monomers[n].blist[0][0]);
      for(j=1; j<=monomers[n].blist[0][0]; j++)
	fprintf(stream, "%d %d %d\n", monomers[n].blist[j][0], monomers[n].blist[j][1], monomers[n].blist[j][2]);
    }
  }
  fclose(stream);
  
  /* add faces  */ 
  f=0;

  for(i=0; i<num_beads; i++) 
    monomers[i].face_id[0]=0;    /* tracks the number of faces the monomer is on */

  for(i=0; i<num_beads; i++) {
    for(j=1; j<=monomers[i].blist[0][0]; j++) {
      if(!check_face(i, monomers[i].blist[j][0], monomers[i].blist[j][1], faces, f)) {
	if(i < Ntype[0]*Nbeads[0])
	  faces[f].DP_id = i/Nbeads[0];
	else
	  faces[f].DP_id = Ntype[0]+(i-Ntype[0]*Nbeads[0])/Nbeads[1];

	faces[f].vertices[0] = i;
	faces[f].vertices[1] = monomers[i].blist[j][0];
	faces[f].vertices[2] = monomers[i].blist[j][1];
	
	for(d=0; d<DIMS; d++) {
	  addface=1;
	  for(n=1; n<=monomers[faces[f].vertices[d]].face_id[0]; n++) 
	    if(monomers[faces[f].vertices[d]].face_id[n] == f)
	      addface = 0;
	  
	  if(addface == 1) {
	    monomers[faces[f].vertices[d]].face_id[0]++;
	    monomers[faces[f].vertices[d]].face_id[monomers[faces[f].vertices[d]].face_id[0]]=f;
	  }	  
	}
	
	f++;
      }
      if(!check_face(i, monomers[i].blist[j][0], monomers[i].blist[j][2], faces, f)) {
	if(i < Ntype[0]*Nbeads[0])
	  faces[f].DP_id = i/Nbeads[0];
	else
	  faces[f].DP_id = Ntype[0]+(i-Ntype[0]*Nbeads[0])/Nbeads[1];

	faces[f].vertices[0] = i;
	faces[f].vertices[1] = monomers[i].blist[j][0];
	faces[f].vertices[2] = monomers[i].blist[j][2];

	for(d=0; d<DIMS; d++) {
	  addface=1;
	  for(n=1; n<=monomers[faces[f].vertices[d]].face_id[0]; n++) 
	    if(monomers[faces[f].vertices[d]].face_id[n] == f)
	      addface = 0;
	  
	  if(addface == 1) {
	    monomers[faces[f].vertices[d]].face_id[0]++;
	    monomers[faces[f].vertices[d]].face_id[monomers[faces[f].vertices[d]].face_id[0]]=f;
	  }	  
	}

	f++;
      }
    }
  }

  printf("%d faces added. %d faces / sphere0  %d faces / sphere1, total faces %d\n", f, Nfaces[0], Nfaces[1], Nfaces[0]*Ntype[0]+Nfaces[1]*Ntype[1]);

  if(f != Nfaces[0]*Ntype[0]+Nfaces[1]*Ntype[1])
    fatal_err("Number of faces does not match", -1);

  
  for(type=0; type<NTYPES; type++) { 
    for(i=0; i<maxNbeads; i++)
      free(v[type][i]);
  }
  
  free(v[0]);
  free(v[1]);

  return 1;

}

/* check if monomer m has been added to the blist of monomer k */
int check_bond(int k, int m, struct monomer *monomers)
{
  int i;
  for (i=1 ; i<=monomers[k].blist[0][0] ; i++)
    if (monomers[k].blist[i][0] == m)
      return 1;
  return 0;
}

/* check if a face already exists */
int check_face(int a, int b, int c, struct face *faces, int nface)
{
  int i,j,d,check;
  int aa, bb, cc;
  for(i=0; i<nface; i++) {
    aa = faces[i].vertices[0];
    bb = faces[i].vertices[1];
    cc = faces[i].vertices[2];
    
    if(aa == a)
      if((bb == b && cc == c) || ((bb == c) && (cc == b)))
	return 1;
    if(aa == b)
      if((bb == a && cc == c) || ((bb == c) && (cc == a)))
	return 1;
    if(aa == c)
      if((bb == a && cc == b) || ((bb == b) && (cc == a)))
	return 1;
  }
  return 0;
}

/* return the midpoint between a and b */
int get_midpoint(int a, int b, struct monomer *monomers, int **midpoints)
{
  int i;

  for (i=1 ; i<=monomers[a].blist[0][0] ; i++)
    if (monomers[a].blist[i][0] == b)  /* if a is bonded to b */
      return midpoints[a][i];
  for (i=1 ; i<=monomers[b].blist[0][0] ; i++)
    if (monomers[b].blist[i][0] == a)  /* if b is bonded to a */
      return midpoints[b][i];

  fprintf(stderr, "get_midpoint: monomer %d is not bonded to monomer %d!\n", a, b);
  exit(1);
  return -1;
}

/* determine vertices position and bonding list */
void setup(double ***v, int blist[NVERTICES][NBONDS+1][3], double radius[NTYPES])
{
  int i;
  double a, b, edgelength;

  for (i=0 ; i<NTYPES ; i++) {
    /* set the coordinates of an icosahedron of the given radius */
    edgelength = 4.*radius[i]/sqrt(10.+2.*sqrt(5.));
    a = edgelength/2.;
    b = edgelength*(1.+sqrt(5.))/4.;

    v[i][0][0] = 0.;    v[i][0][1] =  a;    v[i][0][2] =  b;
    v[i][1][0] = 0.;    v[i][1][1] =  a;    v[i][1][2] = -b;
    v[i][2][0] = 0.;    v[i][2][1] = -a;    v[i][2][2] =  b;
    v[i][3][0] = 0.;    v[i][3][1] = -a;    v[i][3][2] = -b;

    v[i][4][0] =  a;    v[i][4][1] =  b;    v[i][4][2] = 0.;
    v[i][5][0] =  a;    v[i][5][1] = -b;    v[i][5][2] = 0.;
    v[i][6][0] = -a;    v[i][6][1] =  b;    v[i][6][2] = 0.;
    v[i][7][0] = -a;    v[i][7][1] = -b;    v[i][7][2] = 0.;

    v[i][8][0] =  b;    v[i][8][1] = 0.;    v[i][8][2] =  a;
    v[i][9][0] = -b;    v[i][9][1] = 0.;    v[i][9][2] =  a;
    v[i][10][0] =  b;   v[i][10][1] = 0.;   v[i][10][2] = -a;
    v[i][11][0] = -b;   v[i][11][1] = 0.;   v[i][11][2] = -a;
  }

  /* set up the bonding list */
  /* blist[i][j][0] is the j-th neighbor of vertex i */
  /* blist[i][j][1] is the neighbor of vertex i that makes a triangle together with blist[i][j][0] */
  /* blist[i][j][2] is another such neighbor (see below) */
  /*                i
                   /|\
                  / | \
  blist[i][j][1] /  |  \
                 \  |  / blist[i][j][2]
                  \ | /
                   \|/
              blist[i][j][0]     */
  /* it is important that with right hand rule, i-blist[i][j][1]-blist[i][j][0] and i-blist[i][j][0]-blist[i][j][2] both point toward outside of the sphere!! */

  blist[0][1][0] = 2;   blist[0][1][1] = 9;   blist[0][1][2] = 8;
  blist[0][2][0] = 4;   blist[0][2][1] = 8;   blist[0][2][2] = 6;
  blist[0][3][0] = 6;   blist[0][3][1] = 4;   blist[0][3][2] = 9;
  blist[0][4][0] = 8;   blist[0][4][1] = 2;   blist[0][4][2] = 4;
  blist[0][5][0] = 9;   blist[0][5][1] = 6;   blist[0][5][2] = 2;

  blist[1][1][0] = 3;   blist[1][1][1] = 10;  blist[1][1][2] = 11;
  blist[1][2][0] = 4;   blist[1][2][1] = 6;   blist[1][2][2] = 10;
  blist[1][3][0] = 6;   blist[1][3][1] = 11;  blist[1][3][2] = 4;
  blist[1][4][0] = 10;  blist[1][4][1] = 4;   blist[1][4][2] = 3;
  blist[1][5][0] = 11;  blist[1][5][1] = 3;   blist[1][5][2] = 6;

  blist[2][1][0] = 0;   blist[2][1][1] = 8;   blist[2][1][2] = 9;
  blist[2][2][0] = 5;   blist[2][2][1] = 7;   blist[2][2][2] = 8;
  blist[2][3][0] = 7;   blist[2][3][1] = 9;   blist[2][3][2] = 5;
  blist[2][4][0] = 8;   blist[2][4][1] = 5;   blist[2][4][2] = 0;
  blist[2][5][0] = 9;   blist[2][5][1] = 0;   blist[2][5][2] = 7;

  blist[3][1][0] = 1;   blist[3][1][1] = 11;  blist[3][1][2] = 10;
  blist[3][2][0] = 5;   blist[3][2][1] = 10;  blist[3][2][2] = 7;
  blist[3][3][0] = 7;   blist[3][3][1] = 5;   blist[3][3][2] = 11;
  blist[3][4][0] = 10;  blist[3][4][1] = 1;   blist[3][4][2] = 5;
  blist[3][5][0] = 11;  blist[3][5][1] = 7;   blist[3][5][2] = 1;


  blist[4][1][0] = 0;   blist[4][1][1] = 6;   blist[4][1][2] = 8;
  blist[4][2][0] = 1;   blist[4][2][1] = 10;  blist[4][2][2] = 6;
  blist[4][3][0] = 6;   blist[4][3][1] = 1;   blist[4][3][2] = 0;
  blist[4][4][0] = 8;   blist[4][4][1] = 0;   blist[4][4][2] = 10;
  blist[4][5][0] = 10;  blist[4][5][1] = 8;   blist[4][5][2] = 1;

  blist[5][1][0] = 2;   blist[5][1][1] = 8;   blist[5][1][2] = 7;
  blist[5][2][0] = 3;   blist[5][2][1] = 7;   blist[5][2][2] = 10;
  blist[5][3][0] = 7;   blist[5][3][1] = 2;   blist[5][3][2] = 3;
  blist[5][4][0] = 8;   blist[5][4][1] = 10;  blist[5][4][2] = 2;
  blist[5][5][0] = 10;  blist[5][5][1] = 3;   blist[5][5][2] = 8;

  blist[6][1][0] = 0;   blist[6][1][1] = 9;   blist[6][1][2] = 4;
  blist[6][2][0] = 1;   blist[6][2][1] = 4;   blist[6][2][2] = 11;
  blist[6][3][0] = 4;   blist[6][3][1] = 0;   blist[6][3][2] = 1;
  blist[6][4][0] = 9;   blist[6][4][1] = 11;  blist[6][4][2] = 0;
  blist[6][5][0] = 11;  blist[6][5][1] = 1;   blist[6][5][2] = 9;

  blist[7][1][0] = 2;   blist[7][1][1] = 5;   blist[7][1][2] = 9;
  blist[7][2][0] = 3;   blist[7][2][1] = 11;  blist[7][2][2] = 5;
  blist[7][3][0] = 5;   blist[7][3][1] = 3;   blist[7][3][2] = 2;
  blist[7][4][0] = 9;   blist[7][4][1] = 2;   blist[7][4][2] = 11;
  blist[7][5][0] = 11;  blist[7][5][1] = 9;   blist[7][5][2] = 3;


  blist[8][1][0] = 0;   blist[8][1][1] = 4;   blist[8][1][2] = 2;
  blist[8][2][0] = 2;   blist[8][2][1] = 0;   blist[8][2][2] = 5;
  blist[8][3][0] = 4;   blist[8][3][1] = 10;  blist[8][3][2] = 0;
  blist[8][4][0] = 5;   blist[8][4][1] = 2;   blist[8][4][2] = 10;
  blist[8][5][0] = 10;  blist[8][5][1] = 5;   blist[8][5][2] = 4;

  blist[9][1][0] = 0;   blist[9][1][1] = 2;   blist[9][1][2] = 6;
  blist[9][2][0] = 2;   blist[9][2][1] = 7;   blist[9][2][2] = 0;
  blist[9][3][0] = 6;   blist[9][3][1] = 0;   blist[9][3][2] = 11;
  blist[9][4][0] = 7;   blist[9][4][1] = 11;  blist[9][4][2] = 2;
  blist[9][5][0] = 11;  blist[9][5][1] = 6;   blist[9][5][2] = 7;

  blist[10][1][0] = 1;  blist[10][1][1] = 3;  blist[10][1][2] = 4;
  blist[10][2][0] = 3;  blist[10][2][1] = 5;  blist[10][2][2] = 1;
  blist[10][3][0] = 4;  blist[10][3][1] = 1;  blist[10][3][2] = 8;
  blist[10][4][0] = 5;  blist[10][4][1] = 8;  blist[10][4][2] = 3;
  blist[10][5][0] = 8;  blist[10][5][1] = 4;  blist[10][5][2] = 5;

  blist[11][1][0] = 1;  blist[11][1][1] = 6;  blist[11][1][2] = 3;  
  blist[11][2][0] = 3;  blist[11][2][1] = 1;  blist[11][2][2] = 7;
  blist[11][3][0] = 6;  blist[11][3][1] = 9;  blist[11][3][2] = 1;
  blist[11][4][0] = 7;  blist[11][4][1] = 3;  blist[11][4][2] = 9;
  blist[11][5][0] = 9;  blist[11][5][1] = 7;  blist[11][5][2] = 6;
}


/* determine vertices position and bonding list */
void setup_template(double **v, int ***blist, int nvertices, double radius, char *work_dir)
{
  int i,j;
  double a, b, edgelength;
  char filename[MAX_NL];
  FILE *stream;

  /* set the coordinates of an icosahedron of the given radius */
  sprintf(filename, "%s/init/template.init", work_dir);
  stream =fopen(filename, "r");
  
  for(i=0; i<nvertices; i++) {
    fscanf(stream, "%le %le %le\n", &v[i][0], &v[i][1], &v[i][2]);
    //    v[i][0] *= (radius / 4.43);
    //    v[i][1] *= (radius / 4.43);
    //    v[i][2] *= (radius / 4.43);
  }

  for(i=0; i<nvertices; i++) {
    fscanf(stream, "%*s %d\n", &blist[i][0][0]);
    for(j=1; j<=blist[i][0][0]; j++)
      fscanf(stream, "%d %d %d", &blist[i][j][0], &blist[i][j][1], &blist[i][j][2]);
  }

  fclose(stream);
}
