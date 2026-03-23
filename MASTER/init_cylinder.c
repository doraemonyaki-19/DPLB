/* Model cylinder as point particles w/ EV jointed by vertical and perpenddicular springs */
/* linear (angular) forces are not added yet */
/* Coupled with LBE code to add hydrodynamic forces */ 

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

int cyl_init(int **node_map, struct DP_param *cyl_pm, int nDP, int nbeads, struct DP *cylinders, struct monomer *monomers, char *work_dir)
{
  extern int n_proc;
  extern int max_x, max_y, max_z;
  extern int wall_flag;

  int i,j,k,m, pnum;
  int n0 = nbeads;
  int n, nxy;
  int x[DIMS];
  int Ncylinder = 1;   /* number of cylinders */
  int n_circle = cyl_pm->num_beads/(max_x*cyl_pm->NDP);   /* number of beads per circle */
  int num_beads = n0+cyl_pm->num_beads; /* total number beads on cylinder */
  int Ntype[1];
  int n_height[1];
  int N_per_cylinder[1];
  int type, type2, start;
  int offset;
  int *bnodes;
  int *tmp_ip;

  double maxsize[DIMS];
  double bondlength, radius, angle;
  double rperiod=cyl_pm->rperiod;     /* period of the variation of the tube radius */
  double ramp = cyl_pm->ramp;         /* amplitude of radius variation */
  double height[1], maxlen[1];
  double dx, dy, dz, r2, xpos;
  double xm, xp, ym, yp, ymax, ymin;
  double slope;

  char filename[MAX_NL];
  FILE *stream=0;

  int test=TRUE;
  char test_file[MAX_NL] = {"data/init.dat"};
  FILE *file_ptr=0;

  maxsize[0]=(double)max_x;
  maxsize[1]=(double)max_y;
  maxsize[2]=(double)max_z;
  
  /* assign boundary link nodes */
  cyl_pm->nlinks = 2*((max_x+2)*(max_z+2)+(max_x+2)*(max_y+2)+(max_y+2)*(max_z+2))+1;
  printf("Max %d boundary links\n", cyl_pm->nlinks);
  cyl_pm->bnodes = imatrix(cyl_pm->nlinks, Num_Dir);

  for(n=0; n < cyl_pm->nlinks; n++) 
    for(i=0; i<Num_Dir; i++)
      cyl_pm->bnodes[n][i] = 0;

  Ntype[0] = 1;
  n_height[0] = max_x;
  N_per_cylinder[0] = cyl_pm->num_beads;
  radius = cyl_pm->radius;
  angle = 1.0/(double)n_circle*2.*M_PI;

  if (test) {
    if(cyl_pm->initconfig == 1) 
      printf("Adding 1 cylinders of size %dx%d, radius =%lf, %d beads, angle =%le\n", n_circle, n_height[0], radius, N_per_cylinder[0], angle);
    else if(cyl_pm->initconfig == 2) 
      printf("Adding 1 slit of size %dx%dx%d, center length =%lf, height=%d, center width=%lf, %d beads, contraction slope =%le\n", max_x,max_y,max_z,cyl_pm->rperiod, max_z, radius, cyl_pm->num_beads, cyl_pm->slope);
    else if(cyl_pm->initconfig == 3) 
      printf("Adding 1 curved channel of size %dx%dx%d\n", max_x,max_y,max_z);
    else if(cyl_pm->initconfig == 6) 
      printf("Adding 1 trapozoidal channel of size %dx%dx%d\n", max_x,max_y,max_z);

    printf("spring=%d n0=%d num_beads=%d init=%d\n", cyl_pm->spring, n0, num_beads, cyl_pm->initconfig);
  }

  /* Initialize monomer properties */
  for(i=n0; i<num_beads; i++) {
    if(cyl_pm->ev_type == 0)   /* HS */
      monomers[i].radius = 0.5;
    else if(cyl_pm->ev_type == 1)   /* WCA */
      monomers[i].radius = 0.5;
    else {
      printf("wrong EV type!!\n");
      exit(1);
    }

    monomers[i].rho = 1.0;
    monomers[i].dr2 = 0.0;
    monomers[i].blist[0][0] = 4;

    for(j=0; j<DIMS; j++) {
      monomers[i].vel[j]=0.0;
      monomers[i].force[j]=0.0;
      monomers[i].force0[j]=0.0;
      monomers[i].fluid_vel[j]=0.0;
    }
  }

  /* determine bond length (stretched length) */
  if(cyl_pm->spring == 0)  /* FENE */
    // bondlength = cylinder_pm->Q_fene*(monomers[0].radius*2.0)*MAX_EXT;
    bondlength = 0.98; /* no stretch */
  else if(cyl_pm->spring == 2)  /* harmonic */
    bondlength = 1.0;
  else {
    printf("wrong spring type!!\n");
    exit(1);
  }

  /* generate a new config */
  /* start to grow cylinders */
  /* cylinders with sinusoidal variation in diameter */
  if (cyl_pm->initconfig == 1) {
    /* all cylinders are perfect cylinders. the only randomness is in the position and direction of the cylinders */
    for(n=0; n<Ncylinder; n++) {
      if (test)
        printf("Adding cylinder %d\n", n);

      /* define the position of monomers on the rings */
      /* now put in beads of the cylinder */
      for (k=0; k<max_x; k++) {

	if((double)k <= rperiod/2.0)
	  radius = cyl_pm->radius + ramp*cos((double)k/(rperiod)*2.0*M_PI);
	else if((double)k > rperiod/2.0 && (double)k < maxsize[0]-rperiod/2.0)
	  radius = cyl_pm->radius-ramp;
	else
	  radius = cyl_pm->radius + ramp*cos(((double)k-maxsize[0]-rperiod)/(rperiod)*2.0*M_PI);

	for (i=0; i<n_circle; i++)
	  {
	    m=n0+k*n_circle+i;
	    monomers[m].pos_pbc[1] = radius*cos((double)i*angle) + (maxsize[1]-1.0)/2.0;
	    monomers[m].pos_pbc[2] = radius*sin((double)i*angle) + (maxsize[2]-1.0)/2.0;
	    monomers[m].pos_pbc[0] = (double)k;
	  }
	
	for(m=n0; m<num_beads; m++) {
	  for (j=0 ; j<DIMS ; j++) {
	    monomers[m].pos[j] = box(monomers[m].pos_pbc[j],maxsize[j]);
	    monomers[m].pos0[j] = monomers[m].pos_pbc[j];
	    monomers[m].pos_lst[j] = monomers[m].pos[j];
	  }
 	  monomers[m].DP_id = nDP+n ;
	}
      }
    }
    num_beads = m-n0;
  }
  /* slit channels with contracting and expansion segs */
  if (cyl_pm->initconfig == 2) {
    /* define the position of monomers on the rings */
    /* now put in beads of the cylinder */
    /* Match the slit design of Quinn et al., Ann Biomed Eng 2011 */
    /* 2.7  micron high slit => 60 x 60 x 5 Dx^3 */
    /* X:  0 to 15, line with slope Dy/Dx = sqrt(3) for Y=[6 .. 32] and [37 .. 54] */
    /* Z:  0 to 5 */
    /* n_height = # of rings/rectangles along the x-direction */
    /* dx = distance between rings */
    /* contracting seg */

    slope = cyl_pm->slope;
    n_height[0] = cyl_pm->height;
    dx = maxsize[0]/(double)n_height[0];

    xpos = 0.0;
    xm = (max_x-cyl_pm->rperiod)/2;
    xp = (max_x-cyl_pm->rperiod)/2+cyl_pm->rperiod;

    ym = (maxsize[1]-1)/2.0-radius/2-slope*xm;
    yp = (maxsize[1]-1)/2.0+radius/2+slope*xm;

    printf("%d %d %lf %lf %lf %lf %lf %lf\n", n_height[0], cyl_pm->N_per_ring, dx, slope, xm, xp, ym, yp);

    m = n0;
    for (k=0; k<n_height[0]; k++) {
      /* define the y-bounds of the slit */
      if(xpos >=0.0 && xpos < xm) {
	ymin = (ym+slope*xpos);
	ymax = (yp-slope*xpos);
      }
      else if(xpos >=xm && xpos<=xp) {
	ymin = ym+slope*xm;
	ymax = yp-slope*xm;
      }
      else if(xpos > xp && xpos <=(maxsize[0]+1e-6)) {
	ymin = (ym-slope*(xpos-xp-xm));
	ymax = (yp+slope*(xpos-xp-xm));
      }
      
      /* length of y walls */
      //n_circle = (int)(ymax-ymin)+1;
      n_circle = max_y;
      dy = (ymax-ymin)/(double)n_circle;
      for(j=0; j<=n_circle; j++) {
	monomers[m].pos_pbc[1] = ymin+dy*j;
	monomers[m].pos_pbc[2] = 1.5;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }
      //      printf("%d %lf %lf %lf\n", m-1, monomers[m-1].pos_pbc[0], monomers[m-1].pos_pbc[1], monomers[m-1].pos_pbc[2]);
      for(j=1; j<=(max_z-3); j++) {
	monomers[m].pos_pbc[1] = ymax;
	monomers[m].pos_pbc[2] = 1.5+j;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }
      //      printf("%d %lf %lf %lf\n", m-1, monomers[m-1].pos_pbc[0], monomers[m-1].pos_pbc[1], monomers[m-1].pos_pbc[2]);
      for(j=1; j<=n_circle; j++) {
	monomers[m].pos_pbc[1] = ymax-dy*j;
	monomers[m].pos_pbc[2] = max_z-1.5;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }
      //      printf("%d %lf %lf %lf\n", m-1, monomers[m-1].pos_pbc[0], monomers[m-1].pos_pbc[1], monomers[m-1].pos_pbc[2]);
      for(j=1; j<(max_z-3); j++) {
	monomers[m].pos_pbc[1] = ymin;
	monomers[m].pos_pbc[2] = (max_z-1.5)-j;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }

      xpos += dx;
      //    printf("%d %lf %lf %lf\n", m-1, monomers[m-1].pos_pbc[0], monomers[m-1].pos_pbc[1], monomers[m-1].pos_pbc[2]);
    }

    num_beads = m-n0;
    cyl_pm->num_beads = num_beads;

    //  printf("m=%d, %d beads\n", m, num_beads);
    for(i=n0; i<n0+num_beads; i++) {
      for (j=0 ; j<DIMS ; j++) {
	monomers[i].pos[j] = box(monomers[i].pos_pbc[j],maxsize[j]);
	monomers[i].pos0[j] = monomers[i].pos_pbc[j];
	monomers[i].pos_lst[j] = monomers[i].pos[j];
      }
      monomers[i].DP_id = nDP ;
    }
  } 
  /* curved slit channels, cyl_pm->rperiod = radius of curvature; cyl_pm->radius = slit width */
  if (cyl_pm->initconfig == 3) {
    /* define the position of monomers on the rings */
    /* now put in beads of the cylinder */
    
    n_height[0] = cyl_pm->height;
    dx = maxsize[0]/(double)n_height[0];

    xpos = 0.0;
    ym = sqrt(cyl_pm->rperiod*cyl_pm->rperiod-(xpos-(maxsize[0]+1)/2.0)*(xpos-(maxsize[0]+1)/2.0));
    m = n0;
    for (k=0; k<n_height[0]; k++) {
      /* define the y-bounds of the slit */
      ymin = sqrt(cyl_pm->rperiod*cyl_pm->rperiod-(xpos-(maxsize[0]+1)/2.0)*(xpos-(maxsize[0]+1)/2.0))-ym + 2.0;
      ymax = ymin + cyl_pm->radius;
 
      /* length of y walls */
      n_circle = max_y;
      dy = (ymax-ymin)/(double)n_circle;
      for(j=0; j<=n_circle; j++) {
	monomers[m].pos_pbc[1] = ymin+dy*j;
	monomers[m].pos_pbc[2] = 1.5;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }
      //      printf("%d %lf %lf %lf\n", m-1, monomers[m-1].pos_pbc[0], monomers[m-1].pos_pbc[1], monomers[m-1].pos_pbc[2]);
      for(j=1; j<=(max_z-3); j++) {
	monomers[m].pos_pbc[1] = ymax;
	monomers[m].pos_pbc[2] = 1.5+j;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }
      //      printf("%d %lf %lf %lf\n", m-1, monomers[m-1].pos_pbc[0], monomers[m-1].pos_pbc[1], monomers[m-1].pos_pbc[2]);
      for(j=1; j<=n_circle; j++) {
	monomers[m].pos_pbc[1] = ymax-dy*j;
	monomers[m].pos_pbc[2] = max_z-1.5;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }
      //      printf("%d %lf %lf %lf\n", m-1, monomers[m-1].pos_pbc[0], monomers[m-1].pos_pbc[1], monomers[m-1].pos_pbc[2]);
      for(j=1; j<(max_z-3); j++) {
	monomers[m].pos_pbc[1] = ymin;
	monomers[m].pos_pbc[2] = (max_z-1.5)-j;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }

      xpos += dx;
      //    printf("%d %lf %lf %lf\n", m-1, monomers[m-1].pos_pbc[0], monomers[m-1].pos_pbc[1], monomers[m-1].pos_pbc[2]);
    }

    num_beads = m-n0;
    cyl_pm->num_beads = num_beads;

    //  printf("m=%d, %d beads\n", m, num_beads);
    for(i=n0; i<n0+num_beads; i++) {
      for (j=0 ; j<DIMS ; j++) {
	monomers[i].pos[j] = box(monomers[i].pos_pbc[j],maxsize[j]);
	monomers[i].pos0[j] = monomers[i].pos_pbc[j];
	monomers[i].pos_lst[j] = monomers[i].pos[j];
      }
      monomers[i].DP_id = nDP ;
    }
  } 
  /* cosine slit channels, cyl_pm->rperiod = radius of curvature; cyl_pm->radius = slit width */
  if (cyl_pm->initconfig == 5) {
    /* define the position of monomers on the rings */
    /* now put in beads of the cylinder */
    offset = cyl_pm->offset;
    n_height[0] = cyl_pm->height;
    dx = maxsize[0]/(double)n_height[0];

    xpos = 0.0;

    m = n0;
    for (k=0; k<n_height[0]; k++) {
      /* define the y-bounds of the slit */
      if((double)k <= offset) {
	ymax = cyl_pm->radius + ramp + (maxsize[1])/2.0;
	ymin = -(cyl_pm->radius + ramp) + (maxsize[1])/2.0;
      }
      else if((double)k <= (rperiod/2.0 + offset)) {
	ymax = cyl_pm->radius + ramp*cos((double)(k-offset)/(rperiod)*2.0*M_PI) + (maxsize[1])/2.0;
	ymin = -(cyl_pm->radius + ramp*cos((double)(k-offset)/(rperiod)*2.0*M_PI)) + (maxsize[1])/2.0;
      }
      else if((double)k > (rperiod/2.0+offset) && (double)k < (maxsize[0]-rperiod/2.0)) {
	ymax = cyl_pm->radius-ramp + (maxsize[1])/2.0;
	ymin = -(cyl_pm->radius-ramp) + (maxsize[1])/2.0;
      }
      else {
	ymax = cyl_pm->radius + ramp*cos(((double)(k)-maxsize[0]-rperiod)/(rperiod)*2.0*M_PI) + (maxsize[1])/2.0;;
	ymin = -(cyl_pm->radius + ramp*cos(((double)(k)-maxsize[0]-rperiod)/(rperiod)*2.0*M_PI)) + (maxsize[1])/2.0;;
      }
      
      /* length of y walls */
      n_circle = max_y;
      dy = (ymax-ymin)/(double)n_circle;
      for(j=0; j<=n_circle; j++) {
	monomers[m].pos_pbc[1] = ymin+dy*j;
	monomers[m].pos_pbc[2] = 1.5;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }
      //      printf("%d %lf %lf %lf\n", m-1, monomers[m-1].pos_pbc[0], monomers[m-1].pos_pbc[1], monomers[m-1].pos_pbc[2]);
      for(j=1; j<=(max_z-3); j++) {
	monomers[m].pos_pbc[1] = ymax;
	monomers[m].pos_pbc[2] = 1.5+j;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }
      //      printf("%d %lf %lf %lf\n", m-1, monomers[m-1].pos_pbc[0], monomers[m-1].pos_pbc[1], monomers[m-1].pos_pbc[2]);
      for(j=1; j<=n_circle; j++) {
	monomers[m].pos_pbc[1] = ymax-dy*j;
	monomers[m].pos_pbc[2] = max_z-1.5;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }
      //      printf("%d %lf %lf %lf\n", m-1, monomers[m-1].pos_pbc[0], monomers[m-1].pos_pbc[1], monomers[m-1].pos_pbc[2]);
      for(j=1; j<(max_z-3); j++) {
	monomers[m].pos_pbc[1] = ymin;
	monomers[m].pos_pbc[2] = (max_z-1.5)-j;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }

      xpos += dx;
      //    printf("%d %lf %lf %lf\n", m-1, monomers[m-1].pos_pbc[0], monomers[m-1].pos_pbc[1], monomers[m-1].pos_pbc[2]);
    }

    num_beads = m-n0;
    cyl_pm->num_beads = num_beads;

    //  printf("m=%d, %d beads\n", m, num_beads);
    for(i=n0; i<n0+num_beads; i++) {
      for (j=0 ; j<DIMS ; j++) {
	monomers[i].pos[j] = box(monomers[i].pos_pbc[j],maxsize[j]);
	monomers[i].pos0[j] = monomers[i].pos_pbc[j];
	monomers[i].pos_lst[j] = monomers[i].pos[j];
      }
      monomers[i].DP_id = nDP ;
    }
  }
  /* trapezoidal channel  */
  if (cyl_pm->initconfig == 6) {
    /* Channel from Chun Yang's group */
    /* bottom and side walls are flat */
    /* top wall is sloped with offset y0 and slope rad_amp  */
    double zm, zp;
    offset = cyl_pm->offset;
    slope = cyl_pm->ramp;
    n_height[0] = cyl_pm->height;
    dx = maxsize[0]/(double)n_height[0];
    xpos = 0.0;

    zm = 0; 
    zp = maxsize[2]-1.0;

    ym = offset;
    yp = offset+slope*zp;

    m = n0;
    for (k=0; k<n_height[0]; k++) {
      /* define the y-bounds of the slit */
      /*Mark walls on all sides */
      for(j=0; j<=ym; j++) {
	monomers[m].pos_pbc[1] = j;
	monomers[m].pos_pbc[2] = zm;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }
      for(j=zm+1; j<=zp; j++) {
	monomers[m].pos_pbc[1] = offset+slope*j;
	monomers[m].pos_pbc[2] = j;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }
      for(j=1; j<=yp; j++) {
	monomers[m].pos_pbc[1] = offset+slope*zp-j;
	monomers[m].pos_pbc[2] = zp;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }
      if(monomers[m-1].pos_pbc[1] == 0.0) zp=zp-1;
      for(j=zp; j>zm; j--) {
	monomers[m].pos_pbc[1] = 0;
	monomers[m].pos_pbc[2] = j;
	monomers[m].pos_pbc[0] = xpos;
	m++;
      }

      xpos += dx;
      //    printf("%d %lf %lf %lf\n", m-1, monomers[m-1].pos_pbc[0], monomers[m-1].pos_pbc[1], monomers[m-1].pos_pbc[2]);
    }
    
    num_beads = m-n0;
    cyl_pm->num_beads = num_beads;
    cyl_pm->N_per_ring = num_beads / n_height[0];

    printf("m=%d, %d beads %d beads per ring\n", m, num_beads, cyl_pm->N_per_ring);
    for(i=n0; i<n0+num_beads; i++) {
      for (j=0 ; j<DIMS ; j++) {
	monomers[i].pos[j] = box(monomers[i].pos_pbc[j],maxsize[j]);
	monomers[i].pos0[j] = monomers[i].pos_pbc[j];
	monomers[i].pos_lst[j] = monomers[i].pos[j];
      }
      monomers[i].DP_id = nDP ;
    }
  }  
  /* read in the initial cylinder configuration from file */
  else if(cyl_pm->initconfig == 4) {
    sprintf(filename, "%s/init/init.config", work_dir);
    printf("init config file %s\n", filename);
    stream = fopen(filename, "r");      
    fscanf(stream, "%d %d %d %d %d", &cyl_pm->Ntype[0], &n_circle, &n_height[0], &cyl_pm->num_beads, &n0);

    Ncylinder = 1;
    num_beads = n0+cyl_pm->num_beads;

    Ntype[0] = cyl_pm->Ntype[0];
    N_per_cylinder[0] = n_circle*n_height[0];

    fscanf(stream, "%le %le", &cylinders[0].disp2, &cylinders[0].dr2);

    for(n=n0; n<num_beads; n++) {
      fscanf(stream, "%le %le %le %le %le %le %le %le %le %le", &monomers[n].pos[0], &monomers[n].pos[1], &monomers[n].pos[2], &monomers[n].pos_pbc[0], &monomers[n].pos_pbc[1], &monomers[n].pos_pbc[2], &monomers[n].vel[0], &monomers[n].vel[1], &monomers[n].vel[2], &monomers[n].dr2);

      if(n_proc == 0)
	printf("#monomer %d at (%le %le %le) (%le %le %le)\n", n, monomers[n].pos[0], monomers[n].pos[1], monomers[n].pos[2], monomers[n].pos_pbc[0], monomers[n].pos_pbc[1], monomers[n].pos_pbc[2]);
    }
       
    fclose(stream);
  }

  /* initialize monomer properties */

  for(n=n0; n<num_beads+n0; n++) {
    if(n_proc == 0)
      printf("#monomer %d at (%le %le %le) (%le %le %le)\n", n, monomers[n].pos[0], monomers[n].pos[1], monomers[n].pos[2], monomers[n].pos_pbc[0], monomers[n].pos_pbc[1], monomers[n].pos_pbc[2]);

    monomers[n].rho=1.0;
    monomers[n].dr2=0.0;

    for(j=0; j<DIMS; j++) {
      monomers[n].vel[j] = 0.0;
      monomers[n].fluid_vel[j] = 0.0;
      monomers[n].force[j]=0.0;
      monomers[n].force0[j]=0.0;
    }
  }

  /* write config to a file */
  if (test) {
    file_name (test_file, work_dir, -1);
    file_ptr = fopen (test_file, "w");
    if (file_ptr == 0)  fatal_err("Could not open test.dat", -1);
    for (m=0 ; m<num_beads ; m++)
      fprintf(file_ptr, "%le %le %le\n", monomers[m].pos_pbc[0], monomers[m].pos_pbc[1], monomers[m].pos_pbc[2]);
    fclose(file_ptr);
  }

  /* initialize cylinder properties */
  for(j=0; j<DIMS; j++) {
    cylinders[0].com[j]=0.0;
    cylinders[0].com0[j] = cylinders[0].com[j];
    cylinders[0].comold[j] = cylinders[0].com[j];
  }

  /* calculate cylinder area and volume */
  //  cylinder_props(cyl_pm, cylinders, monomers);

  cylinders[0].disp2 = 0.0;
  cylinders[0].dr2 = 0.0;
  for(j=0; j<DIMS; j++){
    cylinders[0].com0[j] = cylinders[0].com[j];
    cylinders[0].comold[j] = cylinders[0].com[j];
  }

  return 1;

}


