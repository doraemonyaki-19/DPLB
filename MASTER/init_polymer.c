/* Following Ahlrichs and Dunweg, modeling polymer segments as point forces in the fluid */
/* Coupled with LBE code to add hydrodynamic forces */
/***********************************************************************
ASLB : hybrid lattice Boltzmann - Brownian dynamics algorithm
Copyright (C) 2009 Yeng-Long Chen

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
#include <time.h>
#define DEBUG 0

int polymer_init(struct DP_param *chain_pm, struct DP_param *cyl_pm, struct DP *chains, struct monomer *monomers, char *work_dir, int n0, int id0)
{
  extern int n_proc;
  extern double tau;
  extern int max_x, max_y, max_z;
  extern int add_noise;
  extern int wall_flag;
  extern unsigned long seed;

  signed long iseed=seed;
  int i,j,i1,j1, type;
  int n, n_mon;
  int n_overlap=0, wall_overlap=0;
  int pnum, pnum1;
  int Ntype[NTYPES], Nbeads[NTYPES];
  int num_beads=0, num_chains=0;
  int flag;

  double maxsize[DIMS];
  double bondlength, bondmag;
  double bond[DIMS];
  double tuberadius, dr, dy, dz;
  char filename[160];
  FILE *stream;

  srand((unsigned)time(0));
  iseed = rand();

  // need to choose type of chain //
  maxsize[0]=(double)max_x;
  maxsize[1]=(double)max_y;
  maxsize[2]=(double)max_z;

  for (type=0 ; type<NTYPES ; type++) {
    Ntype[type] = chain_pm->Ntype[type];
    Nbeads[type] = chain_pm->N_per_DP[type];
    num_beads += Ntype[type]*Nbeads[type];
    num_chains += Ntype[type];
  }

  printf("Adding %d chain0, %d beads/chain, %d chain1, %d beads/chain\n", Ntype[0], Nbeads[0], Ntype[1], Nbeads[1]);

  /* Initialize monomer velocity and forces */
  for(i=n0; i<n0+num_beads; i++) {
    if(chain_pm->ev_type == 1)
      monomers[i].radius = 0.5;
    else if(chain_pm->ev_type == 2)
      monomers[i].radius = chain_pm->Ss;

    monomers[i].rho = 1.0;
    monomers[i].dr2 = 0.0;

    for(j=0; j<DIMS; j++) {
      monomers[i].pos[0]=0.0;
      monomers[i].vel[j]=0.0;
      monomers[i].force[j]=0.0;
      monomers[i].force0[j]=0.0;
      monomers[i].fluid_vel[j]=0.0;
    }
  }

  monomers[0].vel[0]=0.0;

  /* grow a chain from chain end */

  if(chain_pm->initconfig != 3) {
    pnum = n0;
    for(i=0; i<num_chains; i++) {
      /* generate coord. of all chain ends */
      monomers[pnum].DP_id=i+id0;

      for(j=0; j<DIMS; j++) {
	//	monomers[pnum].pos[j]=rand_num(&iseed, monomers[0].radius, (maxsize[j]-1.0-monomers[0].radius));
	monomers[pnum].pos[j]=(double)rand()/((double)(RAND_MAX)+1.0)*(maxsize[j]-1.0-1.0)+1.0;
        monomers[pnum].pos_pbc[j]=monomers[pnum].pos[j];
      }
      
      if(n_proc == 0) {
	printf("#chain %d \n", i);
	printf("#monomer %d at (%le %le %le)\n", pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], monomers[pnum].pos[2]);
      }

      /* check whether bead is inside tube */
      if(cyl_pm->NDP > 0) {
	if(cyl_pm->initconfig == 1) {
	  if(monomers[pnum].pos[0] <= cyl_pm->rperiod/2.0)
	    tuberadius = cyl_pm->radius + cyl_pm->ramp*cos(monomers[pnum].pos[0]/(cyl_pm->rperiod)*2.0*M_PI);
	  else if(monomers[pnum].pos[0] > cyl_pm->rperiod/2.0 && monomers[pnum].pos[0] < maxsize[0]-cyl_pm->rperiod/2.0)
	    tuberadius = cyl_pm->radius-cyl_pm->ramp;
	  else
	    tuberadius = cyl_pm->radius + cyl_pm->ramp*cos((monomers[pnum].pos[0]-maxsize[0]-cyl_pm->rperiod)/(cyl_pm->rperiod)*2.0*M_PI);
	  
	  dy = (monomers[pnum].pos[1] - (maxsize[1]-1)/2.0)*(monomers[pnum].pos[1] - (maxsize[1]-1)/2.0);
	  dz = (monomers[pnum].pos[2] - (maxsize[2]-1)/2.0)*(monomers[pnum].pos[2] - (maxsize[2]-1)/2.0);

	  dy =sqrt(dy);
	  dz = sqrt(dz);
	  dy = dy+1.5;
	  dz = dz+1.5;
	  dr = sqrt(dy*dy+dz*dz);

	  if(dr > tuberadius) {
	    printf("%d Tube radius %le Overlap dr=%le found at monomer %d (%le %le %le)\n", 
		   i, tuberadius, sqrt(dr), pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], monomers[pnum].pos[2]);
	    i--;
	    continue;
	  } 
	}  
	
	if(cyl_pm->initconfig == 6) 
	  if ((monomers[pnum].pos[1] + monomers[pnum].radius*1.0 + 1.05) > (cyl_pm->offset+cyl_pm->ramp*monomers[pnum].pos[2])) {
	    i--;
	    continue;
	  }
      }

      /* check for overlap between monomers for non-overlapping potential */
      if(chain_pm->ev_type < 5) {
	for(pnum1=0; pnum1 < pnum; pnum1++)
	  if(check_monoverlap(monomers[pnum], monomers[pnum1])) {
	    if(n_proc == 0)
	      printf("Chain end overlap found at monomer %d (%le %le %le) and monomer %d (%le %le %le)\n", pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], monomers[pnum].pos[2], pnum1, monomers[pnum1].pos[0], monomers[pnum1].pos[1], monomers[pnum1].pos[2]);
	    i--;
	    break;
	  }
      }

      if(i < Ntype[0])
	n_mon = Nbeads[0];
      else
	n_mon = Nbeads[1];

      pnum += n_mon;
    }
  }

  /* random chain configuration */
  if(chain_pm->initconfig == 1) {
    /* add beads to the chain */
    for(type=0; type < NTYPES; type++) {
      n_mon = Nbeads[type];
      for(i=0; i<Ntype[type]; i++) {
	for(n=1; n<n_mon; n++) {
	  if(type == 0) {
	    pnum = n0+i*Nbeads[0]+n;
	    monomers[pnum].DP_id=id0+i;
	  }
	  else if(type == 1) {
	    pnum = n0+Ntype[0]*Nbeads[0]+i*Nbeads[1]+n;
	    monomers[pnum].DP_id=id0+i+Ntype[0];
	  }

	  monomers[pnum].blist[0][0] = 1;
	  monomers[pnum].blist[1][0] = pnum-1;
	  
	  if(DEBUG == 1)
	    printf("Type %d chain %d monomer %d pnum %d\n", type, i, n, pnum);
	  
	  /* set bond length */
	  if(chain_pm->spring == 0)
	    bondlength = chain_pm->Q_fene[type] * (monomers[0].radius*2.0) * rand_num(&iseed, 0.2, 0.9);
	  else if(chain_pm->spring == 1)
	    //	    bondlength = 1.0+rand_num(&iseed, 0.0,1.0)/10.0;
	    bondlength = chain_pm->nks*chain_pm->sigma_k * rand_num(&iseed, 0.2, 0.5);
	  else if(chain_pm->spring == 2)
	    bondlength = 1.0+rand_num(&iseed, 0.0,1.0)/10.0;
	  
	  bondmag = 0.0;
	  for(j=0; j<DIMS; j++) {
	    bond[j] = rand_num(&iseed, -1.0, 1.0);
	    bondmag += bond[j]*bond[j];
	  }
	  
	  bondmag = sqrt(bondmag);
	  
	  for(j=0; j<DIMS; j++) {
	    bond[j] *= bondlength/bondmag;
	    monomers[pnum].pos_pbc[j]=monomers[pnum-1].pos_pbc[j] + bond[j];
	    monomers[pnum].pos[j]=box(monomers[pnum].pos_pbc[j], maxsize[j]);
	  }
	  
	  /* Check Wall and Tube overlap */
	  if(check_walloverlap(monomers[pnum]) == 1) {
	    if(n_proc == 0)
	      printf("%d Wall Overlap found at monomer %d (%le %le %le)\n", 
		     n, pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], monomers[pnum].pos[2]);
	    n_overlap++;
	    n--;
	    continue;
	  }

	  if (cyl_pm->NDP != 0) {
	    if(cyl_pm->initconfig == 1) { 
	      if(monomers[pnum].pos[0] <= cyl_pm->rperiod/2.0)
		tuberadius = cyl_pm->radius + cyl_pm->ramp*cos(monomers[pnum].pos[0]/(cyl_pm->rperiod)*2.0*M_PI);
	      else if(monomers[pnum].pos[0] > cyl_pm->rperiod/2.0 && monomers[pnum].pos[0] < maxsize[0]-cyl_pm->rperiod/2.0)
		tuberadius = cyl_pm->radius-cyl_pm->ramp;
	      else
		tuberadius = cyl_pm->radius + cyl_pm->ramp*cos((monomers[pnum].pos[0]-maxsize[0]-cyl_pm->rperiod)/(cyl_pm->rperiod)*2.0*M_PI);
	      
	      dy = (monomers[pnum].pos[1] - (maxsize[1]-1.0)/2.0)*(monomers[pnum].pos[1] - (maxsize[1]-1.0)/2.0);
	      dz = (monomers[pnum].pos[2] - (maxsize[2]-1.0)/2.0)*(monomers[pnum].pos[2] - (maxsize[2]-1.0)/2.0);
	      
	      dy =sqrt(dy);
	      dz = sqrt(dz);
	      dy = dy+1.5;
	      dz = dz+1.5;
	      dr = sqrt(dy*dy+dz*dz);
	      
	      if(dr > tuberadius) {
		printf("%d Tube radius %le Overlap dr %le found at monomer %d (%le %le %le)\n", 
		       n, tuberadius, sqrt(dr), pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], monomers[pnum].pos[2]);
		n_overlap++;
		n--;
		continue;
	      } 
	    }
	  	     
	    if(cyl_pm->initconfig == 6) 
	      if ((monomers[pnum].pos[1] + monomers[pnum].radius*1.0 + 1.05) > (cyl_pm->offset+cyl_pm->ramp*monomers[pnum].pos[2])) {
		if(n_proc == 0)
		  printf("%d Tube Overlap found at monomer %d (%le %le %le)\n", 
			 n, pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], monomers[pnum].pos[2]);
		n_overlap++;
		n--;
		continue;
	      }
	  }
	  
	  /* check for overlap  between monomers for non-overlap potential*/
	  if(chain_pm->ev_type < 5) {
	    for(pnum1=0; pnum1 < pnum; pnum1++) {
	      if(check_monoverlap(monomers[pnum], monomers[pnum1])==1) {
		if(n_proc == 0)
		  printf("%d Overlap found at monomer %d (%le %le %le) and monomer %d (%le %le %le)\n", 
			 n, pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], monomers[pnum].pos[2], 
			 pnum1, monomers[pnum1].pos[0], monomers[pnum1].pos[1], monomers[pnum1].pos[2]);
		n_overlap++;
		n--;
		break;
	      }
	    }
	  }
	  
	  if(n_overlap > 1000) {
	    n = 1;
	    n_overlap = 0;
	    continue;
	  }
	  
	  if(n_proc == 0)
	    printf("#monomer %d at (%le %le %le)\n", pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], monomers[pnum].pos[2]);
	}
      }
    }
    return 1;
  }
     
  /* chain init config is stretched */
  /* add beads to the chain */
  if(chain_pm->initconfig == 2) {
    for(type=0; type < NTYPES; type++) {
      n_mon = Nbeads[type];

      for(i=0; i< Ntype[type]; i++) {
	for(n=1; n<n_mon; n++) {
	  if(type == 0) {
	    pnum = n0+i*Nbeads[0]+n;
	    monomers[pnum].DP_id=id0+i;
	  }
	  else if(type == 1) {
	    pnum = n0+Ntype[0]*Nbeads[0]+i*Nbeads[1]+n;
	    monomers[pnum].DP_id=id0+i+Ntype[0];
	  }

	  monomers[pnum].blist[0][0] = 1;
	  monomers[pnum].blist[1][0] = pnum-1;
	  
	  /* set bond length */
	  if(chain_pm->spring == 0)
	    bondlength = chain_pm->Q_fene[type]*(monomers[0].radius*2.0)*MAX_EXT;
	  else if(chain_pm->spring == 1)
	    bondlength = chain_pm->nks*chain_pm->sigma_k*MAX_EXT;
	  else if(chain_pm->spring == 2)
	    bondlength = 1.1;
	  
	  bond[0]=bondlength;
	  bond[1]=0.0;
	  bond[2]=0.0;
	  
	  for(j=0; j<DIMS; j++) {
	    monomers[pnum].pos_pbc[j]=monomers[pnum-1].pos_pbc[j] + bond[j];
	    monomers[pnum].pos[j]=box(monomers[pnum].pos_pbc[j], maxsize[j]);
	  }
	  
	  /* check for tube overlap */
	  if (cyl_pm->NDP != 0) {
	    if(cyl_pm->initconfig == 1) {
	      if(monomers[pnum].pos[0] <= cyl_pm->rperiod/2.0)
		tuberadius = cyl_pm->radius + cyl_pm->ramp*cos(monomers[pnum].pos[0]/(cyl_pm->rperiod)*2.0*M_PI);
	      else if(monomers[pnum].pos[0] > cyl_pm->rperiod/2.0 && monomers[pnum].pos[0] < maxsize[0]-cyl_pm->rperiod/2.0)
		tuberadius = cyl_pm->radius-cyl_pm->ramp;
	      else
		tuberadius = cyl_pm->radius + cyl_pm->ramp*cos((monomers[pnum].pos[0]-maxsize[0]-cyl_pm->rperiod)/(cyl_pm->rperiod)*2.0*M_PI);
	      
	      dr = (monomers[pnum].pos[1] - (maxsize[1]-1)/2.0)*(monomers[pnum].pos[1] - (maxsize[1]-1)/2.0)+(monomers[pnum].pos[2] - (maxsize[2]-1)/2.0)*(monomers[pnum].pos[2] - (maxsize[2]-1)/2.0);
	      
	      if(dr > tuberadius*tuberadius) {
		n_overlap++;
		/* regenerate the first bead */
		pnum = pnum - n;
		monomers[pnum].pos[1]=(double)rand()/((double)(RAND_MAX)+1.0)*(maxsize[1]-1.0-1.0)+1.0;
		monomers[pnum].pos_pbc[1]=monomers[pnum].pos[1];
		n = 1;
		continue;
	      } 
	    }
	      
	    if(cyl_pm->initconfig == 6) 
	      if ((monomers[pnum].pos[1] + monomers[pnum].radius*1.0 + 1.05) > (cyl_pm->offset+cyl_pm->ramp*monomers[pnum].pos[2])) {
		n_overlap++;
		
		/* regenerate the first bead */
		pnum = pnum - n;
		monomers[pnum].pos[1]=(double)rand()/((double)(RAND_MAX)+1.0)*(maxsize[1]-1.0-1.0)+1.0;
		monomers[pnum].pos_pbc[1]=monomers[pnum].pos[1];
		n = 1;
		
		continue;
	      }
	  }
	
	  /* check for overlap  between monomers for non-overlap potential */
	  for(pnum1=0; pnum1 < pnum; pnum1++) {
	    if(check_monoverlap(monomers[pnum], monomers[pnum1])==1) {
	      n_overlap++;
	      /* regenerate the first bead */
	      pnum = pnum - n;
	      monomers[pnum].pos[1]=(double)rand()/((double)(RAND_MAX)+1.0)*(maxsize[1]-1.0-1.0)+1.0;
	      monomers[pnum].pos_pbc[1]=monomers[pnum].pos[1];
	      n = 1;
	      break;
	    }
	  }

	  if(n_proc == 0)
	    printf("#monomer %d at (%le %le %le)\n", pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], monomers[pnum].pos[2]);
	}
      }
    }
    return 2;
  }

  /* initial lattice chain configuration */
  if(chain_pm->initconfig == 3) {
    int px,py,pz;
    px= (int)maxsize[0];
    py= (int)maxsize[1];
    pz= (int)maxsize[2];
      
    for(type = 0; type < NTYPES; type ++) {
      for(i=0; i<Ntype[type]; i++) {
	/* printf("Adding polymer %d\n", i); */
	/* generate coord. of first chain end */
	if(type == 0)
	  monomers[pnum].DP_id=id0+i;
	else if(type == 1)
	  monomers[pnum].DP_id=id0+i+Ntype[0];

	if(i != 0) {
	  monomers[pnum].blist[0][0] = 1;
	  monomers[pnum].blist[1][0] = pnum-1;
	}

	monomers[pnum].pos_pbc[2]=maxsize[2]/2.0-0.2*(chain_pm->nks*chain_pm->sigma_k*(double)(Nbeads[type]-1.0));
	monomers[pnum].pos_pbc[1]=(double)((i)%(py-2))+1.0;
	monomers[pnum].pos_pbc[0]=(double)((i/(py-2))%(px-3))+2.0;
	  
	for(j=0; j<DIMS; j++)
	  monomers[pnum].pos[j] = monomers[pnum].pos_pbc[j];
	  
	if(n_proc == 0) {
	  printf("#chain %d \n", i);
	  printf("#monomer %d at (%le %le %le)\n", pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], monomers[pnum].pos[2]);
	}

	/* add beads to the chain */
	for(n=1; n<n_mon; n++) {
	  if(type == 0) {
	    pnum = n0+i*Nbeads[0]+n;
	    monomers[pnum].DP_id=id0+i;
	  }
	  else if(type == 1) {
	    pnum = n0+Ntype[0]*Nbeads[0]+i*Nbeads[1]+n;
	    monomers[pnum].DP_id=id0+i+Ntype[0];
	  }
	    
	  /* set bond length */
	  if(chain_pm->spring == 0)
	    bondlength = chain_pm->Q_fene[type] * (monomers[0].radius*2.0) * rand_num(&iseed, 0.2, 0.9);
	  else if(chain_pm->spring == 1)
	    bondlength = chain_pm->nks*chain_pm->sigma_k * rand_num(&iseed, 0.2, 0.7);
	  else if(chain_pm->spring == 2)
	    bondlength = 1.0+rand_num(&iseed, -1.0,1.0)/30.0;
	    
	  bondmag = 0.0;
	  bond[0] = rand_num(&iseed, 0.0, 1.0);
	  bond[1] = 0.0;
	  bond[2] = rand_num(&iseed, 0.0, 1.0);
	  for(j=0; j<DIMS; j++) 
	    bondmag += bond[j]*bond[j];
	  
	  bondmag = sqrt(bondmag);
	    
	  for(j=0; j<DIMS; j++) {
	    bond[j] *= bondlength/bondmag;
	    monomers[pnum].pos_pbc[j]=monomers[pnum-1].pos_pbc[j] + bond[j];
	    monomers[pnum].pos[j]=monomers[pnum].pos_pbc[j];
	  }
	  monomers[pnum].pos_pbc[1]=monomers[pnum-1].pos_pbc[1];
	  monomers[pnum].pos[1]=monomers[pnum].pos_pbc[1];
	    
	  if(check_walloverlap(monomers[pnum])==1) {
	    n--;
	    continue;
	  }
	  if(n_proc == 0)
	    printf("#monomer %d at (%le %le %le)\n", pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], monomers[pnum].pos[2]);
	}
      }
    }
    return 3;
  }

  /* read in the initial chain configuration from file */
  if(chain_pm->initconfig == 4) {
    sprintf(filename, "%s/init/initpoly.config", work_dir);
    printf("init config file %s\n", filename);
    stream = fopen(filename, "r");      
    for(type = 0; type < NTYPES; type++)
      fscanf(stream, "%d %d ", &Ntype[type], &Nbeads[type]);
      
    for(i=0; i<num_chains; i++) 
      fscanf(stream, "%le %le", &chains[i].disp2, &chains[i].dr2);

    for(type = 0; type < NTYPES; type++) {
      for(i=0; i<Ntype[type]; i++) {
	if(n_proc == 0)
	  printf("#chain %d \n", i);
	pnum = n0;

	if(i != 0) {
	  monomers[pnum].blist[0][0] = 1;
	  monomers[pnum].blist[1][0] = pnum-1;
	}	  

	for(n=0; n<Nbeads[type]; n++) {
	  if(type == 0)
	    monomers[pnum].DP_id=id0+i;
	  else if(type == 1)
	    monomers[pnum].DP_id=id0+i+Ntype[0];

	  fscanf(stream, "%le %le %le %le %le %le %le %le %le %le", &monomers[pnum].pos[0], &monomers[pnum].pos[1], &monomers[pnum].pos[2], &monomers[pnum].pos_pbc[0], &monomers[pnum].pos_pbc[1], &monomers[pnum].pos_pbc[2], &monomers[pnum].vel[0], &monomers[pnum].vel[1], &monomers[pnum].vel[2], &monomers[pnum].dr2);
	  if(n_proc == 0)
	    printf("#monomer %d at (%le %le %le) (%le %le %le)\n", pnum, monomers[pnum].pos[0], monomers[pnum].pos[1], monomers[pnum].pos[2], monomers[pnum].pos_pbc[0], monomers[pnum].pos_pbc[1], monomers[pnum].pos_pbc[2]);
	  pnum++;
	}
      }
    }
    fclose(stream);
    return 4;
  }

  return 0;
}

int check_monoverlap(struct monomer monomer1, struct monomer monomer2)
{
  int j;
  double maxsize[DIMS];
  double r2;
  double radius2 = (monomer1.radius+monomer2.radius)*(monomer1.radius+monomer2.radius);

  r2=0.0;
  for(j=0; j<DIMS; j++) 
    r2+=(monomer1.pos_pbc[j]-monomer2.pos_pbc[j])*(monomer1.pos_pbc[j]-monomer2.pos_pbc[j]);

  if(r2 < radius2*1.5) 
    return 1;
  else 
    return 0;
}
 
int check_walloverlap(struct monomer monomer1) 
{
  extern int max_x, max_y, max_z;
  extern int wall_flag;
  int j;
  double maxsize[DIMS];
  double r2;

  maxsize[0]=(double)max_x;
  maxsize[1]=(double)max_y;
  maxsize[2]=(double)max_z;

  /* between monomer and wall */
  if(wall_flag == 0)
    return 0;
  else if(wall_flag >= 1)
    if(monomer1.pos_pbc[1] > maxsize[1]-1.0 || monomer1.pos_pbc[1] < 0.0)
      return 1;
    else if(wall_flag >= 2)
      if(monomer1.pos_pbc[2] > maxsize[2]-1.0 || monomer1.pos_pbc[2] < 0.0)
	return 1;
      else if(wall_flag >= 3)
	if(monomer1.pos_pbc[0] > maxsize[0]-1.0 || monomer1.pos_pbc[0] < 0.0)
	  return 1;
	else
	  return 0;
  
  return 0;
}
    
