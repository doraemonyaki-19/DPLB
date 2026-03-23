/*  define inside nodes for deformable particles */
/*	   BNODES_DP    Initialize inside nodes for deformable particles
	   BNODES_DPDEL Delelte nodes no longer in DP

*/
#define DP_DEBUG 0

#include "header.h"

void bnodes_dp(int **node_map, struct DP_param *sphere_pm, struct DP *sph, struct monomer *mon, int nstep)
{
  extern int    max_x, max_y, max_z;
  extern int    num_x;

  int leftbound, rightbound;
  int type, index;
  int Ntype[NTYPES];
  int Nbeads[NTYPES];
  int st_offset[DIMS];
  int maxsize[DIMS];
  int **nx, nxy, xy, xy1, z, z1, px[DIMS];
  int d, m, n, n1,n2, i, j, k, l;
  int **com_shift;
  double dr;
  int ***tmpnode, **tmp_pp, *tmp_p;  /* temporary nodemap */

  maxsize[0] = max_x;
  maxsize[1] = max_y;
  maxsize[2] = max_z;

  for(i=0; i<NTYPES; i++) {
    Ntype[i] = sphere_pm->Ntype[i];
    Nbeads[i] = sphere_pm->N_per_DP[i];
  }

  /* Update the DP com and stretch */
  sphere_props_short(sphere_pm, sph, mon);
  
  nx = (int **) calloc(sphere_pm->NDP, sizeof(*nx));
  if(nx == 0) fatal_err("cannot allocate nx", -1);
  for(i=0; i<sphere_pm->NDP; i++) {
    nx[i] = (int *) calloc(DIMS, sizeof(int));
    if(nx[i] == 0) fatal_err("cannot allocate nx[i]", -1);
  }

  com_shift = (int **) calloc(sphere_pm->NDP, sizeof(*com_shift));
  if(com_shift == 0) fatal_err("cannot allocate com_shift", -1);
  for(i=0; i<sphere_pm->NDP; i++) {
    com_shift[i] = (int *) calloc(DIMS, sizeof(int));
    if(com_shift[i] == 0) fatal_err("cannot allocate com_shift[i]", -1);
  }

  /* allocate memory for temporary node map */
  /* temporary node map spans the stretch dimensions of a DP */
  tmpnode = (int ***) calloc(sphere_pm->NDP, sizeof(*tmpnode));
  if(tmpnode == 0) fatal_err("cannot allocate temp node map", -1);

  for(n=0; n<sphere_pm->NDP; n++) {
    for(j=0; j<DIMS; j++) {
      nx[n][j] = (int)(sph[n].tstretch[j] + 6);
     
      if(nx[n][j] > maxsize[j]) 
	nx[n][j] = maxsize[j];
    }

    nxy = (nx[n][0]+1)*nx[n][1]+(nx[n][1]+1);

    tmp_pp = (int **) calloc(nxy, sizeof(*tmp_pp));
    if(tmp_pp == 0) fatal_err("cannot allocate temp node map 1", -1);
    tmpnode[n] = tmp_pp;

    for(i = 0; i<nxy; i++) {
      tmp_p = (int *) calloc(nx[n][2]+1, sizeof(*tmp_p));
      if(tmp_p == 0) fatal_err("cannot allocate temp node map 2", -1);
      tmpnode[n][i] = tmp_p;
    }

    /* Initialize node map */
    for(i=0; i<nxy; i++)
      for(j=0; j<=nx[n][2]; j++)
	tmpnode[n][i][j]=0;
  }

  /* Test to see if the particle position and shape has changed */
  double move = 0.0;
  int test_move = 0;
  for(n=0; n<sphere_pm->NDP; n++) {
    for(j=0; j<DIMS; j++) {
      st_offset[j] = 1;
      
      move += (sph[n].tcom[j]-sph[n].tcomold[j])*(sph[n].tcom[j]-sph[n].tcomold[j])+(sph[n].tstretch[j]-sph[n].tstretch_old[j])*(sph[n].tstretch[j]-sph[n].tstretch_old[j]);
    }
  }

  test_move = 0;
  if (move > 0.01 || nstep == 0) 
    test_move = 1;

  /*define nodes inside each DP */
  if(test_move == 1) {
#pragma omp parallel for num_threads(NTHREAD)
    for(n=0; n<sphere_pm->NDP; n++) {
      int i, j, k;
      int xy, z, xy1, leftbound, rightbound;
      int n1, n2;
      int tmpposlat[DIMS][2];
      double tpos[DIMS];  /* thread-local temporary position */
      int type = (n<Ntype[0] ? 0 : 1);
      int index = (type==0 ? n*Nbeads[0] : Ntype[0]*Nbeads[0] + (n-Ntype[0])*Nbeads[1]);

      if(DP_DEBUG == 1)
	printf("com[0] = (%lf %lf %lf) comold = (%lf %lf %lf) stretch = (%lf %lf %lf) stretch_old= (%lf %lf %lf)\n", sph[0].tcom[0], sph[0].tcom[1], sph[0].tcom[2], sph[0].tcomold[0], sph[0].tcomold[1], sph[0].tcomold[2], sph[0].tstretch[0], sph[0].tstretch[1], sph[0].tstretch[2], sph[0].tstretch_old[0], sph[0].tstretch_old[1], sph[0].tstretch_old[2]);

      for(d=0; d<DIMS; d++) {
	sph[n].tcomold[d] = sph[n].tcom[d];
	sph[n].tstretch_old[d] = sph[n].tstretch[d];
      }

      /* First, mark nodes around the particle boundary vertices and midpoint vertices */
      for(j=0; j<DIMS; j++)
	com_shift[n][j] = (int)(sph[n].tstretch[j]/2.0+st_offset[j])-(int)(sph[n].tcom[j]);

      for(i=0; i<Nbeads[type]; i++) {
	/*copy to temporary map.*/
	for(j=0; j<DIMS; j++) {
	  tpos[j] = mon[index+i].pos_pbc[j] + com_shift[n][j] + 1.0;
	  tmpposlat[j][0] = (int)(tpos[j]);
	  tmpposlat[j][1] = (int)(tpos[j])+1;

	  if(tmpposlat[j][0] < 0)
	    fatal_err("monomer outside temp nodemap in bnode_dp", -1);
	}

	/*mark nodes surrounding monomers */
	for(j=0; j<2; j++)
	  for(k=0; k<2; k++)
	    for(l=0; l<2; l++) {
	      xy = (tmpposlat[0][j]*nx[n][1])+tmpposlat[1][k];
	      z = tmpposlat[2][l];
	      tmpnode[n][xy][z] = 2;
	    }

	/* also mark the nodes neighboring the midpoint 'phantom' vertices' */
	/*copy to temporary map.*/
	for(m=1; m<=mon[index+i].blist[0][0]; m++) {
	  n2 = mon[index+i].blist[m][0];
	  n1 = index+i;
	  for(j=0; j<DIMS; j++) {
	    tpos[j] = (mon[n2].pos_pbc[j]+mon[n1].pos_pbc[j])/2.0  // phantom vertex at midpoint
	      + com_shift[n][j] + 1.0;
	    tmpposlat[j][0] = (int)(tpos[j]);
	    tmpposlat[j][1] = (int)(tpos[j])+1;
	  }

	  /*mark nodes surrounding monomers */
	  for(j=0; j<2; j++)
	    for(k=0; k<2; k++)
	      for(l=0; l<2; l++) {
		xy = (tmpposlat[0][j]*nx[n][1])+tmpposlat[1][k];
		z = tmpposlat[2][l];
		tmpnode[n][xy][z] = 2;
	      }
	}
      }

      /* scan in all dimensions to find boundaries, then mark all nodes within the nodes */
      for(k=0; k<=nx[n][2]; k++) {
      	/* scan in x */
      	for(j=0; j<=nx[n][1]; j++) {
      	  leftbound = 1;
      	  rightbound = 0;

      	  for(i=0; i<nx[n][0]; i++) {
      	    xy = i*nx[n][1]+j;
      	    xy1 = (i+1)*nx[n][1]+j;
	  
      	    if(tmpnode[n][xy][k] == 0 && tmpnode[n][xy1][k] == 2) {
      	      leftbound = i+1;
	      tmpnode[n][xy1][k] = 0;
	  
      	      break;
      	    }
      	  }
	
      	  for(i=nx[n][0]; i>0; i--) {
      	    xy = i*nx[n][1]+j;
      	    xy1 = (i-1)*nx[n][1]+j;
      	    if(tmpnode[n][xy][k] ==0 && tmpnode[n][xy1][k] == 2) {
      	      rightbound = i-1;
	      tmpnode[n][xy1][k] = 0;
	  
      	      break;
      	    }
      	  }
      
      	  for(i=leftbound; i<rightbound; i++) {
      	    xy = i*nx[n][1]+j;
      	    tmpnode[n][xy][k] = -1;

	    if(DP_DEBUG == 1)
	      printf("ScanX x (%d %d %d) left %d right %d %d \n", i, j, xy, leftbound, rightbound, tmpnode[n][xy][k]);
      	  }
      	}

      	/* scan in y */
      	for(i=0; i<=nx[n][0]; i++) {
      	  leftbound = 1;
      	  rightbound = 0;
 
      	  for(j=0; j<nx[n][1]; j++) {
      	    xy = i*nx[n][1]+j;
      	    xy1 = i*nx[n][1]+j+1;

      	    if(tmpnode[n][xy][k] == 0 && tmpnode[n][xy1][k] != 0) {
	      leftbound = j+1;
	      tmpnode[n][xy1][k] = 0;

      	      break;
      	    }
      	  }
	
      	  for(j=nx[n][1]; j>0; j--) {
      	    xy = i*nx[n][1]+j;
      	    xy1 = i*nx[n][1]+(j-1);
	        
      	    if(tmpnode[n][xy][k] ==0 && tmpnode[n][xy1][k] != 0) {
      	      rightbound = j-1;
	      tmpnode[n][xy1][k] = 0;

      	      break;
      	    }
      	  }
      
      	  for(j=leftbound; j<rightbound; j++) {
      	    xy = i*nx[n][1]+j;
      	    tmpnode[n][xy][k]-=5;

	    if(DP_DEBUG == 1)
	      printf("ScanY j (%d %d %d) left %d right %d %d \n", i, j, xy, leftbound, rightbound, tmpnode[n][xy][k]);
      	  }
      	}
      }

      /* scan in z */
      for(i=0; i<=nx[n][0]; i++) {
      	for(j=0; j<=nx[n][1]; j++) {
      	  xy = i*nx[n][1]+j;
      	  leftbound = 0;
      	  rightbound = 0;

      	  for(k=0; k<nx[n][2]; k++) {
      	    if(tmpnode[n][xy][k] == 0 && tmpnode[n][xy][k+1] != 0) {
      	      leftbound = k+1;
	      tmpnode[n][xy][k+1] = 0;
      	      break;
      	    }
      	  }

      	  for(k=nx[n][2]; k>0; k--) {
      	    if(tmpnode[n][xy][k] ==0 && tmpnode[n][xy][k-1] != 0) {
      	      rightbound = k-1;
      	      tmpnode[n][xy][k-1] = 0;
      	      break;
      	    }
      	  }

      	  for(k=leftbound+1; k<rightbound; k++) {
	    tmpnode[n][xy][k]-=3;
	    if(DP_DEBUG == 1)
	      printf("ScanZ k (%d %d) left %d right %d %d \n", i, j, leftbound, rightbound, tmpnode[n][xy][k]);
      	  }
      	}
      }
    } 
  }

  /*map tmpnode to node_map */
  /* To avoid removing nodes between DP that are right next to each other */
  /* We clear all outside DP nodes first, then mark inside DP nodes */
  /* convert tmpnode coords to node_map coords */
  if(test_move == 1) {
#pragma omp parallel for num_threads(NTHREAD)
    for(n=0; n<sphere_pm->NDP; n++) {
      int i, j, k;
      int xy, z;
      int px[DIMS]; 
      
      for(i=0; i<=nx[n][0]; i++) {
	px[0] = box2((i-1-com_shift[n][0]), maxsize[0])%num_x+1;
	if(px[0] < 0) fatal_err("px[0] outside nodemap in bnode_dp", -1);

	for(j=0; j<=nx[n][1]; j++) {
	  px[1] = box2((j-1-com_shift[n][1]), maxsize[1])+1;
	  if(px[1] < 0) fatal_err("px[1] outside nodemap in bnode_dp", -1);

	  xy = px[0]*(maxsize[1]+2)+px[1];
	  
	  for(k=0; k<=nx[n][2]; k++) {
	    px[2] = box2((k-1-com_shift[n][2]), maxsize[2])+1;
	    if(px[2] < 0) fatal_err("px[2] outside nodemap in bnode_dp", -1);
	    
	    if(DP_DEBUG == 1)
	      printf("%d (%d %d %d) (%lf %lf %lf) (%lf %lf %lf) (%d %d %d) (%d %d %d) %d\n", n, com_shift[n][0], com_shift[n][1], com_shift[n][2], sph[n].tcom[0], sph[n].tcom[1], sph[n].tcom[2], sph[n].tstretch[0], sph[n].tstretch[1], sph[n].tstretch[2], i, j, k, px[0], px[1], px[2], xy);
	    
	    /* compare with node_map to delete overlap nodes outside DP */
	    
	    if(node_map[xy][px[2]] == 2 && tmpnode[n][i*nx[n][1]+j][k] == 0 )
	      node_map[xy][px[2]] = 0;
	  }
	}
      }
    }
  }

  if(test_move == 1) {
    for(n=0; n<sphere_pm->NDP; n++) {
      for(i=0; i<=nx[n][0]; i++) {
	px[0] = box2((i-1-com_shift[n][0]), maxsize[0])%num_x+1;
      
	for(j=0; j<=nx[n][1]; j++) {
	  px[1] = box2((j-1-com_shift[n][1]), maxsize[1])+1;
	  
	  xy = px[0]*(maxsize[1]+2)+px[1];
	  
	  for(k=0; k<=nx[n][2]; k++) {
	    px[2] = box2((k-1-com_shift[n][2]), maxsize[2])+1;
	    
	    /* mark nodes inside DP */
	    if(tmpnode[n][i*nx[n][1]+j][k] == -9) 
	      node_map[xy][px[2]] = 2;

	    if(DP_DEBUG == 1)
	      node_map[xy][px[2]] = tmpnode[n][i*nx[n][1]+j][k];

	  }
	}
      }
    }
  }

  for(i=0; i<sphere_pm->NDP; i++) {
    nxy = (nx[i][0]+1)*nx[i][1]+(nx[i][1]+1);
    for(j=0; j<nxy; j++)
      free(tmpnode[i][j]);
  }
  
  for(i=0; i<sphere_pm->NDP; i++) {
    free(tmpnode[i]);
    free(nx[i]);
    free(com_shift[i]);
  }

  free(tmpnode);
  free(nx);
  free(com_shift);
} 


void sphere_props_short(struct DP_param *sphere_pm, struct DP *spheres, struct monomer *mono)
{
  int i, j, k, d, n;
  int type, index;
  int NDP = sphere_pm->NDP;
  int Ntype[NTYPES];
  int Nbeads[NTYPES];
  int maxsize[DIMS];
  double temp, tempstretch;

  for(i=0; i<NTYPES; i++) {
    Nbeads[i] = sphere_pm->N_per_DP[i];
    Ntype[i] = sphere_pm->Ntype[i];
  }

  /* calculate the sphere center-of-mass */
  for(i=0; i<NDP; i++) {
    type = (i<Ntype[0] ? 0 : 1);
    index = (type==0 ? i*Nbeads[0] : Ntype[0]*Nbeads[0] + (i-Ntype[0])*Nbeads[1]);

    for(d=0; d<DIMS; d++) 
      spheres[i].tcom[d] = 0.0;
/* 
    for(j=0; j<Nbeads[type]; j++)
      for(d=0; d<DIMS; d++) 
        spheres[i].tcom[d] += mono[index+j].pos_pbc[d];

    for(d=0; d<DIMS; d++) 
      spheres[i].tcom[d] = spheres[i].tcom[d]/Nbeads[type];
*/
    /* calculate the sphere  stretch */
    /* determine the maximum stretch in x, y, and z */
    for(d=0; d<DIMS; d++) 
      spheres[i].tstretch[d] = 0.0;

    for(j=0; j<Nbeads[type]-1; j++) {
      for(k=j+1; k<Nbeads[type]; k++) {
        for(d=0; d<DIMS; d++) {
	  tempstretch = 0.0;

          temp = mono[index+k].pos_pbc[d]-mono[index+j].pos_pbc[d];
          tempstretch = sqrt(temp*temp);
	  if(tempstretch > spheres[i].tstretch[d])  {
	    spheres[i].tstretch[d] = tempstretch;
	    spheres[i].tcom[d] = (mono[index+k].pos_pbc[d]+mono[index+j].pos_pbc[d])/2.0;
	  }
        }
      }
    }
  }
}
