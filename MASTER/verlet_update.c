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
#define DEBUG 0

/* void fluc_vel(Float *force); */
void temp_rescale(struct monomer *mon, struct DP_param *sphere_pm);
//void face_force(double sigma_LJ, int n1, int p1, int p2, int p3, double h, double *normal, struct DP_param *, struct monomer *);

void verlet_update(struct monomer *mon, struct face *faces, struct DP *DP, struct DP_param *sphere_pm, struct DP_param *polym_pm, struct DP_param *cyl_pm, Float ***velcs_df, int **node_map, int n_step, VSLStreamStatePtr rngstream)
{
  extern int max_x, max_y, max_z;
  extern int wall_flag;

  double maxsize[3];
  long i,j;
  int d, bead0;
  int num_beads = sphere_pm->num_beads+polym_pm->num_beads+cyl_pm->num_beads;
  int num_sphpo = sphere_pm->num_beads+polym_pm->num_beads;
  int overlap;
  double dt = sphere_pm->dt;
  double h_dt = dt/2.0;
  double h_dt2= h_dt*h_dt;
  double trialpos[DIMS], trialvel[DIMS];
  double dr[DIMS];
  double mon_mass = Rho_Fl*sphere_pm->monmass;
  double x = sphere_pm->fric*dt/sphere_pm->monmass;
  FILE *stream;

  maxsize[0]=(double)max_x;
  maxsize[1]=(double)max_y;
  maxsize[2]=(double)max_z;

  /* Determine DP COM for force calculations */
  for(i=0; i<sphere_pm->Ntype[0]; i++) {
    for(d=0; d<DIMS; d++)
      DP[i].com[d] = 0.0;

    for(j=i*sphere_pm->N_per_DP[0]; j<(i+1)*sphere_pm->N_per_DP[0]; j++) 
      for(d=0; d<DIMS; d++)
	DP[i].com[d] += mon[j].pos_pbc[d];
    
    for(d=0; d<DIMS; d++) 
      DP[i].com[d] /= sphere_pm->N_per_DP[0];
  }

  for(i=sphere_pm->Ntype[0]; i<sphere_pm->Ntype[0]+sphere_pm->Ntype[1]; i++) {
      for(d=0; d<DIMS; d++)
      DP[i].com[d] = 0.0;

      bead0=sphere_pm->Ntype[0]*sphere_pm->N_per_DP[0];
      for(j=bead0+(i-sphere_pm->Ntype[0])*sphere_pm->N_per_DP[1]; j<bead0+(i+1-sphere_pm->Ntype[0])*sphere_pm->N_per_DP[1]; j++) 
      for(d=0; d<DIMS; d++)
	DP[i].com[d] += mon[j].pos_pbc[d];
    
    for(d=0; d<DIMS; d++)
      DP[i].com[d] /= sphere_pm->N_per_DP[1];
  }

  if(sphere_pm->verlet == 0) /* no position update */
    {
      get_forces(sphere_pm, polym_pm, cyl_pm, mon, faces, DP, velcs_df, node_map, n_step, rngstream);

#pragma omp parallel for num_threads(NTHREAD)
      for(i=0;i<num_beads;i++) {
	int d;
	double mon_mass;

	if(i < sphere_pm->num_beads)
	  mon_mass = Rho_Fl*sphere_pm->monmass;
	else if(i < num_sphpo)
	  mon_mass = Rho_Fl*polym_pm->monmass;
	else
	  mon_mass = Rho_Fl*cyl_pm->monmass;

	for(d=0; d<DIMS; d++) {
	  mon[i].vel[d] += mon[i].force[d]/mon_mass*dt;
	  mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
	}
      }
    } 
  else if(sphere_pm->verlet==1)  /* velocity verlet update */
    {
#pragma omp parallel for num_threads(NTHREAD)
      for(i=0;i<num_beads;i++)	{
	int d;
	double mon_mass, dr[DIMS];
	if(i < sphere_pm->num_beads) 
	  mon_mass = Rho_Fl*sphere_pm->monmass;
	else if(i < num_sphpo)
	  mon_mass = Rho_Fl*polym_pm->monmass;
	else
	  mon_mass = Rho_Fl*cyl_pm->monmass;
	
	for(d=0; d<DIMS; d++) {
	  mon[i].pos_tmp2[d]=mon[i].pos_tmp[d];
	  mon[i].pos_tmp[d]=mon[i].pos_pbc[d];
	  dr[d]=dt*mon[i].vel[d]+2.0*h_dt2*mon[i].force0[d]/mon_mass;
	  mon[i].pos_pbc[d]=mon[i].pos_pbc[d]+dr[d];
	}
	
	if(i < num_sphpo) 
	  for(d=0; d<DIMS; d++)
	    mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
      }

      /* Check for bead-DP overlap and perform bead reflection, only for nlist */
      if(sphere_pm->ev_type == 3) 
	check_overlap(mon, faces, DP, sphere_pm, n_step);
	//	if(check_overlap( mon, faces, DP, sphere_pm) == 1) {
	  //	  for(i=0; i<20; i++)
	  //	    finestep_verlet(dt/20.0, mon, faces, DP, sphere_pm, polym_pm, cyl_pm, velcs_df, node_map, n_step, rngstream);	 
	//	}
      //      } 
      
      get_forces(sphere_pm, polym_pm, cyl_pm, mon, faces, DP, velcs_df, node_map, n_step, rngstream);
      //      if(cyl_pm->NDP > 0 && cyl_pm->initconfig != 6)
      //      	get_forces_cyl(sphere_pm, polym_pm, cyl_pm, mon, velcs_df, node_map, n_step, rngstream);

#pragma omp parallel for num_threads(NTHREAD)
      for(i=0;i<num_beads;i++) {
	double mon_mass;
	int d;
	if(i < sphere_pm->num_beads) 
	  mon_mass = Rho_Fl*sphere_pm->monmass;
	else if(i < num_sphpo)
	  mon_mass = Rho_Fl*polym_pm->monmass;
	else
	  mon_mass = Rho_Fl*cyl_pm->monmass;
	
	//	printf("bead %d force = (%le %le %le)\n", i, mon[i].force[0], mon[i].force[1], mon[i].force[2]);
	for(d=0; d<DIMS; d++) {
	  mon[i].vel_tmp[d] = mon[i].vel[d];
	  mon[i].vel[d]+=h_dt*(mon[i].force0[d]+mon[i].force[d])/mon_mass;
	  mon[i].force0[d]=mon[i].force[d];
	}
      }
    }
  else if(sphere_pm->verlet == 2) /* explicit 1st order */
    {
      get_forces(sphere_pm, polym_pm, cyl_pm, mon, faces, DP,  velcs_df, node_map, n_step, rngstream);
      //      if(cyl_pm->NDP > 0 && cyl_pm->initconfig != 6)
      //      	get_forces_cyl(sphere_pm, polym_pm, cyl_pm, mon, velcs_df, node_map, n_step, rngstream);
      
      for(i=0;i<num_beads;i++) {
	if(i < sphere_pm->num_beads) 
	  mon_mass = Rho_Fl*sphere_pm->monmass;
	else if(i < num_sphpo)
	  mon_mass = Rho_Fl*polym_pm->monmass;
	else
	  mon_mass = Rho_Fl*cyl_pm->monmass;

	for(d=0; d<DIMS; d++) {
	  trialvel[d] = (mon[i].vel[d]+dt*mon[i].force[d]/mon_mass);
	  dr[d] = dt*(trialvel[d]+mon[i].vel[d])/2.0;
	  trialpos[d]=mon[i].pos_pbc[d]+dr[d];
	}

	for(d=0; d<DIMS; d++) {
	  mon[i].vel[d]=trialvel[d];
	  mon[i].pos_pbc[d]=trialpos[d];
	}

	if(i < num_sphpo) 
	  for(d=0; d<DIMS; d++)
	    mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);

      }
    }
  else if(sphere_pm->verlet == 3) /* implicit 1st order */
    {
      for(i=0;i<num_beads;i++) {
	if(i < sphere_pm->num_beads) 
	  mon_mass = Rho_Fl*sphere_pm->monmass;
	else if(i < num_sphpo)
	  mon_mass = Rho_Fl*polym_pm->monmass;
	else
	  mon_mass = Rho_Fl*cyl_pm->monmass;

	for(d=0; d<DIMS; d++) {
	  mon[i].pos_tmp2[d]=mon[i].pos_tmp[d];
	  mon[i].pos_tmp[d]=mon[i].pos_pbc[d];
	  mon[i].pos_pbc[d] +=mon[i].vel[d]*dt;
	}

	if(i < num_sphpo) 
	  for(d=0; d<DIMS; d++)
	    mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
      }

      get_forces(sphere_pm, polym_pm, cyl_pm, mon, faces, DP, velcs_df, node_map, n_step, rngstream);
      //      if(cyl_pm->NDP > 0 && cyl_pm->initconfig != 6)
      //	get_forces_cyl(sphere_pm, polym_pm, cyl_pm, mon, velcs_df, node_map, n_step, rngstream);

      for(i=0;i<num_beads;i++) {
	if(i < sphere_pm->num_beads) 
	  mon_mass = Rho_Fl*sphere_pm->monmass;
	else if(i < num_sphpo)
	  mon_mass = Rho_Fl*polym_pm->monmass;
	else
	  mon_mass = Rho_Fl*cyl_pm->monmass;

	for(d=0; d<DIMS; d++) {
	  trialvel[d] = (mon[i].vel[d]+dt*(mon[i].force[d])/mon_mass);
	  dr[d] = dt*(trialvel[d]+mon[i].vel[d])/2.0;
	  mon[i].vel[d]=trialvel[d];
	  mon[i].pos_pbc[d]=mon[i].pos_tmp[d]+dr[d];
	}

	if(i < num_sphpo) 
	  for(d=0; d<DIMS; d++)
	    mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
      }
    }
}

void temp_rescale(struct monomer *mon, struct DP_param *sphere_pm)
{
  int i, d;
  double kT_ob; 
  double scale;
  double avg_vel2[DIMS];
  
  for(d=0; d<DIMS; d++)
    avg_vel2[d] = 0.0;

  for(i=0; i<sphere_pm->num_beads; i++) {
    for(d=0; d<DIMS; d++)
      avg_vel2[d] += mon[i].vel[d]*mon[i].vel[d];
  }

  kT_ob = 0.0;
  for(d=0; d<DIMS; d++)
    kT_ob += avg_vel2[d];

  kT_ob /= sphere_pm->num_beads;

  if(kT_ob < 1e-10)
    scale = 0.0;
  else
    scale = sphere_pm->kT / (kT_ob*Rho_Fl/DIMS);

  sphere_pm->tempscale = scale;
}

/* check for overlap between DP */
/* if a bead penetrate another DP, move it back outside, reflect boundary condition */
int check_overlap(struct monomer *mon, struct face *face, struct DP *DP, struct DP_param *sphere_pm, int n_step)
{
  static int t_counter;
  extern int max_x, max_y, max_z;
  extern int wall_flag;
  const int stack = 6;
  signed long iseed;

  double sigma = 2.0*mon[0].radius;
  double cutoff = 1.5*sigma;
  double cutoff2 = cutoff*cutoff;
  double dt = sphere_pm->dt;
  double h_dt = dt/2.0;
  double h_dt2= h_dt*h_dt;

  int flag;
  int n, i, j, d, n1, n2, nn1, q, k;
  int n0, nf, iface;
  int p1, p2, p3;
  int maxsize[DIMS];
  double p1x, p1y, p1z;
  double p2x, p2y, p2z;
  double p3x, p3y, p3z;
  double det, num1, num2, num3, sum;
  double dr, dr1, dr2, dr1_tmp, dr2_tmp; 
  double mon_mass, drtemp[DIMS];
  double temp, h, h1, h2;
  double com[DIMS], min_com[stack][DIMS], d_com[DIMS], com_tmp[DIMS], DPcom[DIMS];
  double q21[DIMS], q31[DIMS], q32[DIMS];
  double q11[DIMS], q12[DIMS], q13[DIMS], q11mag, q12mag, q13mag;
  double dv11[DIMS], dv21[DIMS], dv31[DIMS];
  double q11_tmp[DIMS], q12_tmp[DIMS], q21_tmp[DIMS], q31_tmp[DIMS], q32_tmp[DIMS], q21q31_tmp[DIMS];
  double q21dv31[DIMS], dv21q31[DIMS], dv21dv31[DIMS];

  int min_face[stack], min_bead[stack];
  double normal[stack][DIMS], norm_tmp[stack][DIMS], pos[stack][DIMS];
  double normalmag, normtmpmag;
  double min_face_dr[stack];
  double hpos[stack][DIMS], tempvel[stack][DIMS], temppos[stack][DIMS], temp_dcom[stack][DIMS], tempforce[stack][DIMS];
  double temptpvel[stack][DIMS], atemp[stack], temp3[DIMS];
  double min_dcom[stack][DIMS], min_d2com[stack][DIMS],  mincom_tmp[stack][DIMS];

  maxsize[0]=max_x;
  maxsize[1]=max_y;
  maxsize[2]=max_z;

  t_counter++;

  /* find the two faces bead n1 is closest to */
  n0 = 0; nf = sphere_pm->num_beads;
  for(n1=n0; n1<nf; n1++) {
    nn1 = n1;

    for(q=0; q<stack; q++) {
      min_bead[q] = sphere_pm->num_beads + 1;
      min_face_dr[q] =100;
      min_face[q] = -1;
    }

    for(n=1; n <= mon[nn1].list[0]; n++) {
      n2 = mon[nn1].list[n];
      if((n2 < sphere_pm->num_beads) && (mon[nn1].DP_id != mon[n2].DP_id)) {
	dr2=0.0;
	for(d=0; d<DIMS; d++) {
	  temp=mon[nn1].pos_pbc[d] - mon[n2].pos_pbc[d];
	  if ((d==0 && wall_flag<3) || (d==2 && wall_flag<2) || (d==1 && wall_flag<1))
	    temp=n_image(temp, maxsize[d]);
	  dr2+=temp*temp;
	}

	/* calculate bead - face_com distances if the bead-bead distance is smaller than cutoff */
	if(dr2 > cutoff2) 
	  continue;
	else {
	  for(i=1; i<=mon[n2].face_id[0]; i++) {
	    /* check if face has already been looped.  If so, skip */
	    flag = 0;
	    iface = mon[n2].face_id[i];

	    for(k=0; k<stack; k++) 
	      if(mon[n2].face_id[i] == min_face[k]) {
		flag = 1;
		break;
	      }
	    if(flag == 1) continue;

	    for(d=0; d<DIMS; d++) {
	      com[d] = 0.0;
	      com_tmp[d] = 0.0;

	      for(j=0; j<DIMS; j++) {
		com[d]+=mon[face[iface].vertices[j]].pos_pbc[d];
		com_tmp[d]+=mon[face[iface].vertices[j]].pos_tmp[d];
	      }
	      com[d]/=(double)DIMS;
	      com_tmp[d]/=(double)DIMS;
	    }
	    
	    dr=0.0;
	    for(d=0; d<DIMS; d++) {
	      temp = mon[nn1].pos_pbc[d] - com[d];
	      if ((d==0 && wall_flag<3) || (d==2 && wall_flag<2) || (d==1 && wall_flag<1))
		d_com[d]=n_image(temp, maxsize[d]);
	      else
		d_com[d]=temp;

	      dr += d_com[d]*d_com[d];
	    }
	    
	    dr = sqrt(dr);
	    
	    if(dr > cutoff) continue;

	    if(dr < min_face_dr[stack-1]) {
	      min_face_dr[stack-1] = dr;
	      min_bead[stack-1] = n2;
	      min_face[stack-1] = mon[n2].face_id[i];
	      
	      for(d=0; d<DIMS; d++) {
		min_com[stack-1][d]=com[d];
		mincom_tmp[stack-1][d]=com_tmp[d];

		min_dcom[stack-1][d]=d_com[d];

		temp = mon[nn1].pos_tmp[d] - mincom_tmp[stack-1][d];
		if ((d==0 && wall_flag<3) || (d==2 && wall_flag<2) || (d==1 && wall_flag<1))
		  min_d2com[stack-1][d]=n_image(temp, maxsize[d]);
		else
		  min_d2com[stack-1][d]=temp;
	      }
	    }
	  
	    for(q=stack-2; q>=0; q--) {
	      if(dr < min_face_dr[q]) {
		min_face_dr[q+1] = min_face_dr[q];
		min_bead[q+1] = min_bead[q];
		min_face[q+1] = min_face[q];

		min_face_dr[q] = dr;
		min_bead[q] = n2;
		min_face[q] = mon[n2].face_id[i];

		for(d=0; d<DIMS; d++) {
		  min_com[q+1][d]=min_com[q][d];
		  mincom_tmp[q+1][d]= mincom_tmp[q][d];
		  min_dcom[q+1][d]=min_dcom[q][d];
		  min_d2com[q+1][d] = min_d2com[q][d];

		  min_com[q][d]=com[d];
		  mincom_tmp[q][d]=com_tmp[d];

		  min_dcom[q][d]=d_com[d];
	
		  temp = mon[nn1].pos_tmp[d] - mincom_tmp[q][d];
		  if ((d==0 && wall_flag<3) || (d==2 && wall_flag<2) || (d==1 && wall_flag<1))
		    min_d2com[q][d]=n_image(temp, maxsize[d]);
		  else
		    min_d2com[q][d]=temp;
		}
	      }
	    }
	  }
	}
      }
    }

    /* if there are no nearby faces, continue to next bead */
    if(min_face_dr[0] == 100 || min_face[0] == -1) {
      flag = 0;
      continue;
    }
    
    for(i=0; i<stack; i++) {
      /* loop over all vertices of the closest two faces */
      if(min_face_dr[i] > cutoff || min_bead[i] == sphere_pm->num_beads+1) 
	continue;

      int ii=1-i;
      
      p1 = face[min_face[i]].vertices[0];
      p2 = face[min_face[i]].vertices[1];
      p3 = face[min_face[i]].vertices[2];
      for (d=0 ; d<DIMS ; d++) {
	q21[d] = mon[p2].pos_pbc[d] - mon[p1].pos_pbc[d];
	q31[d] = mon[p3].pos_pbc[d] - mon[p1].pos_pbc[d];
	q32[d] = mon[p3].pos_pbc[d] - mon[p2].pos_pbc[d];

	q11_tmp[d] = mon[nn1].pos_tmp[d] - mon[p1].pos_tmp[d];
	q21_tmp[d] = mon[p2].pos_tmp[d] - mon[p1].pos_tmp[d];
	q31_tmp[d] = mon[p3].pos_tmp[d] - mon[p1].pos_tmp[d];
	q32_tmp[d] = mon[p3].pos_tmp[d] - mon[p2].pos_tmp[d];

	q11[d] = mon[nn1].pos_pbc[d] - mon[p1].pos_pbc[d];
	q12[d] = mon[nn1].pos_pbc[d] - mon[p2].pos_pbc[d];
	q13[d] = mon[nn1].pos_pbc[d] - mon[p3].pos_pbc[d];
	
	dv21[d] = mon[p2].vel[d] - mon[p1].vel[d];
	dv31[d] = mon[p3].vel[d] - mon[p1].vel[d];
	dv11[d] = mon[nn1].vel[d] - mon[p1].vel[d];
      }
    
      q11mag = 0.0; q12mag = 0.0; q13mag = 0.0;
      for(d=0; d<DIMS; d++) {
	q11mag += q11[d]*q11[d];
	q12mag += q12[d]*q12[d];
	q13mag += q13[d]*q13[d];
      }

      /* if the bead is very close to any one of the vertices, do not reflect and let LJ force repel */
      product(q21, q31, temp3);
      product(q21_tmp, q31_tmp, q21q31_tmp);

      product(q21_tmp, dv31, q21dv31);
      product(dv21, q31_tmp, dv21q31);
      product(dv21, dv31, dv21dv31);

      for(d=0; d<DIMS; d++)
	normal[i][d] = temp3[d];

      normalmag = 0.0;
      normtmpmag = 0.0;
      for(d=0; d<DIMS; d++) {
	normalmag += normal[i][d]*normal[i][d];
	normtmpmag += q21q31_tmp[d]*q21q31_tmp[d];
      }
      normalmag = sqrt(normalmag);
      normtmpmag = sqrt(normtmpmag);
      
      for(d=0; d<DIMS; d++) {
	normal[i][d] /= normalmag;
	norm_tmp[i][d] = q21q31_tmp[d]/normtmpmag;
      }

      /* determine the sign of the normal by comparing with distance of face center+normal from particle center */
      for(d=0; d<DIMS; d++)
	DPcom[d] = DP[mon[min_bead[i]].DP_id].com[d];

      dr1=0.0;
      dr1_tmp = 0.0;
      double normaldotdr1 = 0.0;
      for(d=0; d<DIMS; d++) {
	temp = min_com[i][d]-DPcom[d];
	if ((d==0 && wall_flag<3) || (d==2 && wall_flag<2) || (d==1 && wall_flag<1))
	  temp=n_image(temp, maxsize[d]);
	dr1+=temp*temp;
	normaldotdr1 += temp*normal[q][d];
	
	temp = mincom_tmp[i][d]-DPcom[d];
	if ((d==0 && wall_flag<3) || (d==2 && wall_flag<2) || (d==1 && wall_flag<1))
	  temp=n_image(temp, maxsize[d]);
	dr1_tmp+=temp*temp;
      }

      dr1 = sqrt(dr1);
      dr1_tmp = sqrt(dr1_tmp);
      normaldotdr1 /= dr1;
      if(normaldotdr1 < 0.0) normaldotdr1 = -normaldotdr1;
    
      dr2=0.0;
      dr2_tmp=0.0;
      for(d=0; d<DIMS; d++) {
	temp = min_com[i][d]+dr1/2.0*normal[i][d]-DPcom[d];
	if ((d==0 && wall_flag<3) || (d==2 && wall_flag<2) || (d==1 && wall_flag<1))
	  temp=n_image(temp, maxsize[d]);
	dr2+=temp*temp;

	temp = mincom_tmp[i][d]+dr1_tmp/2.0*norm_tmp[i][d]-DPcom[d];
	if ((d==0 && wall_flag<3) || (d==2 && wall_flag<2) || (d==1 && wall_flag<1))
	  temp=n_image(temp, maxsize[d]);
	dr2_tmp+=temp*temp;
      }
          
      dr2 = sqrt(dr2);
      dr2_tmp = sqrt(dr2_tmp);
      
      /* added condition to check for pinched concave surface */
      /* if cosine between the face normal and the COM-to-face vector < 0.5 */
      /* use the particle of n1 to determine the normal sign*/
      if(normaldotdr1 < 0.5) {
	for(j=0; j<DIMS; j++) {
	  temp = mon[n1].pos_pbc[j]-com_tmp[j];
	  dr1+=temp*temp;
	  normaldotdr1 += temp*normal[q][j];
	}
	dr1=sqrt(dr1);
	normaldotdr1 /= dr1;
	if(normaldotdr1 < 0.0) normaldotdr1 = -normaldotdr1;

	dr2=0.0;
	for(j=0; j<DIMS; j++) {
	  temp = mon[n1].pos_pbc[j]-com_tmp[j]+dr1/2.0*normal[q][j];
	  dr2+=temp*temp;
	}
      
	dr2 = sqrt(dr2);

	if(dr2 > dr1) 
	  for(d=0; d<DIMS; d++)
	    normal[q][d] = -1.0 * normal[q][d];
      }
      else {
	if(dr1 > dr2) 
	  for(d=0; d<DIMS; d++) {
	    normal[i][d] = -1.0 * normal[i][d];
	    norm_tmp[i][d] = -1.0 * norm_tmp[i][d];
	  }
      }
    
      /* Check if bead nn1 overlaps with face 0 */
      h=normal[i][0]*min_dcom[i][0]+normal[i][1]*min_dcom[i][1]+normal[i][2]*min_dcom[i][2]; 
      h1=norm_tmp[i][0]*min_d2com[i][0]+norm_tmp[i][1]*min_d2com[i][1]+norm_tmp[i][2]*min_d2com[i][2]; 

      if(h < 0.0) {
	if(DEBUG == 1)
	  printf("Test 1 step=%d count=%d n1=%d q=%d dr=%le n2=%d face=%d h=%le h1=%le\n", n_step, t_counter, nn1, i, min_face_dr[i], n2, min_face[i], h, h1);

	double a, a0, a1, a2, a3, inner, b1, b2;
	a0 = 0.0;
	a1 = 0.0;
	a2 = 0.0;
	a3 = 0.0;
	
	for(d=0; d<DIMS; d++) {
	  a0 += q11_tmp[d]*q21q31_tmp[d];
	  a1 += q11_tmp[d]*(q21dv31[d]+dv21q31[d]) + dv11[d]*q21q31_tmp[d];
	  a2 += q11_tmp[d]*dv21dv31[d]+dv11[d]*(q21dv31[d]+dv21q31[d]);
	  a3 += dv11[d]*dv21dv31[d];
	}
	
	a = -(a0 / a1) / dt;
	inner = a1*a1-4*a2*a0;
	b1=0.0, b2=0.0;
	if(inner > 0.0) {
	  b1 = (-a1 + sqrt(a1*a1-4*a2*a0))/(2*a2) / dt;
	  b2 = (-a1 - sqrt(a1*a1-4*a2*a0))/(2*a2) / dt;
	}
	
	/* Reflect bead position across face */
	if(a > 1.0 || a < 0.0) {
	  if(DEBUG == 1) 
	    printf("Test 1a step=%d count=%d n1=%d q=%d face=%d a=%le, a0=%le a1=%le b1=%le b2=%le\n", n_step, t_counter, nn1, i, min_face[i], a, a0, a1, b1, b2);
	  continue;
	}
	
	for(d=0; d<DIMS; d++) 
	  pos[i][d] = mon[nn1].pos_tmp[d]+a*mon[nn1].vel[d]*dt;
	
	p1x = mon[p1].pos_pbc[0];   p2x = mon[p2].pos_pbc[0]; p3x = mon[p3].pos_pbc[0];
	p1y = mon[p1].pos_pbc[1];   p2y = mon[p2].pos_pbc[1]; p3y = mon[p3].pos_pbc[1];
	p1z = mon[p1].pos_pbc[2];   p2z = mon[p2].pos_pbc[2]; p3z = mon[p3].pos_pbc[2];
	
	det = p1x*(p2y*p3z-p3y*p2z)+p1y*(p3x*p2z-p2x*p3z)+p1z*(p2x*p3y-p3x*p2y);
	num1 = (pos[i][0]*(p2y*p3z-p3y*p2z)+pos[i][1]*(p2z*p3x-p3z*p2x)+pos[i][2]*(p3y*p2x-p2y*p3x))/det;
	num2 = (pos[i][0]*(p3y*p1z-p3z*p1y)+pos[i][1]*(p3z*p1x-p1z*p3x)+pos[i][2]*(p3x*p1y-p3y*p1x))/det;
	num3 = (pos[i][0]*(p1y*p2z-p1z*p2y)+pos[i][1]*(p1z*p2x-p1x*p2z)+pos[i][2]*(p1x*p2y-p2x*p1y))/det;
      
	/* check if n1 does cross within the face */
	if(num1 < 0.0 || num2 < 0.0 || num3 < 0.0) {
	  if(DEBUG == 1) 
	    printf("Test 1a step=%d count=%d n1=%d q=%d face=%d num1=%le num2=%le num3=%le\n", n_step, t_counter, nn1, i, min_face[i], num1, num2, num3);
	  continue;
	}
	if(num1 > 1.0 || num2 > 1.0 || num3 > 1.0) {
	  if(DEBUG == 1) 
	    printf("Test 1a step=%d count=%d n1=%d q=%d face=%d num1=%le num2=%le num3=%le\n", n_step, t_counter, nn1, i, min_face[i], num1, num2, num3);
	  continue;
	}
	
	sum = num1+num2+num3; 
	num1 /= sum; 
	num2 /= sum; 
	num3 /= sum; 
	
	/* Calculation using the procedure of Fedosov */  
	/* The results of xi and ze are the same as num2 and num3 */
	double p1postmp[DIMS], p2postmp[DIMS], p3postmp[DIMS];
	double ga[DIMS], qa21[DIMS], qa31[DIMS]; 
	double qa21qa31, gqa21, gqa31, qa21mag, qa31mag, xi, ze, denom;
	
	for(d=0; d<DIMS; d++) {
	  p1postmp[d] = (mon[p1].pos_pbc[d] + mon[p1].vel[d]*a*dt );
	  p2postmp[d] = (mon[p2].pos_pbc[d] + mon[p2].vel[d]*a*dt );
	  p3postmp[d] = (mon[p3].pos_pbc[d] + mon[p3].vel[d]*a*dt );
	}
	
	for(d=0; d<DIMS; d++) {
	  ga[d] = pos[i][d] - p1postmp[d];
	  qa21[d] = p2postmp[d] - p1postmp[d];
	  qa31[d] = p3postmp[d] - p1postmp[d];
	}
	
	qa21mag = 0.0; qa31mag = 0.0; gqa21=0.0; gqa31 = 0.0; qa21qa31 = 0.0;
	for(d=0; d<DIMS; d++) {
	  qa21mag += qa21[d]*qa21[d];
	  qa31mag += qa31[d]*qa31[d];
	  gqa21 += ga[d]*qa21[d];
	  gqa31 += ga[d]*qa31[d];
	  qa21qa31 += qa21[d]*qa31[d];
	}
      
	denom = (qa21mag*qa31mag - qa21qa31*qa21qa31);
	xi = (gqa21*qa31mag - gqa31*qa21qa31) / denom;
	ze = (gqa31*qa21mag - gqa21*qa21qa31) / denom;
	
	/* check if n1 does cross within the face */
	if(xi < 0.0 || ze < 0.0 || (xi+ze) > 1.0) {
	  if(DEBUG == 1) 
	    printf("Test 1a step=%d count=%d n1=%d q=%d face=%d xi=%le ze=%le\n", n_step, t_counter, nn1, i, min_face[i], xi, ze);
	  continue;
	}
	
	/*
	  num1 = 1-xi-ze;
	  num2 = xi;
	  num3 = ze;
	*/
	
	/* bounce back */
	mon_mass = Rho_Fl*sphere_pm->monmass;
	for(d=0; d<DIMS; d++) {
	  tempvel[i][d] = -mon[nn1].vel[d] + 2*(num1*mon[p1].vel[d]+num2*mon[p2].vel[d]+num3*mon[p3].vel[d]);
	  temppos[i][d] = pos[i][d] + tempvel[i][d]*(1.0-a)*dt;

	  //tempvel[i][j] = -mon[nn1].vel[j] + 2*((1-xi-ze)*mon[p1].vel[j]+xi*mon[p2].vel[j]+ze*mon[p3].vel[j]);
	  //	  tempforce[i][d] = -mon[nn1].force0[d] + 2*(num1*mon[p1].force0[d]+num2*mon[p2].force0[d]+num3*mon[p3].force0[d]);
	}
	
	for(d=0; d<DIMS; d++)
	  temp_dcom[i][d] = temppos[i][d] - min_com[i][d];
	hpos[i][0]=normal[i][0]*temp_dcom[i][0]+normal[i][1]*temp_dcom[i][1]+normal[i][2]*temp_dcom[i][2];

	if(DEBUG >= 1) {
	  printf("Particle overlap n1=%d i=%d h=%le h1=%le h2=%le\n", nn1, i, h, h1, h2);
	  printf("step=%d t_count=%d face=%d a=%le, b1=%le b2=%le\n", n_step, t_counter, min_face[i], a, b1, b2);
	  printf("n2=%d %d face=%d %d dr=%le %le\n", min_bead[0], min_bead[1], min_face[0], min_face[1], min_face_dr[0], min_face_dr[1]);
	  printf("a=%le a0=%le, a1=%le a2=%le, a3=%le b1=%le b2=%le\n", a, a0, a1, a2, a3, b1, b2);
	  //      printf("qa21mag=%le, qa31mag=%le, qa21qa31=%le, gqa21=%le, gqa31=%le\n", qa21mag, qa31mag, qa21qa31, gqa21, gqa31);
	  printf("xi=%le, ze=%le, num1=%le, num2=%le, num3=%le sum=%le\n", xi, ze, num1, num2, num3, sum);
	  printf("norm= %le %le %le \n", normal[i][0], normal[i][1], normal[i][2]);
	  printf("norm_tmp = %le %le %le \n", norm_tmp[i][0], norm_tmp[i][1], norm_tmp[i][2]);
	  printf("min_dcom = %le %le %le\n", min_dcom[i][0], min_dcom[i][1], min_dcom[i][2]);
	  printf("min_d2com = %le %le %le\n", min_d2com[i][0], min_d2com[i][1], min_d2com[i][2]);
	  printf("DPcom = %le %le %le\n", DPcom[0], DPcom[1], DPcom[2]);
	  printf("com = %le %le %le\n", min_com[i][0], min_com[i][1], min_com[i][2]);
	  printf("com_tmp = %le %le %le\n", mincom_tmp[i][0], mincom_tmp[i][1], mincom_tmp[i][2]);
	  printf("n1 %d n2 %d p1 p2 p3 %d %d %d\n", nn1, min_bead[i], p1, p2, p3);
	  printf("n1pos %lf %lf %lf\n", mon[nn1].pos_pbc[0], mon[nn1].pos_pbc[1], mon[nn1].pos_pbc[2]);
	  printf("n1tmp %lf %lf %lf\n", mon[nn1].pos_tmp[0], mon[nn1].pos_tmp[1], mon[nn1].pos_tmp[2]);
	  printf("n1vel %lf %lf %lf\n", mon[nn1].vel[0], mon[nn1].vel[1], mon[nn1].vel[2]);
	  printf("n1frc %lf %lf %lf\n", mon[nn1].force0[0], mon[nn1].force0[1], mon[nn1].force0[2]);
	  printf("p1pos %lf %lf %lf\n", mon[p1].pos_pbc[0], mon[p1].pos_pbc[1], mon[p1].pos_pbc[2]);
	  printf("p1tmp %lf %lf %lf\n", mon[p1].pos_tmp[0], mon[p1].pos_tmp[1], mon[p1].pos_tmp[2]);
	  printf("p2pos %lf %lf %lf\n", mon[p2].pos_pbc[0], mon[p2].pos_pbc[1], mon[p2].pos_pbc[2]);
	  printf("p2tmp %lf %lf %lf\n", mon[p2].pos_tmp[0], mon[p2].pos_tmp[1], mon[p2].pos_tmp[2]);
	  printf("p3pos %lf %lf %lf\n", mon[p3].pos_pbc[0], mon[p3].pos_pbc[1], mon[p3].pos_pbc[2]);
	  printf("p3tmp %lf %lf %lf\n", mon[p3].pos_tmp[0], mon[p3].pos_tmp[1], mon[p3].pos_tmp[2]);
	  printf("p1_vx=%le p1_vy=%le p1_vz=%le\n", mon[p1].vel[0], mon[p1].vel[1], mon[p1].vel[2]);
	  printf("p2_vx=%le p2_vy=%le p2_vz=%le\n", mon[p2].vel[0], mon[p2].vel[1], mon[p2].vel[2]);
	  printf("p3_vx=%le p3_vy=%le p3_vz=%le\n", mon[p3].vel[0], mon[p3].vel[1], mon[p3].vel[2]);

	  /* There is a scenario for h1 < 0 */
	  /* At previous timestep t', the bead position p' is outside a face f' that is not the same nearby face f as at time t*/
	  /* and the normal distance h1 between p' and f is < 0. */
	  /* However, the bounce back calculation depends only on the crossing face , which is f. */
	  /* So this is resolved. */

	  if(h1 < 0.0) {
	    printf("h1 < 0! i=%d\n", i);
	    printf("n1 %d n2 %d p1 p2 p3 %d %d %d\n", nn1, min_bead[i], p1, p2, p3);
	    printf("n1tmp %lf %lf %lf\n", mon[nn1].pos_tmp[0], mon[nn1].pos_tmp[1], mon[nn1].pos_tmp[2]);
	    int p1_tmp, p2_tmp, p3_tmp;
	    for(q=0; q< stack; q++) {
	      p1_tmp=face[min_face[q]].vertices[0];
	      p2_tmp=face[min_face[q]].vertices[1];
	      p3_tmp=face[min_face[q]].vertices[2];
	      printf("t' minface[%d] %d p1 %d p2 %d p3 %d \n", min_face[q], q, p1_tmp, p2_tmp, p3_tmp); 
	      printf("p1_tmp = %le %le %le\n", mon[p1_tmp].pos_pbc[0], mon[p1_tmp].pos_pbc[1], mon[p1_tmp].pos_pbc[2]);
	      printf("p2_tmp = %le %le %le\n", mon[p2_tmp].pos_pbc[0], mon[p2_tmp].pos_pbc[1], mon[p2_tmp].pos_pbc[2]);
	      printf("p3_tmp = %le %le %le\n", mon[p3_tmp].pos_pbc[0], mon[p3_tmp].pos_pbc[1], mon[p3_tmp].pos_pbc[2]);
	      printf("p1_vx=%le p1_vy=%le p1_vz=%le\n", mon[p1].vel_tmp[0], mon[p1].vel_tmp[1], mon[p1].vel_tmp[2]);
	      printf("p2_vx=%le p2_vy=%le p2_vz=%le\n", mon[p2].vel_tmp[0], mon[p2].vel_tmp[1], mon[p2].vel_tmp[2]);
	      printf("p3_vx=%le p3_vy=%le p3_vz=%le\n", mon[p3].vel_tmp[0], mon[p3].vel_tmp[1], mon[p3].vel_tmp[2]);
	    }
	  }
	}

	for(d=0; d<DIMS; d++) {
	  mon[nn1].vel[d] = tempvel[i][d];
	  mon[nn1].pos_pbc[d] = temppos[i][d];
	  mon[nn1].pos[d]=box(mon[nn1].pos_pbc[d], maxsize[d]);
	}

	if(DEBUG >= 1) {  
	  printf("bb2 n1=%d n2=%d face=%d face_dr=%le h=%le h1=%le\n", nn1, min_bead[i], min_face[i], min_face_dr[i], hpos[i][0], h1);
	  printf("monpos %le %le %le\n", mon[nn1].pos_pbc[0], mon[nn1].pos_pbc[1], mon[nn1].pos_pbc[2]);
	  printf("monvel %le %le %le\n", mon[nn1].vel[0], mon[nn1].vel[1], mon[nn1].vel[2]);
      	}
      }
    }
    
    if(DEBUG == 1) {
      //       if(mon[nn1].force0[0] > 10 || mon[nn1].force0[1] > 10 || mon[nn1].force0[2] > 10 || mon[nn1].force0[0] < -10 || mon[nn1].force0[1] < -10 || mon[nn1].force0[2] < -10) {
      /*
	for(d=0; d<DIMS; d++)
	min_dcom[0][d] = mon[nn1].pos[d] - min_com[0][d];
	h=normal[0][0]*min_dcom[0][0]+normal[0][1]*min_dcom[0][1]+normal[0][2]*min_dcom[0][2];
      */
      //       }
    }
  }
  return 0;
}


/*  Replaced by check_overlap between bead and neighboring DP in verlet_update */
/*  Adding repulsion between a bead and a nearby face */
/* void face_force(double sigma_LJ, int n1, int p1, int p2, int p3, double h, double *normal, struct DP_param *sphere_pm, struct monomer *mon) */
/* { */
/*   int i; */
/*   double eps = sphere_pm->eps*sphere_pm->kT/sigma_LJ; */
/*   double aterm, num1, num2, num3, det, sum; */
/*   double pos[DIMS], force[DIMS]; */
/*   double p1x,p2x,p3x,p1y,p2y,p3y,p1z,p2z,p3z; */

/*   aterm = 0.25*(sigma_LJ/h); */
/*   aterm = aterm*aterm; */
/*   aterm = aterm*aterm*aterm; */
  
/*   p1x = mon[p1].pos_pbc[0];   p2x = mon[p2].pos_pbc[0]; p3x = mon[p3].pos_pbc[0];  */
/*   p1y = mon[p1].pos_pbc[1];   p2y = mon[p2].pos_pbc[1]; p3y = mon[p3].pos_pbc[1];  */
/*   p1z = mon[p1].pos_pbc[2];   p2z = mon[p2].pos_pbc[2]; p3z = mon[p3].pos_pbc[2];  */

/*   /\* Reflect bead position across face *\/ */
/*   /\* Apply face force *\/ */
/*   for(i=0; i<DIMS; i++) { */
/*     pos[i] = mon[n1].pos_pbc[i]-h*normal[i]; */
/*     force[i]= 24.0*(eps*(2.0*(aterm*aterm)-aterm)) * normal[i]; */

/*     det = p3x*p2y*p1z - p2x*p3y*p1z - p3x*p1y*p2z+ p1x*p3y*p2z + p2x*p1y*p3z - p1x*p2y*p3z; */
/*     num1 = (pos[2]*p2y*p1z-pos[1]*p3y*p1z-pos[2]*p1y*p2z+pos[0]*p3y*p2z+pos[1]*p1y*p3z-pos[0]*p2y*p3z)/det; */
/*     num2 = -(pos[2]*p2x*p1z - pos[1]*p3x*p1z - pos[2]*p1x*p2z + pos[0]*p3x*p2z + pos[1]*p1x*p3z - pos[0]*p2x*p3z)/det; */
/*     num3 = (pos[2]*p2x*p1y-pos[1]*p3x*p1y-pos[2]*p1x*p2y+pos[0]*p3x*p2y+pos[1]*p1x*p3y-pos[0]*p2x*p3y)/det; */
/*   } */

/*   for(i=0; i<DIMS; i++) { */
/*     mon[n1].force[i] -= force[i]; */
/*     mon[n1].evforce[i] -= force[i]; */
/*   } */
  
/*   for (i=0 ; i<DIMS ; i++) { */
/*     mon[p1].force[i] += force[i]*num1; */
/*     mon[p2].force[i] += force[i]*num2; */
/*     mon[p3].force[i] += force[i]*num3; */
/*     mon[p1].evforce[i] += force[i]*num1; */
/*     mon[p2].evforce[i] += force[i]*num2; */
/*     mon[p3].evforce[i] += force[i]*num3; */
/*   } */
/* } */

void finestep_verlet(double dt, struct monomer *mon, struct face *faces, struct DP *DP, struct DP_param *sphere_pm, struct DP_param *polym_pm, struct DP_param *cyl_pm, Float ***velcs_df, int **node_map, int n_step, VSLStreamStatePtr rngstream)
{
  extern int max_x, max_y, max_z;
  extern int wall_flag;

  double maxsize[3];
  long i,j;
  int d, bead0;
  int num_beads = sphere_pm->num_beads+polym_pm->num_beads+cyl_pm->num_beads;
  int num_sphpo = sphere_pm->num_beads+polym_pm->num_beads;
  double h_dt = dt/2.0;
  double h_dt2= h_dt*h_dt;
  double trialpos[DIMS], trialvel[DIMS];
  double dr[DIMS];
  double mon_mass = Rho_Fl*sphere_pm->monmass;
  double x = sphere_pm->fric*dt/sphere_pm->monmass;
  FILE *stream;

  maxsize[0]=(double)max_x;
  maxsize[1]=(double)max_y;
  maxsize[2]=(double)max_z;

  /* Determine DP COM for force calculations */
  for(i=0; i<sphere_pm->Ntype[0]; i++) {
    for(d=0; d<DIMS; d++)
      DP[i].com[d] = 0.0;

    for(j=i*sphere_pm->N_per_DP[0]; j<(i+1)*sphere_pm->N_per_DP[0]; j++) 
      for(d=0; d<DIMS; d++)
	DP[i].com[d] += mon[j].pos_pbc[d];
    
    for(d=0; d<DIMS; d++) 
      DP[i].com[d] /= sphere_pm->N_per_DP[0];
  }

  for(i=sphere_pm->Ntype[0]; i<sphere_pm->Ntype[0]+sphere_pm->Ntype[1]; i++) {
    for(d=0; d<DIMS; d++)
      DP[i].com[d] = 0.0;

    bead0=sphere_pm->Ntype[0]*sphere_pm->N_per_DP[0];
    for(j=bead0+(i-sphere_pm->Ntype[0])*sphere_pm->N_per_DP[1]; j<bead0+(i+1-sphere_pm->Ntype[0])*sphere_pm->N_per_DP[1]; j++) 
      for(d=0; d<DIMS; d++)
	DP[i].com[d] += mon[j].pos_pbc[d];
    
    for(d=0; d<DIMS; d++)
      DP[i].com[d] /= sphere_pm->N_per_DP[1];
  }

  if(sphere_pm->verlet == 0) /* no position update */
    {
      get_forces(sphere_pm, polym_pm, cyl_pm, mon, faces, DP, velcs_df, node_map, n_step, rngstream);
      //      if(cyl_pm->NDP > 0 && cyl_pm->initconfig != 6)
      //	get_forces_cyl(sphere_pm, polym_pm, cyl_pm, mon, velcs_df, node_map, n_step, rngstream);

#pragma omp parallel for num_threads(NTHREAD)
      for(i=0;i<num_beads;i++) {
	int d;
	double mon_mass;

	if(i < sphere_pm->num_beads)
	  mon_mass = Rho_Fl*sphere_pm->monmass;
	else if(i < num_sphpo)
	  mon_mass = Rho_Fl*polym_pm->monmass;
	else
	  mon_mass = Rho_Fl*cyl_pm->monmass;

	for(d=0; d<DIMS; d++) {
	  mon[i].vel[d] += mon[i].force[d]/mon_mass*dt;
	  mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
	}
      }
    } 
  else if(sphere_pm->verlet==1)  /* velocity verlet update */
    {
#pragma omp parallel for num_threads(NTHREAD)
      for(i=0;i<num_beads;i++)	{
	int d;
	double mon_mass, dr[DIMS];
	if(i < sphere_pm->num_beads) 
	  mon_mass = Rho_Fl*sphere_pm->monmass;
	else if(i < num_sphpo)
	  mon_mass = Rho_Fl*polym_pm->monmass;
	else
	  mon_mass = Rho_Fl*cyl_pm->monmass;
	
	for(d=0; d<DIMS; d++) {
	  dr[d]=dt*mon[i].vel[d]+2.0*h_dt2*mon[i].force0[d]/mon_mass;
	  mon[i].pos_pbc[d]=mon[i].pos_pbc[d]+dr[d];
	}
	
	if(i < num_sphpo) 
	  for(d=0; d<DIMS; d++)
	    mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
      }
      
      get_forces(sphere_pm, polym_pm, cyl_pm, mon, faces, DP, velcs_df, node_map, n_step, rngstream);
      
#pragma omp parallel for num_threads(NTHREAD)
      for(i=0;i<num_beads;i++) {
	double mon_mass;
	int d;
	if(i < sphere_pm->num_beads) 
	  mon_mass = Rho_Fl*sphere_pm->monmass;
	else if(i < num_sphpo)
	  mon_mass = Rho_Fl*polym_pm->monmass;
	else
	  mon_mass = Rho_Fl*cyl_pm->monmass;

	for(d=0; d<DIMS; d++) {
	  mon[i].vel[d]+=h_dt*(mon[i].force0[d]+mon[i].force[d])/mon_mass;
	  mon[i].force0[d]=mon[i].force[d];
	}
      }
    }
  else if(sphere_pm->verlet == 2) /* explicit 1st order */
    {
      get_forces(sphere_pm, polym_pm, cyl_pm, mon, faces, DP,  velcs_df, node_map, n_step, rngstream);
      //      if(cyl_pm->NDP > 0 && cyl_pm->initconfig != 6)
      //      	get_forces_cyl(sphere_pm, polym_pm, cyl_pm, mon, velcs_df, node_map, n_step, rngstream);
      
      for(i=0;i<num_beads;i++) {
	if(i < sphere_pm->num_beads) 
	  mon_mass = Rho_Fl*sphere_pm->monmass;
	else if(i < num_sphpo)
	  mon_mass = Rho_Fl*polym_pm->monmass;
	else
	  mon_mass = Rho_Fl*cyl_pm->monmass;

	for(d=0; d<DIMS; d++) {
	  trialvel[d] = (mon[i].vel[d]+dt*mon[i].force[d]/mon_mass);
	  dr[d] = dt*(trialvel[d]+mon[i].vel[d])/2.0;
	  trialpos[d]=mon[i].pos_pbc[d]+dr[d];
	}

	for(d=0; d<DIMS; d++) {
	  mon[i].vel[d]=trialvel[d];
	  mon[i].pos_pbc[d]=trialpos[d];
	}

	if(i < num_sphpo) 
	  for(d=0; d<DIMS; d++)
	    mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);

      }
    }
  else if(sphere_pm->verlet == 3) /* implicit 1st order */
    {
      for(i=0;i<num_beads;i++) {
	if(i < sphere_pm->num_beads) 
	  mon_mass = Rho_Fl*sphere_pm->monmass;
	else if(i < num_sphpo)
	  mon_mass = Rho_Fl*polym_pm->monmass;
	else
	  mon_mass = Rho_Fl*cyl_pm->monmass;

	for(d=0; d<DIMS; d++) {
	  mon[i].pos_tmp[d]=mon[i].pos_pbc[d];
	  mon[i].pos_pbc[d] +=mon[i].vel[d]*dt;
	}

	if(i < num_sphpo) 
	  for(d=0; d<DIMS; d++)
	    mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
      }

      get_forces(sphere_pm, polym_pm, cyl_pm, mon, faces, DP, velcs_df, node_map, n_step, rngstream);

      for(i=0;i<num_beads;i++) {
	if(i < sphere_pm->num_beads) 
	  mon_mass = Rho_Fl*sphere_pm->monmass;
	else if(i < num_sphpo)
	  mon_mass = Rho_Fl*polym_pm->monmass;
	else
	  mon_mass = Rho_Fl*cyl_pm->monmass;

	for(d=0; d<DIMS; d++) {
	  trialvel[d] = (mon[i].vel[d]+dt*(mon[i].force[d])/mon_mass);
	  dr[d] = dt*(trialvel[d]+mon[i].vel[d])/2.0;
	  mon[i].vel[d]=trialvel[d];
	  mon[i].pos_pbc[d]=mon[i].pos_tmp[d]+dr[d];
	}

	if(i < num_sphpo) 
	  for(d=0; d<DIMS; d++)
	    mon[i].pos[d]=box(mon[i].pos_pbc[d], maxsize[d]);
      }
    }
}
