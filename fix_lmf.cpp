/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_lmf.h"
#include "atom.h"
#include "input.h"
#include "variable.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "modify.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include "atom.h"
#include "universe.h"
#include "update.h"
#include "lmf_grid.h"
#include <complex>
#include "run.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{XLO=0,XHI=1,YLO=2,YHI=3,ZLO=4,ZHI=5};
enum{NONE=0,EDGE,CONSTANT,VARIABLE};

#define SMALL 0.00001


Fix_LMF::Fix_LMF(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg){
    lmf_sigma=0;
    lmf_write=-1;

    ensemble_flag=0; no_force=0; root_bool=0; slab_flag=0; SOR_flag=0;
    fldflag=0; verbose=0; omega=0; debug=0; DLMF_flag=0;
    lmf_c=0; lmf_tau=0; sc_freq=0; sc_count=0; lmf_it_freq=0;
    
    lmf_pf = force->qqr2e/force->qe2f*4.0*3.14159265359;
    force_flag=0;
    
    n_tot=new int [3];
	grid_spac=new double [3];
	grid_inv_sq=new double [3];
	lmf_lo=new double [3];
	

    for(int dim=0; dim<3; dim++){
        grid_spac[dim]=0; grid_inv_sq[dim]=0; 
        n_tot[dim]=0; lmf_lo[dim]=0;
    }

    int iarg = 3;
    while (iarg < narg) {
        if (strcmp(arg[iarg],"ensemble") == 0){
            ensemble_flag=1;
            iarg+=1;
            }
        else if (strcmp(arg[iarg],"verbose") == 0){
            verbose=1;
            iarg+=1;
            }
        else if (strcmp(arg[iarg],"debug") == 0){
            debug=1;
            iarg+=1;
            }
        else if (strcmp(arg[iarg],"par") == 0){
            if (iarg+3 > narg ) error->all(FLERR,"Illegal fix LMF command");
            lmf_sigma = force->numeric(FLERR,arg[iarg+1]);
            lmf_cutoff = force->numeric(FLERR,arg[iarg+2]); 
            iarg+=3;
        }   
        else if (strcmp(arg[iarg],"dlmf") == 0){
            if (iarg+5 > narg ) error->all(FLERR,"Illegal fix LMF command");
            DLMF_flag=1;
            lmf_c  = force->numeric(FLERR,arg[iarg+1]);
            lmf_tau = force->numeric(FLERR,arg[iarg+2]);
            sc_freq = int(force->numeric(FLERR,arg[iarg+3]));
            lmf_it_freq = int(force->numeric(FLERR,arg[iarg+4]));
            if(lmf_it_freq<=0) error->all(FLERR,"lmf_it_freq must be a positive integer");
            iarg+=5;
        }   
        else if (strcmp(arg[iarg],"slab")==0){
            if (iarg+4 > narg ) error->all(FLERR,"Illegal fix LMF command");
            slab_bool=1;
            grid_spac[0]=force->numeric(FLERR,arg[iarg+1]);
            grid_spac[2]=force->numeric(FLERR,arg[iarg+2]);
            yb=int(force->numeric(FLERR,arg[iarg+3]));
            grid_spac[1]=grid_spac[0];
            if(yb<3) error->all(FLERR,"yb must be greater than 2");
            slab_flag=1;
            iarg+=4;
        }
        else if (strcmp(arg[iarg],"iso")==0){
            if (iarg+2 > narg ) error->all(FLERR,"Illegal fix LMF command");
            double grid_spac_temp=force->numeric(FLERR,arg[iarg+1]);
            for(int dim=0; dim<3; dim++) grid_spac[dim]=grid_spac_temp;
            iarg+=2;
        }
        else if (strcmp(arg[iarg],"write")==0){
            if (iarg+3 > narg ) error->all(FLERR,"Illegal fix LMF command");
            lmf_write=force->numeric(FLERR,arg[iarg+1]);
            char * temp_name = arg[iarg+2];
            lmf_file_name = (char *) malloc(sizeof(char)*(strlen(temp_name)+1));
            strcpy(lmf_file_name,temp_name);
            printf("LMF file name %s \n", lmf_file_name);
            iarg+=3;
        }
        else if (strcmp(arg[iarg],"noforce")==0){
            if (iarg+1 > narg ) error->all(FLERR,"Illegal fix LMF command");
            no_force=1;
            iarg+=1;
        }
        else if (strcmp(arg[iarg],"sor")==0){
            if (iarg+3 > narg ) error->all(FLERR,"Illegal fix LMF command");
            SOR_flag=1;
            omega=force->numeric(FLERR,arg[iarg+1]);
            SOR_tol = force->numeric(FLERR,arg[iarg+2]);
            iarg+=3;
        }
        else if (strcmp(arg[iarg],"force_write") == 0){
            if (iarg+2 > narg ) error->all(FLERR,"Illegal fix analysis command");
            force_file_name = arg[iarg+1];
            force_flag=1;
            int *type=atom->type;
            int nlocal = atom->nlocal;
            for (int i = 0; i < nlocal; i++)
                if (type[i]==1){
                    tag_particle=i;
                    break;
                }
            iarg+=2;
        }
        else {
            printf("arg %s \n", arg[iarg]);
            error->all(FLERR,"Illegal fix LMF command");
        }
    }
    
    if(force_flag){
        if(universe->me==0) {
    	    force_file = fopen(force_file_name,"w");
        }
    } 
    if(DLMF_flag && SOR_flag) error->all(FLERR,"Only one LMF update method allowed");
    if(!DLMF_flag && !SOR_flag) error->all(FLERR,"LMF method not specified");
    for(int dim=0; dim<3; dim++) grid_inv_sq[dim]=1.0/grid_spac[dim]/grid_spac[dim];
    setup_LMF_lattice();
    if(lmf_sigma==0) error->all(FLERR,"Did not specify sigma in Fix LMF");	
    root_bool=0;
    for(int world=0; world<universe->nworlds; world++){
        if(universe->me==universe->root_proc[world]) 
            {
            root_bool=1;
            printf("root %d \n", universe->me);
            }
    }
    
}

Fix_LMF::~Fix_LMF(){
	delete LMF_potential;
    delete smooth_charge; 
    delete LMF_pi;
    delete [] n_tot;
	delete [] grid_spac;
	delete [] grid_inv_sq;
	delete [] lmf_lo;
	free(lmf_file_name);
}


void Fix_LMF::setup_LMF_lattice()
{
    lmf_vol=1.0;
    for(int dim_temp=0; dim_temp<3; dim_temp++){
        lmf_lo[dim_temp]=domain->boxlo[dim_temp];
        double size_temp=domain->boxhi[dim_temp]-lmf_lo[dim_temp];
        if(slab_flag==1 && dim_temp==2){
            lmf_lo[dim_temp]=lmf_lo[dim_temp]-size_temp;
            size_temp=size_temp*float(yb);
        }
        n_tot[dim_temp] = int(ceil(size_temp/grid_spac[dim_temp]))+1;
		grid_spac[dim_temp]= size_temp/double(n_tot[dim_temp]-1);
		lmf_vol*=size_temp;
	}
    
    LMF_potential = new LMF_GRID(n_tot,grid_spac,lmf_lo, ensemble_flag);
	smooth_charge = new LMF_GRID(n_tot,grid_spac,lmf_lo, ensemble_flag);
	
	if(DLMF_flag){
	    LMF_pi = new LMF_GRID(n_tot,grid_spac,lmf_lo, ensemble_flag);
	}
}

int Fix_LMF::modify_param(int narg, char **arg)
{
    if (strcmp(arg[0],"par") == 0){
        if (narg<3 ) error->all(FLERR,"Illegal fix_modify LMF command");
        lmf_sigma = force->numeric(FLERR,arg[1]);
        lmf_cutoff = force->numeric(FLERR,arg[2]); 
        return 3;
    }
    // heat tbd
    else if (strcmp(arg[0],"write") == 0){
        if (narg<2 ) error->all(FLERR,"Illegal fix_modify LMF command");
        lmf_write = force->numeric(FLERR,arg[1]);
        return 2;
    }
    else if (strcmp(arg[0],"sor") == 0){
        if (narg<3 ) error->all(FLERR,"Illegal fix_modify LMF command");
        SOR_flag=1;
        omega=force->numeric(FLERR,arg[1]);
        SOR_tol = force->numeric(FLERR,arg[2]);
        return 3;
    }
    else if (strcmp(arg[0],"dlmf") == 0){
        if (narg<5 ) error->all(FLERR,"Illegal fix_modify LMF command");
        DLMF_flag=1;
        lmf_c=force->numeric(FLERR,arg[1]);
        lmf_tau=force->numeric(FLERR,arg[2]);
        sc_freq = int(force->numeric(FLERR,arg[3]));
        lmf_it_freq = int(force->numeric(FLERR,arg[4]));
        if(lmf_it_freq<=0) error->all(FLERR,"lmf_it_freq must be a positive integer");
        return 5;
    }
    return 0;
}

void Fix_LMF::LMF_update(){
    if(SOR_flag){
        smooth_charge->zero();
        smooth_calc();
        poisson_SOR();
    }
    else if(DLMF_flag){
        if(sc_count==0){
            smooth_charge->zero();
            smooth_calc();
            subtract_av(LMF_potential);
            }
        dlmf_iterate();    
        sc_count=sc_count+1;
        if(sc_count==sc_freq) sc_count=0;
    }
}

int Fix_LMF::poisson_SOR(){
    double Lap=0, base=0; 
	bool EDGE=0;
	double pre_factor = force->qqr2e/force->qe2f*4.0*3.14159;
	double tot_change=0;
	int n_steps;
	bool finish_temp=0;
	for (n_steps=0; n_steps<1000000; n_steps++){
	    tot_change=SOR_iterate(LMF_potential, smooth_charge, pre_factor, n_tot);
	    //if(verbose && debug) printf("Time = %lld n_steps %d tot_change %f omega %f \n", update->ntimestep, n_steps, tot_change, omega);
	    if(tot_change<SOR_tol){
	        finish_temp=1;
	        break;
	        }
	    
	}
	if(!finish_temp)printf("LMF tolerance not reached \n");
	if(verbose) printf("Time = %lld n_steps %d Total difference = %.14f \n", update->ntimestep, n_steps, tot_change);
    return n_steps;
}

void Fix_LMF::dlmf_iterate(){
    double Lap=0, base=0; 
	bool EDGE=0;
	double tot_change=0;
	double lmf_c_s=lmf_c*lmf_c;
	double pi_change;
	double time_step = update->dt/double(lmf_it_freq);
	int * nt = LMF_potential->nt;
	for (int n_steps=0; n_steps<lmf_it_freq; n_steps++){
	    for(int n_0=1; n_0<nt[0]-1; n_0++){
		    for(int n_1=1; n_1<nt[1]-1; n_1++){
			    for(int n_2=1; n_2<nt[2]-1; n_2++){
                    int hi0 = nt[0]-2; int hi1=nt[1]-2;
                    pi_change=lmf_c_s*(LMF_potential->lap_calc(n_0, n_1, n_2)+lmf_pf*smooth_charge->array[n_0][n_1][n_2]);
                    LMF_potential->array[n_0][n_1][n_2]+=time_step*(LMF_pi->array[n_0][n_1][n_2]-LMF_potential->array[n_0][n_1][n_2]/lmf_tau);
                    LMF_pi->array[n_0][n_1][n_2]+=time_step*pi_change;
                    update_boundary(n_0,n_1,n_2,LMF_potential->array, nt);
                    update_boundary(n_0,n_1,n_2,LMF_pi->array, nt);
			    }
		    }
	    }
	}
	if(verbose) printf("Time = %lld Total difference = %.14f \n", update->ntimestep, tot_change);
}

void Fix_LMF::LMF_force(){
    int pos_index[3];

    double **x = atom->x;
    double **f = atom->f;
    double *q = atom->q;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    bigint ntimestep=update->ntimestep;

    double gradient[3];
    if(slab_flag) compute_yb_correction();
    
    pos_index_calc(0,pos_index,x, NULL, 0);
    LMF_potential->gradient_calc(pos_index, gradient);
    if(verbose)printf("yb %f f %f type %d f_long %f pos %f \n", q[0]*yb_correction, f[0][2], type[0], q[0]*force->qe2f*gradient[2], x[0][2]); 
                
    
    if(!no_force){
        for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) {
                double pre_factor = q[i]*force->qe2f; //convert to kcal/mol
                pos_index_calc(i,pos_index,x, NULL, 0);
                LMF_potential->gradient_calc(pos_index, gradient);
                for(int dim=0; dim<3; dim++){
                    f[i][dim]-=pre_factor*gradient[dim];
                }
                if(slab_flag)f[i][2]+=q[i]*yb_correction;
                if(force_flag && universe->me==0 && i==tag_particle)
                    fprintf(force_file,"%d %f %f %f\n",update->ntimestep,-pre_factor*gradient[0], -pre_factor*gradient[1], -pre_factor*gradient[2]);
            } 
    }
    else if(force_flag && universe->me==0){
        double pre_factor = q[tag_particle]*force->qe2f; //convert to kcal/mol
        pos_index_calc(tag_particle,pos_index,x, NULL, 0);
        LMF_potential->gradient_calc(pos_index, gradient);
        fprintf(force_file,"%d %f %f %f %f %f %f \n",update->ntimestep,-pre_factor*gradient[0], -pre_factor*gradient[1], -pre_factor*gradient[2], ew_force[0],ew_force[1],ew_force[2]);
    }        
    
}

void Fix_LMF::smooth_calc(){
    double **x = atom->x;
    int *mask = atom->mask;
    double *q = atom->q;
    int nall = atom->nlocal+atom->nghost;
    double charge;
    int n_steps[3];
    for (int dim=0;dim<3;dim++) n_steps[dim]=ceil(lmf_cutoff/grid_spac[dim]);
    double cut_sq=lmf_cutoff*lmf_cutoff;
    double sigma_inv_sq = 1.0/lmf_sigma/lmf_sigma;

    int pi[3]; double pf = 1.0/pow(sqrt(M_PI),3.0)*sigma_inv_sq/lmf_sigma;
    double offset[3];
    for (int i = 0; i < nall; i++)
        if (mask[i] & groupbit) {
            pos_index_calc(i,pi,x, offset, 1);
            smooth_charge->exp_smooth(pi, n_steps, cut_sq, pf*q[i], sigma_inv_sq, grid_spac, offset);
        }
    
    if(ensemble_flag){
        int num_procs;
        MPI_Comm_size (universe->uworld, &num_procs);
        if(debug)printf("uworld %d \n", num_procs);
        MPI_Comm_size (universe->uorig, &num_procs);
        if(debug)printf("uorig %d \n", num_procs);
        smooth_charge->ens_average(&universe->uorig);
        LMF_potential->ens_average(&universe->uorig);
        if(debug) printf("sc 111 %f sc_g 111 %f \n", smooth_charge->array[1][1][1], smooth_charge->array_g[1][1][1]);
    }

}


void Fix_LMF::compute_yb_correction()
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  //double zprd = domain->zprd;
  int nlocal = atom->nlocal;

  double dipole = 0.0; double qsum=0.0;
  for (int i = 0; i < nlocal; i++) {
    dipole += q[i]*x[i][2];
    qsum+=q[i];
  }
  if(fabs(qsum)>SMALL)printf("System may not be neutral, qsum = %f\n",qsum);

  // sum local contributions to get global dipole moment

  double dipole_all;
  if(ensemble_flag) MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,universe->uorig);
  else MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);
  
  /*
  // need to make non-neutral systems and/or
  //  per-atom energy translationally invariant

  double dipole_r2 = 0.0;
  if (eflag_atom || fabs(qsum) > SMALL) {
    for (int i = 0; i < nlocal; i++)
      dipole_r2 += q[i]*x[i][2]*x[i][2];

    // sum local contributions

    double tmp;
    MPI_Allreduce(&dipole_r2,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    dipole_r2 = tmp;
  }
  */
  // compute corrections
  double pi_temp = 3.14159265358979323846;
  double ffact = force->qqrd2e * (-4.0*pi_temp/lmf_vol);
  yb_correction=ffact*dipole_all;
  //if(verbose)printf("dipole %f dipole_all %f yb_correction %f \n",dipole, dipole_all,yb_correction);
}


void Fix_LMF::pos_index_calc(int i, int * pi, double ** x, double * disp, bool dbool){
    for(int dim=0; dim<3; dim++){
        pi[dim] = int((x[i][dim]-lmf_lo[dim])/grid_spac[dim]);
        if(dbool) disp[dim]=x[i][dim]-pi[dim]*grid_spac[dim]-lmf_lo[dim];
    }
}

int Fix_LMF::setmask(){
  int mask = 0;

  // FLD implicit needs to invoke wall forces before pair style

  if (fldflag) mask |= PRE_FORCE;
  else mask |= POST_FORCE;

  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

void Fix_LMF::setup(int vflag){
 post_force(vflag);
}

void Fix_LMF::pre_force(int vflag){
  post_force(vflag);
}

void Fix_LMF::post_force(int vflag){

  LMF_force();
}

double Fix_LMF::SOR_iterate(LMF_GRID * lmf_pot, LMF_GRID *sc, double pf, int * nt) {
    double t_change=0, sor_change=0;
    double gis_tot=grid_inv_sq[0]+grid_inv_sq[1]+grid_inv_sq[2];
    double pf_sor=omega/gis_tot/2.0;
    int hi_ind[3]; double hi_val=0;
    for(int n_0=1; n_0<nt[0]-1; n_0++){
		for(int n_1=1; n_1<nt[1]-1; n_1++){
			for(int n_2=1; n_2<nt[2]-1; n_2++){
			    int hi0 = nt[0]-2; int hi1=nt[1]-2;
			    /*if(n_0==1 && n_1==1 && n_2==129){
				    printf("sc %f ind %d %d %d hi2 %d \n",pf*sc->array[20][n_1][n_2], n_0,n_1,n_2,nt[2]);
			        printf("pot %f lap %f sc %f sor %f \n",lmf_pot->array[n_0][n_1][n_2], lmf_pot->lap_calc(n_0, n_1, n_2), pf*sc->array[n_0][n_1][n_2], sor_change*sor_change);
			        printf("ne0 %f %f ne1 %f %f ne2 %f %f\n",lmf_pot->array[n_0-1][n_1][n_2], lmf_pot->array[n_0+1][n_1][n_2], lmf_pot->array[n_0][n_1-1][n_2], lmf_pot->array[n_0][n_1+1][n_2], lmf_pot->array[n_0][n_1][n_2-1], lmf_pot->array[n_0][n_1][n_2+1]);
			        printf("ne0 %f %f ne1 %f %f ne2 %f %f\n",pf*sc->array[10][10][129], pf*sc->array[10][20][130], pf*sc->array[15][40][128], pf*sc->array[1][12][129], pf*sc->array[1][18][131], pf*sc->array[20][1][130]);
			    }*/
				sor_change=lmf_pot->lap_calc(n_0, n_1, n_2)+pf*sc->array[n_0][n_1][n_2];
				lmf_pot->array[n_0][n_1][n_2]+=pf_sor*sor_change;
				t_change=t_change+sor_change*sor_change;
				update_boundary(n_0,n_1,n_2,lmf_pot->array, nt);
				/*if(sor_change*sor_change>hi_val){
				    hi_val=sor_change*sor_change;
				    hi_ind[0]=n_0; hi_ind[1]=n_1; hi_ind[2]=n_2;
				}*/
			}
		}
	}
	subtract_av(lmf_pot);
	return (sqrt(t_change/float((nt[0]-2)*(nt[1]-2)*(nt[2]-2))));
}

void Fix_LMF::subtract_av(LMF_GRID * grid_in) {
    double tot_temp=0.0;
    int * nt = grid_in->nt;
    for(int n_0=1; n_0<nt[0]-1; n_0++){
		for(int n_1=1; n_1<nt[1]-1; n_1++){
			for(int n_2=1; n_2<nt[2]-1; n_2++){
			    int hi0 = nt[0]-2; int hi1=nt[1]-2;
			        tot_temp+=grid_in->array[n_0][n_1][n_2];
			}
		}
	}
    tot_temp=tot_temp/double((nt[0]-2)*(nt[1]-2)*(nt[2]-2));
	for(int n_0=0; n_0<nt[0]; n_0++)
		for(int n_1=0; n_1<nt[1]; n_1++)
			for(int n_2=0; n_2<nt[2]; n_2++)
			    grid_in->array[n_0][n_1][n_2]-=tot_temp;
}

void Fix_LMF::update_boundary(int n_0, int n_1, int n_2, double *** array, int * nt){
    int opp_0, opp_1, opp_2;
    double base = array[n_0][n_1][n_2];
    opp_0=return_opposite(n_0,nt[0]);
    opp_1=return_opposite(n_1,nt[1]);
    opp_2=return_opposite(n_2,nt[2]);
    if(opp_2!=n_2){                
        array[n_0][n_1][opp_2]=base;
        array[n_0][opp_1][opp_2]=base;
        array[opp_0][n_1][opp_2]=base;
        array[opp_0][opp_1][opp_2]=base;
    }
    if(opp_1!=n_1){                
        array[n_0][opp_1][n_2]=base;
        array[opp_0][opp_1][n_2]=base;
    }
    if(opp_0!=n_0){                
        array[opp_0][n_1][n_2]=base;
    }
}

int Fix_LMF::return_opposite(int ind,int high){
    if(ind==1) return (high-1);
    else if(ind==(high-2)) return 0;
    else return ind;
}







