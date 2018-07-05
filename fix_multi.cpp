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
#include "fix_multi.h"
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
#include "multi_grid.h"
#include <complex>
#include "run.h"
#include <random>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{XLO=0,XHI=1,YLO=2,YHI=3,ZLO=4,ZHI=5};
enum{NONE=0,EDGE,CONSTANT,VARIABLE};

#define SMALL 0.00001

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


Fix_MULTI::Fix_MULTI(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg){
    
    lmf_sigma=0.0;

    ensemble_flag=0; no_force=0; root_bool=0; slab_flag=0; slab_write=0;
    fldflag=0; verbose=0; debug=0; force_ew_write=0; rad_write=0;
    output_write_freq=100000; error_flag=0;
    
    
    double omega_temp, nr_temp, it_temp;
    
    force_flag=0;

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
        else if (strcmp(arg[iarg],"grid") == 0){
            if (iarg+5 > narg ) error->all(FLERR,"Illegal fix LMF command");
            lmf_sigma = force->numeric(FLERR,arg[iarg+1])/sqrt(2);
            gs_ratio = force->numeric(FLERR,arg[iarg+2]);
            pow_base = force->numeric(FLERR,arg[iarg+3]);
            num_grids = int(force->numeric(FLERR,arg[iarg+4]));
            iarg+=5;
        }  
        else if (strcmp(arg[iarg],"multi")==0){
            if (iarg+5 > narg ) error->all(FLERR,"Illegal fix LMF command");
            ALL_flag=1;
            sigma_in=force->numeric(FLERR,arg[iarg+1])/sqrt(2.0);
            num_rand=int(force->numeric(FLERR,arg[iarg+2]));
            omega_temp=force->numeric(FLERR,arg[iarg+3]);
            it_temp=force->numeric(FLERR,arg[iarg+4]);
            iarg+=5;
        }
        else if (strcmp(arg[iarg],"slab")==0){
            if (iarg+2 > narg ) error->all(FLERR,"Illegal fix LMF command");
            slab_flag=1;
            slab_vol_factor=force->numeric(FLERR,arg[iarg+1]);
            iarg+=2;
        }
        else if (strcmp(arg[iarg],"slab_write")==0){
            if (iarg+4 > narg ) error->all(FLERR,"Illegal fix LMF command");
            slab_flag=1;
            slab_write=1;
            slab_vol_factor=force->numeric(FLERR,arg[iarg+1]);
            output_file_name = arg[iarg+2];
            output_write_freq=int(force->numeric(FLERR,arg[iarg+3]));
            iarg+=4;
        }
        else if (strcmp(arg[iarg],"rad_write")==0){
            if (iarg+5 > narg ) error->all(FLERR,"Illegal fix LMF command");
            rad_write=1;
            max_rad_write=force->numeric(FLERR,arg[iarg+1]);
            rad_write_spacing=force->numeric(FLERR,arg[iarg+2]);
            output_file_name = arg[iarg+3];
            output_write_freq=int(force->numeric(FLERR,arg[iarg+4]));
            iarg+=5;
        }
        else if (strcmp(arg[iarg],"noforce")==0){
            if (iarg+1 > narg ) error->all(FLERR,"Illegal fix LMF command");
            no_force=1;
            iarg+=1;
        }
        else if (strcmp(arg[iarg],"error_calc")==0){
            if (iarg+1 > narg ) error->all(FLERR,"Illegal fix LMF command");
            error_flag=1;
            iarg+=1;
        }
        else if (strcmp(arg[iarg],"force_write") == 0){
            if (iarg+3 > narg ) error->all(FLERR,"Illegal fix analysis command");
            force_write_freq=int(force->numeric(FLERR,arg[iarg+1]));
            force_file_name = arg[iarg+2];
            force_flag=1;
            int *type=atom->type;
            int nlocal = atom->nlocal;
            for (int i = 0; i < nlocal; i++)
                if (type[i]==1){
                    tag_particle=i;
                    break;
                }
            iarg+=3;
        }
        else {
            printf("arg %s \n", arg[iarg]);
            error->all(FLERR,"Illegal fix LMF command");
        }
    }
    if(slab_write || rad_write){
        if(universe->me==0) {
    	    output_file = fopen(output_file_name,"w");
        }
    } 
    if(force_flag){
        if(universe->me==0) {
    	    force_file = fopen(force_file_name,"w");
        }
    } 
    
    if(ALL_flag){
        omega = new double [num_grids];
        it_freq = new int [num_grids];
        for(int i=0; i<num_grids; i++){
            omega[i]=omega_temp;
            it_freq[i]=it_temp;
        }
    }
    else error->all(FLERR,"LMF multi params not specified");
    setup_multi_lattice();
    if(lmf_sigma==0.0) error->all(FLERR,"Did not specify sigma in Fix LMF");	
    root_bool=0;
    for(int world=0; world<universe->nworlds; world++){
        if(universe->me==universe->root_proc[world]) 
            {
            root_bool=1;
            printf("root %d \n", universe->me);
            }
    }
    
    lmf_pf = force->qqr2e/force->qe2f*4.0*M_PI;
    printf("lmf_pf %f \n", lmf_pf);
    gen.seed(time(NULL));
    
}

Fix_MULTI::~Fix_MULTI(){
	
    for(int gn=0; gn<num_grids; gn++){
        delete potential[gn];
	    delete sc[gn];
        delete [] n_tot[gn];
	    delete [] grid_spac[gn];
	    delete [] grid_lo[gn];
	}
	delete [] sc;
	delete [] potential;
	delete [] n_tot;
	delete [] grid_spac;
	delete [] grid_lo;
}

void Fix_MULTI::setup_multi_lattice()
{
    
    n_tot = new int * [num_grids];
    grid_spac = new double * [num_grids];
    grid_lo = new double * [num_grids];
    
    potential = new MULTI_GRID * [num_grids];
    sc = new MULTI_GRID * [num_grids];
    hpf = new double [num_grids];
    
    for(int gn=0; gn<num_grids; gn++){
        n_tot[gn] = new int [3];
        grid_spac[gn] = new double [3];
        grid_lo[gn] = new double [3];
        setup_lmf_lattice(gn);
    }
    
    double * gst = grid_spac[0];
    sc_pf=1.0/double(num_rand)/gst[0]/gst[1]/gst[2];
}

void Fix_MULTI::setup_lmf_lattice(int grid_num)
{
    double * lmf_lo = grid_lo[grid_num];
    double * gs_temp = grid_spac[grid_num];
    int * nt_temp = n_tot[grid_num];
    
    double sigma = pow(pow_base,grid_num)*lmf_sigma;
    double gis_tot=0.0;
    
    for(int dim_temp=0; dim_temp<3; dim_temp++){
        gs_temp[dim_temp]=pow(pow_base,grid_num)*lmf_sigma*sqrt(2)/gs_ratio;
        gis_tot+=1.0/gs_temp[dim_temp]/gs_temp[dim_temp];
        
        lmf_lo[dim_temp]=domain->boxlo[dim_temp];
        double size_temp=domain->boxhi[dim_temp]-lmf_lo[dim_temp];
        if(slab_flag && dim_temp==2){
            lmf_lo[dim_temp]=lmf_lo[dim_temp]-(slab_vol_factor-1.0)*size_temp/2.0;
            size_temp=slab_vol_factor*size_temp;
        }
        if (gs_temp[dim_temp]>size_temp/2.0) error->all(FLERR,"Grid spacing too large");	
        nt_temp[dim_temp] = int(ceil(size_temp/gs_temp[dim_temp]))+2;
		gs_temp[dim_temp] = size_temp/double(nt_temp[dim_temp]-2);	
		printf("grid %d dim %d bounds %f %f gs %f nt %d \n", grid_num, dim_temp, lmf_lo[dim_temp], (nt_temp[dim_temp]-1)*gs_temp[dim_temp]+lmf_lo[dim_temp], gs_temp[dim_temp], nt_temp[dim_temp]);
	}
	
	hpf[grid_num] = omega[grid_num]/gis_tot/2.0;
	printf("hpf %f sigma %f \n",hpf[grid_num], sigma);
	
	bool flags[2];
	double pars[2];

    flags[0]=ensemble_flag; flags[1]=slab_flag; 
	
	sc[grid_num] = new MULTI_GRID(nt_temp, gs_temp, lmf_lo, sigma, flags);
	potential[grid_num] = new MULTI_GRID(nt_temp, gs_temp, lmf_lo, sigma, flags);

}

int Fix_MULTI::modify_param(int narg, char **arg)
{
    if (strcmp(arg[0],"grid") == 0){
        if (narg<5 ) error->all(FLERR,"Illegal fix_modify LMF command");
        lmf_sigma = force->numeric(FLERR,arg[1])/sqrt(2.0);
        gs_ratio = force->numeric(FLERR,arg[2]);
        pow_base = force->numeric(FLERR,arg[3]);
        num_grids = int(force->numeric(FLERR,arg[4]));
        return 5;
    }
    
    else if (strcmp(arg[0],"multi") == 0){
        if (narg<5 ) error->all(FLERR,"Illegal fix_modify LMF command");
        sigma_in = force->numeric(FLERR,arg[1])/sqrt(2.0);
        num_rand=int(force->numeric(FLERR,arg[2]));
        double omega_temp=force->numeric(FLERR,arg[3]);
        double it_temp=force->numeric(FLERR,arg[4]);
        for(int i=0; i<num_grids; i++){
            omega[i]=omega_temp;
            it_freq[i]=it_temp;
        }
        return 5;
    }
    else if (strcmp(arg[0],"write_freq") == 0){
        if (narg<2 ) error->all(FLERR,"Illegal fix_modify LMF command");
        output_write_freq=int(force->numeric(FLERR,arg[1]));
        return 2;
    }
    else if (strcmp(arg[0],"force_write_freq") == 0){
        if (narg<2 ) error->all(FLERR,"Illegal fix_modify LMF command");
        force_write_freq=int(force->numeric(FLERR,arg[1]));
        return 2;
    }
    else if (strcmp(arg[0],"omega_full") == 0){
        if (narg<num_grids+1 ) error->all(FLERR,"Illegal fix_modify LMF command, must specify omega for all grids");
        for(int ng=0; ng<num_grids; ng++){
            double gis_tot=0.0;
            for(int dt=0; dt<3; dt++) gis_tot+=1.0/grid_spac[ng][dt]/grid_spac[ng][dt];
            omega[ng]=force->numeric(FLERR,arg[ng+1]);
            hpf[ng]=omega[ng]/gis_tot/2.0;
        }
        return (num_grids+1);
    }
    else if (strcmp(arg[0],"it_full") == 0){
        if (narg<num_grids+1 ) error->all(FLERR,"Illegal fix_modify LMF command, must specify it for all grids");
        for(int ng=0; ng<num_grids; ng++){
            it_freq[ng]=force->numeric(FLERR,arg[ng+1]);
        }
        return (num_grids+1);
    }
    return 0;
}

void Fix_MULTI::MULTI_update(){
    double st;
    
    bigint ntimestep=update->ntimestep;
    
    sc[0]->zero();
    rand_add();
    st = pow(lmf_sigma,2.0)-pow(sigma_in,2.0);
    sc[0]->spread_charge(st);
    
    for(int gn=0; gn<num_grids-1; gn++){
        sc[gn+1]->zero();
        sc[gn]->add_sc_up(sc[gn+1]);
        st=pow(sc[gn+1]->sigma,2.0)-pow(sc[gn]->sigma,2.0);
        sc[gn+1]->spread_charge(st);
        sc[gn]->add_sc_down(sc[gn+1]);
    }
    for(int gn=0; gn<num_grids; gn++){  
        if(ensemble_flag) sc[gn]->ens_average(&universe->uorig); 
        for (int n_steps=0; n_steps<it_freq[gn]; n_steps++)
            potential[gn]->heat_iterate(sc[gn]->array,hpf[gn],lmf_pf,error_flag);
    }
    
    if(ntimestep%output_write_freq==0 && universe->me==0){
        for(int gn=0; gn<num_grids; gn++){
            if(slab_write) sc[gn]->write_slab(potential[gn]->array,output_file,ntimestep,gn);
            if(rad_write) sc[gn]->write_rad(potential[gn]->array,output_file,ntimestep,gn,max_rad_write,rad_write_spacing*pow(pow_base,gn));
        }
    }
   
    
}

void Fix_MULTI::MULTI_force(){
    int pos_index[3];

    double **x = atom->x;
    double **f = atom->f;
    double *q = atom->q;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    bigint ntimestep=update->ntimestep;

    double gradient[3]; double total_gradient[3]; int dim;

    if(!no_force){
        for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) {
                double pre_factor = q[i]*force->qe2f; //convert to kcal/mol
                for(dim=0; dim<3; dim++) total_gradient[dim]=0.0;
                for(int gn=0; gn<num_grids; gn++){
                    potential[gn]->gradient_calc(x[i], gradient);
                    for(dim=0; dim<3; dim++) total_gradient[dim]+=gradient[dim];
                }
                for(dim=0; dim<3; dim++) f[i][dim]-=pre_factor*total_gradient[dim];
            } 
    }     
    if(force_flag){
        for (int i = 0; i < nlocal; i++){
            if (mask[i] & groupbit){
                if(universe->me==0 && i==tag_particle && update->ntimestep%force_write_freq==0){
                    double pre_factor = q[i]*force->qe2f; //convert to kcal/mol
                    fprintf(force_file,"%d", update->ntimestep);
                    if(force_ew_write) fprintf(force_file," %d %f %f %f",-1,force_ew[0],force_ew[1], force_ew[2]);
                    for(int gn=0; gn<num_grids; gn++){
                        potential[gn]->gradient_calc(x[i], gradient);
                        fprintf(force_file," %d %f %f %f",gn,-pre_factor*gradient[0], -pre_factor*gradient[1], -pre_factor*gradient[2]);
                    }
                    fprintf(force_file,"\n");
                }
            }
        }        
        fflush(force_file);
    }   
             
}

void Fix_MULTI::rand_add(){
    double **x = atom->x;
    int *mask = atom->mask;
    double *q = atom->q;
    int nlocal = atom->nlocal;
    MULTI_GRID * sct = sc[0];
    double pos[3]; double q_add; int ind[3]; int dim;
    std::normal_distribution<double> dist(0.0,sigma_in);
    double *** array = sc[0]->array;
    
    for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
            for(int n=0; n<num_rand; n++){
                q_add=q[i]*sc_pf;
                for(dim=0; dim<3; dim++)
                    pos[dim]=x[i][dim]+dist(gen);
                if(sct->ind_calc(ind, pos)){
                    array[ind[0]][ind[1]][ind[2]] += q_add;
                    sc[0]->update_boundary(ind[0],ind[1],ind[2]);
                }
            }
        } 
}






int Fix_MULTI::setmask(){
  int mask = 0;

  // FLD implicit needs to invoke wall forces before pair style

  if (fldflag) mask |= PRE_FORCE;
  else mask |= POST_FORCE;

  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

void Fix_MULTI::setup(int vflag){
 post_force(vflag);
}

void Fix_MULTI::pre_force(int vflag){
  post_force(vflag);
}

void Fix_MULTI::post_force(int vflag){
  MULTI_force();
}










