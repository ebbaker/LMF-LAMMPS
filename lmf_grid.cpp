/* ----------------------------------------------------------------------
Poisson solver subclass
------------------------------------------------------------------------- */
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "variable.h"
#include "lmf_grid.h"
#include <complex>
#include "math.h"
#include "stdlib.h"
#include "string.h"
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
#include "universe.h"
#include <complex>
#include "run.h"

using namespace LAMMPS_NS;

LMF_GRID::LMF_GRID(int * nt_temp, double * gs_temp, double * lo_temp, bool ens_temp) {
    ens=ens_temp; nt=nt_temp; gs=gs_temp; lo=lo_temp; 
    gis=new double [3];
    for(int dim=0; dim<3; dim++) gis[dim]=1.0/gs[dim]/gs[dim];
    array=create_array(nt);
    if(ens) array_g=create_array(nt);
    else array_g=NULL;
}


LMF_GRID::~LMF_GRID(){
    delete [] array;
    delete [] array_g;
}



void LMF_GRID::gradient_calc(int * pi, double * gr){
    double base, ltemp=0;
    double base_potential = array[pi[0]][pi[1]][pi[2]];
    gr[0]=(array[pi[0]+1][pi[1]][pi[2]]-base_potential)/gs[0];
	gr[1]=(array[pi[0]][pi[1]+1][pi[2]]-base_potential)/gs[1];
	gr[2]=(array[pi[0]][pi[1]][pi[2]+1]-base_potential)/gs[2];
}

void LMF_GRID::exp_smooth(int *ind, int *n_steps, double cut_sq, double pf, double sigma_inv_sq, double * grid_spac, double * offset){
    double delta[3], delt_sq[3], r_sq;
    for(int n_0 = -n_steps[0]; n_0 <= n_steps[0]; n_0++){
        if(n_0+ind[0]<0 || n_0+ind[0]>=nt[0]) continue;
        delta[0]=offset[0]-n_0*grid_spac[0];
        delt_sq[0]=delta[0]*delta[0];
        for(int n_1 = -n_steps[1]; n_1 <= n_steps[1];n_1++){
            if(n_1+ind[1]<0 || n_1+ind[1]>=nt[1]) continue;
            delta[1]=offset[1]-n_1*grid_spac[1];
            delt_sq[1]=delta[1]*delta[1];
            for(int n_2 = -n_steps[2]; n_2 <= n_steps[2]; n_2++){
                if(n_2+ind[2]<0 || n_2+ind[2]>=nt[2]) continue;
                delta[2]=offset[2]-n_2*grid_spac[2];
                delt_sq[2]=delta[2]*delta[2];
                r_sq=delt_sq[0]+delt_sq[1]+delt_sq[2];
                if(r_sq<cut_sq){
                    double charge_temp=pf*exp(-r_sq*sigma_inv_sq);
                    array[n_0+ind[0]][n_1+ind[1]][n_2+ind[2]]+=charge_temp;
                    }
            }
        }
    } 	
}

void LMF_GRID::ens_average(MPI_Comm *uorig){
    if(ens){
        int num_procs;
        double temp_array[nt[2]];
        MPI_Comm_size (*uorig, &num_procs);
        for(int n_0=0; n_0<nt[0]; n_0++){
            for(int n_1=0; n_1<nt[1]; n_1++){
                MPI_Allreduce(array[n_0][n_1],temp_array,nt[2],MPI_DOUBLE,MPI_SUM,*uorig);
                for(int n_2=0; n_2<nt[2]; n_2++){
                    array[n_0][n_1][n_2]=temp_array[n_2]/double(num_procs);
                }
            }
        }
    }
}

double LMF_GRID::lap_calc(int n_0, int n_1, int n_2){
    double temp=0.0;
    double baset = 2.0*array[n_0][n_1][n_2];
	temp += (array[n_0+1][n_1][n_2]+array[n_0-1][n_1][n_2]-baset)*gis[0];
	temp += (array[n_0][n_1+1][n_2]+array[n_0][n_1-1][n_2]-baset)*gis[1];
	temp += (array[n_0][n_1][n_2+1]+array[n_0][n_1][n_2-1]-baset)*gis[2];
	return temp;
}

void LMF_GRID::zero(){
    for(int n_0=0; n_0<nt[0]; n_0++){
        for(int n_1=0; n_1<nt[1]; n_1++){
            for(int n_2=0; n_2<nt[2]; n_2++){
                array[n_0][n_1][n_2]=0;
                if (ens) array_g[n_0][n_1][n_2]=0;
            }	
        }		
    }	
}

void LMF_GRID::subtract_av(double total_field){
	for(int n_0=0; n_0<nt[0]; n_0++){
		for(int n_1=0; n_1<nt[1]; n_1++){
			for(int n_2=0; n_2<nt[2]; n_2++){
				array[n_0][n_1][n_2]-=total_field;
			}	
		}		
	}		
}

double *** LMF_GRID::create_array(int * size){
	double *** ar = new double ** [size[0]];
	for(int n0=0; n0<size[0]; n0++){
		ar[n0]=new double * [size[1]];
		for(int n1=0; n1<size[1]; n1++){
			ar[n0][n1]= new double [size[2]];
			for(int n2=0; n2<size[2]; n2++)
				ar[n0][n1][n2]=0.0;
		}
	}
	return ar;
}