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
#include "fix_analysis.h"
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
#include "pair.h"
#include "neighbor.h"
#include "change_box.h"
#include <complex>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{XLO=0,XHI=1,YLO=2,YHI=3,ZLO=4,ZHI=5};
enum{NONE=0,EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */


analysis::analysis(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg){
  
  verbose=0; 
  fldflag=0;
  wdyn=0;
  grflag=0;
  autflag=0;
  
  gr_oo=gr_oh=gr_hh=gr_oo_u=gr_oh_u=gr_hh_u=NULL;

	aut_count=0; gr_count=0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"gr") == 0){
      if (iarg+5 > narg ) error->all(FLERR,"Illegal fix analysis command");
	  	grflag=1;
	  	gr_file_name = arg[iarg+1];
	  	gr_freq  = force->numeric(FLERR,arg[iarg+2]);
	  	gr_spac = force->numeric(FLERR,arg[iarg+3]);
	  	gr_max = force->numeric(FLERR,arg[iarg+4]);
	  	iarg+=5;
		}
	  else if (strcmp(arg[iarg],"verbose") == 0){
	  	if (iarg+1 > narg ) error->all(FLERR,"Illegal fix analysis command");
	  	verbose=1;
	  	iarg+=1;
		}
	  else if (strcmp(arg[iarg],"aut") == 0){
        if (iarg+5 > narg ) error->all(FLERR,"Illegal fix analysis command");
        aut_file_name = arg[iarg+1];
        aut_write_freq  = force->numeric(FLERR,arg[iarg+2]);
        n_aut = force->numeric(FLERR,arg[iarg+3]);
        num_lags = int(force->numeric(FLERR,arg[iarg+4]));
        autflag=1;
        iarg+=5;
	  }
    else {
  		printf("arg %s \n", arg[iarg]);
  		error->all(FLERR,"Illegal fix analysis command");
  	}
  }
  
  if(grflag){
    gr_tot=int(gr_max/gr_spac);
    gr_oo=create_1d_array(gr_tot);
    gr_oh=create_1d_array(gr_tot);
    gr_hh=create_1d_array(gr_tot);
    gr_oo_u=create_1d_array(gr_tot);
    gr_oh_u=create_1d_array(gr_tot);
    gr_hh_u=create_1d_array(gr_tot);
    if(universe->me==0) {
    	printf("gr_tot=%d gr_spac=%f\n", gr_tot, gr_spac);
    	gr_file = fopen(gr_file_name,"w");
    }
  }
  
  if(autflag){
  	charge_autvals = create_2d_array(num_lags,n_aut);
  	charge_autocorrelation = create_2d_array(num_lags,n_aut);
  	if(universe->me==0) {
    	printf("n_aut=%d num_lags=%d \n", n_aut, num_lags);
    	aut_file = fopen(aut_file_name,"w");
    }
  }

}

analysis::~analysis(){
  if(grflag){
  	delete [] gr_oo;
  	delete [] gr_oh;
  	delete [] gr_hh;
  	delete [] gr_oo_u;
  	delete [] gr_oh_u;
  	delete [] gr_hh_u;
  }
  if(autflag){
  	delete [] charge_autvals;
  	delete [] charge_autocorrelation;
  }
}

double * analysis::create_1d_array(int size){
	double * array = new double [size];
	for(int n=0; n<size; n++)
		array[n]=0.0;
	return array;
}

std::complex<double> ** analysis::create_2d_array(int size1, int size2){
	std::complex<double> ** array = new std::complex<double> * [size1];
	for(int n=0; n<size1; n++){
		array[n]=new std::complex<double> [size2];
		for(int m=0; m<size2; m++)
			array[n][m]=0.0;
	}
	return array;
}

void analysis::zero_1d_array(double * array, int size){
	for(int n=0; n<size; n++)
		array[n]=0.0;
}

void analysis::zero_2d_array(std::complex<double> ** array, int size1, int size2){
	for(int n=0; n<size1; n++)
		for(int m=0; m<size2; m++)
			array[n][m]=0.0;
}


void analysis::calculate_gr(){
  int itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  
  double rsq, dist;
  
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal+atom->nghost;
  
  zero_1d_array(gr_oo,gr_tot);
  zero_1d_array(gr_oh,gr_tot);
  zero_1d_array(gr_hh,gr_tot);
  
  int o_count=0;
  int h_count=0;

  for(int i=0; i<nlocal; i++){
    xtmp=x[i][0];
    ytmp=x[i][1];
    ztmp=x[i][2];
    itype=type[i];
    if(itype==1) o_count++;
    if(itype==2) h_count++;
  	for(int j=0; j<nall; j++){
  		delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
  		jtype = type[j];
      dist = sqrt(rsq);
      int ind = round(dist/gr_spac);
			if (ind<gr_tot){
				if(itype==1){
					if(jtype==1) gr_oo[ind]+=1.0;
					else if(jtype==2) gr_oh[ind]+=1.0; 
					else printf("Types not correct for g(r)\n");
				}
				else if(itype==2){
					if(jtype==2) gr_hh[ind]+=1.0; 
					else if(jtype!=1) printf("Types not correct for g(r)\n");
				}
				else printf("Types not correct for g(r)\n");
			}
		}
	}
 
		
	MPI_Allreduce(gr_oo,gr_oo_u,gr_tot,MPI_DOUBLE,MPI_SUM,universe->uorig);
	MPI_Allreduce(gr_oh,gr_oh_u,gr_tot,MPI_DOUBLE,MPI_SUM,universe->uorig);
	MPI_Allreduce(gr_hh,gr_hh_u,gr_tot,MPI_DOUBLE,MPI_SUM,universe->uorig);

	for(int i=1; i<gr_tot; i++){ 
			double factor = 1.0/4.0/float(i*i)/gr_spac/gr_spac/gr_spac/3.14159265359/.03342703/float(universe->nworlds);
			gr_oo[i]=factor*gr_oo_u[i]/float(o_count);
			gr_oh[i]=factor*gr_oh_u[i]/float(o_count);
			gr_hh[i]=factor*gr_hh_u[i]/float(h_count);
	}
}

void analysis::write_gr(){
  
  if(universe->me==0){
    fprintf(gr_file,"Time %d %d %f %d \n", update->ntimestep, gr_tot, gr_spac, gr_freq);
  	for(int i=1; i<gr_tot; i++)
  		fprintf(gr_file, "%d %f %f %f %f \n", i, i*gr_spac, gr_oo[i], gr_oh[i], gr_hh[i]);
  	fflush(gr_file);
  	}
}

void analysis::write_aut(){
  if(universe->me==0){
    fprintf(aut_file,"Time %d %d %d %f %d \n", update->ntimestep, n_aut, aut_write_freq, domain->boxhi[0]-domain->boxlo[0], num_lags);
  	for(int f_dim=0; f_dim<num_lags-1; f_dim++)
    	for(int n_0=0; n_0<n_aut; n_0++)
	  		fprintf(aut_file, "%d %d %f %f \n", f_dim, n_0, std::abs(charge_autocorrelation[f_dim][n_0]), std::arg(charge_autocorrelation[f_dim][n_0]));
	fflush(aut_file);
	}
}


void analysis::update_aut_vals(){
	double **x = atom->x;
	double *q = atom->q;
	int *mask = atom->mask;
	for(int n_0=0; n_0<n_aut; n_0++){
		charge_autvals[num_lags-1][n_0]=0.0;
	}
	double boxlo = domain->boxlo[0];
	double boxlen = domain->boxhi[0]-boxlo;
	for (int i = 0; i < atom->nlocal; i++)
        if (mask[i] & groupbit) {
          for(int n_t=0; n_t<n_aut; n_t++){
                double pos_val = float(n_t)*((x[i][0]-boxlo)/boxlen);
                charge_autvals[num_lags-1][n_t]=charge_autvals[num_lags-1][n_t]+std::exp(std::complex<double>(2.0*3.1415*1.0i*pos_val))*std::complex<double>(q[i]);
                }
        }
}

void analysis::update_autocorr(){
  int num_temp=aut_count-num_lags;
	for(int f_dim=0; f_dim<num_lags; f_dim++){
    	for(int n_0=0; n_0<n_aut; n_0++){
			charge_autocorrelation[f_dim][n_0]=charge_autvals[num_lags-1][n_0]*std::conj(charge_autvals[num_lags-1-f_dim][n_0])/std::complex<double>(num_temp)+std::complex<double>(float(num_temp-1)/float(num_temp))*charge_autocorrelation[f_dim][n_0];
        }
	}

}

void analysis::lag_back(){
    for(int f_dim=0; f_dim<num_lags-1; f_dim++){
    	for(int n_0=0; n_0<n_aut; n_0++){
		  	charge_autvals[f_dim][n_0]=charge_autvals[f_dim+1][n_0];
		}
	}

}

void analysis::aut_update()
{
  aut_count++;
  if(autflag==1){
      lag_back();
      update_aut_vals();
      if(aut_count>num_lags) update_autocorr();
  }
}

int analysis::setmask(){
  int mask = 0;

  // FLD implicit needs to invoke wall forces before pair style

  if (fldflag) mask |= PRE_FORCE;
  else mask |= POST_FORCE;

  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

int analysis::modify_param(int narg, char **arg){
  if (strcmp(arg[0],"gr") == 0){
	  	if (narg<2 ) error->all(FLERR,"Illegal fix_modify LMF command");
    	grflag=1;
    	return 2;
  	}
  return 0;
}


/* ---------------------------------------------------------------------- */

void analysis::setup(int vflag){
 post_force(vflag);
}

void analysis::pre_force(int vflag){
  post_force(vflag);
}

void analysis::post_force(int vflag){
  if(grflag && gr_count>gr_freq){
  	 calculate_gr();
  	 write_gr();
  	 gr_count=0;
  }
  if(autflag){
  	aut_update();
  	if(aut_count>aut_write_freq){
  		write_aut();
  		zero_2d_array(charge_autvals,num_lags,n_aut);
  		zero_2d_array(charge_autocorrelation,num_lags,n_aut);
  		aut_count=0;
  	}
  }
}


