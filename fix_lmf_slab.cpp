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
#include "fix_lmf_slab.h"
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
#include "change_box.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{XLO=0,XHI=1,YLO=2,YHI=3,ZLO=4,ZHI=5};
enum{NONE=0,EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */


Fix_LMF_slab::Fix_LMF_slab(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg){
  lmf_sigma=0;
  lmf_slab_write=0;
  
  ensemble_flag=0; dens_flag=0; no_force=0; de_flag=0; root_bool=0;
  E_Field_Variable=0; verbose=0; verbose_full=0;
  fldflag=0;
  wdyn=0;
  check_charge=0;

  if(!domain->xperiodic) printf("X dimension must be periodic for LMF slab geometry");
  if(!domain->yperiodic) printf("Y dimension must be periodic for LMF slab geometry");
  if(domain->zperiodic) printf("Z dimension must be fixed for LMF slab geometry");

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
	else if (strcmp(arg[iarg],"verbose_full") == 0){
	  	verbose=1;
	  	verbose_full=1;
	  	iarg+=1;
		}
	else if (strcmp(arg[iarg],"par") == 0){
	  	if (iarg+4 > narg ) error->all(FLERR,"Illegal fix LMF command");
    	lmf_sigma = force->numeric(FLERR,arg[iarg+1]);
    	lmf_cutoff = force->numeric(FLERR,arg[iarg+2]); 
    	bound_len_d = force->numeric(FLERR,arg[iarg+3]); 
    	iarg+=4;
  	}   
  	else if (strcmp(arg[iarg],"gs")==0){
	  	if (iarg+2 > narg ) error->all(FLERR,"Illegal fix LMF command");
	  	grid_spac=force->numeric(FLERR,arg[iarg+1]);
	  	iarg+=2;
  	}
    else if (strcmp(arg[iarg],"write")==0){
        if (iarg+3 > narg ) error->all(FLERR,"Illegal fix LMF command");
        lmf_slab_write=force->numeric(FLERR,arg[iarg+1]);
        char * temp_name = arg[iarg+2];
        lmf_slab_file_name = (char *) malloc(sizeof(char)*(strlen(temp_name)+1));
        strcpy(lmf_slab_file_name,temp_name);
        printf("Slab file name %s \n", lmf_slab_file_name);
        iarg+=3;
    }
    else if (strcmp(arg[iarg],"noforce")==0){
        if (iarg+1 > narg ) error->all(FLERR,"Illegal fix LMF command");
        no_force=1;
        iarg+=1;
    }
    else if (strcmp(arg[iarg],"check_charge")==0){
        if (iarg+1 > narg ) error->all(FLERR,"Illegal fix LMF command");
        check_charge=1;
        iarg+=1;
    }
    else if (strcmp(arg[iarg],"dwrite")==0){
      if (iarg+3>narg) error->all(FLERR,"Illegal fix LMF command");
      dyn_file_name=arg[iarg+1];
      dyn_file = fopen(dyn_file_name, "w");
      dyn_width = force->numeric(FLERR,arg[iarg+2]);
      iarg+=3;
    }
    /*else if (strcmp(arg[iarg],"tsmooth")==0){
      if (iarg+2>narg) error->all(FLERR,"Illegal fix LMF command");
      l_alpha = force->numeric(FLERR,arg[iarg+1]);
      iarg+=2;
    }*/
    else {
        printf("arg %s \n", arg[iarg]);
        error->all(FLERR,"Illegal fix LMF command");
    }
  }
  
  setup_LMF_lattice_slab();
  if(lmf_sigma==0) error->all(FLERR,"Did not specify sigma in Fix LMF");

	grid_inv_sq=1.0/grid_spac/grid_spac;	
	
	root_bool=0;
	for(int world =0 ; world<universe->nworlds; world++){
		if(universe->me==universe->root_proc[world]) 
			{
			root_bool=1;
			printf("root %d \n", universe->me);
			}
		else if(verbose_full==0) verbose=0;
	}

}

Fix_LMF_slab::~Fix_LMF_slab(){
  delete [] LMF_slab;
	delete [] old_LMF_slab;
	delete [] smooth_slab; 
	delete [] charge_dens; 
	delete [] o_dens; 
	delete [] h_dens; 
	delete [] smooth_slab_g;
	delete [] charge_dens_g;
	delete [] o_dens_g;
	delete [] h_dens_g;
	delete [] smooth_slab_wall;
	delete [] smooth_slab_wall_g;
	free(lmf_slab_file_name);
	

}

int Fix_LMF_slab::setmask(){
  int mask = 0;

  // FLD implicit needs to invoke wall forces before pair style

  if (fldflag) mask |= PRE_FORCE;
  else mask |= POST_FORCE;

  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

double * Fix_LMF_slab::create_array_slab(int size){
	double * array = new double [size];
	for(int n=0; n<size; n++)
		array[n]=0.0;
	return array;
}

void Fix_LMF_slab::setup_LMF_lattice_slab(){
	int dim_temp=2;
	lmf_hi=domain->boxhi[dim_temp]+bound_len_d;
	lmf_lo=domain->boxlo[dim_temp]-bound_len_d;
	double size_temp = lmf_hi-lmf_lo;

	n_tot = int(ceil(size_temp/float(grid_spac)))+2;
	grid_spac= size_temp/double(n_tot-2);
	
	if(root_bool)printf("gs %.6f box_hi %.6f box_lo %.6f bound_len_d %.6f lmf_hi_lo %.6f %.6f n_tot %d \n", grid_spac, domain->boxhi[dim_temp], domain->boxlo[dim_temp], bound_len_d, lmf_hi, lmf_lo, n_tot);
	LMF_slab = create_array_slab(n_tot);
	old_LMF_slab = create_array_slab(n_tot);
	smooth_slab = create_array_slab(n_tot);
	smooth_slab_wall = create_array_slab(n_tot);
	charge_dens = create_array_slab(n_tot);
	smooth_slab_g = create_array_slab(n_tot);
	smooth_slab_wall_g =create_array_slab(n_tot);
	charge_dens_g = create_array_slab(n_tot);
	o_dens_g= create_array_slab(n_tot);
	h_dens_g= create_array_slab(n_tot);
	o_dens= create_array_slab(n_tot);
	h_dens= create_array_slab(n_tot);
	//smooth_slab_temp = create_array_slab(n_tot);
	if(verbose)printf("Slab arrays created \n");
}

int Fix_LMF_slab::modify_param(int narg, char **arg){
  if (strcmp(arg[0],"par") == 0){
	  	if (narg<4 ) error->all(FLERR,"Illegal fix_modify LMF command");
    	lmf_sigma = force->numeric(FLERR,arg[1]);
    	lmf_cutoff = force->numeric(FLERR,arg[2]); 
    	bound_len_d = force->numeric(FLERR,arg[3]); 
    	return 4;
  	}
  else if (strcmp(arg[0],"dens_equil") == 0){
	  	if (narg<6 ) error->all(FLERR,"Illegal fix_modify LMF command");
	  	de_start = float(update->ntimestep);
    	de_dist = force->numeric(FLERR,arg[1]);
    	de_scale = force->numeric(FLERR,arg[2]);
    	int num_runs = force->numeric(FLERR,arg[4]);
    	Run run_obj(lmp); char * run_command[1]; run_command[0]=arg[3];
    	de_greater=1;
    	float * dens_vals = new float[num_runs]; float * length_vals = new float[num_runs]; float * scale_vals = new float[num_runs];
    	for(int runs = 0; runs<num_runs; runs++){
    		de_flag=0;
    		run_obj.command(1,run_command);	
    		de_flag=1;
    		run_obj.command(1,run_command);
    		dens_vals[runs]= de_dens/de_time_count; length_vals[runs]=sqrt((domain->boxhi[0]-domain->boxlo[0])*(domain->boxhi[1]-domain->boxlo[1]));
    		scale_vals[runs] = de_scale;
    		dequil();
    		if(universe->me==0)printf("run %d \n", runs);
    	}
    	de_flag=0;
    	if(universe->me==0) {
    		FILE * equil = fopen(arg[5],"w");
    		for(int runs=0; runs<num_runs; runs++) fprintf(equil, "%d %f %f %f  \n", runs, scale_vals[runs], dens_vals[runs], length_vals[runs]);
    	  fprintf(equil, "%d %f %f  \n", num_runs, de_scale, sqrt((domain->boxhi[0]-domain->boxlo[0])*(domain->boxhi[1]-domain->boxlo[1])));
    	  fclose(equil);
    	}
    	return 6;	
  	}
	else if (strcmp(arg[0],"write") == 0){
  		if (narg<2 ) error->all(FLERR,"Illegal fix_modify LMF command");
    	lmf_slab_write = force->numeric(FLERR,arg[1]);
    	return 2;
  	}
  else if (strcmp(arg[0],"dwrite")==0){
    if (narg<2 ) error->all(FLERR,"Illegal fix_modify LMF command"); 
    wdyn=1;
    wdyn_num=force->numeric(FLERR,arg[1]);
    if(universe->me==0)
      fprintf(dyn_file, "Dyn %d %f %d \n", int(update->ntimestep), E_Field_Period, wdyn_num);
    wdyn_count=0;
    return 2;
  	}
  return 0;
}

void Fix_LMF_slab::LMF_update_slab(){
  zero_array_slab(smooth_slab, n_tot);
  zero_array_slab(smooth_slab_wall, n_tot);
  if(verbose)printf("Updating \n");
  slab_calc();
  if(de_flag==1)
	dens_equilibrate();
}



void Fix_LMF_slab::subtract_av(double * array, int size){
	double tot_temp=0;
	for(int n_0=0; n_0<size; n_0++) tot_temp=tot_temp+array[n_0];

	tot_temp=tot_temp/float(size);
	for(int n_0=0; n_0<size; n_0++)array[n_0]=array[n_0]-tot_temp;

}

void Fix_LMF_slab::zero_array_slab(double * array, int size){
	for(int n_0=0; n_0<size; n_0++) array[n_0]=0.0;
}

void Fix_LMF_slab::slab_calc(){
    double **x = atom->x;
    int *mask = atom->mask;
    double *q = atom->q;
    int *type=atom->type;
    double **v=atom->v;
    int nall = atom->nlocal+atom->nghost;
    int *molindex=atom->molindex;
    tagint ** bond_atom=atom->bond_atom;

    double cut_sq=lmf_cutoff*lmf_cutoff;

    double sigma_inv_sq = 1.0/lmf_sigma/lmf_sigma;

    int n_steps=ceil(lmf_cutoff/grid_spac);

    int n_0, ind, m;
    int sdim=2;
    double delta, delt_sq, pos, r_sq, offset;

    int num_procs, my_rank;
    MPI_Comm_rank(universe->uorig, &my_rank);
    int n_worlds = universe->nworlds;

    MPI_Comm_size (universe->uorig, &(num_procs)); 

    double area = double(domain->boxhi[0]-domain->boxlo[0])*double(domain->boxhi[1]-domain->boxlo[1]);

    double pre_factor = 1.0/area/double(n_worlds)/grid_spac;
    double smooth_pre = 1.0/sqrt(M_PI)/lmf_sigma/area/double(n_worlds);
    double int_factor = force->qqr2e/force->qe2f*4.0*M_PI;


    zero_array_slab(charge_dens,n_tot);
    zero_array_slab(h_dens, n_tot);
    zero_array_slab(o_dens, n_tot);

    if(verbose)printf("Calculating Charge %d \n", my_rank);

    bool charge_bool=0;

    for (int i = 0; i < atom->nlocal; i++){
        if(check_charge==1){
                if(type[i]==1 || type[i]==2)
                    charge_bool=1;
                else charge_bool=0;
            }
        else charge_bool=1;
        if (mask[i] & groupbit) {
            double charge_pre = q[i]*smooth_pre;
            ind=int(floor((x[i][sdim]-lmf_lo)/grid_spac));
            offset=x[i][sdim]-float(ind)*grid_spac-lmf_lo;
            
            if (offset>grid_spac*.5) {
                ind+=1;
                offset=offset-grid_spac;	
            }
            charge_dens[ind]+=pre_factor*q[i]; 
            for(n_0 = -n_steps; n_0 <= n_steps; n_0++){
                delta=offset-n_0*grid_spac;
                r_sq=delta*delta;
                if(r_sq<cut_sq){
                    smooth_slab[n_0+ind]+=charge_pre*exp(-r_sq*sigma_inv_sq);
                    if(charge_bool==0) smooth_slab_wall[n_0+ind]+=charge_pre*exp(-r_sq*sigma_inv_sq);
                }
                
            }
            if(type[i]==2) {
                h_dens[ind]+=pre_factor;
            }
            else if (type[i]==1){
                o_dens[ind]+=pre_factor;
            }	
        }
    }
		
	//for(int i=0; i<n_tot; i++)
	//    smooth_slab[i]=smooth_slab[i]*(1.0-l_alpha)+smooth_slab_temp[i]*l_alpha;

	if(verbose)printf("Charge calculated %d \n", my_rank);
	
	if(ensemble_flag){
		MPI_Allreduce(smooth_slab,smooth_slab_g,n_tot,MPI_DOUBLE,MPI_SUM,universe->uorig);
		MPI_Allreduce(smooth_slab_wall,smooth_slab_wall_g,n_tot,MPI_DOUBLE,MPI_SUM,universe->uorig);
		MPI_Allreduce(charge_dens,charge_dens_g,n_tot,MPI_DOUBLE,MPI_SUM,universe->uorig);
		MPI_Allreduce(h_dens,h_dens_g,n_tot,MPI_DOUBLE,MPI_SUM,universe->uorig);
		MPI_Allreduce(o_dens,o_dens_g,n_tot,MPI_DOUBLE,MPI_SUM,universe->uorig);			
	}
	if(verbose)printf("MPI_Allreduce Done %d \n", my_rank);
	for(int n_0=0; n_0<n_tot; n_0++) {
        smooth_slab[n_0]=smooth_slab_g[n_0];
        smooth_slab_wall[n_0]=smooth_slab_wall_g[n_0];
        charge_dens[n_0]=charge_dens_g[n_0];
        h_dens[n_0]=h_dens_g[n_0];
        o_dens[n_0]=o_dens_g[n_0];
        LMF_slab[n_0]=0.0;
    }

	double_int(smooth_slab,LMF_slab,grid_spac,n_tot,int_factor);
	
	if(wdyn==1){
        write_dyn();
        if(wdyn_count>wdyn_num) wdyn=0;
    }

}

void Fix_LMF_slab::double_int(double * array, double * out, double delt, int size, double pre_fact){
	out[0]=0.0;
	double * e_field_temp;
	e_field_temp = new double[size];
	e_field_temp[0]=0.0;
	delt=.5*delt;
	for (int i=0; i<size-1; i++){
		e_field_temp[i+1]=e_field_temp[i]+delt*(array[i]+array[i+1]);
		//printf("i %d %d \n",i, size);
		}
	delt=-pre_fact*delt;
	for (int i=0; i<size-1; i++)
		out[i+1]=out[i]+delt*(e_field_temp[i]+e_field_temp[i+1]);
	
	delete []e_field_temp;
}

void Fix_LMF_slab::LMF_force_slab(){
    double grad;
    int pos_index;

    double **x = atom->x;
    double **f = atom->f;
    double *q = atom->q;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int *type = atom->type;

    bool inside;
    bigint ntimestep=update->ntimestep;
    int sdim=2;

    double temp_force;
    if(!no_force){
        for (int i = 0; i < nlocal; i++)
            if (mask[i] & groupbit) {
                pos_index = int((x[i][sdim]-lmf_lo)/grid_spac);
                if(pos_index<0 || pos_index>(n_tot-2)) error->all(FLERR,"Particle outside LMF box");
                double pre_factor = q[i]*force->qe2f; //convert to kcal/mol
                temp_force=pre_factor*(LMF_slab[pos_index+1]-LMF_slab[pos_index])/grid_spac;
                if(type[i]<3) f[i][sdim] -= temp_force;
            } 
        }
}

void Fix_LMF_slab::dens_equilibrate(){
	double time_dif = float(update->ntimestep)-de_start;
	double dens_temp = 0.0; double middle = (lmf_hi+lmf_lo)/2.0; double count_temp =0.0; 
	for (int i=0; i<n_tot; i++){
		double pos_temp = lmf_lo+i*grid_spac;
		if(pos_temp>(middle-de_dist) && pos_temp<(middle+de_dist))
			{
			dens_temp+=o_dens_g[i];
			count_temp+=1.0;
			}
	}
	de_dens+=dens_temp/count_temp;
	de_time_count+=1.0;
}


void Fix_LMF_slab::dequil(){
    if(root_bool)printf("Equilibrating z_length \n");
		double water_dens = .03342639;
	  char * box_command[8];
	  double sf = 1.0+de_scale; 
		de_dens=de_dens/de_time_count;
		ChangeBox cb(lmp);
		char o1[2]="x"; char sc[6]="scale"; char o2[2]="y"; char buf[10]; char gr[4]="all"; char rem[6]="remap";
		box_command[1]=o1; box_command[4]=o2; box_command[2]=sc; box_command[5]=sc; box_command[7]=rem;
		if(de_dens<water_dens){
			sprintf(buf,"%f", 1.0/sf); 
			if(de_greater) {
				if(root_bool)printf("changing scale %f \n", de_scale/2.0);
				de_scale=de_scale/2.0;
				}
			de_greater=0;
		}
		else{
			sprintf(buf,"%f", sf); 
			if(!de_greater) {
				if(root_bool)printf("changing scale %f \n", de_scale/2.0);
				de_scale=de_scale/2.0;
			}
			de_greater=1;
		}
		box_command[3]=buf; box_command[6]=buf;
		box_command[0]=gr;
		if(root_bool) for (int i=0; i<7; i++) printf("%s \n", box_command[i]);
		if(root_bool) printf("density %f \n", de_dens);
		cb.command(8, box_command);
		de_start=float(update->ntimestep);
		de_dens=0.0;
		de_time_count=0.0;
		if(root_bool)printf("Done equilibrating z_length\n");
}

void Fix_LMF_slab::write_dyn(){
    int middle = int(float(n_tot)/2.0);
    int width = int(dyn_width/grid_spac);
    double slope =  (LMF_slab[middle+width]-LMF_slab[middle-width])/float(2.0*width*grid_spac);
    wdyn_count+=1;
    if(universe->me==0)
      fprintf(dyn_file, "%d %d %f %f \n",int(update->ntimestep), wdyn_count, slope, E_Field_Variable);    
}

/* ---------------------------------------------------------------------- */

void Fix_LMF_slab::setup(int vflag){
 post_force(vflag);
}

void Fix_LMF_slab::pre_force(int vflag){
  post_force(vflag);
}

void Fix_LMF_slab::post_force(int vflag){

  LMF_force_slab();
}


