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
#include "fix_elec_field.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "stdlib.h"
#include "string.h"

#include "fix_lmf_slab.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */


FieldElectric::FieldElectric(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {
  first_check=1;
  fldflag=0; full_slab=0;
  
  int iarg = 3;
  while (iarg < narg) {
  	if (strcmp(arg[iarg],"par") == 0){
	  	if (iarg+3 > narg ) error->all(FLERR,"Illegal fix LMF command");
    	efield = force->numeric(FLERR,arg[iarg+1]);
    	period = force->numeric(FLERR,arg[iarg+2]); 
    	if(period>.00001) freq = 2.0*M_PI/double(period);
    	else freq=0.0;
    	iarg+=3;
  	}   
  	else if (strcmp(arg[iarg],"full_slab") == 0){
  		if (iarg+1 > narg ) error->all(FLERR,"Illegal fix LMF command");
	  	full_slab=1;
    	iarg+=1;
  	} 
  	else {
  		printf("arg %s \n", arg[iarg]);
  		error->all(FLERR,"Illegal fix LMF command");
  	}
  }
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type=atom->type;
  surface_num=0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
    	if(type[i]==3 || type[i]==4)
    		surface_num+=.5;
    }
  double surface_num_g = 0;
  MPI_Allreduce(&surface_num,&surface_num_g,1,MPI_DOUBLE,MPI_SUM,lmp->world);
  surface_num=surface_num_g;
  double area = (domain->boxhi[0]-domain->boxlo[0])*(domain->boxhi[1]-domain->boxlo[1])/surface_num;
  var1=force->qe2f/force->qqr2e/4.0/3.141592653*efield*area; 
  //printf("charge %.6f epsilon %.6f \n", var1, force->qe2f/force->qqr2e/4.0/3.141592653);
}

FieldElectric::~FieldElectric(){


}

int FieldElectric::setmask(){
  int mask = 0;

  // FLD implicit needs to invoke wall forces before pair style

  if (fldflag) mask |= PRE_FORCE;
  else mask |= POST_FORCE;

  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */



void FieldElectric::pre_force(int v_flag)
{
	double area = (domain->boxhi[0]-domain->boxlo[0])*(domain->boxhi[1]-domain->boxlo[1])/surface_num;
  var1=force->qe2f/force->qqr2e/4.0/3.141592653*efield*area; 
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

void FieldElectric::post_force(int v_flag)
{
  double delta,rinv,r2inv,r4inv,r10inv,fwall;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *mask = atom->mask;
  int nall = atom->nmax;
  int *type=atom->type;



  int onflag = 0;
  if(first_check==1){
  	first_time=float(update->ntimestep);
  	first_check=0;
  }
  double area = (domain->boxhi[0]-domain->boxlo[0])*(domain->boxhi[1]-domain->boxlo[1])/surface_num;
  var1=force->qe2f/force->qqr2e/4.0/3.141592653*efield*area; 
  if (fabs(freq)>=0.000001)
  	var2=sin(double(freq)*(double(update->ntimestep)-double(first_time)));
  else var2=1.0;
  int lmf_slab = modify->lmf_slab_check();
  if(lmf_slab!=-1)
  	{
  	Fix_LMF_slab * slab_ptr = dynamic_cast <Fix_LMF_slab *> (modify->fix[lmf_slab]);
  	slab_ptr->E_Field_Variable=efield*var2; slab_ptr->E_Field_Period=period;
  	}
  //if(int(update->ntimestep-first_time)%50<1)printf("var2 %.6f efield %.6f var1 %.6f time_dif %.6f freq %.10f \n", var2, efield, var1, double(update->ntimestep)-double(first_time), freq);
  
  double q_temp=var1*var2; double ftemp = force->qe2f*efield*var2;
  for (int i = 0; i < nall; i++)
    if (mask[i] & groupbit) {
    	if(!full_slab){
    		if(type[i]==3)
    			q[i]=q_temp;
    		else if(type[i]==4)
    			q[i]=-q_temp;
    	}
    	else{
    		f[i][2]-=q[i]*ftemp;
    	}
    }

}

int FieldElectric::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"per") == 0){
	  	if (narg<2 ) error->all(FLERR,"Illegal fix_modify LMF command");
    	period = force->numeric(FLERR,arg[1]); 
    	if(period>.00001) freq = 2.0*M_PI/double(period);
    	else freq=0.0;
    	first_time=float(update->ntimestep);
    	return 2;
  	}
  else if (strcmp(arg[0],"ef") == 0){
	  	if (narg<2 ) error->all(FLERR,"Illegal fix_modify LMF command");
    	efield = force->numeric(FLERR,arg[1]);
    	return 2;
  	}
  return 0;
}

void FieldElectric::setup(int vflag){
 post_force(vflag);
}
