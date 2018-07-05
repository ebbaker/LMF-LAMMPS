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

#ifdef FIX_CLASS

FixStyle(lmf_slab,Fix_LMF_slab)

#else

#ifndef LMP_LMF_SLAB_H
#define LMP_LMF_SLAB_H

#include "fix.h"
#include <complex>

namespace LAMMPS_NS {

//Assumes z direction is fixed slab, X Y periodic

class Fix_LMF_slab : public Fix {
 public:
    Fix_LMF_slab(class LAMMPS *, int, char **);
    ~Fix_LMF_slab();
    void dequil();

    bool ensemble_flag, dens_flag, no_force, de_flag, de_greater, root_bool;

    double de_dens, de_start, de_av, de_time, de_dist, de_time_count, de_scale;

    double *LMF_slab, *old_LMF_slab, *smooth_slab, *smooth_slab_wall, *charge_dens, *o_dens, *h_dens; 
    double *smooth_slab_g, *smooth_slab_wall_g, *charge_dens_g, *o_dens_g, *h_dens_g;
    //double *smooth_slab_temp;
    void setup_LMF_lattice_slab();
    double * create_array_slab(int);

    int modify_param(int , char **);

    void pre_force(int);
    void post_force(int);

    void LMF_force_slab();
    void LMF_update_slab();

    void setup(int);
    int setmask();

    int lmf_slab_write;
    bool check_charge;
    char * lmf_slab_file_name, * dyn_file_name;
    double E_Field_Variable, E_Field_Period;
    double lmf_sigma, lmf_cutoff, grid_spac, grid_inv_sq;
    double lmf_hi, lmf_lo, bound_len_d; //, l_alpha;

    bool verbose, verbose_full;
    bool wdyn;
    int wdyn_num, wdyn_count;
    double dyn_width;
  
    FILE * lmf_slab_file;
    FILE * dyn_file;

    int n_tot;

    void write_dyn();

    protected:

    int fldflag;
    void double_int(double *, double *, double, int, double);
    void dens_equilibrate();
    void subtract_av(double * , int );
    void zero_array_slab(double * , int);
    void slab_calc();
};

}

#endif
#endif

