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

FixStyle(lmf,Fix_LMF)

#else

#ifndef LMP_LMF_H
#define LMP_LMF_H

#include "fix.h"
#include "lmf_grid.h"
#include <complex>

namespace LAMMPS_NS {

class Fix_LMF : public Fix {
    public:
        Fix_LMF(class LAMMPS *, int, char **);
        ~Fix_LMF();
        
        LMF_GRID *LMF_potential, *smooth_charge, *LMF_pi;
        
        bool ensemble_flag, no_force, root_bool, slab_flag, debug, SOR_flag, DLMF_flag, slab_bool;
        
        double lmf_sigma, lmf_cutoff;
        double lmf_c, lmf_tau, lmf_pf;
        double *grid_spac, *grid_inv_sq, *lmf_lo;
        int yb, *n_tot, sc_freq, sc_count, lmf_it_freq;
        double SOR_tol, omega;

        void setup_LMF_lattice();
        void LMF_force();
        void LMF_update();
        void smooth_calc();
        
        char * lmf_file_name;
        FILE * lmf_file;
        int lmf_write;
        
        void compute_yb_correction();
        double yb_correction, lmf_vol;
        
        void pre_force(int);
        void post_force(int);
        void setup(int);
        int setmask();
        int modify_param(int , char **);
        
        
        FILE * force_file;
        char * force_file_name;
        bool force_flag;
        int tag_particle;
        double ew_force[3];
        
        
    private:
        int fldflag;
        bool verbose;
        int poisson_SOR();
        void dlmf_iterate();
        void smooth_calc(double ***);
        void update_boundary(int, int, int, double ***, int *);
        int return_opposite(int,int);
        void subtract_av(LMF_GRID *);
        void pos_index_calc(int, int *, double **, double  *, bool);
        double SOR_iterate(LMF_GRID *, LMF_GRID *, double, int *);
        
};

}

#endif
#endif