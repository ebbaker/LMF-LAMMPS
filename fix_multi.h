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

FixStyle(multi,Fix_MULTI)

#else

#ifndef LMP_MULTI_H
#define LMP_MULTI_H

#include "fix.h"
#include "multi_grid.h"
#include <complex>
#include <random>

namespace LAMMPS_NS {

class Fix_MULTI : public Fix {
    public:
        Fix_MULTI(class LAMMPS *, int, char **);
        ~Fix_MULTI();
        
        MULTI_GRID ** potential, ** sc;
        
        bool ensemble_flag, no_force, root_bool, slab_flag, debug, ALL_flag, slab_write, rad_write, force_ew_write;
        bool error_flag;
        
        double lmf_pf, lmf_sigma, sigma_in, gs_ratio, pow_base, slab_vol_factor, sc_pf;
        double ** grid_spac, ** grid_lo, *omega, *hpf;
        int ** n_tot, * it_freq;
        int num_grids, output_write_freq, force_write_freq, num_rand;
        
        double force_ew[3], max_rad_write, rad_write_spacing;
        
        std::mt19937 gen;
        
        void MULTI_force();
        void MULTI_update();
        void pot_update(int);
        
        void pre_force(int);
        void post_force(int);
        void setup(int);
        int setmask();
        int modify_param(int , char **);
        
        
        FILE * force_file, * output_file;
        char * force_file_name, * output_file_name;
        bool force_flag;
        int tag_particle;
        double ew_force[3];
        
        
    private:
    
        void rand_add();
        int fldflag;
        bool verbose;
        void setup_multi_lattice();
        void setup_lmf_lattice(int);
        
};

}

#endif
#endif