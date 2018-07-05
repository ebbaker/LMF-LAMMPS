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

FixStyle(analysis,analysis)

#else

#ifndef LMP_ANALYSIS_H
#define LMP_ANALYSIS_H

#include "fix.h"
#include <complex>

namespace LAMMPS_NS {

class analysis : public Fix {
 public:
    analysis(class LAMMPS *, int, char **);
    ~analysis();
  
    bool grflag, autflag;

    int gr_freq, gr_tot, gr_count;
    char * gr_file_name;
    double * gr_oo, * gr_oh, *gr_hh;
    double * gr_oo_u, * gr_oh_u, *gr_hh_u;
    double gr_spac, gr_max, l_alpha;

    int n_aut, aut_count, num_lags, aut_write_freq, force_write_freq;
    int tag_particle;
    std::complex<double> ** charge_autvals, ** charge_autocorrelation;
    std::complex<double> ** charge_autvals_t, ** charge_autocorrelation_t;
    char * aut_file_name;

    void update_aut_vals();
    void update_autocorr();
    void aut_update();
    void lag_back();

    void calculate_gr();
    void write_gr();
    void write_aut();

    double * create_1d_array(int);
    std::complex<double> ** create_2d_array(int, int);

    int modify_param(int , char **);

    void pre_force(int);
    void post_force(int);

    void setup(int);
    int setmask();


    bool verbose;
    bool wdyn;
    int wdyn_num, wdyn_count;

    FILE * gr_file;
    FILE * aut_file;
  
    protected:

    int fldflag;
    void zero_1d_array(double * , int );
    void zero_2d_array(std::complex<double> ** , int, int );

};

}

#endif
#endif

