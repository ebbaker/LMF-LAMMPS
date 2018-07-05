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

FixStyle(efield_mine,FieldElectric)

#else

#ifndef LMP_FIELD_ELECTRIC_H
#define LMP_FIELD_ELECTRIC_H

#include "fix_wall.h"

namespace LAMMPS_NS {

class FieldElectric : public Fix {
 public:

  FieldElectric(class LAMMPS *, int, char **);
  ~FieldElectric();
  void pre_force(int);
  void post_force(int);
  bool first_check, full_slab;
  double first_time;
  int modify_param(int , char **);
  double surface_num;
  
  int setmask();
  void setup(int);

 private:
  int fldflag;
  double var1, var2, efield, period, freq;
};

}

#endif
#endif
