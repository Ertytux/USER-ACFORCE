/* -*- c++ -*- ----------------------------------------------------------
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

FixStyle(addacforce, FixAddACForce)
#else

#ifndef LMP_FIX_ADDACFORCE_H
#define LMP_FIX_ADDACFORCE_H

#include <cmath>
#include <cstring>

#include "math_vector.h"

namespace ACFORCE {

struct C_Vector_3 {
  LAMMPS_NS::vector ar;

  C_Vector_3() { LAMMPS_NS::vec_null(ar); }
  C_Vector_3(const C_Vector_3 &rhs) {
    memcpy(ar, rhs.ar, sizeof(LAMMPS_NS::vector));
  }

  C_Vector_3 &operator=(const C_Vector_3 &rhs) {
    if (this != &rhs) {
      memcpy(ar, rhs.ar, sizeof(LAMMPS_NS::vector));
    }
    return *this;
  }

  double &operator[](const int &i) { return ar[i % 3]; }

  C_Vector_3 &operator+=(const C_Vector_3 &prhs) {
    ar[0] += prhs.ar[0];
    ar[1] += prhs.ar[1];
    ar[2] += prhs.ar[2];
    return *this;
  }

  C_Vector_3 operator+(const C_Vector_3 &rhs) const {
    C_Vector_3 C(*this);
    C += rhs;
    return C;
  }

  C_Vector_3 &operator-=(const C_Vector_3 &prhs) {
    ar[0] -= prhs.ar[0];
    ar[1] -= prhs.ar[1];
    ar[2] -= prhs.ar[2];
    return *this;
  }

  C_Vector_3 operator-(const C_Vector_3 &rhs) const {
    C_Vector_3 C(*this);
    C -= rhs;
    return C;
  }

  double operator*(const C_Vector_3 &rhs) const {
    return ar[0] * rhs.ar[0] + ar[1] * rhs.ar[1] + ar[2] * rhs.ar[2];
  }

  C_Vector_3 &operator*=(const double &rhs) {
    ar[0] *= rhs;
    ar[1] *= rhs;
    ar[2] *= rhs;
    return *this;
  }
  C_Vector_3 operator*(const double &rhs) const {
    C_Vector_3 C(*this);
    C *= rhs;
    return C;
  }
  friend C_Vector_3 operator*(const double &lhs, const C_Vector_3 &rhs);
};

inline C_Vector_3 operator*(const double &lhs, const C_Vector_3 &rhs) {
  C_Vector_3 C(rhs);
  C *= lhs;
  return C;
}

}  // namespace ACFORCE

#include "fix.h"

namespace LAMMPS_NS {

class FixAddACForce : public Fix {
 public:
  FixAddACForce(class LAMMPS *, int, char **);
  ~FixAddACForce();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

 private:
  double foriginal[4], foriginal_all[4], time_origin;
  int force_flag, ilevel_respa, maxatom;
  double **sforce;
  double dt, dtv, dtf;
  double xi, densi, sounde, kt, eta, beta, psi;
  ACFORCE::C_Vector_3 k;
};

}  // namespace LAMMPS_NS

#endif
#endif

    /* ERROR/WARNING messages:

    E: Illegal ... command

    Self-explanatory.  Check the input script syntax and compare to the
    documentation for the command.  You can use -echo screen as a
    command-line option when running LAMMPS to see the offending line.

    E: Region ID for fix addforce does not exist

    Self-explanatory.

    E: Variable name for fix addcforce does not exist

    Self-explanatory.

    E: Variable for fix addcforce is invalid style

    Self-explanatory.

    E: Cannot use variable energy with constant force in fix addcforce

    This is because for constant force, LAMMPS can compute the change
    in energy directly.

    E: Must use variable energy with fix addcforce

    Must define an energy vartiable when applyting a dynamic
    force during minimization.

    */
