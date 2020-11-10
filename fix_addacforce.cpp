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

#include "fix_addacforce.h"

#include <stdlib.h>
#include <string.h>

#include "accelerator_kokkos.h"
#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "respa.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

using namespace FixConst;

using namespace ACFORCE;

enum { NONE, CONSTANT, EQUAL, ATOM };

enum { XLO = 0, XHI = 1, YLO = 2, YHI = 3, ZLO = 4, ZHI = 5 };

/* ---------------------------------------------------------------------- */

FixAddACForce::FixAddACForce(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg), sforce(NULL) {
  double dd = 0;
  if (narg < 13) error->all(FLERR, "Illegal fix addacforce command");
  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
  dt = dtf = dtv = 0;
  // Get source power
  if (strcmp(arg[3], "NULL") == 0)
    error->all(FLERR, "Illegal fix addacforce command");
  else {
    xi = force->numeric(FLERR, arg[3]);
  }
  // Get equilibrium density
  if (strcmp(arg[4], "NULL") == 0)
    error->all(FLERR, "Illegal fix addacforce command");
  else {
    densi = force->numeric(FLERR, arg[4]);
  }

  // Get sound speed
  if (strcmp(arg[5], "NULL") == 0)
    error->all(FLERR, "Illegal fix addacforce command");
  else {
    sounde = force->numeric(FLERR, arg[5]);
  }

  // Get wavevector modulus ie \frac{2\pi}{\lambda}
  if (strcmp(arg[6], "NULL") == 0)
    error->all(FLERR, "Illegal fix addacforce command");
  else {
    kt = force->numeric(FLERR, arg[6]);
  }

  // Get wavevector direction
  if (strcmp(arg[7], "NULL") == 0)
    error->all(FLERR, "Illegal fix addacforce command");
  else {
    k.ar[0] = force->numeric(FLERR, arg[7]);
  }
  if (strcmp(arg[8], "NULL") == 0)
    error->all(FLERR, "Illegal fix addacforce command");
  else {
    k.ar[1] = force->numeric(FLERR, arg[8]);
  }
  if (strcmp(arg[9], "NULL") == 0)
    error->all(FLERR, "Illegal fix addacforce command");
  else {
    k.ar[2] = force->numeric(FLERR, arg[9]);
  }

  dd = kt / sqrt(k * k);
  k *= dd;

  if (strcmp(arg[10], "NULL") == 0)
    error->all(FLERR, "Illegal fix addacforce command");
  else {
    eta = force->numeric(FLERR, arg[10]);
  }

  if (strcmp(arg[11], "NULL") == 0)
    error->all(FLERR, "Illegal fix addacforce command");
  else {
    beta = force->numeric(FLERR, arg[11]);
  }

  if (strcmp(arg[12], "NULL") == 0)
    error->all(FLERR, "Illegal fix addacforce command");
  else {
    psi = force->numeric(FLERR, arg[12]);
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  maxatom = 1;
  memory->create(sforce, maxatom, 4, "addacforce:sforce");
  time_origin = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

FixAddACForce::~FixAddACForce() { memory->destroy(sforce); }

/* ---------------------------------------------------------------------- */

int FixAddACForce::setmask() {
  datamask_read = datamask_modify = 0;

  int mask = 0;

  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddACForce::init() {
  dt = update->dt;
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

void FixAddACForce::setup(int vflag) {
  if (strstr(update->integrate_style, "verlet"))
    post_force(vflag);
  else {
    ((Respa *)update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag, ilevel_respa, 0);
    ((Respa *)update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixAddACForce::min_setup(int vflag) { post_force(vflag); }

/* ---------------------------------------------------------------------- */

void FixAddACForce::post_force(int vflag) {
  imageint *image = atom->image;
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  C_Vector_3 unwrap, ff;
  double delta = (update->ntimestep - time_origin) * dt;  // simulation time
  double omega = sounde * kt;
  double rphase = 0, tphase = omega * delta, lmass = 0, xip = sqrt(xi);
  double st1 = sin(tphase), ct1 = cos(tphase), cr1 = 0, sr1 = 0;
  double /*st2 = 2.0 * st1 * ct1, ct2 = ct1 * ct1 - st1 * st1, */ cr2 = 0,
                                                                  sr2 = 0;
  double rho = 0, cr3 = 0, sr3 = 0, cr5 = 0, fm = 0, c02 = sounde * sounde;

  if (atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce, maxatom, 4, "addacforce:sforce");
  }

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;
  if (xi > 0.0) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        domain->unmap(x[i], image[i], unwrap.ar);

        foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];

        if (rmass) {
          lmass = rmass[i];
        } else {
          lmass = mass[type[i]];
        }

        rphase = k * unwrap;

        sr1 = sin(rphase);
        cr1 = cos(rphase);
        sr2 = 2.0 * sr1 * cr1;
        cr2 = cr1 * cr1 - sr1 * sr1;
        cr3 = cr2 * cr1 - sr2 * sr1;
        sr3 = cr2 * sr1 + sr2 * cr1;
        cr5 = cr3 * cr2 - sr3 * sr2;

        //(1 + x Cos[t w] Sin[Kr]) Subscript[\[Rho], 0]
        rho = densi * (1.0 + xip * ct1 * sr1);
        // 1/16  x^2  Sin[2 Kr] Subscript[\[Rho], 0] (c0^2 (-8 + x^2 Sin[Kr]^2)
        // + 8 psi Subscript[\[Rho], 0])
        fm = 0.0625 * xi * sr2 * densi *
             (c02 * (-8.0 + xi * sr1 * sr1) + 8.0 * psi * densi);

        // 1/128 c0^2  x^3 ((-16 + x^2) Cos[3 Kr] -
        // x^2 Cos[5 Kr]) Subscript[\[Rho], 0] Cos[ wt]
        fm += 0.0078125 * c02 * xi * xip * ((-16 + xi) * cr3 - xi * cr5) *
              densi * ct1;

        //-(1/2) (1 + beta) c0 eta k Kv x Cos[Kr] Sin[wt]

        fm -= 0.5 * (1.0 + beta) * eta * omega * cr1 *
              st1;  // Viscosity correction

        fm *= lmass / densi;

        ff = k * fm;

        foriginal[0] += 0.5 * lmass * sounde * sounde * (rho - densi);
        f[i][0] += ff.ar[0];
        f[i][1] += ff.ar[1];
        f[i][2] += ff.ar[2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixAddACForce::post_force_respa(int vflag, int ilevel, int iloop) {
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddACForce::min_post_force(int vflag) { post_force(vflag); }

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixAddACForce::compute_scalar() {
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal, foriginal_all, 4, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixAddACForce::compute_vector(int n) {
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal, foriginal_all, 4, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return foriginal_all[n + 1];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixAddACForce::memory_usage() {
  return (double)(maxatom * 4 * sizeof(double));
}

/* ---------------------------------------------------------------------- */
