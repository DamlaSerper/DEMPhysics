/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <cmath>
#include <string.h>
#include <stdlib.h>
#include "fix_centrifugal.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

/* Centrifugal or centripetal force (magnitudes equal) found as particle
   mass times angular velocity of the rotating system times radial distance
   from the axis of the centrifuge. Assuming here that the centrifuge axis
   is vertical (aligned with the y axis) and passes through the origin.

   This is not enabled for multisphere particles at present and a warning
   will be issued if using this fix with multisphere particles is
   attempted.*/

FixCentrifugal::FixCentrifugal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal fix centrifugal command");
  omega = force->numeric(FLERR,arg[3]); // Angular velocity in SI units
  omegasq = omega*omega;

  // No centrifugal force applied to particles within the inlet pipe
  inletpiperad = force->numeric(FLERR,arg[4]);
  inletpipey = force->numeric(FLERR,arg[5]);

  // No centrifugal force applied to particles outside the basket
  basketr = force->numeric(FLERR,arg[6]);

  if (basketr < 0.0 || inletpiperad < 0.0) 
    error->all(FLERR,"Inappropriate parameters in fix centrifugal command");
}

/* ---------------------------------------------------------------------- */

int FixCentrifugal::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCentrifugal::init()
{
  /* Check for the presence of multisphere particles and issue an error
     if present */
  int nms = modify->n_fixes_style("multisphere");
  if(nms > 0)
    error->all(FLERR,"Support for fix multisphere not implemented in fix centrifugal");
}

/* ---------------------------------------------------------------------- */

void FixCentrifugal::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    int nlevels_respa = ((Respa *) update->integrate)->nlevels;
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixCentrifugal::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCentrifugal::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  double r; // Distance from centre in x-z plane
  double fmag; // Magnitude of centrifugal force

  /* Note that f[i][0] += fmag*x[i][0]/r is equivalent to
     f[i][0] += rmass[i]*omegasq*x[i][0]*/

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	r = sqrt(x[i][0]*x[i][0] + x[i][2]*x[i][2]);
	//fmag = rmass[i]*omegasq*r;
	//f[i][0] += fmag*x[i][0]/r;
	//f[i][2] += fmag*x[i][2]/r;

	/* Apply the centrifugal force only to particles outside the inlet
	   pipe but within the basket.*/
	if (r < basketr && (r > inletpiperad || x[i][1] < inletpipey)) {	
	  f[i][0] += rmass[i]*omegasq*x[i][0];
	  f[i][2] += rmass[i]*omegasq*x[i][2];
	}
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	r = sqrt(x[i][0]*x[i][0] + x[i][2]*x[i][2]);
	//fmag = mass[type[i]]*omegasq*r;
	//f[i][0] += fmag*x[i][0]/r;
	//f[i][2] += fmag*x[i][2]/r;

	/* Apply the centrifugal force only to particles outside the inlet
	   pipe but within the basket.*/
	if (r < basketr && (r > inletpiperad || x[i][1] < inletpipey)) {
	  f[i][0] += mass[type[i]]*omegasq*x[i][0];
	  f[i][2] += mass[type[i]]*omegasq*x[i][2];
	}
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixCentrifugal::post_force_respa(int vflag, int ilevel, int iloop)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCentrifugal::min_post_force(int vflag)
{
  post_force(vflag);
}
