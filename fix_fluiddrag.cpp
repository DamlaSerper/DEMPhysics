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
#include "fix_fluiddrag.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "math_const.h"
#include "error.h"
#include "force.h"
#include "comm.h" // Add this library to print to log file

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

/* The drag force exerted by a fluid on a particle is given by 0.5 times the
   drag coefficient times the cross-sectional area of the particle times the
   fluid density times the square of relative velocity between the particle
   and fluid. Assuming here that the centrifuge axis is vertical (aligned
   with the y axis) and passes through the origin. Also assuming that the
   inlet pipe is vertical and the top circle of this cylinder is at (0,0,0).

   This is not enabled for multisphere particles at present and a warning
   will be issued if using this fix with multisphere particles is
   attempted.*/

FixFluiddrag::FixFluiddrag(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 11) error->all(FLERR,"Illegal fix fluiddrag command");
  omega = force->numeric(FLERR,arg[3]); // Angular velocity in SI units
  fluiddensity = force->numeric(FLERR,arg[4]);
  fluidvisc = force->numeric(FLERR,arg[5]);
  // No drag force applied to particles within the inlet pipe
  inletpiperad = force->numeric(FLERR,arg[6]);
  inletpipey = force->numeric(FLERR,arg[7]);
  //size multiplier
  size_mult = force->numeric(FLERR,arg[8]);
  // basket radius
  basketrad = force->numeric(FLERR,arg[9]);
  // fluid velocity at inlet pipe
  vfluidver_0 = force->numeric(FLERR,arg[10]);
  
  
  if (omega < 0.0) 
    error->all(FLERR,"Currently limited to positive omega in fix fluiddrag command to make the maths easier");

  if (fluiddensity < 0.0 || fluidvisc < 0.0 || inletpiperad < 0.0) 
    error->all(FLERR,"Inappropriate parameters in fix fluiddrag command");
}

/* ---------------------------------------------------------------------- */

int FixFluiddrag::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFluiddrag::init()
{
  /* Check for the presence of multisphere particles and issue an error
     if present */
  int nms = modify->n_fixes_style("multisphere");
  if(nms > 0)
    error->all(FLERR,"Support for fix multisphere not implemented in fix fluiddrag");
}

/* ---------------------------------------------------------------------- */

void FixFluiddrag::setup(int vflag)
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

void FixFluiddrag::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFluiddrag::post_force(int vflag)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double r; // Distance of particle from centre in x-z plane
  double vfluid,vfluid_ver,umag,umag_ver; // Magnitudes of the fluid and relative velocity vectors
  double vxfluid,vzfluid; // Signed components of fluid velocity vector
  double relvx,relvy,relvz; // Relative velocity components between fluid and particle
  double re,re_r,re_ver,re_r_ver,cd,cd_r,cd_ver,cd_r_ver,cd_mult,cd_mult_ver; // Reynolds number, drag coefficient
  double aparticle,aparticle_r; // Cross-sectional area of particle
  double fdmag,fdmag_r,fdmag_ver,fdmag_r_ver; // Magnitude of drag force
  double compoundreterm = 2.0*fluiddensity/fluidvisc;
  double ratio, ratio_ratio, ratio_ver, ratio_ratio_ver; // ratio of fdmag/fdmag_r) and (fdmag/fdmag_r)/pow(size_mult*3);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      r = sqrt(x[i][0]*x[i][0] + x[i][2]*x[i][2]);
      vfluid = omega*r;

      // Apply the drag force only to particles outside the inlet pipe
      if (r > inletpiperad || x[i][1] < inletpipey) {
	/* Now find the signed components of fluid velocity. This is a bit
	   confusing because the x component of fluid velocity at x = 0 is a
	   maximum. The signs also need to be calculated carefully.*/
	vxfluid = vfluid*fabs(x[i][2])/r;
	if (x[i][2] > 0.0) vxfluid *= -1.0;

	vzfluid = vfluid*fabs(x[i][0])/r;
	if (x[i][0] < 0.0) vzfluid *= -1.0;

	relvx = vxfluid - v[i][0];
	relvz = vzfluid - v[i][2];

	umag = sqrt(relvx*relvx + relvz*relvz);

	// _r stands for non-scaled (real) particle sizes
	re = compoundreterm*umag*radius[i];
	re_r = compoundreterm*umag*radius[i]/size_mult;
	
	if (re < 1.0e-6) {
		continue; // In case umag is very, very small which would lead to problems
	} else if (re < 0.2) {
		cd = 24/re; 
	} else if (re < 1000) {
		cd = 24/re + 6.3/pow(re,0.313);
	} else if (re < 2.0e5) {
		cd = 0.44;
	} else {
		cd = 0.1;
	}
	
	if (re_r < 1.0e-6) {
		continue; // In case umag is very, very small which would lead to problems
	} else if (re_r < 0.2) {
		cd_r = 24/re_r; 
	} else if (re_r < 1000) {
		cd_r = 24/re_r + 6.3/pow(re_r,0.313);
	} else if (re_r < 2.0e5) {
		cd_r = 0.44;
	} else {
		cd_r = 0.1;
	}
	
	cd_mult = size_mult*(cd_r/cd);

	aparticle = MY_PI*radius[i]*radius[i];
	aparticle_r = MY_PI*(radius[i]/size_mult)*(radius[i]/size_mult);
	fdmag = cd_mult*0.5*cd*aparticle*fluiddensity*umag*umag;
	fdmag_r = 0.5*cd_r*aparticle_r*fluiddensity*umag*umag;
	

	/* Apply this drag force in the appropriate direction (oppose motion
	   between particle and fluid, e.g., a particle in a stagnant fluid
	   will be brought to a stop while a particle in a moving fluid
	   will approach the fluid's velocity).*/
	f[i][0] += fdmag*relvx/umag;
	f[i][2] += fdmag*relvz/umag;

	ratio = (fdmag/fdmag_r);
	ratio_ratio = ratio / pow(size_mult,3);
	
	
	// Here everything is for vertical direction, above it was horizontal
	vfluid_ver = (inletpiperad*inletpiperad*vfluidver_0)/(basketrad*basketrad);
	relvy = vfluid_ver - v[i][1];
	umag_ver = fabs(relvy);
	re_ver = compoundreterm*umag_ver*radius[i];
	re_r_ver = compoundreterm*umag_ver*radius[i]/size_mult;
	if (re_ver < 1.0e-6) {
		continue; // In case umag is very, very small which would lead to problems
	} else if (re_ver < 0.2) {
		cd_ver = 24/re_ver; 
	} else if (re_ver < 1000) {
		cd_ver = 24/re_ver + 6.3/pow(re_ver,0.313);
	} else if (re_ver < 2.0e5) {
		cd_ver = 0.44;
	} else {
		cd_ver = 0.1;
	}
	if (re_r_ver < 1.0e-6) {
		continue; // In case umag is very, very small which would lead to problems
	} else if (re_r_ver < 0.2) {
		cd_r_ver = 24/re_r_ver; 
	} else if (re_r_ver < 1000) {
		cd_r_ver = 24/re_r_ver + 6.3/pow(re_r_ver,0.313);
	} else if (re_r_ver < 2.0e5) {
		cd_r_ver = 0.44;
	} else {
		cd_r_ver = 0.1;
	}
	cd_mult_ver = size_mult*(cd_r_ver/cd_ver);
	fdmag_ver = cd_mult_ver*0.5*cd_ver*aparticle*fluiddensity*umag_ver*umag_ver;
	fdmag_r_ver = 0.5*cd_r_ver*aparticle_r*fluiddensity*umag_ver*umag_ver;
	f[i][1] += fdmag_ver*relvy/umag_ver;
	ratio_ver = (fdmag_ver/fdmag_r_ver);
	ratio_ratio_ver = ratio_ver / pow(size_mult,3);
	
	// For printing to the log file
	if (comm->me == 0 && logfile) {
		fprintf(logfile, "Particle Reynolds Number - Scaled, Horizontal (re) %1.10e \n Particle Reynolds Number - Non-Scaled, Horizontal (re_r) %1.10e \n Drag Coefficient - Scaled, Horizontal (cd) %1.1e \n Drag Coefficient - Non-Scaled, Horizontal (cd_r) %1.1e \n Drag Coefficient Multiplier - Horizontal (cd_mult) %1.1e \n Drag Force - Scaled, Horizontal (fdmag) %1.10e \n Drag Force - Non-Scaled, Horizontal (fdmag_r) %1.10e \n Ratio of Drag Force Scaled to Drag Force Non-Scaled - Horizontal (ratio) %1.10e \n Ratio of Drag Force Ratio to Size Multiplier's Cube - Horizontal (ratio_ratio) %1.10e \n", re, re_r, cd, cd_r, cd_mult, fdmag, fdmag_r,ratio, ratio_ratio); 
		fprintf(logfile, "Particle Reynolds Number - Scaled, Vertical (re_ver) %1.10e \n Particle Reynolds Number - Non-Scaled, Vertical (re_r_ver) %1.10e \n Drag Coefficient - Scaled, Vertical (cd_ver) %1.1e \n Drag Coefficient - Non-Scaled, Vertical (cd_r_ver) %1.1e \n Drag Coefficient Multiplier - Vertical (cd_mult_ver) %1.1e \n Drag Force - Scaled, Vertical (fdmag_ver) %1.10e \n Drag Force - Non-Scaled, Vertical (fdmag_r_ver) %1.10e \n Ratio of Drag Force Scaled to Drag Force Non-Scaled - Vertical (ratio) %1.10e \n Ratio of Drag Force Ratio to Size Multiplier's Cube - Vertical (ratio_ratio_ver) %1.10e \n", re_ver, re_r_ver, cd_ver, cd_r_ver, cd_mult_ver, fdmag_ver, fdmag_r_ver, ratio_ver, ratio_ratio_ver); 
	}
	
	// For printing to the screen
	//fprintf(logfile, "Particle Reynolds Number - Scaled, Horizontal (re) %1.10e \n Particle Reynolds Number - Non-Scaled, Horizontal (re_r) %1.10e \n Drag Coefficient - Scaled, Horizontal (cd) %1.1e \n Drag Coefficient - Non-Scaled, Horizontal (cd_r) %1.1e \n Drag Coefficient Multiplier - Horizontal (cd_mult) %1.1e \n Drag Force - Scaled, Horizontal (fdmag) %1.10e \n Drag Force - Non-Scaled, Horizontal (fdmag_r) %1.10e \n Ratio of Drag Force Scaled to Drag Force Non-Scaled - Horizontal (ratio) %1.10e \n Ratio of Drag Force Ratio to Size Multiplier's Cube - Horizontal (ratio_ratio) %1.10e \n", re, re_r, cd, cd_r, cd_mult, fdmag, fdmag_r,ratio, ratio_ratio); 
	
	//fprintf(screen, "Particle Reynolds Number - Scaled, Vertical (re_ver) %1.10e \n Particle Reynolds Number - Non-Scaled, Vertical (re_r_ver) %1.10e \n Drag Coefficient - Scaled, Vertical (cd_ver) %1.1e \n Drag Coefficient - Non-Scaled, Vertical (cd_r_ver) %1.1e \n Drag Coefficient Multiplier - Vertical (cd_mult_ver) %1.1e \n Drag Force - Scaled, Vertical (fdmag_ver) %1.10e \n Drag Force - Non-Scaled, Vertical (fdmag_r_ver) %1.10e \n Ratio of Drag Force Scaled to Drag Force Non-Scaled - Vertical (ratio) %1.10e \n Ratio of Drag Force Ratio to Size Multiplier's Cube - Vertical (ratio_ratio_ver) %1.10e \n", re_ver, re_r_ver, cd_ver, cd_r_ver, cd_mult_ver, fdmag_ver, fdmag_r_ver, ratio_ver, ratio_ratio_ver);  %1.10e \n", re, re_r, cd, cd_r, cd_mult, fdmag, fdmag_r,ratio, ratio_ratio); 
     }
   }
 }
/* ---------------------------------------------------------------------- */

void FixFluiddrag::post_force_respa(int vflag, int ilevel, int iloop)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFluiddrag::min_post_force(int vflag)
{
  post_force(vflag);
}
