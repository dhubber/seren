! NBODY_SPH_SETUP.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Performs initialisation routines to prepare for hybrid N-body/SPH simulation.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_sph_setup
  use filename_module, only : restart
  use particle_module
  use time_module
  use diagnostics_module
  implicit none

  integer :: p

  debug1("Setting up the simulation [nbody_sph_setup.F90]")

  if (ptot == 0) stop 'No SPH particles present'

! Setting up for different particle types
  call types

! Set-up HEALPix rays
#if defined(HEALPIX)
  call initialize_HP_sources
#endif

! Copy all sink data to 'star' structure
  call copy_sinks_to_stars

! Initialize certain variables before first force calculation
  call initialize_nbody_sph_variables_1

! Build and stock trees for first time
  call tree_update(nbuild,nstock)

! Make initial guesses of h either using tree or a global estimate
#if defined(BH_TREE) && !defined(CONSTANT_H)
  if (.not. restart .and. ptot > 0) call BHhydro_hguess
#elif !defined(CONSTANT_H)
  if (.not. restart .and. ptot > 0) call h_guess
#endif

! Calculating initial SPH quantities (h, neibs, rho, etc.)
  call sph_update

! Calculate HEALPix quantities if activated.
#if defined(HEALPIX)
  call HP_update
#endif

! Initialize all thermal properties depending on options used
#if defined(HYDRO)
  call initialize_thermal_properties
#endif

! Updating tree properties
#if defined(BH_TREE) && defined(SELF_GRAVITY)
  call BHgrav_stock
#endif

! Calculate hydro forces on all SPH particles
#if defined(HYDRO)
  call sph_hydro_forces
#endif

! Calculate gravitational forces on all SPH particles
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  call sph_grav_forces
#endif

! Calculate gravitational forces on all sink particles
#if defined(FORCE_SPLITTING)
  call nbody_sph_star_split_forces
#else
  call nbody_sph_star_forces
#endif

! Now calculate higher (2nd and 3rd) time derivatives of stars (in order to 
! calculate initial timestep in N-body integration).
! (ref : Aarseth 200?)
#if defined(NBODY_HERMITE4)
  call nbody_hermite4_extra_terms
#endif

! Add wind contribution to acceleration here (for now)
#if defined(STELLAR_WIND)
  do p=1,ptot
     if (sph(p)%accdo) sph(p)%a(1:NDIM) = &
          &sph(p)%a(1:NDIM) + sph(p)%a_wind(1:NDIM)
  end do
#endif

! Initialize other key variables (after initial force calculation)
  call initialize_nbody_sph_variables_2
  call copy_stars_to_sinks

! Calculate initial diagnostic quantities
#if defined(DEBUG_DIAGNOSTICS)
  call diagnostics
  etot0 = etot
#endif

  return
END SUBROUTINE nbody_sph_setup
