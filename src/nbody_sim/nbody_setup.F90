! NBODY_SETUP.F90
! D. A. Hubber - 11/2/2008
! Sets up N-body simulation from existing gas and sink particle data.
! Performs the following procedures
! 1 - Identify all particles which are bound to one sink only (e.g. such as 
!     those in disks) and accretes them to sinks.  
! 2 - Copies all sink data to star arrays
! 3 - Calculate acceleration, 1st, 2nd and 3rd time acceleration derivatives 
!     of all stars (Makino & Aarseth 1992) 
! 4 - Initialize some arrays values and time variables
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_setup
  use particle_module
  use sink_module
  use kernel_module
  use time_module
  use Nbody_module
  use diagnostics_module
  use type_module
  implicit none

  integer :: s                     ! Star counter

  debug1("Set-up N-body simulation [nbody_setup.F90]")

! Accrete SPH particles bound to sinks (e.g. from long-lived disks)
#if defined(SPH_SIMULATION) || defined(NBODY_SPH_SIMULATION)
  call nbody_accrete_bound_particles
#endif

! Copy all sink data to 'star' structure
  call copy_sinks_to_stars

! Calculate the gravitational acceleration and 'jerk' of stars
  star(1:SMAX)%accdo = .true.
  call nbody_grav_forces

! Now calculate higher (2nd and 3rd) time derivatives of stars (in order to 
! calculate initial timestep in N-body integration; ref : Aarseth 200?)
#if defined(NBODY_HERMITE4)
  call nbody_hermite4_extra_terms
#endif

! Initialise star particle info not contained in IC file
  do s=1,stot
     star(s)%rold(1:NDIM)  = star(s)%r(1:NDIM)
     star(s)%vold(1:NDIM)  = star(s)%v(1:NDIM)
     star(s)%a0(1:NDIM)    = star(s)%a(1:NDIM)
     star(s)%adot0(1:NDIM) = star(s)%adot(1:NDIM)
#if defined(BINARY_COM_MOTION)
     star(s)%binid = -1
#endif
  end do

! Set ptot to zero so SPH particles are ignored in diagnostics
  ptot = 0

#if defined(DEBUG_DIAGNOSTICS)
  call diagnostics
  etot0 = etot
#endif

! Initialize time variables
  nbody_lastsnap = time
  ntempnext      = nsteps + ntempstep
  nextsnap       = time + snaptime
  lastsnap       = time
  ntempnext      = nsteps + ntempstep
  ndiagnext      = nsteps + ndiagstep
  n              = 0
  nresync        = n

  do s=1,stot
     star(s)%nlast = n
     star(s)%nlevel = 0
  end do

  return
END SUBROUTINE nbody_setup
