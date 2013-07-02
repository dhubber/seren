! NBODY_INTEGRATE.F90
! D. A. Hubber - 13/9/2010
! Integration step
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_integrate
  use Nbody_module
  use time_module
  implicit none

  integer :: iteration            ! No. of correction iterations

  debug2("Performing next N-body integration step [nbody_integrate.F90]")

! Calculate new timesteps for all N-body particles
  call nbody_timesteps

! Advance time and integer times
  n      = n + 1
  nsteps = nsteps + 1
  time   = time + timestep

! Integration scheme to advance particle positions and velocities
  call nbody_advance

! Predict-evaluate-correct, i.e. P(EC)^n iteration loop
! ----------------------------------------------------------------------------
  do iteration=1,npec

     ! Compute force and force derivatives at beginning of timestep
     call nbody_grav_forces
     
     ! Calculate higher-order derivatives and correction forces 
     ! for previous steps
     call nbody_correction_terms

  end do
! ----------------------------------------------------------------------------

! Once iterations are finished, set all star properties for end of step
  call nbody_end_step

! Output diagnostics and/or snapshot files
  call nbody_output

  return
END SUBROUTINE nbody_integrate
