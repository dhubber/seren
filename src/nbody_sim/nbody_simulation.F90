! NBODY_SIMULATION.F90
! D. A. Hubber - 23/7/2010
! Main control routine for N-body simulation.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_simulation
  use definitions
  use seren_sim_module
  use filename_module
  use particle_module, only : ptot
  use neighbour_module, only : pp_gather
  use time_module, only : time, endtime
#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  use sink_module, only : stot, sink_frac
  use Nbody_module, only : nbody_frac, nbody_endtime
#endif
  implicit none

  debug1("Main N-body simulation [nbody_simulation.F90]")

! Set flags for N-body simulation
  sph_sim       = .false.
  nbody_sph_sim = .false.
  nbody_sim     = .true.


! If there are more than 2 sinks, continue integration of sink orbits 
! using N-body integrator
! ============================================================================
  if (stot > 1 .and. time < nbody_endtime) then

     ! Set-up all N-body arrays, variables etc..
     call nbody_setup

     ! Main N-body integration loop
     ! -----------------------------------------------------------------------
     do 

        ! Exit main loop if simulation has reached endtime
        if (time >= nbody_endtime) exit

        ! Perform hermite integration step
        call nbody_integrate

     end do
     ! -----------------------------------------------------------------------

     debug1("N-body simulation completed")

  end if
! ============================================================================


  return
END SUBROUTINE nbody_simulation
