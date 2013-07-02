! SINK_ADVANCE.F90
! C. P. Batty & D. A. Hubber - 23/8/2007
! Advance positions and velocities of sink particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sink_advance
  implicit none

  debug2("Advancing positions of all sink particles [sink_advance.F90]")
  debug_timing("ADVANCE")

#if defined(EULER)
  call advance_sink_euler
#elif defined(RUNGE_KUTTA2)
  call advance_sink_RK
#elif defined(LEAPFROG_KDK)
  call advance_sink_leapfrog_kdk
#elif defined(LEAPFROG_DKD)
  call advance_sink_leapfrog_dkd
#endif

#if defined(USE_MPI)
  call sink_share
#endif

  return
END SUBROUTINE sink_advance
