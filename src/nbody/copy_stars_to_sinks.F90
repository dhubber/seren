! COPY_STARS_TO_SINKS.F90
! D. A. Hubber - 19/6/2008
! Copy star data array to sink data array for ease of outputting to file.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE copy_stars_to_sinks
  use sink_module
  use Nbody_module
  implicit none

  integer :: s            ! sink/star counter

  debug2("Copying sink data to star arrays [copy_sinks_to_stars.F90]")

  do s=1,stot  
     sink(s)%r(1:NDIM)    = real(star(s)%r(1:NDIM),PR)
     sink(s)%v(1:NDIM)    = real(star(s)%v(1:NDIM),PR)
     sink(s)%m            = real(star(s)%m,PR)
     sink(s)%h            = real(star(s)%h,PR)
     sink(s)%radius       = real(star(s)%radius,PR)
     sink(s)%a(1:NDIM)    = real(star(s)%a(1:NDIM),PR)
     sink(s)%rold(1:NDIM) = real(star(s)%vold(1:NDIM),PR)
     sink(s)%vold(1:NDIM) = real(star(s)%vold(1:NDIM),PR)
     sink(s)%angmom(1:3)  = real(star(s)%angmom(1:3),DP)
     sink(s)%gpe          = real(star(s)%gpe,PR)
     sink(s)%dmdt         = 0.0_DP
#if defined(DEBUG_FORCES)
     sink(s)%agrav(1:NDIM) = real(star(s)%a(1:NDIM),PR)
#endif
  end do

  return
END SUBROUTINE copy_stars_to_sinks
