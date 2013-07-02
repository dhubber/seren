! COPY_SINKS_TO_STARS.F90
! D. A. Hubber - 19/6/2008
! Copy sink data array to star data array either for determining binary 
! properties, or switching from the SPH to the N-body simulation. 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE copy_sinks_to_stars
  use sink_module
  use Nbody_module
  implicit none

  integer :: s            ! sink/star counter

  debug2("Copying sink data to star arrays [copy_sinks_to_stars.F90]")

  do s=1,stot  
     star(s)%ncreate      = sink(s)%ncreate
     star(s)%tcreate      = sink(s)%tcreate
     star(s)%r(1:NDIM)    = real(sink(s)%r(1:NDIM),DP)
     star(s)%v(1:NDIM)    = real(sink(s)%v(1:NDIM),DP)
     star(s)%m            = real(sink(s)%m,DP)
     star(s)%h            = real(sink(s)%h,DP)
     star(s)%radius       = real(sink(s)%radius,DP)
     star(s)%a(1:NDIM)    = real(sink(s)%a(1:NDIM),DP)
     star(s)%rold(1:NDIM) = real(sink(s)%vold(1:NDIM),DP)
     star(s)%vold(1:NDIM) = real(sink(s)%vold(1:NDIM),DP)
     star(s)%angmom(1:3)  = real(sink(s)%angmom(1:3),DP)
     star(s)%gpe          = real(sink(s)%gpe,DP)
     star(s)%dmdt         = 0.0_DP
  end do

  return
END SUBROUTINE copy_sinks_to_stars
