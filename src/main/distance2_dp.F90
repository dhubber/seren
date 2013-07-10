! DISTANCE2_DP.F90
! C. P. Batty & D. A. Hubber - 11/12/2006
! Calculates the relative position vector between two particles (pp-p), for 
! any dimensionality and when periodic boundary conditions are employed. 
! Also returns the relative distance squared.  In this subroutine, the 
! position of the first particle p, i.e. rp(1:NDIM), and the id of the 
! second particle, pp, are passed.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE distance2_dp(rp,pp,dr,drsqd)
  use definitions
  use particle_module, only : sph
#if defined(PERIODIC)
  use periodic_module, only : periodic_half, periodic_size
#endif
  implicit none

  integer, intent(in) :: pp                  ! id of second particle, pp
  real(kind=DP), intent(in)  :: rp(1:NDIM)   ! Position of particle p
  real(kind=DP), intent(out) :: dr(1:NDIM)   ! Vector displacements (pp-p)
  real(kind=DP), intent(out) :: drsqd        ! Separation squared

! First compute x-direction (always needed)
  dr(1) = real(sph(pp)%r(1),DP) - rp(1)
#if defined(PERIODIC_X) && !defined(GHOST_PARTICLES)
  if (dr(1) >  real(periodic_half(1),DP)) &
       &dr(1) = dr(1) - real(periodic_size(1),DP)
  if (dr(1) < -real(periodic_half(1),DP)) &
       &dr(1) = dr(1) + real(periodic_size(1),DP)
#endif

! Compute y-dimension (for 2-D and 3-D cases)
#if NDIM==2 || NDIM==3
  dr(2) = real(sph(pp)%r(2),DP) - rp(2)
#if defined(PERIODIC_Y) && !defined(GHOST_PARTICLES)
  if (dr(2) >  real(periodic_half(2),DP)) &
       &dr(2) = dr(2) - real(periodic_size(2),DP)
  if (dr(2) < -real(periodic_half(2),DP)) &
       &dr(2) = dr(2) + real(periodic_size(2),DP)
#endif
#endif

! Compute z-dimension (for 3-D case only)
#if NDIM==3
  dr(3) = real(sph(pp)%r(3),DP) - rp(3)
#if defined(PERIODIC_Z) && !defined(GHOST_PARTICLES)
  if (dr(3) >  real(periodic_half(3),DP)) &
       &dr(3) = dr(3) - real(periodic_size(3),DP)
  if (dr(3) < real(-periodic_half(3),DP)) &
       &dr(3) = dr(3) + real(periodic_size(3),DP)
#endif
#endif

! Compute distance squared
#if NDIM==1
  drsqd = dr(1)*dr(1)
#elif NDIM==2
  drsqd = dr(1)*dr(1) + dr(2)*dr(2)
#elif NDIM==3
  drsqd = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
#endif

  return
END SUBROUTINE distance2_dp
