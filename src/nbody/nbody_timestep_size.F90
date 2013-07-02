! NBODY_TIMESTEP_SIZE.F90
! D. A. Hubber - 13/9/2010
! Calculate appropriate N-body timestep for star s.  
! Currently uses Aarseth timestep for 4th-order Hermite scheme.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_timestep_size(s,dt)
  use Nbody_module
  implicit none

  integer, intent(in) :: s             ! Star id
  real(kind=DP), intent(out) :: dt     ! Physical timestep

  real(kind=DP) :: asqd                ! Acceleration squared
  real(kind=DP) :: a1sqd               ! jerk squared
  real(kind=DP) :: a2sqd               ! 2nd derivative squared
  real(kind=DP) :: a3sqd               ! 3rd derivative squared
  real(kind=DP) :: dt2                 ! dt*dt


! Timestep for 4-th order Hermite scheme
! ----------------------------------------------------------------------------
#if defined(NBODY_HERMITE4)
  asqd  = real(dot_product(star(s)%a,star(s)%a),DP)
  a1sqd = real(dot_product(star(s)%adot,star(s)%adot),DP)
  a2sqd = real(dot_product(star(s)%a2dot,star(s)%a2dot),DP)
  a3sqd = real(dot_product(star(s)%a3dot,star(s)%a3dot),DP)
  dt2   = (sqrt(asqd*a2sqd) + a1sqd) / (sqrt(a1sqd*a3sqd) + a2sqd)
  dt    = nbody_timemult*sqrt(dt2)

#if defined(DEBUG_TIMESTEP_SIZE)
  write(6,*) "Timestep for star ",s,dt,dt2
  write(6,*) "asqd : ",asqd,a1sqd,a2sqd,a3sqd
#endif

! ----------------------------------------------------------------------------
#elif defined(NBODY_HERMITE6)

#endif
! ----------------------------------------------------------------------------

  return
END SUBROUTINE nbody_timestep_size
