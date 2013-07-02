! DIRECT_SPH_HERMITE4_GRAVITY.F90
! D. A. Hubber - 26/01/2012
! Calculate the gravitational acceleration of star s due to all other 
! stars in the simulation.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE direct_sph_hermite4_gravity(s,hs,rs,vs,agravs,adots,potp)
  use interface_module, only : gravity_hermite4
  use definitions
  use particle_module
  use Nbody_module
  implicit none

  integer, intent(in) :: s                     ! Star i.d.
  real(kind=DP), intent(in) :: hs              ! Smoothing length of star s
  real(kind=DP), intent(in) :: rs(1:NDIM)      ! Position of star s
  real(kind=DP), intent(in) :: vs(1:NDIM)      ! Velocity of star s
  real(kind=DP), intent(out) :: agravs(1:NDIM) ! Grav. accel of star s
  real(kind=DP), intent(out) :: adots(1:NDIM)  ! Jerk of star s
  real(kind=DP), intent(out) :: potp           ! Potential of star s
  
  integer :: p                        ! Particle counter
  real(kind=DP) :: adottemp(1:NDIM)   ! Aux. variable for calculating jerk
  real(kind=DP) :: atemp(1:NDIM)      ! Auxilary accel variable 
  real(kind=DP) :: dpotp              ! Aux. star pot variable

  debug3("[nbody_gravity.F90] : ",s)

! Zero arrays
  adots(1:NDIM)  = 0.0_DP
  agravs(1:NDIM) = 0.0_DP
  potp           = 0.0_DP


! Loop over all other stars and calculate net gravitational acceleration
! ----------------------------------------------------------------------------
  do p=1,ptot
#if defined(NBODY_HERMITE4) && defined(MEANH_GRAVITY)
     call gravity_hermite4_meanh(0.5_DP*(hs + real(sph(p)%h,DP)),&
          &real(sph(p)%m,DP),rs(1:NDIM),real(sph(p)%r(1:NDIM),DP),vs(1:NDIM),&
          &real(sph(p)%v(1:NDIM),DP),atemp(1:NDIM),adottemp(1:NDIM),dpotp)
#elif defined(NBODY_HERMITE4)
     call gravity_hermite4(1.0_DP/hs,real(sph(p)%h,DP),real(sph(p)%m,DP),&
          &rs(1:NDIM),real(sph(p)%r(1:NDIM),DP),vs(1:NDIM),&
          &real(sph(p)%v(1:NDIM),DP),atemp(1:NDIM),adottemp(1:NDIM),dpotp)
#endif
     agravs(1:NDIM) = agravs(1:NDIM) + atemp(1:NDIM)
     adots(1:NDIM) = adots(1:NDIM) + adottemp(1:NDIM)
     potp = potp + dpotp
     
#if defined(DEBUG_HERMITE4)
     write(6,*) "star :",s,p,atemp(1:NDIM),adots(1:NDIM),potp
#endif

  end do
! ----------------------------------------------------------------------------

  return
END SUBROUTINE direct_sph_hermite4_gravity
