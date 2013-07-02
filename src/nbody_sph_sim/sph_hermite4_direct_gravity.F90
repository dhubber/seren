! SPH_HERMITE4_DIRECT_GRAVITY.F90
! D. A. Hubber - ..
! Calculate gravitational forces on star particles due to all SPH particles 
! using the mean-h Hermite scheme.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_hermite4_direct_gravity(s,hs,rs,vs,agrav,adot,pot)
  use interface_module, only : gravity_hermite4_meanh
  use definitions
  use type_module, only : pgravitystart, pgravityend
  use particle_module, only : ptot,sph
  use sink_module, only : stot
  implicit none

  integer, intent(in) :: s                      ! Id of current particle
  real(kind=DP), intent(in) :: hs               ! Smoothing length of s
  real(kind=DP), intent(in) :: rs(1:NDIM)       ! Position of star s
  real(kind=DP), intent(in) :: vs(1:NDIM)       ! Velocity of star s
  real(kind=DP), intent(out) :: agrav(1:NDIM)   ! Gravitational acceleration
  real(kind=DP), intent(out) :: adot(1:NDIM)    ! ..
  real(kind=DP), intent(out) :: pot             ! Gravitational potential

  integer       :: pp                 ! particle identifier
  real(kind=DP) :: adottemp(1:NDIM)   ! Aux. variable for calculating jerk
  real(kind=DP) :: atemp(1:NDIM)      ! Auxilary accel variable 
  real(kind=DP) :: dpotp              ! Aux. star pot variable
  real(kind=DP) :: hmean              ! Mean h

  debug3("Calculating gravitational force [nbody_sph_star_direct_gravity.F90] for ",s)

! Zero arrays
  adot(1:NDIM)  = 0.0_DP
  agrav(1:NDIM) = 0.0_DP
  pot           = 0.0_DP

! Loop over all other stars and calculate net gravitational acceleration
! ----------------------------------------------------------------------------
  do pp=1,ptot
#if defined(NBODY_HERMITE4) && defined(GRAD_H_SPH) && defined(MEANH_GRAVITY)
     call gravity_hermite4_gradh_meanh(hs,real(sph(pp)%h,DP),&
          &real(sph(pp)%m,DP),rs(1:NDIM),&
          &real(sph(pp)%r(1:NDIM),DP),vs(1:NDIM),&
          &real(sph(pp)%v(1:NDIM),DP),0.0_PR,real(sph(pp)%zo,DP),&
          &atemp(1:NDIM),adottemp(1:NDIM),dpotp)
#elif defined(NBODY_HERMITE4) && defined(MEANH_GRAVITY)
     call gravity_hermite4_meanh(0.5_DP*(hs + real(sph(pp)%h,DP)),&
          &real(sph(pp)%m,DP),rs(1:NDIM),&
          &real(sph(pp)%r(1:NDIM),DP),vs(1:NDIM),&
          &real(sph(pp)%v(1:NDIM),DP),atemp(1:NDIM),&
          &adottemp(1:NDIM),dpotp)
#elif defined(NBODY_HERMITE4)
     call gravity_hermite4(1.0_DP/hs,real(sph(pp)%h,DP),&
          &real(sph(pp)%m,DP),rs(1:NDIM),&
          &real(sph(pp)%r(1:NDIM),DP),vs(1:NDIM),&
          &real(sph(pp)%v(1:NDIM),DP),atemp(1:NDIM),&
          &adottemp(1:NDIM),dpotp)
#endif
     agrav(1:NDIM) = agrav(1:NDIM) + atemp(1:NDIM)
     !adot(1:NDIM) = adot(1:NDIM) + adottemp(1:NDIM)
     pot = pot + dpotp
     
#if defined(DEBUG_HERMITE4)
     write(6,*) "star :",s,pp,atemp(1:NDIM),adottemp(1:NDIM),dpotp
#endif

  end do
! ----------------------------------------------------------------------------


  return
END SUBROUTINE sph_hermite4_direct_gravity
