! GRAVITY_GRADH_MEANH.F90
! C. P. Batty & D. A. Hubber - 24/6/2007
! Calculates gravitational force between masses mp and mpp at positions 
! rp and rpp.  Passes positions, masses and smoothing lengths rather 
! than particle or sink ids. 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE gravity_gradh_meanh(hp,hpp,mpp,rp,rpp,zo_p,zo_pp,atemp,dpotp)
  use interface_module, only : distance3,ewald_force,wgrav,wpot,w2
  use definitions
  use kernel_module
  implicit none

  real(kind=PR), intent(in)  :: hp             ! Smoothing length of p
  real(kind=PR), intent(in)  :: hpp            ! Smoothing length of pp
  real(kind=PR), intent(in)  :: mpp            ! Mass of pp
  real(kind=PR), intent(in)  :: rp(1:NDIM)     ! Position of p
  real(kind=PR), intent(in)  :: rpp(1:NDIM)    ! Position of pp
  real(kind=PR), intent(in)  :: zo_p           ! Zeta/Omega for p
  real(kind=PR), intent(in)  :: zo_pp          ! Zeta/Omega for pp
  real(kind=PR), intent(out) :: atemp(1:NDIM)  ! Acceleration vector
  real(kind=PR), intent(out) :: dpotp          ! Gravitational potential

  real(kind=PR) :: dr(1:NDIM)        ! Relative position vector
  real(kind=PR) :: drmag             ! Distance
  real(kind=PR) :: drsqd             ! Distance squared
  real(kind=PR) :: grav              ! aux variable
  real(kind=PR) :: hmean             ! Arithmetic mean smoothing length
  real(kind=PR) :: invdrmag          ! 1.0 / drmag
  real(kind=PR) :: invh              ! 1.0 / h
#if defined(EWALD)
  real(kind=PR) :: eaccel(1:NDIM)    ! Ewald grav. acceleration
#endif

! Calculate relative displacement vector between bodies
#if defined(EWALD)
  call distance3(rp(1:NDIM),rpp(1:NDIM),dr(1:NDIM),drsqd)
#else
  dr(1:NDIM) = rpp(1:NDIM) - rp(1:NDIM)
  drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif
  drmag = sqrt(drsqd) + SMALL_NUMBER
  invdrmag = 1.0_PR / drmag
  hmean = 0.5_PR*(hp + hpp)

! Calculate contribution of gravity from particle p
  if (drmag < KERNRANGE*hmean) then
     invh = 1.0_PR / hmean
     grav = invh*invh*wgrav(drmag*invh)
     dpotp = invh*wpot(drmag*invh)
  else
     grav = invdrmag*invdrmag
     dpotp = invdrmag
  end if

! Calculate meanh grad-h correction term for grav. force for p
  if (drmag < KERNRANGE*hp) then
     invh = 1.0_PR / hp
     grav = grav + 0.5_PR*zo_p*w2(drmag*invh)*invh**4
  end if

! Calculate meanh grad-h correction term for grav. force for pp
  if (drmag < KERNRANGE*hpp) then
     invh = 1.0_PR / hpp
     grav = grav + 0.5_PR*zo_pp*w2(drmag*invh)*invh**4
  end if

! Record acceleration in output vector      
  atemp(1:NDIM) = mpp*invdrmag*grav*dr(1:NDIM)
  dpotp = mpp*dpotp
  
! Add Ewald periodic correction if required
#if defined(EWALD)
  call ewald_force(dr(1:NDIM),mpp,eaccel(1:NDIM))
  atemp(1:NDIM) = atemp(1:NDIM) + eaccel(1:NDIM)
#endif

  return
END SUBROUTINE gravity_gradh_meanh
