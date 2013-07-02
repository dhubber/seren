! GRAVITY_NBODY.F90
! C. P. Batty & D. A. Hubber - 24/6/2007
! Calculates gravitational force between masses mp and mpp at positions 
! rp and rpp.  Passes positions, masses and smoothing lengths rather 
! than particle or sink ids. 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE gravity_nbody(mpp,rp,rpp,atemp,dpotp)
  use interface_module, only : distance3,ewald_force
  use definitions
  use kernel_module
  implicit none

  real(kind=PR), intent(in)  :: mpp            ! Mass of pp
  real(kind=PR), intent(in)  :: rp(1:NDIM)     ! Position of p
  real(kind=PR), intent(in)  :: rpp(1:NDIM)    ! Position of pp
  real(kind=PR), intent(out) :: atemp(1:NDIM)  ! Acceleration vector
  real(kind=PR), intent(out) :: dpotp          ! Gravitational potential

  real(kind=PR) :: dr(1:NDIM)         ! Relative position vector
  real(kind=PR) :: drmag              ! Distance
  real(kind=PR) :: drsqd              ! Distance squared
  real(kind=PR) :: grav               ! aux variable
  real(kind=PR) :: invdrmag           ! 1.0 / drmag
#if defined(EWALD)
  real(kind=PR) :: eaccel(1:NDIM)     ! Ewald grav. acceleration
#endif

! Calculate relative displacement vector between bodies
#if defined(EWALD)
  call distance3(rp(1:NDIM),rpp(1:NDIM),dr(1:NDIM),drsqd)
#else
  dr(1:NDIM) = rpp(1:NDIM) - rp(1:NDIM)
  drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif
  drmag    = sqrt(drsqd) + SMALL_NUMBER
  invdrmag = 1.0_PR / drmag
  grav     = mpp*invdrmag**3
  dpotp    = mpp*invdrmag
  atemp(1:NDIM) = grav*dr(1:NDIM)
  
! Add Ewald periodic correction if required
#if defined(EWALD)
  call ewald_force(dr(1:NDIM),mpp,eaccel(1:NDIM))
  atemp(1:NDIM) = atemp(1:NDIM) + eaccel(1:NDIM)
#endif

  return
END SUBROUTINE gravity_nbody
