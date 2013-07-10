! ADD_EXTERNAL_GRAVITATIONAL_FORCE.F90
! D. A. Hubber - 12/06/2010
! Add a specified external gravitational field to the simulation.  
! Field is selected by the Makefile variable EXTERNAL_FORCE.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE add_external_gravitational_force(rp,vp,atemp,adottemp,pottemp)
  use particle_module
  implicit none

  real(kind=PR), intent(in) :: rp(1:NDIM)         ! Position
  real(kind=PR), intent(in) :: vp(1:NDIM)         ! Velocity
  real(kind=DP), intent(out) :: atemp(1:NDIM)     ! Grav. acceleration
  real(kind=DP), intent(out) :: adottemp(1:NDIM)  ! Jerk (1st deriv of accel.)
  real(kind=DP), intent(out) :: pottemp           ! Grav. potential

#if defined(UDS_POTENTIAL) || defined(NFW1996_POTENTIAL) || defined(PLUMMER_POTENTIAL)
  real(kind=PR) :: drmag                          ! Distance from origin
  real(kind=PR) :: drsqd                          ! Distance squared
  real(kind=PR) :: invdrmag                       ! 1 / drmag
#endif
#if defined(PLUMMER_POTENTIAL)
  real(kind=PR), parameter :: mplummer = 833.0_PR ! Plummer mass
  real(kind=PR), parameter :: rplummer = 1.0_PR   ! Plummer radius
#endif


! Plummer potential field
! ----------------------------------------------------------------------------
#if defined(PLUMMER_POTENTIAL)
  drsqd = dot_product(rp(1:NDIM),rp(1:NDIM))
  atemp(1:NDIM) = -mplummer*rp(1:NDIM)*(drsqd + rplummer**2)**(-1.5_PR)
  adottemp(1:NDIM) = 3.0_PR*mplummer*(drsqd + rplummer**2)**(-2.5_PR)*&
       &dot_product(rp(1:NDIM),vp(1:NDIM))*rp(1:NDIM) - &
       &mplummer*(drsqd + rplummer**2)**(-1.5_PR)*vp(1:NDIM)
  pottemp = 2.0_PR*mplummer*(drsqd + rplummer**2)**(-0.5_PR)
#endif
! ----------------------------------------------------------------------------


! Uniform-density sphere potential field
! ----------------------------------------------------------------------------
#if defined(UDS_POTENTIAL)
  drsqd = dot_product(rp(1:NDIM),rp(1:NDIM))
  if (drsqd <= 1.0_PR) then
     atemp(1:NDIM) = -rp(1:NDIM)
     pottemp = 0.0_PR
  else
     drmag = sqrt(drsqd) + SMALL_NUMBER
     invdrmag = 1.0_PR / drmag
     atemp(1:NDIM) = -rp(1:NDIM)*(invdrmag**3)
     pottemp = -2.0_PR*invdrmag
  end if
#endif
! ----------------------------------------------------------------------------


! Navarro, Frenk & White (1996) dark-matter potential
! ----------------------------------------------------------------------------
#if defined(NFW1996_POTENTIAL)
  drsqd = dot_product(rp(1:NDIM),rp(1:NDIM))
  drmag = sqrt(drsqd) + SMALL_NUMBER
  invdrmag = 1.0_PR / drmag
  if (drmag <= 0.05_PR) then
     atemp(1:NDIM) = -0.9_PR*(0.5_PR - drmag*(2.0_PR/3.0_PR &
          & - drmag*(0.75_PR - 0.8_PR*drmag)))*rp(1:NDIM)*invdrmag
     pottemp = 0.0_PR
  else
     atemp(1:NDIM) = -0.9_PR*(log(1.0_PR + drmag) + 1.0_PR/(1.0_PR + drmag) &
          & - 1.0_PR)*rp(1:NDIM)*(invdrmag**3)
     pottemp = 0.0_PR
  end if
#endif
! ----------------------------------------------------------------------------


  return
END SUBROUTINE add_external_gravitational_force
