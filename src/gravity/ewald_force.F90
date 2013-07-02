! EWALD_FORCE.F90
! A. McLeod & D. A. Hubber - 21/01/2008
! Calculates periodic gravity force using Ewald method.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE ewald_force(dr,m,eaccel)
  use definitions
  use ewald_module
  implicit none

  real(kind=PR), intent(in)  :: dr(1:NDIM)     ! Relative displacement vector
  real(kind=PR), intent(in)  :: m              ! Mass of particle pp
  real(kind=PR), intent(out) :: eaccel(1:NDIM) ! Ewald grav acceleration

  integer :: i                 ! Loop variables and grid cell identifiers
  integer :: j
  integer :: k                 
  integer :: signx(1:NDIM)     ! Octant identifier
  real(kind=PR) :: dx(1:NDIM)  ! Position within a cell 0 > dx > 1
#if NDIM==1
  real(kind=PR) :: fint(1:2)   ! Linear interpolation weighting factor
#elif NDIM==2
  real(kind=PR) :: fint(1:4)   ! Bilinear interpolation weighting factor
#elif NDIM==3
  real(kind=PR) :: fint(1:8)   ! Trilinear interpolation weighting factor
#endif

! Transform dr into units of grid cell size
  dx(1:NDIM) = dr(1:NDIM)*ewsizeil(1:NDIM)      
  signx(1:NDIM) = -1

! Find octant of dr, and make dx positive
  do i=1,NDIM
     if (dx(i) < 0.0_PR) then
        dx(i) = -dx(i)                            
        signx(i) = 1
     end if
  end do

! Find grid cell (grid cells start at 1)
  i = min(int(dx(1)) + 1,ewsize(1) - 1)
#if NDIM==2 || NDIM==3
  j = min(int(dx(2)) + 1,ewsize(2) - 1)
#endif
#if NDIM==3
  k = min(int(dx(3)) + 1,ewsize(3) - 1)
#endif

! ----------------------------------------------------------------------------
#if NDIM==3
! Subtract number of grid cells to leave remainder
  dx = dx - real((/i,j,k/),PR) + 1.0_PR     

! Interpolate to get more accurate forces
  fint(1) = (1.0_PR - dx(1))*(1.0_PR - dx(2))*(1.0_PR - dx(3))
  fint(2) = (1.0_PR - dx(1))*(1.0_PR - dx(2))*(dx(3))
  fint(3) = (1.0_PR - dx(1))*(dx(2))*(1.0_PR - dx(3))
  fint(4) = (1.0_PR - dx(1))*(dx(2))*(dx(3))
  fint(5) = (dx(1))*(1.0_PR - dx(2))*(1.0_PR - dx(3))
  fint(6) = (dx(1))*(1.0_PR - dx(2))*(dx(3))
  fint(7) = (dx(1))*(dx(2))*(1.0_PR - dx(3))
  fint(8) = (dx(1))*(dx(2))*(dx(3))

! Add correction terms with interpolation
  eaccel(1:3) = m*signx(1:3)*&
       &(fint(1)*fcorr(1:3,i,j,k)     + fint(2)*fcorr(1:3,i,j,k+1)   + &
       & fint(3)*fcorr(1:3,i,j+1,k)   + fint(4)*fcorr(1:3,i,j+1,k+1) + &
       & fint(5)*fcorr(1:3,i+1,j,k)   + fint(6)*fcorr(1:3,i+1,j,k+1) + &
       & fint(7)*fcorr(1:3,i+1,j+1,k) + fint(8)*fcorr(1:3,i+1,j+1,k+1))

! ----------------------------------------------------------------------------
#elif NDIM==2
  dx = dx - real((/i,j/),PR)

! Interpolate to get more accurate forces
  fint(1) = (1.0_PR - dx(1))*(1.0_PR - dx(2))
  fint(2) = (1.0_PR - dx(1))*(dx(2))
  fint(3) = (dx(1))*(1.0_PR - dx(2))
  fint(4) = (dx(1))*(dx(2))

  eaccel(1:2) = m*signx(1:2)*&
       &(fint(1)*fcorr(1:2,i,j,1)   + fint(2)*fcorr(1:2,i,j+1,1) + &
       & fint(3)*fcorr(1:2,i+1,j,1) + fint(4)*fcorr(1:2,i+1,j+1,1))
#endif
! ----------------------------------------------------------------------------

  return
END SUBROUTINE ewald_force
