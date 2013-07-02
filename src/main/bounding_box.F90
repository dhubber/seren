! BOUNDING_BOX.F90
! C. P. Batty & D. A. Hubber - 15/11/2010
! Finds bounding box for particles between pstart and pend in main arrays.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE bounding_box(nptcls,r,r_max,r_min)
  use particle_module, only : PR
  use type_module
! #if defined(H_RHO)
!   use particle_module, only : rextent
! #endif
  implicit none

  integer, intent(in) :: nptcls                    ! no. of particles
  real(kind=PR), intent(in) :: r(1:NDIM,1:nptcls)  ! positions of particles
  real(kind=PR), intent(out) :: r_max(1:NDIM)      ! maximum r value
  real(kind=PR), intent(out) :: r_min(1:NDIM)      ! minimum r value

  integer :: k                                     ! dimension counter
  integer :: p                                     ! particle counter

  debug2("Finding bounding box for selected particles [bounding_box.F90]")

  r_max(1:NDIM) = -BIG_NUMBER
  r_min(1:NDIM) =  BIG_NUMBER

  do p=1,nptcls
     r_max = max(r_max,r(1:NDIM,p))
     r_min = min(r_min,r(1:NDIM,p))
  end do

  ! The following is a dirty side effect!
! ! Calculate domain extent for bisection method in grad-h method
! #if defined(H_RHO)
!   rextent = maxval(r_max(1:NDIM) - r_min(1:NDIM))
! #endif

  return
END SUBROUTINE bounding_box
