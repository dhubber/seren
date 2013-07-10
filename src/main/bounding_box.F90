! BOUNDING_BOX.F90
! C. P. Batty & D. A. Hubber - 15/11/2010
! Finds bounding box for particles between pstart and pend in main arrays.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE bounding_box(nptcls,r,r_max,r_min)
  use particle_module, only : PR
  use type_module
  implicit none

  integer, intent(in) :: nptcls                    ! No. of particles
  real(kind=PR), intent(in) :: r(1:NDIM,1:nptcls)  ! Positions of particles
  real(kind=PR), intent(out) :: r_max(1:NDIM)      ! Maximum r value
  real(kind=PR), intent(out) :: r_min(1:NDIM)      ! Minimum r value

  integer :: k                                     ! Dimension counter
  integer :: p                                     ! Particle counter

  debug2("Finding bounding box for selected particles [bounding_box.F90]")

  r_max(1:NDIM) = -BIG_NUMBER
  r_min(1:NDIM) =  BIG_NUMBER

  do p=1,nptcls
     r_max = max(r_max,r(1:NDIM,p))
     r_min = min(r_min,r(1:NDIM,p))
  end do

  return
END SUBROUTINE bounding_box
