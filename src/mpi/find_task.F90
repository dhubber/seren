! FIND_TASK.F90
! A. McLeod - 23/10/2008
! Subroutine to rapidly find which task a point is in
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine find_task(rp,task)
    use definitions
    use mpi_decomposition_module

    real(kind=PR), intent(in) :: rp(1:NDIM)  ! Particle position
    integer, intent(out) :: task     ! Task number

    integer :: i                     ! ..
    integer :: l                     ! Where we are in the tree
    real(kind=PR) :: minr(1:NDIM)    ! Minimum of left-hand branch
    real(kind=PR) :: maxr(1:NDIM)    ! Maximum of left-hand branch

    i = 0

    ! ------------------------------------------------------------------------
    do l = 1,MPItreedepth

      ! Take d to next level
      i = i * 2

      ! Min and max of left hand branch
      minr(1:NDIM) = MPItreemin(l)%data(i)%xyz   
      maxr(1:NDIM) = MPItreemax(l)%data(i)%xyz

      if (any(rp < minr) .OR. any(rp > maxr)) then
        ! If any of these happen, particle must be in the right hand branch - i=i+1
        i = i + 1
        cycle
      end if

      ! Else particle must be in the left hand branch - no change to i (???)

    end do
    ! ------------------------------------------------------------------------

! Min and max of left hand branch
!     minr(1:NDIM) = MPItreemin(1:NDIM,MPItreedepth,i)   
!     maxr(1:NDIM) = MPItreemax(1:NDIM,MPItreedepth,i)

    ! (???)
    task = MPIrevgeometry(i)

    return
  end subroutine find_task
