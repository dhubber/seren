! INSERTION_SORT.F90
! C. P. Batty & D. A. Hubber - 16/5/2006
! Orders a list of numbers in ascending order using insertion sort.  
! Should only be used for ordering small lists (N < 20). 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE insertion_sort_real(nsort,iarray,rarray)
  use definitions
  implicit none

  integer, intent(in) :: nsort                     ! No. of items to sort
  integer, intent(inout) :: iarray(1:nsort)        ! Array of ids
  real(kind=PR), intent(inout) :: rarray(1:nsort)  ! Array of reals to sort

  integer :: i                     ! Counter in sort
  integer :: iaux                  ! Original sort index
  integer :: j                     ! Counter in array
  real(kind=PR) :: raux            ! Aux. variable for sorting

! Loop through list and reorder integer ids  
  do j=2,nsort
     iaux = iarray(j)
     raux = rarray(j)
     do i=j-1,1,-1
        if (rarray(i) <= raux) exit
        rarray(i+1) = rarray(i)
        iarray(i+1) = iarray(i)
     end do
     rarray(i+1) = raux
     iarray(i+1) = iaux
  end do

  return
END SUBROUTINE insertion_sort_real



! ============================================================================
SUBROUTINE insertion_sort_dp(nsort,iarray,darray)
  use definitions
  implicit none

  integer, intent(in) :: nsort                     ! No. of items to sort  
  integer, intent(inout) :: iarray(1:nsort)        ! Array of ids
  real(kind=DP), intent(inout) :: darray(1:nsort)  ! Array of doubles to sort

  integer :: i                     ! Counter in sort
  integer :: iaux                  ! Original sort index
  integer :: j                     ! Counter in array
  real(kind=DP) :: daux            ! Aux. variable for sorting

! Loop through list and reorder integer ids 
  do j=2,nsort
     iaux = iarray(j)
     daux = darray(j) 
     do i=j-1,1,-1    
        if (darray(i) <= daux) exit
        darray(i+1) = darray(i)
        iarray(i+1) = iarray(i)
     end do
     darray(i+1) = daux 
     iarray(i+1) = iaux
  end do

  return
END SUBROUTINE insertion_sort_dp



! ============================================================================
SUBROUTINE insertion_sort_int(nsort,iarray)
  implicit none

  integer, intent(in)    :: nsort            ! No. of elements to be sorted
  integer, intent(inout) :: iarray(1:nsort)  ! List of integers

  integer :: i                  ! Counter in sort
  integer :: iaux               ! Aux. variable for sorting
  integer :: j                  ! Counter in list

  do j=2,nsort
     iaux = iarray(j)
     do i=j-1,1,-1
        if (iarray(i) <= iaux) exit
        iarray(i+1) = iarray(i)
     end do
     iarray(i+1) = iaux
  end do

  return
END SUBROUTINE insertion_sort_int
