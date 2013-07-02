! HEAPSORT.F90
! D. A. Hubber & A. P. Whitworth - 26/6/2007
! Sorts a list of real variables into ascending order using heapsort.  
! Accepts the list of real variables to be sorted (rarray) and the list 
! of corresponding ids (iarray).  Returns the list of ids in order, but 
! leaves the array of reals unchanged.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE heapsort_real(nsort,rarray,iarray)
  use definitions
  implicit none

  integer, intent(in)  :: nsort                    ! No. of values to be sorted
  integer, intent(inout) :: iarray(1:nsort)        ! Ids to be sorted
  real(kind=PR), intent(inout) :: rarray(1:nsort)  ! Values to be sorted

  integer :: ibuf             ! ID-buffer for swapping
  integer :: s                ! location in binary heap
  integer :: snow             ! location being looked at now
  integer :: scom             ! location with which it is being compared
  integer :: itemp(1:nsort)   ! IDs of sorted values

  debug2("Sorting list of real variables using heapsort [heapsort_real.F90]")

  do s=1,nsort
     itemp(s) = nsort + 1 - s
  end do

! To build the binary heap
! ----------------------------------------------------------------------------
  do s=2,nsort
     snow = s
1    if (snow == 1) cycle
     scom = snow/2
     if (rarray(itemp(scom)) >= rarray(itemp(snow))) cycle
     ibuf = itemp(snow)
     itemp(snow) = itemp(scom)
     itemp(scom) = ibuf
     snow = scom
     goto 1
  end do
     
! To invert the binary heap
! ----------------------------------------------------------------------------
  do s=nsort,2,-1
     ibuf  = itemp(s)
     itemp(s) = itemp(1)
     itemp(1) = ibuf
     snow  = 1
2    scom  = 2*snow
     if (scom >= s) cycle
     if ((rarray(itemp(scom+1)) > rarray(itemp(scom))) .and. &
          & (scom + 1 < s)) scom = scom + 1
     if (rarray(itemp(scom)) <= rarray(itemp(snow))) cycle
     ibuf = itemp(snow)
     itemp(snow) = itemp(scom)
     itemp(scom) = ibuf
     snow = scom
     goto 2
  end do

! Change aux ids to particle ids
  do s=1,nsort
     itemp(s) = iarray(itemp(s))
  end do
  do s=1,nsort
     iarray(s) = itemp(s)
  end do
  

  return
END SUBROUTINE heapsort_real
