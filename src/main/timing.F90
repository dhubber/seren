! TIMING.F90
! D. A. Hubber - 14/2/2008
! Calculates time between current and last call to this routine. 
! Used for general purpose timing of code segments, as opposed to 
! using external programs such as gprof.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE timing(current_mark)
  use definitions
  use timing_module
#if defined(USE_MPI)
  use mpi
#endif
  implicit none

  character(len=*), intent(in) :: current_mark   ! Current mark string

  integer :: i                                   ! Aux. counter
  integer :: imark                               ! Current mark id
  integer (kind=ILP) :: current_itime            ! Current integer time
  integer (kind=ILP) :: dti                      ! Integer time interval
  real(kind=DP) :: current_rtime                 ! Current real time
  real(kind=DP) :: dtr                           ! Real time interval

! Loop through ids and find which marker this is.  
  imark = 0
  do i=1,mark_tot
     if (trim(current_mark) == trim(marker_id(i))) then
        imark = i
        exit
     end if
  end do

! If it doesn't exist, create a new spot in the timing tables
  if (imark == 0) then
     mark_tot = mark_tot + 1
     imark = mark_tot 
     marker_id(mark_tot) = trim(current_mark)
  end if

! Get current times using intrinsic Fortran functions
  call system_clock(current_itime)
#if defined(USE_MPI)
  current_rtime = MPI_WTIME()
#else
  call cpu_time(current_rtime)
#endif

! Calculate timings and store in arrays
! ----------------------------------------------------------------------------
  if (last_id > 0) then
 
     ! Calculate time since last timing
     dtr = current_rtime - last_rtime
     dti = current_itime - last_itime
     if (dtr < 0.0_DP) dtr = 0.0_DP
     if (dti < 0) dti = 0

     ! Record new accumulated time
     iblock(last_id) = iblock(last_id) + dti
     rblock(last_id) = rblock(last_id) + dtr

     ! Advance total time
     rtime = rtime + dtr
     itime = itime + dti

  end if
! ----------------------------------------------------------------------------
  
! Record last mark in order to record timings in next call
  last_id = imark
  last_itime = current_itime
  last_rtime = current_rtime

  return
END SUBROUTINE timing
