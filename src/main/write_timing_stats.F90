! WRITE_TIMING_STATS.F90
! D. A. Hubber - 14/2/2008
! Writes out all timing statistics to file for quick analysis of code profile.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_timing_stats
   use definitions
   use timing_module
   use time_module
   use filename_module, only : run_dir,run_id
#if defined(USE_MPI)
   use mpi_communication_module
   use mpi
#endif
   implicit none

   integer :: i                    ! Aux. loop counters
   integer :: j                    ! Aux. loop counters
   integer, parameter :: iunit=1   ! Unit for writing output
   real(kind=DP) :: ipercentage    ! percentage of total integer time
   real(kind=DP) :: rpercentage    ! percentage of total real time
   character(len=256) :: filename  ! filename extension for data output
   character(len=256) :: fs                                ! Format string
#if defined(USE_MPI)
   integer            :: d                                 ! Domain counter
   integer            :: ierr                              ! MPI error value
   integer            :: imark                             ! Current mark id
   integer(kind=ILP)  :: tot_itime                         ! Total integer time taken
   real(kind=DP)      :: tot_rtime                         ! Total real time taken
   integer(kind=ILP)  :: temp_itime                        ! Integer time taken count
   real(kind=DP)      :: temp_rtime                        ! Real time taken count
   integer(kind=ILP), allocatable :: remote_nhydrocomp(:)  ! Number of hydro calculations
#if defined(GRAVITY)
   integer(kind=ILP), allocatable :: remote_ngravcomp(:)   ! Number of hydro calculations
#endif
   integer, allocatable :: remote_mark_tot(:)              ! Total number of markers
   integer(kind=ILP), allocatable :: remote_iblock(:,:)    ! Integer time taken by each block
   real(kind=DP), allocatable :: remote_rblock(:,:)        ! Real time taken by each block
   character :: gather_marker_id(50)     ! Character array for gathering marker strings
   character :: recv_marker_id(50,0:lastrank) ! Character array for receiving marker strings
   character(len=50), allocatable :: remote_marker_id(:,:) ! Array of marker strings
#endif

   debug2("Write timing information to file [write_timing_stats.F90]")

#if defined(USE_MPI)

   if (rank==0) then
      ! Collect all timing statistics together
      allocate(remote_nhydrocomp(-1:lastrank))
#if defined(GRAVITY)
      allocate(remote_ngravcomp(-1:lastrank))
#endif
      allocate(remote_mark_tot(-1:lastrank))
      allocate(remote_iblock(1:NBLOCKS,-1:lastrank))
      allocate(remote_rblock(1:NBLOCKS,-1:lastrank))
      allocate(remote_marker_id(1:NBLOCKS,-1:lastrank))
      remote_mark_tot(-1) = 0
      remote_iblock(1:NBLOCKS,-1) = 0_ILP
      remote_rblock(1:NBLOCKS,-1) = 0._DP
      WAIT_TIME_MACRO
      call MPI_GATHER(nhydrocomp,1,MPI_INTEGER_ILP,remote_nhydrocomp(0),&
         &1,MPI_INTEGER_ILP,0,MPI_COMM_WORLD,ierr)
#if defined(GRAVITY)
      call MPI_GATHER(ngravcomp,1,MPI_INTEGER_ILP,remote_nhydrocomp(0),&
         &1,MPI_INTEGER_ILP,0,MPI_COMM_WORLD,ierr)
#endif
      call MPI_GATHER(mark_tot,1,MPI_INTEGER,remote_mark_tot(0),&
         &1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(iblock,NBLOCKS,MPI_INTEGER_ILP,remote_iblock(1,0),&
         &NBLOCKS,MPI_INTEGER_ILP,0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(rblock,NBLOCKS,MPI_DOUBLE_PRECISION,remote_rblock(1,0),&
         &NBLOCKS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      do i=1,NBLOCKS
         do j=1,50
            gather_marker_id(j) = marker_id(i)(j:j)
         end do
         call MPI_GATHER(gather_marker_id,50,MPI_CHARACTER,recv_marker_id(1,0),&
            &50,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
         do d=0,lastrank
            do j=1,50
               remote_marker_id(i,d)(j:j) = recv_marker_id(j,d)
            end do
         end do
      end do
      tot_rtime = rtime
      call MPI_REDUCE(MPI_IN_PLACE, tot_rtime, 1, MPI_DOUBLE_PRECISION,&
         &MPI_SUM, 0, MPI_COMM_WORLD,ierr)
      tot_itime = itime
      call MPI_REDUCE(MPI_IN_PLACE, tot_itime, 1, MPI_INTEGER,&
         &MPI_SUM, 0, MPI_COMM_WORLD,ierr)
      CALC_TIME_MACRO
   else
      ! Participate to help the root collect timings
      WAIT_TIME_MACRO
      call MPI_GATHER(nhydrocomp,1,MPI_INTEGER_ILP,0,0,0,0,MPI_COMM_WORLD,ierr)
#if defined(GRAVITY)
      call MPI_GATHER(ngravcomp,1,MPI_INTEGER_ILP,0,0,0,0,MPI_COMM_WORLD,ierr)
#endif
      call MPI_GATHER(mark_tot,1,MPI_INTEGER,0,0,0,0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(iblock,NBLOCKS,MPI_INTEGER_ILP,0,0,0,0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(rblock,NBLOCKS,MPI_DOUBLE_PRECISION,0,0,0,0,MPI_COMM_WORLD,ierr)
      do i=1,NBLOCKS
         do j=1,50
            gather_marker_id(j) = marker_id(i)(j:j)
         end do
         call MPI_GATHER(gather_marker_id,50,MPI_CHARACTER,0,0,0,0,MPI_COMM_WORLD,ierr)
      end do
      tot_rtime = rtime
      call MPI_REDUCE(tot_rtime, 0, 1, MPI_DOUBLE_PRECISION,&
         &MPI_SUM, 0, MPI_COMM_WORLD,ierr)
      tot_itime = itime
      call MPI_REDUCE(tot_itime, 0, 1, MPI_INTEGER,&
         &MPI_SUM, 0, MPI_COMM_WORLD,ierr)
      CALC_TIME_MACRO
   end if

   if (rank/=0) return
  
   remote_nhydrocomp(-1) = sum(remote_nhydrocomp(0:lastrank))
#if defined(GRAVITY)
   remote_ngravcomp(-1) = sum(remote_ngravcomp(0:lastrank))
#endif
   do d=0,lastrank
      do j=1,remote_mark_tot(d)
         imark = 0
         do i=1,remote_mark_tot(-1)
            if (trim(remote_marker_id(j,d)) == trim(remote_marker_id(i,-1))) then
               imark = i
               exit
            end if
         end do
         if (imark == 0) then
            remote_mark_tot(-1) = remote_mark_tot(-1) + 1
            imark = remote_mark_tot(-1)
            remote_marker_id(imark,-1) = trim(remote_marker_id(j,d))
         end if
         remote_iblock(imark,-1) = remote_iblock(imark,-1) + remote_iblock(j,d)
         remote_rblock(imark,-1) = remote_rblock(imark,-1) + remote_rblock(j,d)
      end do
   end do

! Loop over all blocks and print statistics
! ----------------------------------------------------------------------------
  if (sum(remote_mark_tot) > 0) then

     call reorder_markers(remote_mark_tot(-1), remote_iblock(1:mark_tot,-1),&
        &remote_rblock(1:mark_tot,-1), remote_marker_id(1:mark_tot,-1))
           
     filename = trim(adjustl(run_dir))//trim(adjustl(run_id))//".timing"
     open(unit=iunit,file=filename,status="unknown")

     fs = "(1X,a40,2X,i14,2X,1F8.4,2X,1G14.6,2X,1F8.4)"
     write(iunit,*) "SEREN timing statistics"
     write(iunit,'(X,A)') "----------------------------------------&
          &------------------------&
          &----------------------------"
     write(iunit,'(X,A)') "Marker id                                      itime            &
          &i%       rtime            r%"
     write(iunit,'(X,A)') "----------------------------------------&
          &------------------------&
          &----------------------------"
     
     do i=1,remote_mark_tot(-1)
        ipercentage = 100.0_DP * real(iblock(i),DP) / real(itime,DP)
        rpercentage = 100.0_DP * rblock(i) / rtime
        ipercentage = 100.d0 * real(remote_iblock(i,-1),DP) / &
           & real(tot_itime,DP) ! divide by total itime per task
        rpercentage = 100.d0 * remote_rblock(i,-1) / &
           & tot_rtime          ! divide by total rtime per task
        write(iunit,fs) remote_marker_id(i,-1),remote_iblock(i,-1),&
           &ipercentage,remote_rblock(i,-1),rpercentage
     end do
     write(iunit,'(X,A)') "----------------------------------------&
          &------------------------&
          &----------------------------"
     write(iunit,fs) 'TOTAL               ',itime,100.0,rtime,100.0
     write(iunit,'(X,A)') "----------------------------------------&
          &------------------------&
          &----------------------------"
     write(iunit,*)

     write(iunit,*) "SEREN computational rate statistics"
     write(iunit,'(X,A)') "----------------------------------------&
          &------------------------&
          &----------------------------"

     do i=1,remote_mark_tot(-1)
        if (remote_rblock(i,-1) <= SMALL_NUMBER) cycle
#if defined(GRAVITY)
        if (trim(adjustl(remote_marker_id(i,-1))) == "GRAVITY_FORCES") write(iunit,*) &
             &"Rate of gravity calcs : ",real(ngravcomp,DP)/remote_rblock(i,-1)
#endif
        if (trim(adjustl(remote_marker_id(i,-1))) == "HYDRO_FORCES") write(iunit,*) &
             &"Rate of hydro calcs   : ",real(nhydrocomp,DP)/remote_rblock(i,-1)
     end do
     write(iunit,'(X,A)') "----------------------------------------&
          &------------------------&
          &----------------------------"
     
     close(iunit)

  end if
! ----------------------------------------------------------------------------

#else

! Loop over all blocks and print statistics
! ----------------------------------------------------------------------------
  if (mark_tot > 0) then

     filename = trim(adjustl(run_dir))//trim(adjustl(run_id))//".timing"
     open(unit=iunit,file=filename,status="unknown")

     fs = "(1X,a20,2X,i10,2X,1F8.4,2X,1G14.6,2X,1F8.4)"

     write(iunit,*) "SEREN timing statistics"
     write(iunit,*) "----------------------------------------&
          &----------------------------"
     write(iunit,*) "Marker id                  itime        &
          &i%       rtime            r%"
     write(iunit,*) "----------------------------------------&
          &----------------------------"
     
     do i=1,mark_tot
        ipercentage = 100.0_DP * real(iblock(i),DP) / real(itime,DP)
        rpercentage = 100.0_DP * rblock(i) / rtime
        write(iunit,fs) marker_id(i),iblock(i),ipercentage,rblock(i),rpercentage
     end do
     write(iunit,*) "----------------------------------------&
          &----------------------------"
     write(iunit,fs) 'TOTAL               ',itime,100.0,rtime,100.0
     write(iunit,*) "----------------------------------------&
          &----------------------------"
     write(iunit,*)

     write(iunit,*) "SEREN computational rate statistics"
     write(iunit,*) "----------------------------------------&
          &----------------------------"
     write(iunit,*) "No. of steps          : ",nsteps
     do i=1,mark_tot
        if (rblock(i) <= SMALL_NUMBER) cycle
#if defined(GRAVITY)
        if (trim(adjustl(marker_id(i))) == "GRAVITY_FORCES") then
           write(iunit,*) "No. of gravity calcs  : ",ngravcomp
           write(iunit,*) "Rate of gravity calcs : ",real(ngravcomp,DP)/rblock(i)
        end if
#endif
        if (trim(adjustl(marker_id(i))) == "HYDRO_FORCES") then
           write(iunit,*) "No. of hydro calc     : ",nhydrocomp
           write(iunit,*) "Rate of hydro calcs   : ",real(nhydrocomp,DP)/rblock(i)
        end if
     end do
     write(iunit,*) "----------------------------------------&
          &----------------------------"
     
     close(iunit)

  end if
! ----------------------------------------------------------------------------
#endif

  return
END SUBROUTINE write_timing_stats

#if defined(USE_MPI)
SUBROUTINE reorder_markers(nlist, ilist, rlist, clist)
  use definitions
  implicit none

  integer, intent(in) :: nlist                   ! No. of pot. neighbours
  integer (kind=ILP), intent(inout) :: ilist(1:nlist)       ! list of ids
  real(kind=DP), intent(inout) :: rlist(1:nlist) ! list of reals to sort
  character(len=50), intent(inout) :: clist(1:nlist) ! list of chars to sort

  integer :: i                        ! counter in sort
  integer :: j                        ! counter in list
  integer(kind=ILP) :: itemp                    ! Aux. variable for sorting
  real(kind=DP) :: rtemp              ! Aux. variable for sorting
  character(len=50) :: ctemp          ! Aux. variable for sorting

  do j=2,nlist
     itemp = ilist(j)
     rtemp = rlist(j)
     ctemp = clist(j)
     do i=j-1,1,-1
        if (rlist(i) >= rtemp) exit
        rlist(i+1) = rlist(i)
        ilist(i+1) = ilist(i)
        clist(i+1) = clist(i)
     end do
     rlist(i+1) = rtemp
     ilist(i+1) = itemp
     clist(i+1) = ctemp
  end do

  return
END SUBROUTINE reorder_markers
#endif
