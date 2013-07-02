! BROADCASTBOUNDINGBOXES.F90
! A. Mcleod - 03/07/08
! Broadcasts the domain bounding boxes to all domains,
! and receives those bounding boxes from other domains.
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine broadcastboundingboxes(boxmin,boxmax)
  use mpi
  use mpi_communication_module
  use particle_module, only : ptot
  use domain_comparisons
  implicit none

  real(kind=PR),intent(inout) :: boxmin(1:NDIM,0:lastrank)  ! boxmin to share
  real(kind=PR),intent(inout) :: boxmax(1:NDIM,0:lastrank)  ! boxmax to share
  real(kind=PR) :: boxminmax(1:2*NDIM,0:lastrank) ! Temp. array of boxmin and boxmax
  integer :: ierr                 ! Return code
#ifdef DEBUG_BOXES
  integer :: d                    ! Loop counter over domains
#endif

  debug2("Broadcasting bounding boxes [broadcastboundingboxes.F90]")

! Pack min and max bounding box into single array
  boxminmax(1:NDIM,rank) = boxmin(1:NDIM,rank)
  boxminmax(NDIM+1:2*NDIM,rank) = boxmax(1:NDIM,rank)

! 'Gather' the bounding boxes from each task into a single array at the root
  WAIT_TIME_MACRO
  if (rank==0) then
     call MPI_GATHER(MPI_IN_PLACE,0,MPI_REAL_PR,boxminmax,2*NDIM,&
          &MPI_REAL_PR,0,MPI_COMM_WORLD,ierr)
     !CHOICE SENDBUF,INTEGER SENDCOUNT,INTEGER SENDTYPE,CHOICE RECVBUF,INTEGER RECVCOUNT,
     !INTEGER RECVTYPE,INTEGER ROOT,INTEGER COMM,INTEGER IERROR)
  else
     call MPI_GATHER(boxminmax(1,rank),2*NDIM,MPI_REAL_PR,&
          &0,0,0,0,MPI_COMM_WORLD,ierr)
     !CHOICE SENDBUF,INTEGER SENDCOUNT,INTEGER SENDTYPE,CHOICE RECVBUF,INTEGER RECVCOUNT,
     !INTEGER RECVTYPE,INTEGER ROOT,INTEGER COMM,INTEGER IERROR)
  end if

! Data is now gathered at the root ready to be broadcasted out
  call MPI_BCAST(boxminmax,2*NDIM*numtasks,MPI_REAL_PR,0,MPI_COMM_WORLD,ierr)
  CALC_TIME_MACRO

! Collect bounding box information back into arrays
  boxmin(1:NDIM,0:lastrank) = boxminmax(1:NDIM,0:lastrank)
  boxmax(1:NDIM,0:lastrank) = boxminmax(NDIM+1:2*NDIM,0:lastrank)

! Check that bounding boxes have been broadcast correctly.
#ifdef DEBUG_BOXES
  write(6,*) "BROADCAST BOUNDING BOXES Writing bounding boxes to screen"
  do d=0,lastrank
     if (rank==0) write(6,'(A,I0,A,3(F0.7,X))') "boxmin(rank ",d,") : ",boxmin(1:NDIM,d)
     if (rank==0) write(6,'(A,I0,A,3(F0.7,X))') "boxmax(rank ",d,") : ",boxmax(1:NDIM,d)
  end do
#endif

  return
end subroutine broadcastboundingboxes
