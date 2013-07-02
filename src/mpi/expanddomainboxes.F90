! EXPANDDOMAINBOXES.F90
! A. Mcleod - 09/07/09
! Expands the domain box to include particles which have moved outside it,
! then broadcasts these updated boxes by first gathering them at the root
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine expanddomainboxes()
   ! Expand and broadcast domain boxes to all tasks, by gathering all at the root
   ! and send broadcasting out to all other tasks
   use mpi
   use mpi_communication_module
   use particle_module
   use periodic_module
   use domain_comparisons
   implicit none

   real(kind=PR) :: dminmax(1:2*NDIM,0:lastrank) ! Temporary array of bbmin and bbmax
   integer                      :: p     ! Loop counter over particles
   integer                      :: ierr  ! Return code
#ifdef DEBUG_BOXES
   integer                      :: d     ! Loop counter over domains
#endif
   real(kind=PR) :: r_max(1:NDIM), r_min(1:NDIM) ! minimum and maximum positions
   real(kind=PR), allocatable :: r(:,:) ! Temporary list of particle positions

   debug2("Expanding domain boxes [expanddomainboxes.F90]")

   debug_timing("EXPAND_DOMAIN_BOXES_PREP")

! Expand domain boundaries if particles have moved

! no idea why the following bit had smoothing kerles in?

!   call maxmin_hydro(r_max,r_min,ptot,SMOO,parray(1:SMOO,1:ptot))
!   r_max(1:NDIM) = -BIG_NUMBER
!   r_min(1:NDIM) =  BIG_NUMBER
!   do p=1,ptot
!      r_max = max(r_max, sph(p)%r + KERNRANGE*sph(p)%h)
!      r_min = min(r_min, sph(p)%r - KERNRANGE*sph(p)%h)
!   end do
!   r_min = r_min - 2._PR*spacing(r_max)
!   r_max = r_max + 2._PR*spacing(r_max)
!   domain_bbmin(1:NDIM,rank) = min(domain_bbmin(1:NDIM,rank),r_min)
!   domain_bbmax(1:NDIM,rank) = max(domain_bbmax(1:NDIM,rank),r_max)

   allocate(r(1:NDIM,1:ptot))
   do p=1,ptot
      r(1:NDIM,p) = sph(p)%r(1:NDIM)
   end do
   call bounding_box(ptot,r,r_max(1:NDIM),r_min(1:NDIM))
   deallocate(r)
   r_min = r_min - 2._PR*spacing(r_max)
   r_max = r_max + 2._PR*spacing(r_max)
   domain_bbmin(1:NDIM,rank) = min(domain_bbmin(1:NDIM,rank),r_min)
   domain_bbmax(1:NDIM,rank) = max(domain_bbmax(1:NDIM,rank),r_max)
   
   debug_timing("EXPAND_DOMAIN_BOXES")

! Pack min and max bounding box into single array
   dminmax(1:NDIM,rank) = domain_bbmin(1:NDIM,rank)
   dminmax(NDIM+1:2*NDIM,rank) = domain_bbmax(1:NDIM,rank)

   ! 'Gather' the bounding boxes from each task into a single array at the root
   WAIT_TIME_MACRO
   if (rank==0) then
      call MPI_GATHER(MPI_IN_PLACE,0,MPI_REAL_PR,dminmax,2*NDIM,&
      &MPI_REAL_PR,0,MPI_COMM_WORLD,ierr)
      !CHOICE SENDBUF,INTEGER SENDCOUNT,INTEGER SENDTYPE,CHOICE RECVBUF,INTEGER RECVCOUNT,
      !INTEGER RECVTYPE,INTEGER ROOT,INTEGER COMM,INTEGER IERROR)
   else
      call MPI_GATHER(dminmax(1,rank),2*NDIM,MPI_REAL_PR,0,0,0,0,MPI_COMM_WORLD,ierr)
      !CHOICE SENDBUF,INTEGER SENDCOUNT,INTEGER SENDTYPE,CHOICE RECVBUF,INTEGER RECVCOUNT,
      !INTEGER RECVTYPE,INTEGER ROOT,INTEGER COMM,INTEGER IERROR)
   end if

! Data is now gathered at the root ready to be broadcasted out
   call MPI_BCAST(dminmax, 6*numtasks, MPI_REAL_PR, 0, MPI_COMM_WORLD, ierr)
   CALC_TIME_MACRO

! Collect bounding box information back into arrays
   domain_bbmin(1:NDIM,0:lastrank) = dminmax(1:NDIM,0:lastrank)
   domain_bbmax(1:NDIM,0:lastrank) = dminmax(NDIM+1:2*NDIM,0:lastrank)

! Check that bounding boxes have been broadcast correctly.
#ifdef DEBUG_BOXES
   write (6,*) "domain_bbmin after = ", domain_bbmin(1:NDIM,rank)
   write (6,*) "domain_bbmax after = ", domain_bbmax(1:NDIM,rank)
   if (rank==0) then
      write(6,*) "EXPAND DOMAIN BOXES Writing domain boxes to screen"
      do d=0,lastrank
         write(6,'(A,I0,A,3(F0.12,X))') "Rank : ",d," domain min = ", domain_bbmin(1:NDIM,d)
         write(6,'(A,I0,A,3(F0.12,X))') "Rank : ",d," domain max = ", domain_bbmax(1:NDIM,d)
      end do
   end if
#endif

   return
end subroutine expanddomainboxes
