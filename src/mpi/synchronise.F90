! SYNCHRONISE.F90
! A. Mcleod - 03/07/08
! Shares important variables between threads at the start of each timestep
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine synchronise(step_min,nmax)
   use particle_module, only : ptot
   use mpi
   use mpi_communication_module
   implicit none

   real(kind=DP), intent(inout) :: step_min     ! Minimum step size
   integer(kind=ILP), intent(inout) :: nmax     ! Maximum timestep level

   integer :: ierr                              ! MPI error variable
   type (synchronisetype) :: sync(0:lastrank)   ! Synchronisation type

   debug_timing("SYNCHRONISE")

   ! Set variables to be gathered
   sync(rank)%totalptot = ptot
   sync(rank)%step_min = step_min
   sync(rank)%loadbalance = do_load_balance
   sync(rank)%nmax = nmax

   WAIT_TIME_MACRO
   ! 'Gather' the totalptots from each task into a single array at the root
   if (rank==0) then
      call MPI_GATHER(MPI_IN_PLACE,1,MPI_SYNCHRONISE,sync(0),1,&
      &MPI_SYNCHRONISE,0,MPI_COMM_WORLD,ierr)
      !CHOICE SENDBUF,INTEGER SENDCOUNT,INTEGER SENDTYPE,CHOICE RECVBUF,INTEGER RECVCOUNT,
      !INTEGER RECVTYPE,INTEGER ROOT,INTEGER COMM,INTEGER IERROR)
      sync(0)%totalptot = sum(sync(0:lastrank)%totalptot)
      sync(0)%step_min = minval(sync(0:lastrank)%step_min)
      sync(0)%loadbalance = any(sync(0:lastrank)%loadbalance)
      sync(0)%nmax = maxval(sync(0:lastrank)%nmax)
   else
      call MPI_GATHER(sync(rank),1,MPI_SYNCHRONISE,0,0,0,0,MPI_COMM_WORLD,ierr)
      !CHOICE SENDBUF,INTEGER SENDCOUNT,INTEGER SENDTYPE,CHOICE RECVBUF,INTEGER RECVCOUNT,
      !INTEGER RECVTYPE,INTEGER ROOT,INTEGER COMM,INTEGER IERROR)
   end if

   ! Data is now gathered at the root ready to be broadcasted out
   call MPI_BCAST(sync(0), 1, MPI_SYNCHRONISE, 0, MPI_COMM_WORLD, ierr)
   CALC_TIME_MACRO

   ! Assign data based on results from root
   totalptot = sync(0)%totalptot
   step_min = sync(0)%step_min
   do_load_balance = sync(0)%loadbalance
   nmax = sync(0)%nmax

   return
end subroutine synchronise
