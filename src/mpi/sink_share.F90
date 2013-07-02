! SINK_SHARE.F90
! A. McLeod - 5/9/12
! Collects and broadcasts sink information
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sink_share
   use particle_module
   use mpi_communication_module
   use mpi
   use sink_module
   use domain_comparisons
   implicit none

   integer             :: ierr                ! MPI error value
   integer             :: d                   ! Domain counter
   integer             :: s                   ! Sink counter
   integer             :: i                   ! Counter
   integer             :: slot                ! Slot for sinks
   integer             :: n_local             ! Number of local sinks
   integer             :: n_sinks(0:lastrank) ! Number of remote sinks  
   integer :: recv(0:endranklist)             ! Array of receive requests
   type(sink_node)     :: local_sinks(1:stot) ! Local sinks
   type(sink_node), allocatable :: remote_sinks(:) ! Remote sinks
   
   debug2("Sharing sink data to all tasks")

   ! This is all unnecessary if we are the only task
   
   if (numtasks == 1) return
   
   n_local = 0
   n_sinks(0:lastrank) = 0
   
   ! Identify local sinks and load into local_sinks if not the root
   if (rank /= 0) then
      do s=1,stot
         if (sink(s)%domain == rank) then
            n_local = n_local + 1
            local_sinks(n_local) = sink(s)
            local_sinks(n_local)%slot = s
         end if
      end do
   end if
   WAIT_TIME_MACRO
   call MPI_GATHER(n_local, 1, MPI_INTEGER, n_sinks, 1, MPI_INTEGER, &
                  &0, MPI_COMM_WORLD, ierr)
   CALC_TIME_MACRO
   
   ! Do the transmission of sinks
   if (rank==0) then
      ! Receive sinks
      i = 0
      slot = 1
      allocate(remote_sinks(1:sum(n_sinks)+1)) ! +1 to prevent out-of-range
      do d=1,lastrank
         if (n_sinks(d) > 0) then
            call MPI_IRECV(remote_sinks(slot), n_sinks(d), MPI_SINK_NODE, d, &
                          &SINK_TAG, MPI_COMM_WORLD, recv(i), ierr)
         else
            recv(i) = MPI_REQUEST_NULL
         end if
         slot = slot + n_sinks(d)
         i = i + 1
      end do
   else
      ! Send local sinks
      WAIT_TIME_MACRO
      if (n_local > 0) call MPI_SEND(local_sinks, n_local, MPI_SINK_NODE, 0,&
                            &SINK_TAG, MPI_COMM_WORLD, ierr)
      CALC_TIME_MACRO
   end if
   
   ! Load sinks into correct places
   if (rank==0) then
      WAIT_TIME_MACRO
      call MPI_WAITALL(lastrank, recv, MPI_STATUSES_IGNORE, ierr)
      CALC_TIME_MACRO
      do i=1,sum(n_sinks)
         s = remote_sinks(i)%slot
         sink(s) = remote_sinks(i)
      end do
   end if
   
   ! Work out if sinks should switch tasks
   ! Redistribute sinks to domains
   if (rank==0) then
      sinkloop: do s=1,stot
         if (.NOT. inside_domain_box(sink(s)%r,sink(s)%domain)) then
            ! Sink has moved outside of it's original domain,
            ! find it a new domain
            do d=0,lastrank
               ! Test if the sink is inside domain d
               if (inside_domain_box(sink(s)%r,d)) then
                  ! Sink is inside domain d
                  sink(s)%domain = d
                  cycle sinkloop
               end if
            end do
            write (6,*) "Sink ", s, " outside domain boundary but not sent anywhere!"
            write (6,*) "Sink is at position ", sink(s)%r
            stop
         end if
      end do sinkloop
   end if
   
   ! Distribute sink data from root to all tasks
   WAIT_TIME_MACRO
   call MPI_BCAST(sink, stot, MPI_SINK_NODE, 0, MPI_COMM_WORLD, ierr)
   CALC_TIME_MACRO

   return
END SUBROUTINE sink_share
