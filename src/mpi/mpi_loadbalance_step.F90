! MPI_LOADBALANCE_STEP
! A. McLeod - 17/03/2011
! Subroutine for organising load-balances
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine mpi_loadbalance_step
   use mpi_communication_module
   use time_module
   use debug_tests
   implicit none
   
   ! Expand domain boxes to contain each task's particles
   if (.NOT.do_load_balance) call expanddomainboxes
   
   ! See if we should be doing a load balance, 
   ! or whether to do one next integration
   ! -------------------------------------------------------------------------
   if (numtasks > 1) then
      if (do_load_balance) then
         call loadbalance
         
         ! Send and receive particles that have left or entered our domain
         call transfer_particles
         call expanddomainboxes
#ifdef DEBUG2
         call check_particles_in_boxes
#endif
         call loadbalance_time
      else if (demand_load_balance) then
         ! If a task really, really wants a load balance, set do_load_balance
         do_load_balance = .TRUE.
      else if (rank==0 .and. &
         & (nsteps - nloadbalance + 1_ILP >= loadbalance_nsteps)) then
         ! If it is time for a loadbalance anyway, set do_load_balance
         ! Do a load balance next step
         do_load_balance = .TRUE.
      end if
   end if
   demand_load_balance = .FALSE.
   ! -------------------------------------------------------------------------

   return
end subroutine mpi_loadbalance_step
