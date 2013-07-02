! LOADBALANCE_TIME.F90
! A. McLeod - 31/07/08
!
! -------------------------------------------------------------------------

#include "macros.h"

! -------------------------------------------------------------------------
subroutine loadbalance_time
   use particle_module
   use mpi_communication_module
   use time_module
   use mpi
   use type_module, only : phydrostart
   implicit none

   integer                  :: ierr            ! MPI error value

   integer(kind=ILP)        :: last_loadbalance_nsteps ! Last loadbalance_nsteps

   real(kind=DP)            :: accdo_per_nstep ! Estimated accdo's per nstep
   real(kind=DP)            :: time_per_accdo  ! Time per accdo last step
   real(kind=DP)            :: fudge_factor    ! Estimate of last time's misprediction

   last_loadbalance_nsteps = loadbalance_nsteps

   ! Calculate the sum of 1/nstep(p) for all particles; this is used
   ! to calculate a predict_acctot across all tasks
   accdo_per_nstep = sum(1._DP/real(2**(level_step-sph(phydrostart:ptot)%nlevel),DP))
   
   if (rank==0) then
      call MPI_REDUCE(MPI_IN_PLACE,accdo_per_nstep,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   else
      call MPI_REDUCE(accdo_per_nstep,0,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   end if

! Here adjust loadbalance_nsteps if it is too quick or too slow
   if (rank==0) then
!       write (6,*) "old loadbalance_nsteps = ", loadbalance_nsteps
!       write (6,*) "aim_time_load_balance = ", aim_time_load_balance
!       write (6,*) "last_loadbalance_time = ", last_loadbalance_time

      ! Multiply by nstepsize to get accdo per nstep
      accdo_per_nstep = accdo_per_nstep * real(nstepsize,DP)
!       write (6,*) "accdo_per_nstep = ", accdo_per_nstep

      ! Damped misprediction factor
      fudge_factor = sqrt(aim_time_load_balance / last_loadbalance_time)
!       write (6,*) "fudge_factor = ", fudge_factor

      ! This all relies on predict_acctot being reasonably accurate...

      time_per_accdo = last_loadbalance_time / real(sum_acctot,DP)
!       write (6,*) "time_per_accdo = ", time_per_accdo

      ! Calculate our desired loadbalance_nsteps
      loadbalance_nsteps = nint( &
         & (aim_time_load_balance / time_per_accdo) / & ! Number of accdo's required
         & accdo_per_nstep ,ILP)  ! Number of accdo's per nstep
!       write (6,*) "new loadbalance_nsteps = ", loadbalance_nsteps

      ! Fudge our value
      loadbalance_nsteps = nint(fudge_factor * real(loadbalance_nsteps,DP),ILP)
!       write (6,*) "fudged loadbalance_nsteps = ", loadbalance_nsteps

      ! Make an even number
      if (mod(loadbalance_nsteps,2) /= 0) &
         loadbalance_nsteps = loadbalance_nsteps + 1
!       write (6,*) "made-even loadbalance_nsteps = ", loadbalance_nsteps

      ! Finally limit to 2 times last time's value
      loadbalance_nsteps = min(2_ILP*last_loadbalance_nsteps,loadbalance_nsteps)
!       write (6,*) "loadbalance_nsteps = ", loadbalance_nsteps

      ! Also limit to min_steps_load_balance
      loadbalance_nsteps = max(min_steps_load_balance,loadbalance_nsteps)
!       write (6,*) "loadbalance_nsteps = ", loadbalance_nsteps

      !and max_steps_load_balance
      loadbalance_nsteps = min(max_steps_load_balance,loadbalance_nsteps)
!       write (6,*) "loadbalance_nsteps = ", loadbalance_nsteps

   end if

   call MPI_BCAST(loadbalance_nsteps,1,MPI_INTEGER_ILP,0,MPI_COMM_WORLD,ierr)

   ! Now we reset the sum_acctot counter
   sum_acctot = 0

   ! Predict what sum_acctot will be locally, and we can compare it later...
   predict_acctot = sum(1._DP/real(2**(level_step-sph(phydrostart:ptot)%nlevel),DP))
   predict_acctot = predict_acctot * real(loadbalance_nsteps,PR) * real(nstepsize,PR)

   return
end subroutine loadbalance_time
