! START.F90
! A. McLeod - 31/07/08
! Subroutine for MPI things that have to happen when Seren finishes
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine mpi_finish(quiet)
   use mpi_communication_module
   use mpi
   implicit none

   logical, intent(in)  :: quiet     ! ..
   integer :: ierr                   ! MPI error value
   
   if (.not. quiet) then
      if (rank == 0) write(6,'(A,F0.3,A)') &
           &"Walltime of ", MPI_WTIME() - walltime, " seconds"
      if (rank == 0) write(6,'(A,F0.3,A)') &
           &"Without starting time ", MPI_WTIME() - non_start_time, " seconds"
   end if
   call mpi_finalize(ierr)

   return
end subroutine mpi_finish
