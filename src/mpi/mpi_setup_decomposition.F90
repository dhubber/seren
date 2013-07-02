! MPI_SETUP_DECOMPOSITION.F90
! A. McLeod - 31/07/08
! ..
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine mpi_setup_decomposition
   use filename_module
   use interface_module, only : read_data
   use mpi_communication_module
   use mpi
   implicit none
   
   integer :: ierr

   ! -------------------------------------------------------------------------
   if (rank==0) then

      ! Reading input snapshot data file
      call read_data(in_file,in_file_form,.TRUE.)
  
      ! Setting up variables for different particle types
      ! Call types_1 as this does not set ptypes for particles, which we
      ! don't need yet
      call types_1
      
      ! Converting all important physical variables to dimensionless code units
      ! This does not rescale parameters, just physical quantities
      call convert_to_code_units_1(.TRUE.)
      
      ! Perform domain decomposition
      call decomposition
      
      ! Clean up in order to re-initialise as separate domain 0
      call clean_up
   end if
   ! -------------------------------------------------------------------------
   
   ! Ensure that the master thread is finished with decomposition
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)

   ! Share data from the master to the slave threads
   call mpi_share_data

   ! Read in domain data file for each domain, produced by decomposition
   call read_data(out_temp,"seren_unform")

   ! Delete temporary decomposition files
   open(1, file=out_temp, status="unknown", form="unformatted")
   close(1,status="delete")
  
   return
end subroutine mpi_setup_decomposition
