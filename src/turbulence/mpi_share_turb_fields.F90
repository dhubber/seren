! MPI_SHARE_TURB_FIELDS.F90
! A.McLeod - 24/11/2012
! Share turbulent acceleration fields to non-root tasks
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE mpi_share_turb_fields
   use definitions
   use turbulence_module
   use mpi_communication_module
   use mpi
   implicit none
   
   integer, parameter :: message_length=NDIM*TURB_GS**3
                                               ! Length of flattened arrays
   integer            :: ierr                  ! MPI error variable
   
   ! NOTE - this is all a bit inefficient as we don't usually need to copy
   ! afield_last (it is just the old afield_next) except at the start
   
   ! Share afield_last and afield_next
   WAIT_TIME_MACRO
   call MPI_BCAST(afield_last, message_length, MPI_COMPLEX_PR, 0, &
                  & MPI_COMM_WORLD, ierr)
   
   call MPI_BCAST(afield_next, message_length, MPI_COMPLEX_PR, 0, &
                  & MPI_COMM_WORLD, ierr)
   
   call MPI_BCAST(turb_last_time, 1, MPI_DOUBLE_PRECISION, 0, &
                  & MPI_COMM_WORLD, ierr)
   
   call MPI_BCAST(turb_next_time, 1, MPI_DOUBLE_PRECISION, 0, &
                  & MPI_COMM_WORLD, ierr)
   CALC_TIME_MACRO
   
   return
   
END SUBROUTINE mpi_share_turb_fields