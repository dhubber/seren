! MPI_SHARE_DATA.F90
! A. McLeod 16/06/08
! Shares initial data between the threads
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine mpi_share_data
  use mpi
  use mpi_communication_module
  use mpi_decomposition_module, only : MPItreedepth
  use particle_module, only : ptot
  use periodic_module
  use definitions
  implicit none

  integer               :: ierr                ! Return code
  
  debug2("Initial MPI data sharing [mpi_share_data.F90]")

  ! Broadcast domain boundaries to all tasks
  call MPI_BCAST(domain_bbmin(1,0),NDIM*numtasks,MPI_REAL_PR,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(domain_bbmax(1,0),NDIM*numtasks,MPI_REAL_PR,0,MPI_COMM_WORLD,ierr)
  local_min = domain_bbmin(1:NDIM,rank)
  local_max = domain_bbmax(1:NDIM,rank)

  ! Broadcast total number of particles to all tasks
  call MPI_BCAST(totalptot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  ! Broadcast MPItreedepth
  call MPI_BCAST(MPItreedepth,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  ! Allocate fromeach
  allocate(grav_fromeach(1:2,0:endranklist))

  return
end subroutine mpi_share_data
