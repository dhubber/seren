! MPI_FIND_DENSEST_PARTICLE.F90
! A. M. McLeod - 4/9/2012
! Identifies task with densest particle

subroutine mpi_find_densest_particle(rho,task)
   use mpi_communication_module
   use mpi
   implicit none
   real(kind=PR), intent(in)  :: rho          ! Density of local particle
   integer, intent(out)       :: task         ! Task with densest particle
   real(kind=PR), allocatable :: rhotasks(:)  ! rhomean values from tasks
   integer                    :: ierr         ! MPI error value

   if (rank==0) then
      allocate(rhotasks(1:numtasks))
      call MPI_GATHER(rho,1,MPI_REAL_PR,rhotasks(1),1,MPI_REAL_PR,0,MPI_COMM_WORLD,ierr)
      task = maxloc(rhotasks,1) - 1 ! Need to subtract one because maxloc
                                    ! returns index from 1 (to numtasks),
                                    ! not 0 (to lastrank)
      deallocate(rhotasks)
   else
      call MPI_GATHER(rho,1,MPI_REAL_PR,0,0,0,0,MPI_COMM_WORLD,ierr)
   end if
   call MPI_BCAST(task, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   return
end subroutine mpi_find_densest_particle
