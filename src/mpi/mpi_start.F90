! MPI_START.F90
! A. McLeod - 31/07/08
! Subroutine for MPI things that have to happen when Seren loads
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine mpi_start_seren
   use mpi_communication_module
   use mpi
#if defined(OPENMP)
   use omp_lib
#endif
   implicit none

   integer :: ierr                      ! MPI error value
#if defined(OPENMP)
   integer :: threadmodel               ! MPI thread model
#endif

   ! Start MPI processes - this splits the program into separate threads
#if defined(OPENMP)
   call mpi_init_thread(MPI_THREAD_MULTIPLE,threadmodel,ierr)
   if (threadmodel /= MPI_THREAD_MULTIPLE) then
      write (6,*) "MPI_THREAD_MULTIPLE not provided!"
      if (threadmodel == MPI_THREAD_SINGLE) &
         & write (6,*) "We only have thread model MPI_THREAD_SINGLE"
      if (threadmodel == MPI_THREAD_FUNNELED) &
         & write (6,*) "We only have thread model MPI_THREAD_FUNNELED"
      if (threadmodel == MPI_THREAD_SERIALIZED) &
         & write (6,*) "We only have thread model MPI_THREAD_SERIALIZED"
      !stop ! Why not give it a shot anyway? Probably fine...
   end if
#else
   call mpi_init(ierr)
#endif

   ! Quit program if there's an error starting MPI
   if (ierr .ne. mpi_success) then
      write (6,*) "Error starting MPI."
      call mpi_abort(mpi_comm_world, 1, ierr)
      stop
   end if

   ! Set rank variables for multiple MPI tasks
   call mpi_comm_rank(mpi_comm_world, rank, ierr)
   call mpi_comm_size(mpi_comm_world, numtasks, ierr)
   lastrank = numtasks - 1
   lastrank2 = max(1,lastrank)
   endranklist = max(0,lastrank-1) ! For lists 0:lastrank-1
   write (MPI_ext,"(I0)") rank

   ! -------------------------------------------------------------------------
   if (rank == 0) then
      write(6,*) "**********************************************************************"
      write(6,*) "*            ****     ******    *****     ******   *     *           *"
      write(6,*) "*           *    *    *         *    *    *        **    *           *"
      write(6,*) "*           *         *         *    *    *        * *   *           *"
      write(6,*) "*            ****     *****     *****     ******   *  *  *           *"
      write(6,*) "*                *    *         *    *    *        *   * *           *"
      write(6,*) "*           *    *    *         *    *    *        *    **           *"
      write(6,*) "*            ****     ******    *    *    ******   *     *           *"
      write(6,*) "*                                                                    *"
      write(6,*) "* |     | ==== =====    ===== ====  ===== ===== =====  ====  |     | *"
      write(6,*) "* |\   /| |   |  |      |     |   \   |     |     |   |    | |\    | *"
      write(6,*) "* | \ / | |   |  |      |     |    |  |     |     |   |    | | \   | *"
      write(6,*) "* |  |  | |===   |      |===  |    |  |     |     |   |    | |  \  | *"
      write(6,*) "* |     | |      |      |     |    |  |     |     |   |    | |   \ | *"
      write(6,*) "* |     | |      |      |     |   /   |     |     |   |    | |    \| *"
      write(6,*) "* |     | |    =====    ===== ====  =====   |   =====  ====  |     | *"
      write(6,*) "*                                                                    *"
      write(6,*) "*                            Version " // SEREN_VERSION // &
        & "                           *"
      write(6,*) "*                              15/01/2013                            *"
      write(6,*) "*                                                                    *"
      write(6,*) "*        Coders : David Hubber, Chris Batty & Andrew McLeod          *"
      write(6,*) "*                 Thomas Bisbas, Krisada Rawiraswattana,             *"
      write(6,*) "*                 Dimitrios Stamatellos, Stefanie Walch,             *"
      write(6,*) "*                 Anthony Whitworth                                  *"
      write(6,*) "*                                                                    *"
      write(6,*) "*               http://www.astro.group.shef.ac.uk/seren              *"
      write(6,*) "**********************************************************************"

      write(6,*) "Number of MPI processes: ", numtasks
#if defined(OPENMP)
      !$OMP PARALLEL
      !$OMP MASTER
      write (6,*) "Number of OPEN MP threads per task: ", OMP_GET_NUM_THREADS()
      !$OMP END MASTER
      !$OMP END PARALLEL
#endif
#ifdef DEBUG1
      write (6,*) "Timer resolution: ", MPI_WTICK()
#endif
   end if
   ! -------------------------------------------------------------------------

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!  Set up MPI derived types
   call mpi_setup_types

   ! Allocate variables containing domain bounding boxes.  Can't go in 
   ! allocate_memory because they will get cleaned up in clean_up
   ! and we need the initial values
   allocate(domain_bbmin(1:NDIM,0:lastrank))
   allocate(domain_bbmax(1:NDIM,0:lastrank))
   allocate(bbmin(1:NDIM,0:lastrank))
   allocate(bbmax(1:NDIM,0:lastrank))
   allocate(activemin(1:NDIM,0:lastrank))
   allocate(activemax(1:NDIM,0:lastrank))

   lastwrite = MPI_WTIME()
   walltime = MPI_WTIME()
   since_loadbalance_time = MPI_WTIME()
   calctime = 0.0_DP
   diag_calctime = 0.0_DP
   do_load_balance = .FALSE.
   demand_load_balance = .FALSE.
   loadbalance_nsteps = min_steps_load_balance

   return
end subroutine mpi_start_seren
