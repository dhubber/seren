! FORCES.F90
! C. P. Batty & D. A. Hubber - 11/1/2007
! Calculates accelerations (hydrodynamic and gravitational) for all particles
! - controls calls to relevant subroutines
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE gravity_export_return(send_in)
  use type_module
  use particle_module
#ifdef SINKS
  use sink_module
#endif
  use mpi_communication_module
  use hydro_module
  use mpi
#ifdef DEBUG_EXPORT
  use filename_module, only : run_id, MPI_ext
#endif
  implicit none

  integer, intent(out) :: send_in(0:endranklist)
                 ! Send requests from original sending of exported particles
  integer :: send(0:endranklist)   ! Array of send requests (one for each OTHER task)
  integer :: recv(0:endranklist)   ! Array of recv requests (one for each OTHER task)
  integer :: d                     ! domain counter
  integer :: slot
  integer :: i, j                  ! loop counters
  integer :: p                     ! particle counter
#ifdef SINKS
  integer :: s                     ! sink counter
#endif
  integer :: stat(MPI_STATUS_SIZE) ! MPI Status variable
  integer :: ierr                  ! Return code
  integer :: pcount                ! Number of particles received per receive
  real (kind=DP) :: agravp(1:NDIM) ! Particle gravitational
                                   ! acceleration from export calculation

#ifdef DEBUG_EXPORT
  character(len=100)    :: out_file ! filename extension for data output
#endif

  debug2("Returning gravity results from exported particles [gravity_export_return.F90]")
  
  debug_timing("FORCES GRAV WAITALL SEND")
! Wait until we finish sending exported particles
  WAIT_TIME_MACRO
  call MPI_WAITALL(lastrank,send_in,MPI_STATUSES_IGNORE,ierr)
  CALC_TIME_MACRO

  debug_timing("GRAV_EXP_RET SEND")

! Can't do this until all sends completed
  deallocate(sendgravity)

! Return particles to their original domains
  slot = 1
  do i=0,lastrank-1
    d = grav_fromeach(1,i)        ! Domain of the current set of particles
#ifdef DEBUG_EXPORT
    write (6,*) "Rank ",rank," sending results of exported of gravity particles to task ", d
    write (6,*) "Message length ", grav_fromeach(2,i), ", from slot ", slot
#endif
    call MPI_ISEND(returngravity(slot), grav_fromeach(2,i), MPI_RETURNGRAVITY, &
         & d, GRAV_RETURN_TAG, MPI_COMM_WORLD,send(i),ierr)
    slot = slot + grav_fromeach(2,i)
  end do

  debug_timing("GRAV_EXP_RET RECV")

! Receive exported particles
  allocate(receivereturngravity(0:endranklist))
  pexportgrav = 0

  do i=0,lastrank-1
    ! Test for any message - wait until one arrives (blocking probe)
    debug_timing("GRAV_EXP_RET RECV (PROBE RECV)")
    WAIT_TIME_MACRO
    call MPI_PROBE(MPI_ANY_SOURCE, GRAV_RETURN_TAG, MPI_COMM_WORLD, stat, ierr)
    CALC_TIME_MACRO
    debug_timing("GRAV_EXP_RET RECV")
    d = stat(MPI_SOURCE)
    call MPI_GET_COUNT(stat,MPI_RETURNGRAVITY,pcount,ierr)

    allocate(receivereturngravity(i)%data(1:pcount+1))

! Post the receive - until this is posted the receive does not actually occur
! and may or may not be buffered
    call MPI_IRECV(receivereturngravity(i)%data(1),pcount,MPI_RETURNGRAVITY,&
         &d,GRAV_RETURN_TAG,MPI_COMM_WORLD,recv(i),ierr)
    pexportgrav = pexportgrav + pcount

    grav_fromeach(2,i) = pcount

  end do

  debug_timing("GRV_EX_RET RECV(WAIT ALL RECV)")
  WAIT_TIME_MACRO
  call MPI_WAITALL(lastrank,recv,MPI_STATUSES_IGNORE,ierr)
  CALC_TIME_MACRO
  debug_timing("GRAV_EXP_RET RECV")
  ! All our receives are now completed

  debug_timing("GRAV_EXP_RET MERGE RESULTS")
! Merge the results of the exportation with the local calculation
  do i=0,lastrank-1
    do j=1,grav_fromeach(2,i)
      p = receivereturngravity(i)%data(j)%p
#ifdef SINKS
      if (p > 0) then
#endif
        agravp = receivereturngravity(i)%data(j)%acc
        sph(p)%a = sph(p)%a + agravp
        sph(p)%gpot = sph(p)%gpot + receivereturngravity(i)%data(j)%pot
#ifdef SINKS
      else
        ! Sink particle
        s = -p
        agravp = receivereturngravity(i)%data(j)%acc
        sink(s)%remote_agrav = sink(s)%remote_agrav + agravp
        sink(s)%remote_gpot = sink(s)%remote_gpot + receivereturngravity(i)%data(j)%pot
      end if
#endif
    end do
  end do

  do i=0,lastrank-1
    deallocate(receivereturngravity(i)%data)
  end do
  deallocate(receivereturngravity)
  
  debug_timing("FORCES GRAV_TREES WAITALL RECV")
! Wait for pruned trees to finish sending if they have not already done so
  WAIT_TIME_MACRO
  call MPI_WAITALL(lastrank,treesend,MPI_STATUSES_IGNORE,ierr)
  CALC_TIME_MACRO

! Wait for all sends to complete
  debug_timing("GRAV_EXP_RET RECV(WAITALL SEND)")
  WAIT_TIME_MACRO
  call MPI_WAITALL(lastrank,send,MPI_STATUSES_IGNORE,ierr)
  CALC_TIME_MACRO

! Can't do this until all sends completed
  deallocate(returngravity)

  return
END SUBROUTINE gravity_export_return
