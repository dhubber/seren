! DO_EXPORT_GRAVITY.F90
! A. McLeod - ????
! ..
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE do_export_gravity(recv)
  use type_module
  use particle_module
#ifdef SINKS
  use sink_module
#endif
  use mpi_communication_module
  use hydro_module
  use time_module
  use mpi
#ifdef DEBUG_EXPORT
  use filename_module, only : run_id, MPI_ext
#endif
  implicit none

  integer, intent(inout) :: recv(0:endranklist)
                 ! Array of recv requests (one for each OTHER task)
  integer :: p                     ! particle counter
  integer :: omp_p                 ! particle counter
  integer :: i, j                  ! Loop counters
  integer :: t                     ! Task number that has finished receiving
  integer :: ierr                  ! MPI error variable
  integer :: slot(0:endranklist)   ! Slots to load particles to
  real (kind=DP) :: agravp(1:NDIM) ! Particle gravitational
                                   ! acceleration from export calculation
  real (kind=DP) :: potp           ! Particle gravitational
                                   ! potential from export calculation
  real (kind=DP) :: sphgpotp       ! Purely a dummy variable in this context
  real (kind=PR) :: rp(1:NDIM)     ! Particle position
  real (kind=PR) :: invhp          ! Inverse smoothing length
  real (kind=PR) :: zo_p           ! Zeta over omega
  real (kind=PR) :: agravmag_p     ! agravmag_p
  real (kind=PR) :: gpot_p         ! gpot_p

#ifdef DEBUG_EXPORT
  character(len=100)    :: out_file ! filename extension for data output
#endif

  debug2("Calculating remote gravity on exported particles [do_export_gravity.F90]")
  debug_timing("DO_EXPORT_GRAVITY")

  allocate(returngravity(1:pexportgrav+1)) ! +1 to avoid out-of-range on MPI_ISEND
                                     ! if zero length send

  ! Calculate where particles will go (in original IRECV order)
  ! This is because we may WAITANY in a different order, but
  ! gravity_export_return expects the order of grav_fromeach
  slot(0) = 1
  do i=1,lastrank-1
    slot(i) = slot(i-1) + grav_fromeach(2,i-1)
  end do

#ifdef DEBUG_EXPORT
  out_file = trim(adjustl(run_id))//".exportgravreturn."//trim(adjustl(MPI_ext))
  OPEN (1,file=out_file,status='unknown',form='formatted')
#endif

! Now calculate gravity forces for the exported particles
  do i=0,lastrank-1
    debug_timing("GRAV EXPORT WAITANY RECV")
    ! Wait until we finish receiving exported particles
    WAIT_TIME_MACRO
    call MPI_WAITANY(lastrank,recv(0:),t,MPI_STATUS_IGNORE,ierr)
    CALC_TIME_MACRO
    debug_timing("DO_EXPORT_GRAVITY")
    t = t - 1
    p = slot(t)
#ifdef DEBUG_EXPORT
    if (ierr /= 0) stop "ierr /= 0!"
#endif
    EXPORT_START_TIME_MACRO
!$OMP PARALLEL DO DEFAULT(SHARED)&
!$OMP& PRIVATE(rp,zo_p,agravmag_p,gpot_p,invhp,sphgpotp,agravp,potp,omp_p)
    do j=1,grav_fromeach(2,t)
    ! Make local copies of properties of particle p
      rp(1:NDIM) = receivegravity(t)%data(j)%r(1:NDIM)
#ifdef GRAD_H_SPH
      zo_p = receivegravity(t)%data(j)%zo
#else
      zo_p = 0._PR
#endif
#if defined(BH_TREE) && !defined(GEOMETRIC_MAC)
      agravmag_p = receivegravity(t)%data(j)%agravmag
#else
      agravmag_p = 0._PR
#endif
      invhp = real(1,DP)/receivegravity(t)%data(j)%h
      ! FIX ME
#if defined(BH_TREE)
      !call BHtreegravity(0,invhp,rp,agravp,potp,sphgpotp,zo_p,agravmag_p,gpot_p)
      ! should work apart from needing agravmag_p
      call BHgrav_accel(0,invhp,zo_p,agravmag_p,rp,agravp,potp)
#else
      !call direct_gravity(0,invhp,rp,agravp,potp,zo_p)
      call direct_sph_gravity(p,invhp,zo_p,rp,agravp,potp)
#endif
!$OMP CRITICAL
      omp_p = p
      p = p + 1
!$OMP END CRITICAL
      returngravity(omp_p)%p = receivegravity(t)%data(j)%p
#ifdef DEBUG_EXPORT
      returngravity(omp_p)%d = receivegravity(t)%data(j)%d
      returngravity(omp_p)%calc_in = rank
      write (1,*) "Particle in slot ", t, p, " from domain ", returngravity(p)%d
#endif
      returngravity(omp_p)%acc(1:NDIM) = agravp(1:NDIM)
      returngravity(omp_p)%pot = potp
#ifdef DEBUG_EXPORT
      write (1,*) returngravity(omp_p)
#endif
    end do
!$OMP END PARALLEL DO
    deallocate(receivegravity(t)%data)
    EXPORT_STOP_TIME_MACRO
  end do

  deallocate(receivegravity)

#ifdef DEBUG_EXPORT
  CLOSE (1)
#endif

  return
END SUBROUTINE do_export_gravity
