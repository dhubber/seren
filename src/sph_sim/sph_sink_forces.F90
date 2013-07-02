! SPH_SINK_FORCES.F90
! C. P. Batty & D. A. Hubber - 2/8/2010
! Calculates gravitational accelerations for all sink particles
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_sink_forces
  use interface_module, only : add_external_gravitational_force,&
       &BHgrav_accel,direct_sink_gravity,direct_sph_gravity
  use type_module
  use particle_module
  use sink_module
  use time_module
#if defined(OPENMP)
  use omp_lib
#endif
#if defined(CELL_WALK)
  use tree_module
#endif
#if defined(USE_MPI)
  use mpi_communication_module
#endif
  implicit none

#if defined(SINKS)
  integer :: s                        ! sink counter
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  real(kind=DP) :: agravs(1:NDIM)     ! Grav. acceleration of sink s
  real(kind=DP) :: atemp(1:NDIM)      ! Aux. grav acceleration vector
  real(kind=DP) :: pots               ! Grav. potential of sink s
#endif
#if defined(EXTERNAL_FORCE)
  real(kind=DP) :: adottemp(1:NDIM)   ! 'Jerk' aux. variable
#endif

  debug2("Calculating grav. forces for all sinks [sink_grav_forces.F90]")


! Zero acceleration arrays here for now
  if (accdo_sinks .and. stot > 0) then
     do s=1,stot
#if defined(USE_MPI)
        if (sink(s)%domain /= rank) cycle
#endif
        sink(s)%a(1:NDIM) = 0.0_PR
        sink(s)%gpot = 0.0_PR
#if defined(DEBUG_FORCES)
        sink(s)%ahydro(1:NDIM) = 0.0_PR
        sink(s)%agrav(1:NDIM) = 0.0_PR
#endif
     end do
  end if


! Loop over all sink particles and calculate grav. forces from SPH particles
! ============================================================================
#if defined(GRAVITY) && defined(SELF_GRAVITY)
  if (accdo_sinks .and. stot > 0) then
     debug2("Calculating forces on all sink particles")
     debug_timing("SINK_FORCES")

     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(agravs,atemp,pots)
     do s=1,stot
#if defined(USE_MPI)
        if (sink(s)%domain /= rank) cycle
#endif

        ! Contribution from remote domain (particles only)
        ! --------------------------------------------------------------------
#if defined(USE_MPI)
        agravs(1:NDIM) = sink(s)%remote_agrav
        sink(s)%gpot = sink(s)%remote_gpot
#else
        agravs(1:NDIM) = 0.0_DP
        atemp(1:NDIM) = 0.0_DP
#endif

        ! SPH self-gravity
        ! --------------------------------------------------------------------
#if defined(SELF_GRAVITY) && defined(BH_TREE)
        call BHgrav_accel(-s,1.0_PR/sink(s)%h,0.0_PR,&
             &sink(s)%agravmag,sink(s)%r(1:NDIM),atemp(1:NDIM),pots)
#elif defined(SELF_GRAVITY) && defined(BINARY_TREE)
        call binary_gravacc(-s,1.0_PR/sink(s)%h,&
             &sink(s)%r(1:NDIM),atemp(1:NDIM),pots)
#elif defined(SELF_GRAVITY)   
        call direct_sph_gravity(-s,1.0_PR/sink(s)%h,0.0_PR,&
             &sink(s)%r(1:NDIM),atemp(1:NDIM),pots)
#endif
#if defined(USE_MPI)
        agravs(1:NDIM) = agravs(1:NDIM) + atemp(1:NDIM)
        sink(s)%gpot = sink(s)%gpot + real(pots,PR)
#else
        agravs(1:NDIM) = atemp(1:NDIM)
        sink(s)%gpot = real(pots,PR)
#endif

        ! Sink gravity contribution
        ! --------------------------------------------------------------------
#if defined(SINKS)
        call direct_sink_gravity(-s,sink(s)%h,&
             &sink(s)%r(1:NDIM),atemp(1:NDIM),pots)
        agravs(1:NDIM) = agravs(1:NDIM) + real(atemp(1:NDIM),PR)
        sink(s)%gpot = sink(s)%gpot + real(pots,PR)
#endif

        ! Compute values for tree MAC
        ! --------------------------------------------------------------------
#if defined(BH_TREE) & defined(EIGEN_MAC)
        sink(s)%agravmag = sink(s)%gpot
#elif defined(BH_TREE) & !defined(GEOMETRIC_MAC)
        sink(s)%agravmag = real(sqrt(dot_product(agravs(1:NDIM),&
             &agravs(1:NDIM))),PR)
#endif

        ! External gravitational force contribution
        ! --------------------------------------------------------------------
#if defined(EXTERNAL_FORCE)
        call add_external_gravitational_force(real(sink(s)%r(1:NDIM),DP),&
             &real(sink(s)%v(1:NDIM),DP),atemp(1:NDIM),adottemp(1:NDIM),pots)
        agravs(1:NDIM) = agravs(1:NDIM) + real(atemp(1:NDIM),PR)
        sink(s)%gpot = sink(s)%gpot + real(pots,PR)
#endif

        sink(s)%a(1:NDIM) = sink(s)%a(1:NDIM) + real(agravs(1:NDIM),PR)
#if defined(DEBUG_FORCES)
        sink(s)%agrav(1:NDIM) = real(agravs(1:NDIM),PR)
#endif

     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------

  end if
  ! --------------------------------------------------------------------------

#endif
! ============================================================================

#endif

  return
END SUBROUTINE sph_sink_forces
