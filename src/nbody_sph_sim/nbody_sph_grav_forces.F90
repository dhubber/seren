! NBODY_SPH_GRAV_FORCES.F90
! C. P. Batty & D. A. Hubber - 11/1/2007
! Calculates gravitational accelerations for all particles 
! (including sink particles) - controls calls to relevant subroutines.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_sph_grav_forces
  use interface_module, only : add_external_gravitational_force,BHgrav_accel,&
       &direct_sph_gravity,nbody_hermite4_direct_gravity,direct_sink_gravity
  use type_module
  use particle_module
  use time_module
  use timing_module
#if defined(OPENMP)
  use omp_lib
#endif
  implicit none

#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  integer :: acctot                   ! No. of particles on acc. step
  integer :: i                        ! Aux. particle counter
  integer :: p                        ! particle counter
  integer, allocatable :: acclist(:)  ! List of particles on acc. step
  real(kind=DP) :: adot(1:NDIM)       ! 'jerk' of star s
  real(kind=DP) :: agravp(1:NDIM)     ! Grav. acceleration of particle p
  real(kind=DP) :: atemp(1:NDIM)      ! Temp. grav. accel vector
  real(kind=DP) :: potp               ! Grav. potential of particle p
#if defined(OPENMP)
  integer :: chunksize                ! Data packet size for dynamic OpenMP
#endif

  debug2("Calculating grav. forces for all particles [nbody_sph_grav_forces.F90]")

  allocate(acclist(1:ptot))


! Loop over all self-gravitating particles and calculate forces
! ============================================================================
#if defined(GRAVITY)

! For multiple particle timesteps, first make a list of all self-gravitating
! SPH particles on an acceleration step, and then parallelize over that list.
  acctot = 0
  do p=pgravitystart,pgravityend
     if (sph(p)%accdo) then
        acctot = acctot + 1
        acclist(acctot) = p
     end if
  end do
#if defined(OPENMP)
  chunksize = int(CHUNKFRAC*real(acctot,PR)/real(omp_get_max_threads(),PR)) + 1
#endif
#if defined(TIMING)
  ngravcomp = ngravcomp + int(acctot,ILP)
#endif

  ! --------------------------------------------------------------------------
  if (acctot > 0) then
     debug_timing("GRAVITY_FORCES")

     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) &
     !$OMP PRIVATE(adot,agravp,atemp,p,potp)
     do i=1,acctot
        p = acclist(i)
        agravp(1:NDIM) = 0.0_DP
        atemp(1:NDIM) = 0.0_DP

        ! SPH self-gravity
        ! --------------------------------------------------------------------
#if defined(BH_TREE)
        call BHgrav_accel_jerk(p,1.0_PR/sph(p)%h,sph(p)%r(1:NDIM),&
             &sph(p)%v(1:NDIM),agravp(1:NDIM),adot(1:NDIM),potp)
#else    
        call sph_hermite4_direct_gravity(p,1.0_DP/real(sph(p)%h,DP),&
             &real(sph(p)%r(1:NDIM),DP),real(sph(p)%v(1:NDIM),DP),&
             &agravp(1:NDIM),adot(1:NDIM),potp)
#endif
        agravp(1:NDIM) = atemp(1:NDIM)
        sph(p)%gpot = potp

        ! Star gravity contribution
        ! --------------------------------------------------------------------
#if defined(NBODY_HERMITE4)
        call nbody_hermite4_direct_gravity(-p,1.0_DP/real(sph(p)%h,DP),&
             &real(sph(p)%r(1:NDIM),DP),real(sph(p)%v(1:NDIM),DP),&
             &agravp(1:NDIM),adot(1:NDIM),potp)
        agravp(1:NDIM) = agravp(1:NDIM) + atemp(1:NDIM)
        sph(p)%gpot = sph(p)%gpot + real(potp,PR)
!#elif defined(NBODY_HERMITE6)
#endif

        ! Compute values for tree MAC
        ! --------------------------------------------------------------------
#if defined(GRAVITY) && defined(BH_TREE) & defined(EIGEN_MAC)
        sph(p)%agravmag = sph(p)%gpot
#elif defined(GRAVITY) && defined(BH_TREE) & !defined(GEOMETRIC_MAC)
        sph(p)%agravmag = real(sqrt(dot_product(agravp(1:NDIM,p),&
             &agravp(1:NDIM,p))),PR)
#endif

        ! External gravitational force contribution
        ! --------------------------------------------------------------------
#if defined(EXTERNAL_FORCE)
        call add_external_gravitational_force(sph(p)%r(1:NDIM),&
             &atemp(1:NDIM),potp)
        agravp(1:NDIM) = agravp(1:NDIM) + atemp(1:NDIM)
        sph(p)%gpot = sph(p)%gpot + real(potp,PR)
#endif

        sph(p)%a(1:NDIM) = sph(p)%a(1:NDIM) + real(agravp(1:NDIM),PR)
#if defined(DEBUG_FORCES)
        sph(p)%sph(p)%a_grav(1:NDIM) = real(agravp(1:NDIM),PR)
#endif
     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------

  end if
  ! --------------------------------------------------------------------------

#endif
! ============================================================================

  deallocate(acclist)

#endif

  return
END SUBROUTINE nbody_sph_grav_forces
