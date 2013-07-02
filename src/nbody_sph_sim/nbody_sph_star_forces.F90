! NBODY_SPH_STAR_FORCES.F90
! D. A. Hubber - 17/09/2010
! ..
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_sph_star_forces
  use interface_module, only : BHgrav_accel,BHgrav_accel_jerk,&
       &direct_sph_gravity,direct_sph_hermite4_gravity
  use sink_module, only : stot
  use Nbody_module
  implicit none

#if defined(GRAVITY)
  integer :: s                        ! star counter
  real(kind=DP) :: adot(1:NDIM)       ! 'jerk' of star s
  real(kind=DP) :: agrav(1:NDIM)      ! grav. accel of star s
  real(kind=DP) :: pot                ! potential of star s
#if defined(NBODY_HERMITE6)
  real(kind=DP) :: adot2(1:NDIM)      ! 'snap' of star s
#endif

  debug2("Calculating grav. acclerations of stars [nbody_sph_star_forces.F90]")
  debug_timing("STAR_GRAVITY")

! First, calculate gravitational accel. due to all SPH particles
! ----------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(adot,agrav,pot) IF (stot > 100)
  do s=1,stot
     if (.not. star(s)%accdo) cycle
!#if defined(BH_TREE)
!     call BHgrav_accel(-s,1.0_PR/real(star(s)%h,PR),0.0_PR,&
!          &real(star(s)%agravmag,PR),real(star(s)%r(1:NDIM),PR),&
!          &agrav(1:NDIM),pot)
!     adot(1:NDIM) = 0.0_DP
!#else
!     call direct_sph_gravity(-s,1.0_PR/real(star(s)%h,PR),&
!          &0.0_PR,real(star(s)%r(1:NDIM),PR),agrav(1:NDIM),pot)
!     adot(1:NDIM) = 0.0_DP
!#endif
!#if defined(BH_TREE)
!     call BHgrav_accel_jerk(s,star(s)%h,star(s)%r(1:NDIM),&
!          &star(s)%v(1:NDIM),agrav(1:NDIM),adot(1:NDIM),pot)
!#else
     call direct_sph_hermite4_gravity(s,star(s)%h,star(s)%r(1:NDIM),&
          &star(s)%v(1:NDIM),agrav(1:NDIM),adot(1:NDIM),pot)
!#endif
     star(s)%a(1:NDIM)    = star(s)%a(1:NDIM) + agrav(1:NDIM)
     star(s)%adot(1:NDIM) = star(s)%adot(1:NDIM) + adot(1:NDIM)
     star(s)%gpot         = star(s)%gpot + pot
     star(s)%gpe          = star(s)%gpe + star(s)%m*pot
  end do
!$OMP END PARALLEL DO


! Calculate all terms for 4th-order Hermite scheme
! ----------------------------------------------------------------------------
#if defined(NBODY_HERMITE4) || defined(NBODY_LEAPFROG_KDK)
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(adot,agrav,pot) IF (stot > 100)
  do s=1,stot
     if (.not. star(s)%accdo) cycle
     call nbody_hermite4_direct_gravity(s,star(s)%h,star(s)%r(1:NDIM),&
          &star(s)%v(1:NDIM),agrav(1:NDIM),adot(1:NDIM),pot)
     star(s)%a(1:NDIM)    = star(s)%a(1:NDIM) + agrav(1:NDIM)
     star(s)%adot(1:NDIM) = star(s)%adot(1:NDIM) + adot(1:NDIM)
     star(s)%gpot         = star(s)%gpot + pot
     star(s)%gpe          = star(s)%gpe + star(s)%m*pot
  end do
!$OMP END PARALLEL DO


! Calculate all terms for 6th-order Hermite scheme
! ----------------------------------------------------------------------------
#elif defined(NBODY_HERMITE6)

#endif
! ----------------------------------------------------------------------------


! Calculate agravmag for all sinks
! ----------------------------------------------------------------------------
  do s=1,stot
#if defined(GRAVITY) && defined(BH_TREE) & defined(EIGEN_MAC)
     star(s)%agravmag = star(s)%gpot
#elif defined(GRAVITY) && defined(BH_TREE) && !defined(GEOMETRIC_MAC)
     star(s)%agravmag = &
          & sqrt(dot_product(star(s)%a(1:NDIM),star(s)%a(1:NDIM)))
#endif
  end do
! ----------------------------------------------------------------------------


! External gravitational force contribution
! --------------------------------------------------------------------
#if defined(EXTERNAL_FORCE)
  do s=1,stot
     call add_external_gravitational_force(star(s)%r(1:NDIM),&
          &star(s)%v(1:NDIM),agrav(1:NDIM),adot(1:NDIM),pot)
     star(s)%a(1:NDIM)    = star(s)%a(1:NDIM) + agrav(1:NDIM)
     star(s)%adot(1:NDIM) = star(s)%adot(1:NDIM) + adot(1:NDIM)
     star(s)%gpot         = star(s)%gpot + pot
     star(s)%gpe          = star(s)%gpe + star(s)%m*pot
  end do
#endif
! ----------------------------------------------------------------------------
#endif

  return
END SUBROUTINE nbody_sph_star_forces
