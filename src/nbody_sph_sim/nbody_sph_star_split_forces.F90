! NBODY_SPH_STAR_SPLIT_FORCES.F90
! D. A. Hubber - 17/09/2010
! ..
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_sph_star_split_forces
  use Nbody_module
  use time_module, only : time
  use sink_module, only : stot
  implicit none

  integer :: s                        ! star counter
  real(kind=DP) :: adot(1:NDIM)       ! 'jerk' of star s
  real(kind=DP) :: agrav(1:NDIM)      ! grav. accel of star s
  real(kind=DP) :: pot                ! potential of star s
#if defined(NBODY_HERMITE6)
  real(kind=DP) :: adot2(1:NDIM)      ! 'snap' of star s
#endif

  debug2("Calculating grav acclerations of stars [nbody_sph_star_forces.F90]")
  debug_timing("STAR_GRAVITY")


! First, calculate gravitational force due to all SPH particles
! ----------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(adot,agrav,pot) IF (stot > 100)
  do s=1,stot
     if (.not. star(s)%accsph) cycle
!     if (.not. star(s)%accdo) cycle
!     if ((.not. star(s)%accsph) .and. (.not. star(s)%accdo)) cycle
#if defined(BH_TREE)
     call BHgrav_accel(-s,1.0_PR/star(s)%h,&
          &star(s)%r(1:NDIM),agrav(1:NDIM),pot)
#else
     call direct_sph_gravity(-s,1.0_PR/star(s)%h,&
          &star(s)%r(1:NDIM),agrav(1:NDIM),pot)
#endif
     star(s)%asph(1:NDIM)    = agrav(1:NDIM)
     star(s)%adotsph(1:NDIM) = 0.0_PR !adot(1:NDIM)
     star(s)%gpotsph         = pot
     star(s)%lastsph         = time
     star(s)%accsph          = .false.
  end do
!$OMP END PARALLEL DO


! Calculate all terms for 4th-order Hermite scheme
! ----------------------------------------------------------------------------
#if defined(NBODY_HERMITE4)
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(adot,agrav,pot) IF (stot > 100)
  do s=1,stot
     if (.not. star(s)%accdo) cycle
     call nbody_hermite4_direct_gravity(s,star(s)%h,star(s)%r(1:NDIM),&
          &star(s)%v(1:NDIM),agrav(1:NDIM),adot(1:NDIM),pot)
     star(s)%astar(1:NDIM)    = agrav(1:NDIM)
     star(s)%adotstar(1:NDIM) = adot(1:NDIM)
     star(s)%gpotstar         = pot
  end do
!$OMP END PARALLEL DO


! Calculate all terms for 6th-order Hermite scheme
! ----------------------------------------------------------------------------
#elif defined(NBODY_HERMITE6)

#endif
! ----------------------------------------------------------------------------


! Add together star and sph contributions to forces
! Also if required, calculate agravmag for all stars
! ----------------------------------------------------------------------------
  do s=1,stot
     star(s)%a(1:NDIM) = star(s)%astar(1:NDIM) + star(s)%asph(1:NDIM) &
          & !+ star(s)%adotsph(1:NDIM)*(time - star(s)%lastsph)
     star(s)%adot(1:NDIM) = star(s)%adotstar(1:NDIM) + star(s)%adotsph(1:NDIM)
     star(s)%gpot = star(s)%gpotstar + star(s)%gpotsph
     star(s)%gpe = star(s)%m*star(s)%gpot
#if defined(GRAVITY) && defined(BH_TREE) & defined(EIGEN_MAC)
     star(s)%agravmag = star(s)%gpot
#elif defined(GRAVITY) && defined(BH_TREE) & !defined(GEOMETRIC_MAC)
     star(s)%agravmag = &
          & sqrt(dot_product(star(s)%a(1:NDIM),star(s)%a(1:NDIM)))
#endif
  end do
! ----------------------------------------------------------------------------


  return
END SUBROUTINE nbody_sph_star_split_forces
