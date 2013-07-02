! NBODY_GRAV_FORCES.F90
! D. A. Hubber - 12/09/2010
! ..
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_grav_forces
  use sink_module, only : stot
  use Nbody_module
  implicit none

  integer :: s                        ! star counter
  real(kind=DP) :: adot(1:NDIM)       ! 'jerk' of star s
  real(kind=DP) :: agrav(1:NDIM)      ! grav. accel of star s
  real(kind=DP) :: pot                ! potential of star s
#if defined(NBODY_HERMITE6)
  real(kind=DP) :: adot2(1:NDIM)      ! 'snap' of star s
#endif

  debug2("Calculating grav acclerations of sinks [nbody_grav_forces.F90]")
  debug_timing("NBODY_GRAV_FORCES")

! Calculate all terms for 4th-order Hermite scheme
! ----------------------------------------------------------------------------
#if defined(NBODY_HERMITE4)
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(adot,agrav,pot) IF (stot > 100)
  do s=1,stot
     if (.not. star(s)%accdo) cycle
     call nbody_hermite4_direct_gravity(s,star(s)%h,star(s)%r(1:NDIM),&
          &star(s)%v(1:NDIM),agrav(1:NDIM),adot(1:NDIM),pot)
     star(s)%a(1:NDIM)    = agrav(1:NDIM)
     star(s)%adot(1:NDIM) = adot(1:NDIM)
     star(s)%gpe          = star(s)%m*pot

#if defined(EXTERNAL_FORCE)
     call add_external_gravitational_force(star(s)%r(1:NDIM),&
          &star(s)%v(1:NDIM),agrav(1:NDIM),adot(1:NDIM),pot)
     star(s)%a(1:NDIM)    = star(s)%a(1:NDIM) + agrav(1:NDIM)
     star(s)%adot(1:NDIM) = star(s)%adot(1:NDIM) + adot(1:NDIM)
     star(s)%gpe          = star(s)%gpe + star(s)%m*pot
#endif

  end do
!$OMP END PARALLEL DO

! Calculate all terms for 6th-order Hermite scheme
! ----------------------------------------------------------------------------
#elif defined(NBODY_HERMITE6)

#endif
! ----------------------------------------------------------------------------

  return
END SUBROUTINE nbody_grav_forces
