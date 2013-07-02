! NBODY_CORRECTION_TERMS.F90
! D. A. Hubber - 13/9/2010
! Calculate correction terms for selected N-body integration scheme at the 
! end of the star's current timestep.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_correction_terms
  use sink_module, only : stot
  use definitions
  use time_module
  use Nbody_module
  implicit none

  integer :: s                      ! Star counter
  real(kind=DP) :: dt               ! Physical timestep
  real(kind=DP) :: dt2              ! dt*dt
  real(kind=DP) :: dt3              ! dt*dt*dt

  debug2("Calculating correction terms [nbody_correction_terms.F90]")
  debug_timing("NBODY_CORRECTION")

! 4th-order Hermite correction terms
! ----------------------------------------------------------------------------
#if defined(NBODY_HERMITE4) || defined(NBODY_LEAPFROG_KDK)
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(dt,dt2,dt3) IF (stot > 500)
  do s=1,stot
     if (.not. star(s)%accdo) cycle
#if defined(BINARY_COM_MOTION)
     if (star(s)%binid > 0) cycle
#endif

     dt  = star(s)%laststep
     dt2 = dt*dt
     dt3 = dt2*dt

     star(s)%a2dot0(1:NDIM) = (-6.0_DP*(star(s)%a0(1:NDIM)-star(s)%a(1:NDIM))&
          & - dt*(4.0_DP*star(s)%adot0(1:NDIM) + 2.0_DP*star(s)%adot(1:NDIM)))&
          & / dt2
     star(s)%a3dot(1:NDIM) = (12.0_DP*(star(s)%a0(1:NDIM) - star(s)%a(1:NDIM))&
          & + 6.0_DP*dt*(star(s)%adot0(1:NDIM) + star(s)%adot(1:NDIM))) / dt3
#if defined(FORCE_SPLITTING)
     star(s)%a2dot0star(1:NDIM) = &
          & (-6.0_DP*(star(s)%a0star(1:NDIM) - star(s)%astar(1:NDIM))&
          & - dt*(4.0_DP*star(s)%adot0star(1:NDIM) &
          & + 2.0_DP*star(s)%adotstar(1:NDIM))) / dt2
     star(s)%a3dotstar(1:NDIM) = (12.0_DP*(star(s)%a0star(1:NDIM) - &
          & star(s)%astar(1:NDIM)) + 6.0_DP*dt*(star(s)%adot0star(1:NDIM) + &
          & star(s)%adotstar(1:NDIM))) / dt3
     star(s)%a2dotstar(1:NDIM) = &
          & star(s)%a2dot0star(1:NDIM) + star(s)%a3dotstar(1:NDIM)*dt
#endif

#if defined(NBODY_LEAPFROG_KDK)
     star(s)%r(1:NDIM) = star(s)%rold(1:NDIM) &
          & + star(s)%vold(1:NDIM)*dt + 0.5_DP*star(s)%a0(1:NDIM)*dt2
     star(s)%v(1:NDIM) = star(s)%vold(1:NDIM) &
          & + 0.5_DP*(star(s)%a0(1:NDIM) + star(s)%a(1:NDIM))*dt
#elif defined(NBODY_HERMITE4) && defined(NBODY_HERMITE4_TS)
     star(s)%r(1:NDIM) = star(s)%rold(1:NDIM) &
          & + 0.5_DP*(star(s)%v(1:NDIM) + star(s)%vold(1:NDIM))*dt &
          & - 0.1_DP*(star(s)%a(1:NDIM) - star(s)%a0(1:NDIM))*dt2 &
          & + (star(s)%adot(1:NDIM) + star(s)%adot0(1:NDIM))*dt3/120.0_DP
     star(s)%v(1:NDIM) = star(s)%vold(1:VDIM) &
          & + 0.5_DP*(star(s)%a(1:VDIM) + star(s)%a0(1:VDIM))*dt &
          & - (star(s)%adot(1:VDIM) - star(s)%adot0(1:VDIM))*dt2/12.0_DP
#elif defined(NBODY_HERMITE4)
     star(s)%r(1:NDIM) = star(s)%rold(1:NDIM) + star(s)%vold(1:NDIM)*dt &
          & + 0.5_DP*star(s)%a0(1:NDIM)*dt2 + ONESIXTH_DP*star(s)%adot0*dt3 &
          & + star(s)%a2dot0(1:NDIM)*dt2*dt2/24.0_DP &
          & + star(s)%a3dot(1:NDIM)*dt3*dt2/120.0_DP
     star(s)%v(1:NDIM) = star(s)%vold(1:VDIM) + star(s)%a0(1:VDIM)*dt &
          & + 0.5_DP*star(s)%adot0(1:VDIM)*dt2 &
          & + star(s)%a2dot0(1:NDIM)*dt3/6.0_DP &
          & + star(s)%a3dot(1:NDIM)*dt2*dt2/24.0_DP
#endif
     star(s)%a2dot(1:NDIM) = star(s)%a2dot0(1:NDIM) + star(s)%a3dot(1:NDIM)*dt

#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(star(s)%r(1:NDIM),star(s)%v(1:VDIM))
#endif

#if defined(DEBUG_HERMITE_CORRECTION_TERMS)
     write(6,*) "CORRECTION TERMS : ",s,dt
     write(6,*) "a0     : ",star(s)%a0(1:NDIM)
     write(6,*) "a      : ",star(s)%a(1:NDIM)
     write(6,*) "adot0  : ",star(s)%adot0(1:NDIM)
     write(6,*) "adot   : ",star(s)%adot(1:NDIM)
     write(6,*) "a2dot0 : ",star(s)%a2dot0(1:NDIM)
     write(6,*) "a2dot  : ",star(s)%a2dot(1:NDIM)
     write(6,*) "a3dot  : ",star(s)%a3dot(1:NDIM)
#endif

  end do
!$OMP END PARALLEL DO


! 6th-order Hermite correction terms
! ----------------------------------------------------------------------------
#elif defined(NBODY_HERMITE6)


#endif
! ----------------------------------------------------------------------------


  return
END SUBROUTINE nbody_correction_terms
