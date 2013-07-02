! SPH_ADVANCE.F90
! C. P. Batty & D. A. Hubber - 23/8/2007
! Calls routines to advance positions and velocities of all particles 
! depending the integration scheme employed.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_advance
  use interface_module, only : advance_boundary_particle,advance_euler,&
       &advance_leapfrog_dkd,advance_leapfrog_kdk,&
       &advance_predictor_corrector,advance_runge_kutta
  use particle_module, only : ptot
  use type_module, only : pboundary,phydrostart
  implicit none

  integer :: p          ! Particle counter

  debug_timing("ADVANCE")


! If using a non-standard energy integrator, call it here.
! ----------------------------------------------------------------------------
#if defined(EXPONENTIAL_COOLING_HEATING)
  debug2("[sph_advance.F90]")

!$OMP PARALLEL DO DEFAULT(SHARED)
  do p=phydrostart,ptot
     call cooling_exponential_integration(p)
  end do
!$OMP END PARALLEL DO
#endif


! First, 'frog-march' boundary particles
! ----------------------------------------------------------------------------
  if (pboundary > 0) then
     debug2("Advancing positions of boundary particles [sph_advance.F90]")
     do p=1,pboundary
        call advance_boundary_particle(p)
     end do
  end if
! ----------------------------------------------------------------------------


! Next, advance hydro particles depending on choice of integration scheme
! ----------------------------------------------------------------------------
  debug2("Advancing positions of all SPH particles [sph_advance.F90]")

!$OMP PARALLEL DO DEFAULT(SHARED)
  do p=phydrostart,ptot
#if defined(EULER)
    call advance_euler(p)
#elif defined(RUNGE_KUTTA2)
    call advance_runge_kutta(p)
#elif defined(LEAPFROG_KDK)
    call advance_leapfrog_kdk(p)
#elif defined(LEAPFROG_DKD)
    call advance_leapfrog_dkd(p)
#endif
  end do
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------------


  return
END SUBROUTINE sph_advance
