! ADVANCE_BOUNDARY_PARTICLE.F90
! D. A. Hubber - 2/2/2009
! Advances the position and velocity of boundary particle p.
! Currently hard-wired for static boundary particles.  
! (N.B. Any frog-marching code should be put here.)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE advance_boundary_particle(p)
  use particle_module
  use time_module
  implicit none

  integer, intent(in) :: p     ! Boundary particle id

  integer(kind=ILP) :: dn      ! Integer timestep since beginning of timestep
  integer(kind=ILP) :: nhalf   ! Half the full integer timestep
  integer(kind=ILP) :: nfull   ! Full integer timestep
#if defined(BOUNDARY_FROG_MARCH)
  real(kind=PR) :: dt          ! Physical time since beginning of timestep
#endif

! Work out integer and real time intervals from beginning of step 
  sph(p)%accdo = .false.
  dn       = n - sph(p)%nlast
  nfull    = 2**(level_step - sph(p)%nlevel)
  nhalf    = nfull / 2_ILP
#if defined(BOUNDARY_FROG_MARCH)
  dt       = real(timestep,PR)*real(dn,PR)
#endif

! Only calculate forces at the half timestep
#if defined(LEAPFROG_KDK) || defined(RUNGE_KUTTA2)
  if (dn == nhalf) sph(p)%accdo = .true.
#endif

! If end of timestep, record as 'old' values
#if defined(RUNGE_KUTTA2) || defined(EULER)
  if (dn == nfull) then
     sph(p)%accdo    = .true.
     sph(p)%nlast    = n
     sph(p)%laststep = timestep*real(nfull,DP)
  end if
#endif

  return
END SUBROUTINE advance_boundary_particle
