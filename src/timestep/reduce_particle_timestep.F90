! REDUCE_PARTICLE_TIMESTEP.F90
! D. A. Hubber - 9/9/2010
! Reduce the timestep of SPH particle p due to too large a difference in 
! timesteps between neighbours.  Sets to timestep level, lnew
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE reduce_particle_timestep(p,lnew)
  use particle_module
  use hydro_module
  use time_module

  integer, intent(in) :: p                ! Particle id
  integer(kind=ILP), intent(in) :: lnew   ! New timestep level

#if defined(DEBUG_REDUCE_TIMESTEP)
  write(6,*) "Reducing timestep for particle ",p
  write(6,*) "nminneib :",sph(p)%nminneib,level_step
  write(6,*) "lold :",sph(p)%nlevel,"    lnew :",lnew
  write(6,*) "nlast :",sph(p)%nlast,"    n :",n
#endif

#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
  sph(p)%nminneib = lnew
#endif
  if (n == sph(p)%nlast) return
  if (n /= sph(p)%nlast) sph(p)%laststep = timestep*real(n - sph(p)%nlast,DP)
#if defined(LEAPFROG_KDK)
  sph(p)%a_old(1:VDIM) = sph(p)%a(1:VDIM)
  if (n /= sph(p)%nlast) sph(p)%r(1:NDIM) = &
       &sph(p)%r_old(1:NDIM) + sph(p)%v_old(1:NDIM)*real(sph(p)%laststep,PR) &
       & + 0.5_PR*sph(p)%a(1:NDIM)*real(sph(p)%laststep,PR)**2
#endif
  sph(p)%r_old(1:NDIM) = sph(p)%r(1:NDIM)
  sph(p)%v_old(1:VDIM) = sph(p)%v(1:VDIM)
#if defined(HYDRO) && defined(ENTROPIC_FUNCTION) && defined(ENERGY_EQN)
  sph(p)%Aold = sph(p)%Aent
#elif defined(HYDRO) && defined(ENERGY_EQN)
  sph(p)%u_old = sph(p)%u
#endif
#if defined(HYDRO) && defined(VISC_TD)
  sph(p)%talpha_old = sph(p)%talpha
#endif

  if (n /= sph(p)%nlast) sph(p)%nlevel = lnew
  if (n /= sph(p)%nlast) sph(p)%nlast = n

! Calculate new accelerations for all particles to ensure we do not keep 
! old and invalid accelerations.
  sph(p)%accdo = .true.

  return
END SUBROUTINE reduce_particle_timestep
