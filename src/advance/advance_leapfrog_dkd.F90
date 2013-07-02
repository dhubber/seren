! ADVANCE_LEAPFROG_DKD.F90
! C. P. Batty & D. A. Hubber - 19/3/2007
! Second order leapfrog drift-kick-drift integration scheme.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE advance_leapfrog_dkd(p)
  use interface_module, only : check_boundary_conditions
  use particle_module
  use time_module
  use hydro_module
  implicit none

  integer, intent(in) :: p      ! Particle id

  integer :: dn                 ! Integer timestep since beginning of timestep
  integer :: nfull              ! Full integer timestep
  integer :: nhalf              ! Half the full integer timestep
  real(kind=PR) :: rp(1:NDIM)   ! Local copy of position 
  real(kind=PR) :: vp(1:VDIM)   ! Local copy of velocity
  real(kind=PR) :: dt           ! Physical time since beginning of timestep

  debug3("Advancing particle ",p)

! Work out integer and real time intervals from beginning of step 
  sph(p)%accdo = .false.
  dn       = n - sph(p)%nlast
  nfull    = 2**(level_step - sph(p)%nlevel)
  nhalf    = nfull / 2
  dt       = real(timestep,PR)*real(dn,PR)


! Advance particles that are below the full timestep
! ----------------------------------------------------------------------------
  if (dn < nfull) then

     rp(1:NDIM) = sph(p)%r_old(1:NDIM) + sph(p)%v_old(1:NDIM)*dt
     vp(1:VDIM) = sph(p)%v_old(1:VDIM) + sph(p)%a(1:VDIM)*dt

     ! Check if particle has left periodic sphere/box
#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(rp(1:NDIM),vp(1:VDIM))
#endif

     ! Only calculate forces at the half timestep
     if (dn == nhalf) sph(p)%accdo = .true.


! Else advance particles that are at or beyond the half-step
! ----------------------------------------------------------------------------
  else

     vp(1:VDIM) = sph(p)%v_old(1:VDIM) + sph(p)%a(1:VDIM)*dt
     rp(1:NDIM) = sph(p)%r_old(1:NDIM) + &
          &0.5_PR*(sph(p)%v_old(1:NDIM) + vp(1:NDIM))*dt

     ! Check if particle has left periodic sphere/box
#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(rp(1:NDIM),vp(1:VDIM))
#endif

  end if
! ----------------------------------------------------------------------------

! Record positions, velocities and other quantities in arrays for all cases 
  sph(p)%r(1:NDIM) = rp(1:NDIM)
  sph(p)%v(1:VDIM) = vp(1:VDIM)
#if defined(ENTROPIC_FUNCTION) && defined(ENTROPY_EQN)
  sph(p)%Aent = sph(p)%Aold + sph(p)%dAdt*dt
  sph(p)%u = (sph(p)%Aent*sph(p)%rho**(gamma - 1.0_PR))/(gamma - 1.0_PR)
#elif defined(EXPLICIT_ENERGY_EQN) && defined(ENERGY_EQN) && !defined(RAD_WS)
  sph(p)%u = sph(p)%u_old + sph(p)%dudt*dt
  if (sph(p)%u < SMALL_NUMBER) sph(p)%u = sph(p)%u_old*exp(-sph(p)%u_old/sph(p)%dudt)
#endif
#if defined(VISC_TD)
  sph(p)%talpha = sph(p)%talpha_old + sph(p)%dalpha_dt*dt
  if (sph(p)%talpha < alpha_min) sph(p)%talpha = alpha_min
#endif

! If end of timestep, record as 'old' values
  if (dn == nfull) then
     sph(p)%nlast        = n
     sph(p)%laststep     = real(nfull,DP)*timestep
     sph(p)%r_old(1:NDIM) = rp(1:NDIM)
     sph(p)%v_old(1:VDIM) = vp(1:VDIM)
#if defined(ENTROPIC_FUNCTION) && defined(ENTROPY_EQN)
     sph(p)%Aold         = sph(p)%Aent
#elif defined(ENERGY_EQN)
     sph(p)%u_old        = sph(p)%u
#endif
#if defined(VISC_TD)
     sph(p)%talpha_old   = sph(p)%talpha
#endif
  end if

  return
END SUBROUTINE advance_leapfrog_dkd
