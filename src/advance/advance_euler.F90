! ADVANCE_EULER.F90
! C. P. Batty & D. A. Hubber - 19/3/2007
! First order Euler integration scheme.  (N.B. Disabled in sanitycheck.F90; 
! only retained for 'historical completeness')
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE advance_euler(p)
  use interface_module, only : check_boundary_conditions
  use particle_module
  use time_module
  use hydro_module
  implicit none

  integer, intent(in) :: p      ! Particle id

  integer :: dn                 ! Integer timestep since beginning of timestep
  integer :: nfull              ! Full integer timestep
  real(kind=PR) :: dt           ! Physical time since beginning of timestep
  real(kind=PR) :: rp(1:NDIM)   ! Local copy of position 
  real(kind=PR) :: vp(1:VDIM)   ! Local copy of velocity

  debug3("Advancing particle ",p)

! Work out integer and real time intervals from beginning of step 
  sph(p)%accdo = .false.
  dn       = n - sph(p)%nlast
  nfull    = 2**(level_step - sph(p)%nlevel)
  dt       = real(timestep,PR)*real(dn,PR)

! Full Euler integration step
  rp(1:NDIM) = sph(p)%r_old(1:NDIM) + sph(p)%v_old(1:NDIM)*dt
  vp(1:VDIM) = sph(p)%v_old(1:VDIM) + sph(p)%a(1:VDIM)*dt

  ! Check if particle has left periodic sphere/box
#if defined(BOUNDARY_CONDITIONS)
  call check_boundary_conditions(rp(1:NDIM),vp(1:VDIM))
#endif

! Record positions, velocities and other quantities in arrays for all cases 
#if !defined(STATIC_PARTICLES)
  sph(p)%r(1:NDIM) = rp(1:NDIM)
  sph(p)%v(1:VDIM) = vp(1:VDIM)
#endif
#if defined(VISC_TD)
  sph(p)%talpha = sph(p)%talpha_old + sph(p)%dalpha_dt*dt
  if (sph(p)%talpha < alpha_min) sph(p)%talpha = alpha_min
#endif

! If end of timestep, record as 'old' values
  if (dn == nfull) then
     sph(p)%accdo         = .true.
     sph(p)%nlast         = n
     sph(p)%laststep      = timestep*real(nfull,DP)
     sph(p)%r_old(1:NDIM) = rp(1:NDIM)
     sph(p)%v_old(1:VDIM) = vp(1:VDIM)
#if defined(VISC_TD)
     sph(p)%talpha_old    = sph(p)%talpha
#endif
  end if


! ----------------------------------------------------------------------------
#if defined(ENERGY_EQN)
  if (eos == "energy_eqn" .and. energy_eqn == "explicit") then 
     sph(p)%u = sph(p)%u_old + sph(p)%dudt*dt
     if (sph(p)%u < SMALL_NUMBER) &
          &sph(p)%u = sph(p)%u_old*exp(-sph(p)%u_old/sph(p)%dudt)
     if (dn == nfull) sph(p)%dudt_old = sph(p)%dudt
  end if
#endif
#if defined(ENERGY_EQN)
  if (dn == nfull) sph(p)%u_old = sph(p)%u
#endif


! Integrate entropy equation if selected
! ----------------------------------------------------------------------------
#if defined(ENTROPY_EQN) && defined(ENTROPIC_FUNCTION)
  if (eos == "entropy_eqn") then
     sph(p)%Aent = sph(p)%Aold + sph(p)%dAdt*dt
     sph(p)%u = (sph(p)%Aent*sph(p)%rho**(gamma - 1.0_PR))/(gamma - 1.0_PR)
  end if
#endif
#if defined(ENTROPIC_FUNCTION)
  if (dn == nfull) sph(p)%Aold = sph(p)%Aent
#endif

! Integrate magnetic fields
! ----------------------------------------------------------------------------
#if defined(MHD)
  sph(p)%B = sph(p)%B_old + sph(p)%dBdt*dt
  sph(p)%phi = sph(p)%phi_old + sph(p)%dphi_dt*dt
  if (dn == nfull) then
     sph(p)%B_old = sph(p)%B
     sph(p)%phi_old = sph(p)%phi
  end if
#endif


  return
END SUBROUTINE advance_euler
