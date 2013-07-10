! ADVANCE_RUNGE_KUTTA.F90
! C. P. Batty & D. A. Hubber - 19/3/2007
! Second order Runge-Kutta integration scheme.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE advance_runge_kutta(p)
  use particle_module
  use time_module
  use hydro_module
  use type_module
  implicit none

  integer, intent(in) :: p     ! Particle id

  character(len=256) :: eos    ! Equation of state for particle p
  integer :: dn                ! Integer timestep since beginning of timestep
  integer :: nfull             ! Full integer timestep
  real(kind=DP) :: dt          ! Physical time since beginning of timestep
  real(kind=PR) :: rp(1:NDIM)  ! Local copy of position 
  real(kind=PR) :: vp(1:VDIM)  ! Local copy of velocity

  debug3("Advancing particle ",p)

! Work out integer and real time intervals from beginning of step
  sph(p)%accdo = .false.
  dn       = n - sph(p)%nlast
  nfull    = 2**(level_step - sph(p)%nlevel)
  dt       = timestep*real(dn,DP)
  eos      = typeinfo(sph(p)%ptype)%eos


! Advance particles that are below or at the half timestep
! ----------------------------------------------------------------------------
  if (dn <= nfull/2) then

     ! First half of the Runge-Kutta integration step
     rp(1:NDIM) = sph(p)%r_old(1:NDIM) + sph(p)%v_old(1:NDIM)*dt
     vp(1:VDIM) = sph(p)%v_old(1:VDIM) + sph(p)%a(1:VDIM)*dt

     ! Check if particle has left periodic sphere/box
#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(rp(1:NDIM),vp(1:VDIM))
#endif

     ! Record velocities if at half timestep
     if (dn == nfull/2) then
        sph(p)%accdo = .true.
        sph(p)%v_half = vp(1:VDIM)
     endif

! Else advance particles that are beyond the half timestep
! ----------------------------------------------------------------------------
  else

     ! Second half of the Runge-Kutta integration step
     rp(1:NDIM) = sph(p)%r_old(1:NDIM) + sph(p)%v_half(1:NDIM)*dt
     vp(1:VDIM) = sph(p)%v_old(1:VDIM) + sph(p)%a(1:VDIM)*dt

     ! Check if particle has left periodic sphere/box
#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(rp(1:NDIM),vp(1:VDIM))
#endif

  end if
! ----------------------------------------------------------------------------


! Record positions, velocities and other quantities in arrays for all cases 
#if !defined(STATIC_PARTICLES)
  sph(p)%r(1:NDIM) = rp(1:NDIM)
  sph(p)%v(1:VDIM) = vp(1:VDIM)
#endif
#if defined(HYDRO) && defined(VISC_TD)
  if (p <= phydroend) sph(p)%talpha = sph(p)%talpha_old + sph(p)%dalpha_dt*dt
  if (p <= phydroend .and. sph(p)%talpha < alpha_min) sph(p)%talpha = alpha_min
#endif

! If end of timestep, record as 'old' values
  if (dn == nfull) then
     sph(p)%accdo        = .true.
     sph(p)%nlast        = n
     sph(p)%laststep     = timestep*real(nfull,DP)
     sph(p)%r_old(1:NDIM) = rp(1:NDIM)
     sph(p)%v_old(1:VDIM) = vp(1:VDIM)
#if defined(HYDRO) && defined(VISC_TD)
     if (p <= phydroend) sph(p)%talpha_old = sph(p)%talpha
#endif
  end if


! Integrate energy equation if selected
! ----------------------------------------------------------------------------
#if defined(EXPLICIT_ENERGY_EQN) && defined(ENERGY_EQN)
  if (eos == "energy_eqn") then 
     sph(p)%u = sph(p)%u_old + sph(p)%dudt*dt
  end if
#endif
#if defined(ENERGY_EQN)
  if (eos == "energy_eqn" .and. dn == nfull) sph(p)%u_old = sph(p)%u
#endif


! Integrate entropy equation if selected
! ----------------------------------------------------------------------------
#if defined(HYDRO) && defined(ENTROPY_EQN) && defined(ENTROPIC_FUNCTION)
  if (eos == "entropy_eqn") then
     sph(p)%Aent = sph(p)%Aold + sph(p)%dAdt*dt
  end if
#endif
#if defined(HYDRO) && defined(ENTROPIC_FUNCTION)
  if (dn == nfull) sph(p)%Aold = sph(p)%Aent
#endif

 
  return
END SUBROUTINE advance_runge_kutta
