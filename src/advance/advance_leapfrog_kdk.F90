! ADVANCE_LEAPFROG_KDK.F90
! C. P. Batty & D. A. Hubber - 19/3/2007
! Second order leapfrog kick-drift-kick integration scheme.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE advance_leapfrog_kdk(p)
  use interface_module, only : check_boundary_conditions
  use particle_module, only : sph
  use hydro_module
  use time_module
  use type_module
  implicit none

  integer, intent(in) :: p      ! Particle id

  character(len=256) :: eos     ! Equation of state for particle p
  integer(kind=ILP) :: dn       ! Integer timestep since beginning of timestep
  integer(kind=ILP) :: nfull    ! Full integer timestep
  integer(kind=ILP) :: nhalf    ! Half integer timestep
  real(kind=PR) :: dt           ! Physical time since beginning of timestep

  debug3("Advancing particle ",p)

! Work out integer and real time intervals from beginning of step 
  sph(p)%accdo = .false.
  dn      = n - sph(p)%nlast
  nfull   = 2_ILP**(level_step - sph(p)%nlevel)
  nhalf   = nfull / 2_ILP
  dt      = real(timestep,PR)*real(dn,PR)
  eos     = typeinfo(sph(p)%ptype)%eos


! Advance particle positions and velocities
! ----------------------------------------------------------------------------
#if !defined(STATIC_PARTICLES)
  sph(p)%r(1:NDIM) = sph(p)%r_old(1:NDIM) + &
       &sph(p)%v_old(1:NDIM)*dt + 0.5_PR*sph(p)%a(1:NDIM)*dt*dt
  sph(p)%v(1:VDIM) = sph(p)%v_old(1:VDIM) + sph(p)%a(1:VDIM)*dt
#endif
  
! Check if particle has left periodic sphere/box
#if defined(BOUNDARY_CONDITIONS)
  call check_boundary_conditions(sph(p)%r(1:NDIM),sph(p)%v(1:VDIM))
#endif
  
#if defined(HYDRO) && defined(VISC_TD)
  if (p <= phydroend) sph(p)%talpha = sph(p)%talpha_old + sph(p)%dalpha_dt*dt
  if (p <= phydroend .and. sph(p)%talpha < alpha_min) sph(p)%talpha = alpha_min
#endif

! If end of timestep, record as 'old' values
  if (dn == nfull) then
     sph(p)%accdo         = .true.
     sph(p)%nlast         = n
     sph(p)%laststep      = timestep*real(nfull,DP)
     sph(p)%r_old(1:NDIM) = sph(p)%r(1:NDIM)
     sph(p)%v_old(1:VDIM) = sph(p)%v(1:VDIM)
#if defined(HYDRO) && defined(VISC_TD)
     if (p <= phydroend) sph(p)%talpha_old = sph(p)%talpha
#endif
  end if


! Integrate energy equation if selected
! ----------------------------------------------------------------------------
#if defined(ENERGY_EQN) && defined(EXPLICIT_ENERGY_EQN)
  if (eos == "energy_eqn") then 
     sph(p)%u = sph(p)%u_old + sph(p)%dudt*dt
     if (sph(p)%u < SMALL_NUMBER) &
          &sph(p)%u = sph(p)%u_old*exp(-sph(p)%u_old/sph(p)%dudt)
     if (dn == nfull) sph(p)%dudt_old = sph(p)%dudt
  end if
#endif
! Energy equation is integrated IMPLICITLY in the chemistry module
! ----------------------------------------------------------------------------
#if defined(ENERGY_EQN) && defined(CHEMCOOL)
  if (eos == "chemcool_eos" ) then
     ! The pdv and shocks are stored in sph(p)%dudt and calculated above. 
     if (dn == nfull) then
        call do_chemcool_step(sph(p)%rho, sph(p)%u, sph(p)%abundances, sph(p)%dudt, dt)
        sph(p)%dudt_old = sph(p)%dudt
     end if
  end if
#endif
#if defined(ENERGY_EQN)
  if (dn == nfull) sph(p)%u_old = sph(p)%u
#endif

! Integrate entropy equation if selected
! ----------------------------------------------------------------------------
#if defined(ENTROPY_EQN) && defined(ENTROPIC_FUNCTION)
  if (eos == "entropy_eqn" .and. dn <= nhalf) then
     sph(p)%Aent = sph(p)%Aold + sph(p)%dAdt*dt
     if (dn == nfull/2_ILP) sph(p)%Ahalf = sph(p)%Aent
  else if (eos == "entropy_eqn") then
     sph(p)%Aent = sph(p)%Ahalf + sph(p)%dAdt*dt_half
  end if
#endif
#if defined(ENTROPIC_FUNCTION)
  if (dn == nfull) sph(p)%Aold = sph(p)%Aent
#endif

  return
END SUBROUTINE advance_leapfrog_kdk
