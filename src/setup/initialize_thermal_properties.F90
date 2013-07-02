! INITIALIZE_THERMAL_PROPERTIES.F90
! D. A. Hubber - 24/06/2009
! Initialize any thermal properties depending on whether the energy equation
! is solved, if Dragon format is used to read in the ICs, and on the 
! equation of state chosen.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_thermal_properties
  use interface_module
  use filename_module
  use particle_module
  use hydro_module
  use type_module
  use Eos_module
  implicit none

#if defined(HYDRO)
  integer :: p                  ! Particle counter

  debug2("Calculate all initial thermal variables [initialize_thermal_properties.F90]")


! Initialise internal energies based on Radiative cooling algorithm
! ----------------------------------------------------------------------------
#if defined(RAD_WS)
!$OMP PARALLEL DO DEFAULT(SHARED)
  do p=1,ptot
!  do p=pgravitystart,ptot
     call find_idens(sph(p)%rho,sph(p)%idens)
     call find_itemp(sph(p)%temp,sph(p)%itemp)
     sph(p)%u = eosenergy(sph(p)%rho,sph(p)%temp,sph(p)%idens,sph(p)%itemp)
     sph(p)%u_old = sph(p)%u
     sph(p)%dudt = 0.0_PR
#if defined(DIFFUSION)
     sph(p)%du_dt_diff = 0.0_PR
     sph(p)%k_cond = 0.0_PR
#endif
     sph(p)%press = Pconst2*sph(p)%temp*sph(p)%rho / &
          &eosmu(sph(p)%rho,sph(p)%temp,sph(p)%idens,sph(p)%itemp)
     sph(p)%sound = sqrt(sph(p)%press/sph(p)%rho)
     call find_temp_from_energy(sph(p)%idens,sph(p)%u,sph(p)%itemp,sph(p)%temp)
  end do
!$OMP END PARALLEL DO


! Energy equation
! ----------------------------------------------------------------------------
#elif defined(ENERGY_EQN) && !defined(U_IMPLICIT_SOLVER)
  do p=1,ptot
     sph(p)%u_old = sph(p)%u
     sph(p)%dudt = 0.0_PR
#if defined(LEAPFROG_KDK)
     sph(p)%dudt_old = 0.0_PR
#endif
#if defined(ENTROPIC_FUNCTION)
     sph(p)%Aent = (gamma - 1.0_PR)*sph(p)%u/sph(p)%rho**(gamma - 1.0_PR)
#endif
  end do


! Entropic equation
! ----------------------------------------------------------------------------
#elif defined(ENTROPY_EQN)
  do p=1,ptot
     sph(p)%Aent = (gamma - 1.0_PR)*sph(p)%u/sph(p)%rho**(gamma - 1.0_PR)
     sph(p)%Aold = sph(p)%Aent
     sph(p)%dAdt = 0.0_PR
  end do


! Otherwise, if using some other equation of state (e.g. barotropic), 
! calculate the temperatures now we have the densities.
! ----------------------------------------------------------------------------
#else
  sph(1:ptot)%temp = 0.0_PR

#endif
! ----------------------------------------------------------------------------

! Update all other thermal properties
  call update_thermal_properties
#endif


  return
END SUBROUTINE initialize_thermal_properties
