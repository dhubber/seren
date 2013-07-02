! EFFECTIVE_GAMMA.F90
! D. A. Hubber - 9/12/2009
! Calculate the effective adiabatic index (ratio of specific heats) 
! for two intereacting particles.  Calculates exact values for isothermal, 
! polytropic and barotropic equations of state.  Currently calculates an 
! approximate value for gamma_eff for the RAD_WS option by asusming the gas 
! is locally polytropic with exponent gamma.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE effective_gamma(p1,p2,gamma_eff)
  use particle_module
  use hydro_module
  implicit none

  integer, intent(in) :: p1                 ! Particle 1 id
  integer, intent(in) :: p2                 ! Particle 2 id
  real(kind=PR), intent(out) :: gamma_eff   ! Effective gamma

#if defined(RAD_WS)
  real(kind=PR) :: dlogP                    ! delta log P
  real(kind=PR) :: dlogrho                  ! delta log rho
#endif

  
! Isothermal equation of state
! ----------------------------------------------------------------------------
#if defined(ISOTHERMAL) 
  gamma_eff = 1.0_PR
#endif
! ---------------------------------------------------------------------------- 


! Barotropic equation of state
! ----------------------------------------------------------------------------
#if defined(BAROTROPIC)
  gamma_eff = 0.5_PR*((1.0_PR + gamma*((sph(p1)%rho/rhobary)**(gammaone))) &
       & / (1.0_PR + ((sph(p1)%rho/rhobary)**(gammaone))) + &
       & (1.0_PR + gamma*((sph(p2)%rho/rhobary)**(gammaone))) / &
       & (1.0_PR + ((sph(p2)%rho/rhobary)**(gammaone))))
#endif
! ---------------------------------------------------------------------------- 

 
! Polytropic equation of state
! ----------------------------------------------------------------------------
#if defined(POLYTROPIC)
  gamma_eff = gamma
#endif
! ----------------------------------------------------------------------------

  
! Stiff equation of state (Monaghan shear flow)
! ----------------------------------------------------------------------------
#if defined(STIFF)
  stop 'Riemann solver not working for STIFF yet'
#endif
! ----------------------------------------------------------------------------
  

! Ideal-gas equation of state when using the internal energy equation
! ----------------------------------------------------------------------------
#if defined(ENERGY_EQN) && !defined(RAD_WS)
  gamma_eff = gamma
#endif
! ----------------------------------------------------------------------------
  
  
! Polytropic-cooling approximation
! ----------------------------------------------------------------------------
#if defined(RAD_WS)
  dlogP = (log10(sph(p1)%press) - log10(sph(p2)%press))
  dlogrho = log10(sph(p1)%rho) - log10(sph(p2)%rho)
  if (abs(dlogrho) > 0.00001_PR) then 
     gamma_eff = dlogP / dlogrho
  else
     gamma_eff = gamma
  end if
#endif
! ----------------------------------------------------------------------------
  
  
! Heating from sinks (stars)
! ----------------------------------------------------------------------------
#if defined(STELLAR_HEAT)
  stop 'Riemann solver not working for STELLAR_HEAT yet'
#endif
! ----------------------------------------------------------------------------
  
  
! Other thermal quantities due to ionizing radiation
! ----------------------------------------------------------------------------
#if defined(IONIZING_UV_RADIATION)
  stop 'Riemann solver not working for UV radiation yet'
#endif
! ----------------------------------------------------------------------------
  

! Ensure (for now) that we do not model a very soft gas
  if (gamma_eff < 1.0_PR) gamma_eff = 1.0_PR


  return
END SUBROUTINE effective_gamma
