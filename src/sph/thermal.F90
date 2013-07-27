! THERMAL.F90
! C. P. Batty & D. A. Hubber - 19/1/2007
! Calculates temperature, pressure, sound speed of all particles depending
! on the chosen equation of state.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE thermal(p)
  use interface_module, only : distance3,eosenergy,eosmu,&
       &find_idens,find_itemp,find_temp_from_energy
  use particle_module
  use hydro_module
  use scaling_module
  use time_module
  use type_module
#if defined(STELLAR_HEAT)
  use sink_module
#endif
#if defined(RAD_WS)
  use Eos_module
#endif
#if defined(IONIZING_UV_RADIATION)
  use HP_module
#endif
#if defined(CHEMCOOL)
  use chemistry_constants
#endif
  implicit none

  integer, intent(in) :: p     ! particle id

#if defined(HYDRO)
  character(len=256) :: eos    ! Equation of state for particle p
#if defined(RAD_WS)
  real(kind=PR) :: mu_bar_p    ! Mean gas particle mass for particle p
#endif
#if defined(IONIZING_UV_RADIATION) && !defined(RAD_WS)
  real(kind=PR) :: invmu       ! 1 / mu
#endif
#if defined(STELLAR_HEAT)
  integer       :: j           ! Aux. loop counter
  real(kind=PR) :: dr(1:NDIM)  ! Relative displacement
  real(kind=PR) :: drsqd       ! Distance squared 
#endif
#if defined(CHEMCOOL)
  double precision :: numdens, numdens_tot
#endif

  eos = typeinfo(sph(p)%ptype)%eos


! Isothermal equation of state
! ----------------------------------------------------------------------------
  if (eos == "isothermal") then
     sph(p)%temp  = isotemp
     sph(p)%press = Pconst*sph(p)%temp*sph(p)%rho
     sph(p)%sound = sqrt(sph(p)%press/sph(p)%rho)
     sph(p)%u     = sph(p)%press/sph(p)%rho 


! Barotropic equation of state
! ----------------------------------------------------------------------------
  else if (eos == "barotropic") then
     sph(p)%temp  = isotemp*(1.0_PR + (sph(p)%rho/rhobary)**(gamma - 1.0_PR))
     sph(p)%press = Pconst*sph(p)%temp*sph(p)%rho
     sph(p)%sound = sqrt(sph(p)%press/sph(p)%rho)
     sph(p)%u     = sph(p)%press/sph(p)%rho 
  
  
! Polytropic equation of state
! ----------------------------------------------------------------------------
  else if (eos == "polytropic") then
     sph(p)%press = Kpoly*(sph(p)%rho**(gamma))
     sph(p)%temp  = Kpoly*(sph(p)%rho**(gamma - 1.0_PR))/Pconst
     sph(p)%sound = sqrt(gamma*sph(p)%press/sph(p)%rho)
     sph(p)%u     = sph(p)%press/sph(p)%rho 


! Ideal-gas equation of state when using the internal energy equation
! ----------------------------------------------------------------------------
#if defined(ENERGY_EQN)
  else if (eos == "energy_eqn") then
     sph(p)%temp = (gamma - 1.0_PR)*sph(p)%u / Pconst
     sph(p)%press = (gamma - 1.0_PR)*sph(p)%rho*sph(p)%u
     sph(p)%sound = sqrt(gamma*sph(p)%press/sph(p)%rho)
#if defined(ENTROPIC_FUNCTION) && defined(ENERGY_EQN)
     sph(p)%Aent = (gamma - 1.0_PR)*sph(p)%u/sph(p)%rho**(gamma - 1.0_PR)
#endif
#endif


! Ideal-gas equation of state when using the internal energy equation
! ----------------------------------------------------------------------------
#if defined(ENTROPY_EQN) && defined(ENTROPIC_FUNCTION)
  else if (eos == "entropy_eqn") then
     sph(p)%u = (sph(p)%Aent*sph(p)%rho**(gamma - 1.0_PR))/(gamma - 1.0_PR)
     sph(p)%press = (gamma - 1.0_PR)*sph(p)%rho*sph(p)%u
     sph(p)%temp  = sph(p)%press / (Pconst*sph(p)%rho)
     sph(p)%sound = sqrt(gamma*sph(p)%press/sph(p)%rho)
#endif
  
  
! Polytropic-cooling approximation
! ----------------------------------------------------------------------------
#if defined(RAD_WS)
  else if (eos == "rad_ws") then
#if defined(IONIZING_UV_RADIATION)
     if (sph(p)%temp < sph(p)%tempmin) then
        sph(p)%temp = sph(p)%tempmin
        call find_idens(sph(p)%rho,sph(p)%idens)
        call find_itemp(sph(p)%tempmin,sph(p)%itemp)
        sph(p)%u = eosenergy(sph(p)%rho,&
             &sph(p)%tempmin,sph(p)%idens,sph(p)%itemp)
     else 
        call find_idens(sph(p)%rho,sph(p)%idens)
        call find_temp_from_energy(sph(p)%idens,sph(p)%u,&
             &sph(p)%itemp,sph(p)%temp)
     end if
#else
     call find_idens(sph(p)%rho,sph(p)%idens)
     call find_temp_from_energy(sph(p)%idens,sph(p)%u,sph(p)%itemp,sph(p)%temp)
#endif
     sph(p)%temp   = sph(p)%temp
     mu_bar_p      = eosmu(sph(p)%rho,sph(p)%temp,sph(p)%idens,sph(p)%itemp)
     sph(p)%press  = Pconst2*sph(p)%temp*sph(p)%rho/mu_bar_p
     sph(p)%sound  = sqrt(sph(p)%press/sph(p)%rho)
#endif
  

! Ideal-gas equation of state when using the internal energy equation
! ----------------------------------------------------------------------------
#if defined(ENERGY_EQN)
  else if (eos == "energy_eqn") then
     sph(p)%temp = (gamma - 1.0_PR)*sph(p)%u / Pconst
     sph(p)%press = (gamma - 1.0_PR)*sph(p)%rho*sph(p)%u
     sph(p)%sound = sqrt(gamma*sph(p)%press/sph(p)%rho)
#if defined(ENTROPIC_FUNCTION) && defined(ENERGY_EQN)
     sph(p)%Aent = (gamma - 1.0_PR)*sph(p)%u/sph(p)%rho**(gamma - 1.0_PR)
#endif
#endif

! Ideal-gas equation of state when using the internal energy equation
! ----------------------------------------------------------------------------
#if defined(ENERGY_EQN) && defined(CHEMCOOL)
  else if (eos == "chemcool_eos") then
     gamma = GAMMA_GAS
     sph(p)%press = (gamma - 1.0_PR)*sph(p)%rho*sph(p)%u
     numdens = sph(p)%rho * rhoscale * rhocgs /  ((1.0 + 4.0 * ABHE) * PROTONMASS)
     numdens_tot = (1.0D0 + ABHE - sph(p)%abundances(1) + sph(p)%abundances(2)) * numdens
     sph(p)%temp = sph(p)%press/ ( numdens_tot * KBOLTZ )
     sph(p)%sound = sqrt(gamma*sph(p)%press/sph(p)%rho)
#endif


! Ideal-gas equation of state when using the internal energy equation
! ----------------------------------------------------------------------------
#if defined(ENTROPY_EQN) && defined(ENTROPIC_FUNCTION)
  else if (eos == "entropy_eqn") then
     sph(p)%u = (sph(p)%Aent*sph(p)%rho**(gamma - 1.0_PR))/(gamma - 1.0_PR)
     sph(p)%press = (gamma - 1.0_PR)*sph(p)%rho*sph(p)%u
     sph(p)%temp  = sph(p)%press / (Pconst*sph(p)%rho)
     sph(p)%sound = sqrt(gamma*sph(p)%press/sph(p)%rho)
#endif
  

! Heating from sinks (stars)
! ----------------------------------------------------------------------------
#if defined(STELLAR_HEAT)
  else if (eos == "stellar_heat") then
     sph(p)%temp = 0.0_PR
     do j=1,STARS
        call distance3(sink(j)%r(1:NDIM),sph(p)%r(1:NDIM),dr(1:NDIM),drsqd)
        sph(p)%temp = sph(p)%temp + (isotemp**4) / drsqd
     end do
     sph(p)%temp  = sph(p)%temp + 10000.0_PR ! (10K)**4
     sph(p)%temp  = sph(p)%temp**0.25_PR
     sph(p)%press = Pconst*sph(p)%rho*sph(p)%temp
     sph(p)%sound = sqrt(sph(p)%press/sph(p)%rho)
#endif


! Quit if invalid EOS is chosen (or not selected in Makefile)
! ----------------------------------------------------------------------------
  else 
     write(6,*) "Invalid choice of EOS for particle ",p,eos
     stop

  end if
! ----------------------------------------------------------------------------
  
  
! Other thermal quantities due to ionizing radiation
! ----------------------------------------------------------------------------
#if defined(IONIZING_UV_RADIATION) && !defined(RAD_WS)
  if (sph(p)%temp < sph(p)%tempmin) then
     sph(p)%temp = max(sph(p)%temp,sph(p)%tempmin)
     invmu = ((sph(p)%temp - Tneut)/mu_ion + &
          &(Tion - sph(p)%temp)/mu_bar)/(Tion - Tneut)
     if (sph(p)%temp .ge. Tion) invmu = 1.0_PR/mu_ion
     if (sph(p)%temp .le. Tneut) invmu = 1.0_PR/mu_bar
     sph(p)%press = invmu*Pconst2*sph(p)%temp*sph(p)%rho
     sph(p)%sound = sqrt(sph(p)%press/sph(p)%rho)
#if defined(ENERGY_EQN)
     if (eos == "energy_eqn") then
        sph(p)%u = invmu*Pconst2*sph(p)%temp / (gamma - 1.0_PR)
        sph(p)%u_old = sph(p)%u
        sph(p)%dudt = 0.0_PR
     end if
#endif
  end if
#endif
  
  
! Calculate Balsara switch factor
! ----------------------------------------------------------------------------
#if defined(VISC_BALSARA)
  sph(p)%balsara = min(1.0_PR,(abs(sph(p)%div_v) / &
       &(sph(p)%balsara + abs(sph(p)%div_v) + &
       &BAL_DENOM*sph(p)%sound/sph(p)%h)))
#endif
#endif
  
  return
END SUBROUTINE thermal
