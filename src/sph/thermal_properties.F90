! THERMAL_PROPERTIES.F90
! D. A. Hubber - 20/03/2011
! ..
! =============================================================================


#include "macros.h"


#if defined(HYDRO)
! =============================================================================
! PRESSURE
! Returns the ideal gas pressure
! =============================================================================
FUNCTION pressure(p)
  use interface_module, only : eosmu,specific_internal_energy
  use particle_module, only : sph
  use type_module
  use hydro_module
  implicit none

  integer, intent(in) :: p                     ! particle id
  real(kind=PR) :: pressure                    ! ..
  real(kind=PR) :: mu_bar_p                    ! ..
 
  if (typeinfo(sph(p)%ptype)%eos == "isothermal") then
     pressure = Pconst*sph(p)%temp*sph(p)%rho
  else if (typeinfo(sph(p)%ptype)%eos == "barotropic") then
     pressure = Pconst*sph(p)%temp*sph(p)%rho
  else if (typeinfo(sph(p)%ptype)%eos == "polytropic") then
     pressure = Kpoly*(sph(p)%rho**gamma)
  else if (typeinfo(sph(p)%ptype)%eos == "energy_eqn") then
     pressure = gammaone*sph(p)%rho*specific_internal_energy(p)
#if defined(RAD_WS)
  else if (typeinfo(sph(p)%ptype)%eos == "rad_ws") then
     mu_bar_p = eosmu(sph(p)%rho,sph(p)%temp,sph(p)%idens,sph(p)%itemp)
     pressure = Pconst2*sph(p)%temp*sph(p)%rho/mu_bar_p
#endif
  else if (typeinfo(sph(p)%ptype)%eos == "stiff") then
     if (sph(p)%rho > 1.0_PR) then
        pressure = (50.0_PR/7.0_PR)*((sph(p)%rho**7) - 1.0_PR)
     else
        pressure = (50.0_PR/7.0_PR)
     end if
  else if (typeinfo(sph(p)%ptype)%eos == "stellar_heat") then
     pressure = Pconst*sph(p)%rho*sph(p)%temp
  end if

END FUNCTION pressure




! =============================================================================
! SOUND_SPEED
! Returns the sound speed
! =============================================================================
FUNCTION sound_speed(p)
  use particle_module, only : sph
  use hydro_module
  use type_module
  implicit none

  integer, intent(in) :: p                     ! particle id
  real(kind=PR) :: sound_speed
 
  if (typeinfo(sph(p)%ptype)%eos == "isothermal") then
     sound_speed = sqrt(sph(p)%press/sph(p)%rho)
  else if (typeinfo(sph(p)%ptype)%eos == "barotropic") then
     sound_speed = sqrt(sph(p)%press/sph(p)%rho)
  else if (typeinfo(sph(p)%ptype)%eos == "polytropic") then
     sound_speed = sqrt(gamma*sph(p)%press/sph(p)%rho)
  else if (typeinfo(sph(p)%ptype)%eos == "energy_eqn") then
     sound_speed = sqrt(gamma*sph(p)%press/sph(p)%rho)
  else if (typeinfo(sph(p)%ptype)%eos == "rad_ws") then
     sound_speed = sqrt(gamma*sph(p)%press/sph(p)%rho)
  else if (typeinfo(sph(p)%ptype)%eos == "stiff") then
     sound_speed = sqrt(gamma*sph(p)%press/sph(p)%rho)
  else if (typeinfo(sph(p)%ptype)%eos == "stellar_heat") then
     sound_speed = sqrt(gamma*sph(p)%press/sph(p)%rho)
  end if

END FUNCTION sound_speed



! =============================================================================
! SPECIFIC_INTERNAL_ENERGY
! Returns the specific internal energy
! =============================================================================
FUNCTION specific_internal_energy(p)
  use particle_module, only : sph
  use hydro_module
  use type_module
  implicit none

  integer, intent(in) :: p                     ! particle id
  real(kind=PR) :: specific_internal_energy    ! ..
 
  if (typeinfo(sph(p)%ptype)%eos == "isothermal") then
     specific_internal_energy = Pconst*isotemp/gammaone
  else if (typeinfo(sph(p)%ptype)%eos == "barotropic") then
     specific_internal_energy = Pconst*sph(p)%temp/gammaone
  else if (typeinfo(sph(p)%ptype)%eos == "polytropic") then
     specific_internal_energy = Kpoly*(sph(p)%rho**gammaone)/gammaone
  else if (typeinfo(sph(p)%ptype)%eos == "energy_eqn") then
     specific_internal_energy = sph(p)%u
  else if (typeinfo(sph(p)%ptype)%eos == "rad_ws") then
     specific_internal_energy = sph(p)%u
  else if (typeinfo(sph(p)%ptype)%eos == "stiff") then
     specific_internal_energy = Pconst*sph(p)%temp/gammaone
  else if (typeinfo(sph(p)%ptype)%eos == "stellar_heat") then
     specific_internal_energy = Pconst*sph(p)%temp/gammaone
  end if

END FUNCTION specific_internal_energy




! =============================================================================
! TEMPERATURE
! Returns the gas temperature
! =============================================================================
FUNCTION temperature(p)
  use interface_module, only : eosmu
  use particle_module, only : sph
  use hydro_module
  use type_module
  implicit none

  integer, intent(in) :: p                     ! particle id
  real(kind=PR) :: temperature                 ! temperature

  real(kind=PR) :: mu_bar_p                    ! mean gas particle mass for p
 
  if (typeinfo(sph(p)%ptype)%eos == "isothermal") then
     temperature = isotemp
  else if (typeinfo(sph(p)%ptype)%eos == "barotropic") then
     temperature = isotemp*(1.0_PR + (sph(p)%rho/rhobary)**(gammaone))
  else if (typeinfo(sph(p)%ptype)%eos == "polytropic") then
     temperature = Kpoly*(sph(p)%rho**gammaone)/Pconst
  else if (typeinfo(sph(p)%ptype)%eos == "energy_eqn") then
     temperature = sph(p)%press/Pconst/sph(p)%rho
#if defined(RAD_WS)
  else if (typeinfo(sph(p)%ptype)%eos == "rad_ws") then
     mu_bar_p = eosmu(sph(p)%rho,sph(p)%temp,sph(p)%idens,sph(p)%itemp)
     temperature = Pconst2*sph(p)%temp*sph(p)%rho/mu_bar_p
#endif
  else if (typeinfo(sph(p)%ptype)%eos == "stiff") then
     if (sph(p)%rho > 1.0_PR) then
        temperature = (50.0_PR/7.0_PR)*((sph(p)%rho**7) - 1.0_PR)
     else
        temperature = (50.0_PR/7.0_PR)
     end if
  else if (typeinfo(sph(p)%ptype)%eos == "stellar_heat") then
     temperature = Pconst*sph(p)%rho*sph(p)%temp
  end if

END FUNCTION temperature
#endif
