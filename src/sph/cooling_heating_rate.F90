! COOLING_HEATING_RATE.F90
! D. A. Hubber - 29/08/2011
! ============================================================================

#include "macros.h"

! ============================================================================


! ============================================================================
! COOLING_RATE
! Calculate cooling/heating rate for particle p
! ============================================================================
FUNCTION cooling_rate(p)
  use particle_module
  use hydro_module
  use type_module
  use scaling_module
  use constant_module
  use time_module
  implicit none

  integer, intent(in) :: p                        ! Particle id
  real(kind=PR) :: cooling_rate                   ! Cooling rate (dudt)

  character(len=256) :: eos                       ! Eqn. of state for p
  real(kind=PR), parameter :: rconst = 256.0_PR   ! Cooling rate constant
  real(kind=PR), parameter :: u_eq = 1.5_PR       ! Equilibrium temp
  real(kind=DP) :: auxscale                       ! Aux. scaling variable

  eos = typeinfo(sph(p)%ptype)%eos


! Simple linear cooling law
! ----------------------------------------------------------------------------
  if (eos == "energy_eqn" .and. cooling_law == "linear1") then

     cooling_rate = rconst*gammaone*(u_eq - sph(p)%u)

! Simple temperature 'switch' for cooling
! ----------------------------------------------------------------------------
  else if (eos == "energy_eqn" .and. cooling_law == "SD93") then

     if (sph(p)%temp < 1.0E4_PR) then
        cooling_rate = 0.0_PR
     else
        auxscale = (Escale*Ecgs*(rscale*rcgs)**3)/(tscale*tcgs)
        cooling_rate = -2.0E-23_PR*sph(p)%rho*&
             &(mscale*m_SI/m_hydrogen)**2/auxscale
     end if
   
! ----------------------------------------------------------------------------
  else if (eos == "energy_eqn" .and. cooling_law /= "none") then

     write(6,*) "Invalid cooling function options selected"
     stop

  end if
! ----------------------------------------------------------------------------


END FUNCTION cooling_rate



! ============================================================================
! COOLING_EXPONENTIAL_INTEGRATION
! Calculate cooling/heating integration for particle p
! ============================================================================
SUBROUTINE cooling_exponential_integration(p)
  use particle_module
  use hydro_module
  use type_module
  use scaling_module
  use constant_module
  use time_module
  implicit none

  integer, intent(in) :: p                       ! Particle id
  real(kind=PR) :: cooling_rate                  ! Cooling rate (dudt)

  character(len=256) :: eos                      ! Eqn. of state for p
  real(kind=PR), parameter :: rconst = 256.0_PR  ! Cooling rate constant
  real(kind=PR), parameter :: u_eq = 1.5_PR      ! Equilibrium temp
  real(kind=DP) :: auxscale                      ! Aux. scaling variable
  integer(kind=ILP) :: dn       ! Integer timestep since beginning of timestep
  real(kind=PR) :: dt           ! Physical time since beginning of timestep

  eos = typeinfo(sph(p)%ptype)%eos

  dn = n - sph(p)%nlast
  dt = real(timestep,PR)*real(dn,PR)

! Simple linear cooling law
! ----------------------------------------------------------------------------
  if (eos == "energy_eqn" .and. cooling_law == "linear1") then

     sph(p)%u = sph(p)%u_old*exp(-rconst*gammaone*dt) + &
          &(u_eq + sph(p)%dudt/rconst/gammaone)*&
          &(1.0_PR - exp(-rconst*gammaone*dt))

! Error message if invalid cooling option is selected
! ----------------------------------------------------------------------------
  else if (eos == "energy_eqn" .and. cooling_law /= "none") then

     write(6,*) "Invalid cooling function options selected"
     stop

  end if
! ----------------------------------------------------------------------------


  return
END SUBROUTINE cooling_exponential_integration
