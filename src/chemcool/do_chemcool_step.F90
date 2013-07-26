subroutine do_chemcool_step(rho, u, abundances, dudt, timestep)

!
! PCC 18-07-2013
!
! This subroutine is the interface between the SPH code and 
! the chemistry solver. It is called in place of doing an 
! explicit update of the energy, and updates the particle
! energy and aboundances, as integrated over the timestep, dt.
! The timestep may be too large for the chemistry and cooling/heating
! so the code is allowed to sub-cycle
!

use scaling_module
use chemistry_constants

implicit none

! Phyiscal parameters from SPH code (code units)
double precision :: abundances(2)
double precision :: rho, u, dudt
double precision :: timestep
double precision :: initial_u, final_u

! The chemistry needs cgs version of the above parameters
double precision :: density, energy, deng_dt

! time loop variables 
double precision :: tmax, dt, dttry, tcurrent, tleft, time_passed
integer :: itime, itime_max
parameter ( itime_max = 100000000 ) ! 1e8 attemps.

! sub-cylcing variables
integer :: istep_try, istep_try_max
parameter ( istep_try_max = 1000 )
integer :: iteration_sucess

!
! Initialise the variables for going into the chemistry solver.
! Everything needs to be converted to cgs, and specific energy 
! (and the rates) needs to be converted to energy density
!

initial_u = u

dt = timestep * tcgs * tscale
tmax = dt

density = (rho * rhocgs * rhoscale) 
energy = (u * ucgs * uscale) * density
deng_dt = (dudt * dudtcgs * dudtscale) * density

!
! The main loop!
!

tcurrent = 0
do itime = 1, itime_max 
   !
   ! reduce timestep if we're going to overshoot the end point
   !
   tleft = tmax - tcurrent
   if ( dt.gt.tleft ) dt = tleft
   !
   ! Do a chem_cool step
   !
   dttry = dt
   do istep_try = 1, istep_try_max
      !
      ! start with the standard dt as a step
      !
      call solve_chem_timestep(dttry, iteration_sucess, density, energy, deng_dt, abundances)
      !
      ! See if it was successful, if not, decrease timestep and try again
      !
      if ( iteration_sucess.gt.0 ) then
         exit
      else
         dttry = dttry/2.0D0
      end if
   end do
   !
   ! keep track of total time advanced.
   !
   tcurrent = tcurrent + dttry
   if ( tcurrent.ge.tmax ) exit
end do

!
! Do some sanity checks here before going back to main code
!



!
! Convert energy (and the rate) back to specific before going back to the main code.
!

u = energy / density / (ucgs * uscale)
final_u = u
! the rate is updated here too, as it will be needed when trying to predict the
! new timestep
dudt =  (final_u - initial_u) / timestep

end subroutine do_chemcool_step
