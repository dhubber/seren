program one_zone_chemcool
!
! This is the test program from the simple H2 chemistry
!

use chemistry_constants

implicit none

! Phyiscal parameters of gas
double precision :: abundances(2)
double precision :: density, temperature

! time loop variables 
double precision :: tmax, dt, dttry, tcurrent, tleft, time_passed
integer :: itime, itime_max
parameter ( itime_max = 100000000 ) ! 1e8 attemps.

! sub-cylcing variables
integer :: istep_try, istep_try_max
parameter ( istep_try_max = 1000 )
integer :: iteration_sucess

!
! Initialise the chemisty
!

print *, "Initailising the chemistry... "
tmax = 100D6 * YEAR  ! 5 Myr
dt = 1000 * YEAR
density = 1.0D0*((1.0 + 4.0 * ABHE) * PROTONMASS)
temperature = 10000.
abundances(1) = 0.1  ! Zero H2 fraction intially
abundances(2) = 1.0D-3  ! Small ionisation fraction initially

print *, "Done!"

!
! The main loop!
!

print *, "Entering the time loop"
tcurrent = 0
time_passed = 0
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
      call solve_chem_timestep(dttry, iteration_sucess, density, temperature, abundances)
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
   ! Output information 
   !
   time_passed = time_passed + dttry
   if ( time_passed.gt.(1e5*YEAR) ) then
      write(*, 1000), tcurrent/YEAR/1e6, abundances(1), abundances(2), temperature, iteration_sucess
      time_passed = 0
   end if
   !
   ! keep track of total time advanced.
   !
   tcurrent = tcurrent + dttry
   if ( tcurrent.ge.tmax ) exit
end do
1000 format(4(1X, F13.5), I3)

666 print *, "Reached the maximum time... Done!"

end program one_zone_chemcool 
