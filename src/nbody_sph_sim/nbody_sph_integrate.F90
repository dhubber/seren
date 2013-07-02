! SPH_INTEGRATE.F90
! C. P. Batty & D. A. Hubber - 10/1/2007
! Main SPH simulation integration subroutine.  
! Controls calls to all routines inside the main SPH integration loop.  
! 1.  - Builds or updates the tree
! 2.  - Computes new smoothing lengths, neighbour lists and densities 
!       for all particles
! 3.  - Calculates thermal properties of all SPH particles
! 4.  - Calculates all hydro forces on all SPH particles
! 5.  - Calculates all gravitational forces on all SPH particles
! 6.  - Calculates all gravitational forces on all sink particles
! 7.  - Writes output files 
! 8.  - Calculates new timesteps for all SPH particles
! 9.  - Checks all neighbour timesteps are comparable
! 10. - Updates thermal properties due to polytropic cooling (if required)
! 11. - Advances the positions and velocities of all particles
! 12. - Reduce the timestep of particles with low-timestep neighbours
! 13. - Searches for new sinks, and accretes particles to existing sinks
! 14. - Removes any escaping, inconsequential particles
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_sph_integrate
  use particle_module
  use time_module
  use Nbody_module
  use sink_module
  use type_module
  implicit none

  integer :: p
  integer :: s

  debug2("Performing next N-body/SPH integration step [nbody_sph_integrate.F90]")

! Calculate new timesteps for particles
  call nbody_sph_timesteps

! Check and flag if any active particles neighbour particles with 
! relatively long timesteps.  
#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
  if (ptot > 0) call check_neighbour_timesteps
#endif

! Update radiative cooling terms
#if defined(HYDRO) && defined(RAD_WS)
  call rad_ws_update
#endif

! Integration scheme to advance particle properties
  call sph_advance
  call nbody_advance
  call copy_stars_to_sinks

! Reduce the timesteps of any particles with short-timestep neighbours 
! to new 'safe' value when the timesteps are properly synchronised.
#if defined(CHECK_NEIGHBOUR_TIMESTEPS) && defined(IMMEDIATE_TIMESTEP_REDUCTION)
  call reduce_timesteps
#endif

! Updating tree properties
  call tree_update(nbuild,nstock)

! Calculate new smoothing lengths, densities and other SPH properties
#if defined(STATIC_PARTICLES)
  sph(1:ptot)%accdo = .false.
#endif
  call sph_update

! Walk HEALPix rays and determine ionisation temperatures, wind accels, etc..
#if defined(HYDRO) && defined(HEALPIX)
  call HP_update
#endif

! Calculating thermal properties for all particles
#if defined(HYDRO)
  call update_thermal_properties
#endif

! Updating tree properties
#if defined(BH_TREE) && defined(SELF_GRAVITY)
  if (nstock == nsteps) call BHgrav_stock
#endif
#if defined(BH_TREE)
#if defined(EULER) || defined(RUNGE_KUTTA2)
  if (nstock == nsteps) nstock = nstock + 1
#else
  if (nstock == nsteps) nstock = nstock + 2
#endif
  if (nbuild == nsteps) nbuild = nbuild + nbuildstep
#endif

! Zero acceleration array of all active particles here (for now)
  do p=1,ptot
     if (sph(p)%accdo) sph(p)%a(1:VDIM) = 0.0_PR
#if defined(GRAVITY)
     if (sph(p)%accdo) sph(p)%gpot = 0.0_PR
#endif
#if defined(DEBUG_FORCES) && defined(GRAVITY)
     if (sph(p)%accdo) sph(p)%a_grav(1:NDIM) = 0.0_PR
#endif
#if defined(DEBUG_FORCES) && defined(HYDRO)
     if (sph(p)%accdo) sph(p)%a_hydro(1:NDIM) = 0.0_PR
#endif
  end do

! Zero relevant star arrays
  do s=1,stot
     if (.not. star(s)%accdo) cycle
     star(s)%a(1:NDIM) = 0.0_DP
     star(s)%adot(1:NDIM) = 0.0_DP
     star(s)%gpot = 0.0_DP
     star(s)%gpe = 0.0_DP
  end do

!#if defined(FORCE_SPLITTING)
!  do s=1,stot
!     if (star(s)%accsph .or. star(s)%accdo) then
!        star(s)%asph(1:NDIM) = 0.0_PR
!        star(s)%adotsph(1:NDIM) = 0.0_PR
!        star(s)%g
!     end if
!     if (star(s)%accdo) then
!        star(s)%astar(1:NDIM) = 0.0_PR
!        star(s)%adotstars(1:NDIM) = 0.0_PR
!     end if
!  end do
!#endif


! Calculate hydro forces on all SPH particles
#if defined(HYDRO)
  call sph_hydro_forces
#endif

! Calculate gravitational forces on all SPH particles
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  call sph_grav_forces
#endif

! Calculate gravitational forces on all sink particles
#if defined(FORCE_SPLITTING) && defined(GRAVITY)
  call nbody_sph_star_split_forces
#elif defined(GRAVITY)
  call nbody_sph_star_forces
#endif

! Add wind contribution to acceleration here (for now)
#if defined(STELLAR_WIND)
  do p=1,ptot
     if (sph(p)%accdo) sph(p)%a(1:NDIM) = &
          &sph(p)%a(1:NDIM) + sph(p)%a_wind(1:NDIM)
  end do
#endif

! Apply correction to integrated quantities if using Leapfrog-KDK scheme
#if defined(LEAPFROG_KDK)
  call leapfrog_kdk_correction_terms
#endif

! Calculate higher-order derivatives and correction forces for previous steps
#if defined(GRAVITY)
  call nbody_correction_terms

! Once iterations are finished, set all star properties for end of step
  call nbody_end_step
#endif

! Output diagnostics and/or snapshot files
  call nbody_sph_output

  return
END SUBROUTINE nbody_sph_integrate
