! SPH_INTEGRATE.F90
! C. P. Batty & D. A. Hubber - 10/1/2007
! Main SPH simulation integration subroutine.  
! Controls calls to all routines inside the main SPH integration loop.  
! 1.  - Builds or updates the tree
! 2.  - If required, creates ghost particles and builds ghost tree
! 3.  - Computes new smoothing lengths, neighbour lists and densities 
!       for all particles
! 4.  - Calculates thermal properties of all SPH particles
! 5.  - Calculates all hydro forces on all SPH particles
! 6.  - Calculates all gravitational forces on all SPH particles
! 7.  - Calculates all gravitational forces on all sink particles
! 8.  - Writes output files 
! 9.  - Calculates new timesteps for all SPH particles
! 10. - Checks all neighbour timesteps are comparable
! 11. - Updates thermal properties due to radiative cooling (if required)
! 12. - Advances the positions and velocities of all particles
! 13. - Reduce the timestep of particles with low-timestep neighbours
! 14. - Searches for new sinks, and accretes particles to existing sinks
! 15. - Removes any escaping, inconsequential particles
! 16. - Adds new particles (e.g. jets)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_integrate
  use particle_module
  use sink_module
  use time_module
  use tree_module
  use type_module
  use neighbour_module
#if defined(USE_MPI)
  use mpi
  use mpi_communication_module
#endif
  use periodic_module
  implicit none

  integer :: p                      ! Particle counter
#if defined(SINKS)
  integer :: s                      ! Sink counter
#endif
  logical :: done       ! Are we done?

  debug2("Performing next SPH integration step [sph_integrate.F90]")

#if defined(USE_MPI)
  start_calctime = calctime
  calcstart = MPI_WTIME()
#endif

! Walk HEALPix rays and determine ionisation temperatures, wind accels, etc..
#if defined(HYDRO) && defined(HEALPIX)
  call HP_update
#endif

! Updating tree properties
  call tree_update(nbuild,nstock)

! If using ghost particles with periodic boundaries, create ghost particles 
! and set-up the tree for further searching.
! ----------------------------------------------------------------------------
#if defined(GHOST_PARTICLES) && defined(PERIODIC) && !defined(USE_MPI)
  debug_timing("GHOSTS")
  do
     pghost = 0
     call search_ghost_particles
#if defined(BH_TREE)
     if (pghost > 0) then
        call BHghost_build
        call BHhydro_stock(cmax_ghost,ctot_ghost,&
             &ltot_ghost,first_cell_ghost,last_cell_ghost,BHghost)
     end if
#endif

     call sph_update
     
     ! If the no. of ghost particles is acceptable, exit loop
     done = (pghost >= 0 .and. pghost <= pghostmax)
     if (done) exit

     ! Check if ghost boundary zones are valid
     !call check_ghost_boundaries

  end do
! ----------------------------------------------------------------------------
#else
#if defined(USE_MPI)
  call mpi_sph_update
#else
  call sph_update
#endif

#endif
! ----------------------------------------------------------------------------

! Calculate all sink-stellar properties here for now
#if defined(STELLAR_FEEDBACK)
  call calculate_sink_properties
#endif

! Calculating thermal properties for all particles
#if defined(HYDRO)
  call update_thermal_properties
#endif

! Copy all particle data to corresponding ghosts now we have thermal info
#if defined(GHOST_PARTICLES) && defined(PERIODIC)
  if (pperiodic > 0) call copy_data_to_ghosts
#endif

! Send gather and scatter ghosts for calculating hydro forces etc.
#if defined(USE_MPI)
  call mpi_scatter_ghosts
#endif

! Updating tree properties
#if defined(BH_TREE) && defined(SELF_GRAVITY)
  if (nstock <= nsteps) call BHgrav_stock
#endif
#if defined(BH_TREE)
#if defined(EULER) || defined(RUNGE_KUTTA2)
  if (nstock <= nsteps) nstock = nsteps + 1
#else
  if (nstock <= nsteps) nstock = nsteps + 2
#endif
  if (nbuild <= nsteps) nbuild = nsteps + nbuildstep
#endif

#if defined(USE_MPI) && defined(SELF_GRAVITY)
! Slice off top levels of gravity tree
  call BHgrav_slice
! Share pruned gravity trees
  call share_pruned_gravtrees
#endif

! Zero acceleration array of all active particles here (for now)
  do p=1,ptot
     if (sph(p)%accdo) then
        sph(p)%a(1:VDIM) = 0.0_PR
#if defined(GRAVITY)
        sph(p)%gpot = 0.0_PR
#endif
! #if defined(DEBUG_FORCES) && defined(GRAVITY)
!         sph(p)%a_grav(1:NDIM) = 0.0_PR
! #endif
#if defined(DEBUG_FORCES) && defined(HYDRO)
        sph(p)%a_hydro(1:NDIM) = 0.0_PR
#endif
     end if
  end do

! Zero relevant sink gravity variables
!#if defined(SINKS)
!  do s=1,stot
!     if (accdo_sinks) sink(s)%a(1:NDIM) = 0.0_DP
!     if (accdo_sinks) sink(s)%gpot = 0.0_PR
!  end do
!#endif

! Calculate gravitational forces on all SPH particles
! Note that in MPI, 'remote' gravity forces are added to sph()%a
! Gravity must therefore be calculated first if it is necessary to know what
! fraction of the total acceleration is gravitational; this is necessary for
! DEBUG_FORCES (a_grav) and MACs other that use agravmag
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  call sph_grav_forces
#endif

! Calculate gravitational forces on all sink particles
#if defined(SINKS)
  call sph_sink_forces
#endif

! Calculate hydro forces on all SPH particles
#if defined(HYDRO)
  call sph_hydro_forces
#endif

! Add wind contribution to acceleration here (for now)
#if defined(STELLAR_WIND)
  do p=1,ptot
     if (sph(p)%accdo) sph(p)%a(1:NDIM) = &
          &sph(p)%a(1:NDIM) + sph(p)%a_wind(1:NDIM)
  end do
#endif

! Add turbulent forcing
#if defined(TURBULENT_FORCING)
  call turb_forces
#endif

! Apply correction to integrated quantities if using Leapfrog-KDK scheme
#if defined(LEAPFROG_KDK)
  call leapfrog_kdk_correction_terms
#endif

! Output diagnostics and/or snapshot files
  call sph_output

#if defined(GHOST_PARTICLES)
  ! Any remaining ghosts are now wrong!
  pghost = 0
  pperiodic = 0
#endif

! Calculate new timesteps for particles
  call sph_timesteps

! Check and flag if any active particles neighbour particles with 
! relatively long timesteps.
#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
  call check_neighbour_timesteps
#endif

! Either synchronise all particle timesteps if flag is activated, or reduce 
! the timesteps of any particles with short-timestep neibs to new 'safe' value
  call reduce_timesteps

! Search for new sinks and accrete particles into existing sinks
#if defined(SINKS) && defined(GRAVITY) && defined(HYDRO) && !defined(NO_ACCRETION)
  call sink_update
#endif

! Advance all time variables
  n      = n + 1
  nsteps = nsteps + 1
  time   = time + timestep

! Update radiative cooling terms
#if defined(HYDRO) && defined(RAD_WS)
  call rad_ws_update
#endif

! Integration scheme to advance particle properties
  call sph_advance
#if defined(SINKS)
  call sink_advance
#endif

#if defined(USE_MPI)
  call mpi_loadbalance_step
#endif

! Remove any outlying particles depending on chosen criteria
#if defined(REMOVE_OUTLIERS)
  call remove_outlying_particles
#endif

! Add any new SPH particles in this routine
#if defined(STUFF)
  call add_new_particles
#endif

! Check if we need to calculate next turbulent forcing term
#if defined(TURBULENT_FORCING)
  call turb_check_time
#endif

#if defined(USE_MPI)
  calctime = calctime + MPI_WTIME() - calcstart
  diag_calctime = calctime - start_calctime
#endif

  return
END SUBROUTINE sph_integrate
