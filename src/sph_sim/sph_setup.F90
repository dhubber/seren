! SPH_SETUP.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Performs initialisation routines to prepare for main SPH simulation.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_setup
#if defined(BH_TREE)
  use interface_module, only : BHhydro_stock,BHhydro_update_hmax
#endif
  use filename_module, only : restart
  use particle_module
  use diagnostics_module
  use time_module
  use type_module
  use tree_module
  implicit none

  integer :: p          ! Particle counter
  logical :: done       ! Are we done?
  
  debug1("Setting up the SPH simulation [sph_setup.F90]")

! Setting up array-limits for different particle types
  call types

! Set-up HEALPix rays
#if defined(HEALPIX)
  call initialize_HP_sources
#endif

! Initialize certain variables before first force calculation
  call initialize_sph_variables_1

! Build and stock trees for first time
  call tree_update(nbuild,nstock)

! Make initial guesses of h either using tree or a global estimate
#if defined(BH_TREE) && !defined(CONSTANT_H)
  if (.not. restart) call BHhydro_hguess
  call BHhydro_update_hmax(cmax_hydro,ctot_hydro,&
       &ltot_hydro,first_cell_hydro,last_cell_hydro,BHhydro)
#elif !defined(CONSTANT_H)
  if (.not. restart) call h_guess
#endif

! Calculating initial SPH quantities (h, neibs, rho, etc.)
#if defined(USE_MPI)
  call mpi_sph_update
#else
  call sph_update
#endif

! Calculate HEALPix quantities if activated.
#if defined(HEALPIX)
  call HP_update
  call tree_update(nbuild,nstock)
#endif

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

  end do
! ----------------------------------------------------------------------------
#endif

! Initialize all thermal properties depending on options used
#if defined(HYDRO)
  call initialize_thermal_properties
#endif

! Copy all particle data to corresponding ghosts now we have thermal info.
#if defined(GHOST_PARTICLES) && defined(PERIODIC)
  if (pperiodic > 0) call copy_data_to_ghosts
#endif

! Send gather and scatter ghosts for calculating hydro forces etc.
#if defined(USE_MPI)
  call mpi_scatter_ghosts
#endif

! Updating tree properties
#if defined(BH_TREE) && defined(SELF_GRAVITY)
  call BHgrav_stock
#endif

#if defined(USE_MPI) && defined(SELF_GRAVITY)
! Slice off top levels of gravity tree
  call BHgrav_slice
! Share pruned gravity trees
  call share_pruned_gravtrees
#endif

! If using OSPH or RTSPH, need to iterate at least once to obtain consistent 
! density with the chosen weighting quantity (e.g. Aent,u,temp)
#if defined(OSPH) || defined(RTSPH)
  call sph_update
  call initialize_thermal_properties
  call sph_update
  call initialize_thermal_properties
#endif

! Calculate gravitational forces on all SPH particles
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  call sph_grav_forces
#endif

! Calculate gravitational forces on all sink particles
#if defined(SINKS) && defined(GRAVITY)
  call sph_sink_forces
#endif

! Calculate hydro forces on all SPH particles
#if defined(HYDRO)
  call sph_hydro_forces
#endif

! Add wind contribution to acceleration here
#if defined(STELLAR_WIND)
  do p=1,ptot
     sph(p)%a(1:NDIM) = sph(p)%a(1:NDIM) + sph(p)%a_wind(1:NDIM)
  end do
#endif

! Add turbulent forcing
#if defined(TURBULENT_FORCING)
  call turb_forces
#endif

! Compute conductivity terms here (for now)
#if defined(RAD_WS) && defined(FLUX_LIMITED_DIFFUSION)
!$OMP PARALLEL DO DEFAULT(SHARED)
  do p=1,ptot
     call find_idens(sph(p)%rho,sph(p)%idens)
     call find_itemp(sph(p)%temp,sph(p)%itemp)
     call conductivity(p)
  end do
!$OMP END PARALLEL DO
#endif

! Calculate initial diagnostic quantities
#if defined(DEBUG_DIAGNOSTICS)
  call diagnostics
  etot0 = etot
#endif

! Initialize other key variables (after initial force calculation)
  call initialize_sph_variables_2

  return
END SUBROUTINE sph_setup
