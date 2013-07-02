! SPH_SIMULATION.F90
! C. P. Batty & D. A. Hubber - 8/12/2006 
! Main control subroutine for SPH simulation.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_simulation
  use interface_module, only : write_data,write_data_debug
  use definitions
  use seren_sim_module
  use filename_module
  use sink_module
  use scaling_module
  use hydro_module
  use particle_module
  use neighbour_module, only : pp_gather
  use time_module, only : time, sph_endtime
#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  use sink_module, only : stot, sink_frac
  use Nbody_module, only : nbody_frac, nbody_endtime
#endif
#if defined(IONIZING_UV_RADIATION)
  use type_module, only : pionized
#endif
#ifdef USE_MPI
  use mpi_communication_module
#endif
  implicit none

#if defined(DEBUG_PLOT_DATA) || defined(DEBUG_GRID_RESULTS)
  character(len=256) :: out_debug    ! data snapshot output
#endif

  debug1("Main SPH simulation [sph_simulation.F90]")

! Set flags for SPH simulation
  sph_sim       = .true.
  nbody_sph_sim = .false.
  nbody_sim     = .false.

! If we do not have any SPH particles, exit setup to switch to N-body
  if (ptot == 0) return

! Setting up SPH simulation 
  call sph_setup

! Write initial snapshot 
  if (.not. inifile) call write_data(out_init,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  out_debug = trim(adjustl(run_dir))//trim(adjustl(run_id))//".debug.ini"
#if defined(USE_MPI)
  out_debug = trim(out_debug)//".MPI."//trim(adjustl(MPI_ext))
#endif
  if (.not. inifile) call write_data_debug(out_debug,rzero(1:NDIM))
#endif
#if defined(DEBUG_GRID_RESULTS)
  out_debug = trim(adjustl(run_dir))//trim(adjustl(run_id))//".grid.ini"
#if defined(USE_MPI)
  out_debug = trim(out_debug)//".MPI."//trim(adjustl(MPI_ext))
#endif
  if (.not. inifile) call write_data_grid_results(out_debug)
#endif
  inifile = .true.

! Main SPH integration loop
! ============================================================================
  do

     ! Exit main loop if simulation has reached endtime
     if (time >= sph_endtime) exit     

     ! Exit main loop if too few particles remain (e.g. from accretion)
#if defined(USE_MPI)
     if (totalptot <= (2*pp_gather*numtasks)) then
        write (6,*) "ptot = ", ptot
        write (6,*) "pp_gather * 2 = ", pp_gather*2
        write (6,*) "numtasks = ", numtasks
        stop "Too few particles, graceful exit not yet implemented"
     end if
#else
     if (ptot <= pp_gather) exit
#endif

     ! Exit loop to switch to N-body integration
#if defined(NBODY_SIMULATION) || defined(NBODY_SIMULATION)
     if (sink_frac >= real(nbody_frac,PR)) exit
#endif

     ! Exit loop to switch to hybrid simulation
#if defined(NBODY_SPH_SIMULATION)
     if (sink_frac >= 1.0_PR) exit
#endif

     ! Exit loop if all particles become ionized
#if defined(IONIZING_UV_RADIATION)
     if (pionized == ptot) exit
#endif

     ! Exit if maximum density has been reached for Polytropic cooling test
#if defined(DEBUG_RAD)
     if (maxval(sph%rho) > (1.0E-2_PR/(rhoscale*rhocgs))) exit
#endif

     ! Exit if 'singularity' has been reached for freefall collapse test
#if defined(FREEFALL_TEST)
     if (maxval(sph(1:ptot)%rho) > 1.E7_PR) exit
#endif

     ! Main SPH simulation control subroutine
     call sph_integrate

  end do
! ============================================================================


! Write final snapshot
  call write_data(out_final,out_file_form)
  call diagnostics
#if defined(DEBUG_PLOT_DATA)
  out_debug = trim(adjustl(run_dir))//trim(adjustl(run_id))//".debug.fin"
#if defined(USE_MPI)
  out_debug = trim(out_debug)//".MPI."//trim(adjustl(MPI_ext))
#endif
  call write_data_debug(out_debug,rzero(1:NDIM))
#endif
#if defined(DEBUG_GRID_RESULTS)
  out_debug = trim(adjustl(run_dir))//trim(adjustl(run_id))//".grid.fin"
#if defined(USE_MPI)
  out_debug = trim(out_debug)//".MPI."//trim(adjustl(MPI_ext))
#endif
  call write_data_grid_results(out_debug)
#endif
#if defined(TIMING)
  call write_timing_stats
#endif

  debug1("SPH simulation completed")

  return
END SUBROUTINE sph_simulation
