! NBODY_SPH_SIMULATION.F90
! D. A. Hubber - 23/7/2010
! Main control routine for hybrid N-body/SPH simulation
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_sph_simulation
  use definitions
  use seren_sim_module
  use filename_module
  use time_module
!  use particle_module, only : ptot
!  use neighbour_module, only : pp_gather
!  use Nbody_module, only : nbody_frac
#if defined(NBODY_SIMULATION)
  use sink_module, only : sink_frac
#endif
#if defined(IONIZING_UV_RADIATION)
  use type_module, only : pionized
#endif
  implicit none

#if defined(DEBUG_PLOT_DATA) && defined(NBODY_SPH_SIMULATION)
  character(len=256) :: out_debug    ! data snapshot output
#endif

  debug1("Hybrid N-body/SPH simluation [nbody_sph_simulation.F90]")

! Set flags for hybrid N-body/SPH simulation
  sph_sim       = .false.
  nbody_sph_sim = .true.
  nbody_sim     = .false.

! Setting up hybrid simulation
  call nbody_sph_setup

! Write initial snapshot 
  if (.not. inifile) call write_data(out_init,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  out_debug = trim(adjustl(run_dir))//trim(adjustl(run_id))//".debug.ini"
  if (.not. inifile) call write_data_debug(out_debug,rzero(1:NDIM))
#endif
  inifile = .true.


! Main N-body/SPH integration loop
! ============================================================================
  do

     ! Exit main loop if simulation has reached endtime
     if (time >= nbody_sph_endtime) exit     

     ! Exit main loop if too few particles remain (e.g. from accretion)
     !if (ptot <= pp_gather) exit

     ! Exit loop to switch to N-body integration
#if defined(NBODY_SIMULATION)
     if (sink_frac >= real(nbody_frac,PR)) exit
#endif

     ! Exit loop if all particles become ionized
#if defined(IONIZING_UV_RADIATION)
     !if (pionized == ptot) exit
#endif

     ! Main SPH simulation control subroutine
     call nbody_sph_integrate

  end do
! ============================================================================


  return
END SUBROUTINE nbody_sph_simulation
