! NBODY_OUTPUT.F90
! D. A. Hubber - 28/01/2008
! Write all output from N-body integrator to screen and to files
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_output
  use definitions
  use sink_module
  use scaling_module
  use time_module
  use Nbody_module
  use filename_module
  implicit none

  character(len=11)  :: file_ext     ! Filename extension for data output
  character(len=6)   :: file_numb    ! Filename extension for data output
  character(len=256) :: out_data     ! Data snapshot filename
  character(len=256) :: out_debug    ! Debug data snapshot filename
  integer :: s                       ! Star counter
  real(kind=DP) :: mtot              ! Total mass
  real(kind=DP) :: rcom(1:NDIM)      ! Position of centre of mass
  real(kind=DP) :: vcom(1:NDIM)      ! Velocity of centre of mass

  debug2("Writing N-body output to screen and file [nbody_output.F90]")
  debug_timing("NBODY_OUTPUT")


! Output periodic snapshot file
! ----------------------------------------------------------------------------
  if (time >= nextsnap .and. snapshot < 100000) then
     lastsnap = nextsnap
     nextsnap = nextsnap + snaptime
     snapshot = snapshot + 1

     call copy_stars_to_sinks

     write(file_numb,"(I5.5)") snapshot
     file_ext = "."//file_numb

     out_data = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
          &trim(adjustl(fileform_ext))//trim(adjustl(file_ext))
     call write_data(out_data,out_file_form) 

     ! Write name of snapshot file to log for potential restart
     open(1,file=restart_log,status="unknown",form="formatted")
     write(1,'(a)') out_data
     write(1,'(a)') out_file_form
     close(1)

     ! Calculate and write binary properties to file
#if defined(BINARY_STATS)
     call binary_search
#endif
#if defined(LAGRANGIAN_RADII)
     call write_lagrangian_radii(rzero)
#endif

  end if
! ----------------------------------------------------------------------------


! Write to screen diagnostic information 
! ----------------------------------------------------------------------------
#if defined(DEBUG_DIAGNOSTICS)
  if (n == nresync .and. nmaxsteps == ndiagnext) then
     ndiagnext = ndiagnext + ndiagstep
     call diagnostics
#if defined(TIMING)
     call write_timing_stats
#endif
  end if
#endif
! ----------------------------------------------------------------------------


! Write star information to file
! ----------------------------------------------------------------------------
#if defined(DEBUG_OUTPUT_STAR_DATA)
  if (nsteps == nsinknext) then
     nsinknext = nsinknext + nsinkstep
     call write_star_data
  end if
#endif


! Write regular time information to screen each step
! ----------------------------------------------------------------------------
#if defined(DEBUG1)
  if (nsteps == noutput) then
100  format (1X,'Step # : ',I10,6X,'snapshot : ',I5)
     write (6,100) nsteps,snapshot
200  format (1X,'Time : ',G14.8,1X,a,4X,'dt : ',G14.8,1X,a)
     write (6,200) time*tscale,adjustl(trim(tunit)),&
          &timestep*tscale,adjustl(trim(tunit))
     noutput = noutput + noutputstep
  end if
#endif


  return
END SUBROUTINE nbody_output
