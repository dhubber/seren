! WRITE_RAD_WS_TEST_DATA.F90
! D. A. Hubber - 8/2/2009
! Write out data relevant for RAD_WS testing.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_rad_ws_test_data
  use particle_module, only : sph,ptot
  use hydro_module
  use scaling_module
  use type_module
  use Eos_module
  use filename_module, only : run_id, run_dir
  use time_module, only : nsteps, time
#if defined(USE_MPI)
  use mpi_communication_module
  use mpi
#endif
  implicit none

#if defined(DEBUG_RAD)
  character(len=256) :: out_file   ! output file name
  integer :: p                     ! Particle counter
  real(kind=PR) :: kappamax        ! kappa of most dense particle
  real(kind=PR) :: radmax          ! Radial distance of most dense particle
  real(kind=PR) :: rhomax          ! Maximum density
  real(kind=PR) :: taumax          ! Optical depth of most dense particle
  real(kind=PR) :: Tmax            ! Temperature of most dense particle
#if defined(USE_MPI)
  integer       :: ierr            ! MPI error value
  integer       :: max_rho_rank(1) ! rank of task with densest particle
  real(kind=PR) :: data_send(1:5)  ! data of densest particle for each task
  real(kind=PR) :: data_recv(1:5, 0:lastrank) ! receive array for above data
#endif

  debug2("Writing data for polytropic cooling test [write_rad_ws_test_data.F90]")

! Find most dense particle and record attributes
  rhomax = 0.0_PR
  do p=pgravitystart,ptot
     if (sph(p)%rho > rhomax) then
        rhomax   = sph(p)%rho
        Tmax     = sph(p)%temp
        kappamax = rad_info(2,p)
        taumax   = rad_info(4,p)
        radmax   = rad_info(8,p)
     end if
  end do
  
#if defined(USE_MPI)
  ! Send the data from each task to the root task
  data_send = (/rhomax, Tmax, kappamax, taumax, radmax/)
  call MPI_GATHER(data_send, 5, MPI_REAL_PR, data_recv, &
     &5, MPI_REAL_PR, 0, MPI_COMM_WORLD, ierr)
  
  if (rank/=0) return
  
  max_rho_rank = maxloc(data_recv(1, :)) - 1
    ! minus 1 to account for maxloc counting from 1
  
  rhomax = data_recv(1, max_rho_rank(1))
  Tmax = data_recv(2, max_rho_rank(1))
  kappamax = data_recv(3, max_rho_rank(1))
  taumax = data_recv(4, max_rho_rank(1))
  radmax = data_recv(5, max_rho_rank(1))
#endif
  
  out_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".rad"

! If first step, overwrite to create new file. Else, append to existing file
  if (nsteps <= 1) then
     open(1,file=out_file,status="unknown",form="formatted",position="rewind")
  else
     open(1,file=out_file,status="unknown",form="formatted",position="append")
  end if

! Write data and close file
  write(1,'(11E15.7)') time*tscale,rhomax*rhoscale*rhocgs,Tmax,&
       &kappamax*kappascale*kappacgs,taumax,&
       &abs(radmax)*Escale*Ecgs/tscale/tcgs/mscale/mcgs,&
       &log10(rhomax*rhoscale*rhocgs),log10(Tmax),&
       &log10(kappamax*kappascale*kappacgs),log10(taumax),&
       &log10(abs(radmax)*Escale*Ecgs/tscale/tcgs/mscale/mcgs)
  close(1)
#endif

  return
END SUBROUTINE write_rad_ws_test_data

