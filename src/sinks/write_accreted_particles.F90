! WRITE_ACCRETED_PARTICLES.F90
! D. A. Hubber - 16/4/2009
! Write particle ids, and other information, of accreted particles to file.
! A separate file of the form "run_id.accretionX" is generated where X
! is the sink id.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_accreted_particles(s,pp_tot,pp_templist)
  use filename_module
  use particle_module
  use sink_module
  use Nbody_module
  use scaling_module
  use time_module
  implicit none

  integer, intent(in) :: s                      ! sink id
  integer, intent(in) :: pp_tot                 ! no. of accreted particles
  integer, intent(in) :: pp_templist(1:pp_tot)  ! list of accreted particles

  character(len=8)  :: file_ext    ! filename extension for data output
  character(len=256) :: out_file   ! filename extension for data output
  logical :: ex                    ! Does file exist already?
  integer :: i                     ! Auxilary counter
  integer :: p                     ! Particle id
  integer :: idummy(1:2)           ! Dummy integer array
  real(kind=PR) :: rdummy(1:5)     ! Dummy real array

  debug2("Writing sink information to files [write_sink_data.F90]")

! Make filename from runid and sink number
  if (s >= 100) then
     write(file_ext,"(I3)") s
  else if (s>=10) then
     write(file_ext,"(I2)") s
  else
     write(file_ext,"(I1)") s
  end if
  out_file = trim(adjustl(run_dir))//trim(adjustl(run_id))&
       &//".accretion"//trim(adjustl(file_ext))

! Check if file exists
  inquire(file=out_file,exist=ex)

! Read through file to synchronise files with time of simulation.
  if (ex) then
     open(1,file=out_file,status="unknown",&
          &form="formatted",position="append")
     do
        read(1,'(I8,21E18.10)',end=10,err=10) idummy(1:2),rdummy(1:5)
        if (rdummy(1) > time*tscale) exit
     end do
10   backspace (1,err=20)
  else
     open(1,file=out_file,status="unknown",form="formatted")
  end if

! Write all useful sink/star information to file
20 do i=1,pp_tot
     p = pp_templist(i)
     write(1,'(2I8,5E18.10)') nsteps,sph(p)%porig,time*tscale,&
          &sph(p)%m*mscale,sph(p)%r(1:NDIM)
  end do

  close(1)

  return
END SUBROUTINE write_accreted_particles
