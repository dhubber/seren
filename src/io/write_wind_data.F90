! WRITE_WIND_DATA.F90
! D. A. Hubber - 25/03/2011
! ..
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_wind_data
  use particle_module
  use filename_module
  use time_module
  use HP_module
  use scaling_module
  implicit none

  character(len=256) :: out_file            ! ..
  logical :: ex
  integer :: i
  integer :: p                              ! ..
  integer, allocatable :: plist(:)          ! ..
  real(kind=PR) :: raverage                 ! ..
  real(kind=PR) :: rdummy(1:3)              ! ..
  real(kind=PR), allocatable :: rholist(:)  ! ..

  debug2("[write_wind_data.F90]")

  allocate(plist(1:ptot))
  allocate(rholist(1:ptot))

  do p=1,ptot
     plist(p) = p
!     rholist(p) = sph(p)%rho
     rholist(p) = dot_product(sph(p)%v(1:NDIM),sph(p)%v(1:NDIM))
  end do

  ! Use heapsort to order list in terms of increasing density
  call heapsort_real(ptot,rholist,plist)

  ! Now compute average position of 100 most dense particles
  raverage = 0.0_PR
  do i=1,100
     p = plist(ptot - 100 + i)
     raverage = raverage + sqrt(dot_product(sph(p)%r(1:NDIM),sph(p)%r(1:NDIM)))
  end do
  raverage = raverage / 100.0_PR

  write(6,*) "Wind bubble : ",raverage*rscale,time*tscale

  ! Determine filename for data
  out_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".wind"

  ! Check if file exists
  inquire(file=out_file,exist=ex)
  
  ! Read through file to synchronise files with time of simulation.
!  if (ex) then
!     open(1,file=out_file,status="unknown",&
!          &form="formatted",position="rewind")
!     do
!        read(1,'(I8,21E18.10)',end=10,err=10) rdummy(1:3)
!        if (rdummy(1) > time*tscale) exit
!     end do
!10   backspace (1,err=20)
!  else
     open(1,file=out_file,status="unknown",form="formatted",position="append")
!  end if

20 write(1,'(3E18.10)') time*tscale,raverage*rscale,&
       &(1.5_PR*M_loss*v_wind)**0.25_PR*sqrt(time)*rscale

  deallocate(rholist)
  deallocate(plist)

  return
END SUBROUTINE write_wind_data
