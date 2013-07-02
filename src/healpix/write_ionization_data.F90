! WRITE_IONIZATION_DATA.F90
! D. A. Hubber - 19/10/2009
! Write out various ionization quantites to file.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_ionization_data
  use HP_module
  use hydro_module
  use particle_module
  use filename_module
  use scaling_module
  use time_module
  use type_module
  implicit none

  character(len=256) :: out_file    ! Name of output file
  integer :: p                      ! Particle counter
  real(kind=PR) :: raverage         ! Average radius of shell particles
#if defined(DEBUG_HP_OUTPUT)
  integer :: pfront                 ! No. of particles at centre of IF
  integer :: shellparticles         ! No. of shell particles in IF
#endif

  debug2("[write_ionization_data.F90]")

! Record no. of ionized particles
  pionized = 0
  raverage = 0.0_PR
  do p=1,ptot
     if (sph(p)%temp > 0.5_PR*(Tneut + Tion) .or. &
          &sph(p)%tempmin > 0.5_PR*(Tneut + Tion)) pionized = pionized + 1
  end do

! Diagnostic output
! ----------------------------------------------------------------------------
#if defined(DEBUG_HP_OUTPUT)
#ifndef NOTEMPSMOOTH
  out_file = trim(adjustl(run_id))//".IFposition"
  open(unit=602,file=out_file,position="append")
  do p=1,ptot
     if (sph(p)%temp > 0.1_PR*Tion .and.sph(p)%temp < 0.9_PR*Tion .or. &
          sph(p)%tempmin > 0.1_PR*Tion .and. sph(p)%tempmin < 0.9_PR*Tion) &
          &write(602,'(2E15.7)') time*tscale,&
          &sqrt(sph(p)%r(1)**2 + sph(p)%r(2)**2 + sph(p)%r(3)**2)
  end do
  close(602)

  shellparticles = 0
  do p=1,ptot
     if (sph(p)%temp > Tneut .and. sph(p)%temp < Tion) then
        shellparticles = shellparticles + 1
     end if
  end do

  pfront = 0
  do p=1,ptot
     if (sph(p)%temp > 0.25_PR*Tion .and. sph(p)%temp < 0.75_PR*Tion) then
        pfront = pfront + 1
        raverage = raverage + &
             &sqrt(dot_product(sph(p)%r(1:NDIM),sph(p)%r(1:NDIM)))
     end if
  end do
  if (pfront > 0) raverage = raverage / real(pfront,PR)

  out_file = trim(adjustl(run_id))//".shellparticles"
  open(unit=621,file=out_file,position="append")
  write(621,*) nsteps,time*tscale,shellparticles
  close(621)
#endif

  out_file = trim(adjustl(run_id))//".maxres"
  open(unit=622,file=out_file,position="append")
  write(622,*) nsteps,time*tscale,HPmaxres
  close(622)
#endif
! ----------------------------------------------------------------------------

  out_file = trim(adjustl(run_id))//".particles"
  open (unit=601,file=out_file,position="append")
  write(601,*) nsteps,time*tscale,pionized,ptot-pionized,raverage
  close(601)

  !write(6,*) "No. ionized particles :",pionized

  return
END SUBROUTINE write_ionization_data

