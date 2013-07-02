! IC_RANDOM_CUBE.F90
! D. A. Hubber - 21/9/2007
! Generates a random distribution of particles depending on dimensionality. 
! For NDIM = 1, a line of equally spaced particles 
!          = 2, a square sheet of randomally placed particles
!          = 3, a cube of randomally placed particles 
! The particles are given equal masses and are scaled such that the 
! total mass is unity.  The line/sheet/cube is placed in the positive 
! x-axis/quadrant/octant such that one corner is at the origin. 
! ptot           : No. of particles in cube
! length         : Side length of cube
! out_file       : Name of output file
! out_file_form  : Output file format
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_random_cube
  use interface_module, only : write_data
  use filename_module, only : out_file_form
  use scaling_module, only : rscale,mscale
  use particle_module
  use type_module
  implicit none

  character(len=256) :: out_file  ! Name of output file
  integer :: k                    ! Dimension counter
  integer :: p                    ! Particle counter
  real(kind=PR) :: length         ! Side-length of line/square/cube
  real(kind=PR) :: randnumb       ! Random number

#ifndef DIMENSIONLESS
  write(6,*) "Compiler flag error : Only works with DIMENSIONLESS flag on"
  stop
#endif

! Set parameters to default values
  call default_parameters

  write(6,*) "------------------------------"
  write(6,*) "       ic_random_cube         "
  write(6,*) "------------------------------"
  write(6,*) 
  write(6,*) "Reading in input parameters"
  write(6,*) "Number of particles : "
  read(5,*) ptot
  write(6,*) "Side length : "
  read(5,*) length
  write(6,*) "Name of output file : "
  read(5,*) out_file
  write(6,*) "Format of output file : "
  read(5,*) out_file_form

  out_file      = trim(adjustl(out_file))
  out_file_form = trim(adjustl(out_file_form))

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types
  call units

! Allocate memory for particle data
  call allocate_memory(.FALSE.)
  rscale    = 1.0_PR
  mscale    = 1.0_PR
  pgas      = ptot
  picm      = 0
  pboundary = 0
  call types

! If NDIM == 1, generate line of equally spaced particles
! ---------------------------------------------------------------------------- 
#if NDIM==1
  do p=1,ptot
     sph(p)%r(1) = length * (real(p) - 0.5_PR) / real(ptot,PR)
     sph(p)%m = 1.0_PR / real(ptot,PR)
  end do

! Else if NDIM ==2 or NDIM == 3, generate random sheet or cube
! ----------------------------------------------------------------------------
#else
  call random_seed
  do p=1,ptot
     do k=1,NDIM
        call random_number(randnumb)
        sph(p)%r(k) = length*randnumb
     end do
     sph(p)%m = 1.0_PR / real(ptot,PR)
     sph(p)%porig = p
  end do
#endif
! ----------------------------------------------------------------------------

! Write data to file
  call write_data(out_file,out_file_form)

! Clean up memory
  call clean_up

  stop
END PROGRAM ic_random_cube

