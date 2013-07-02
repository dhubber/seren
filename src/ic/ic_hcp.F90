! IC_HCP.F90
! D. A. Hubber - 21/9/2007
! Generates a regularly spaced hexagonal close packed array (hcp) of 
! particles.  Currently only works in two dimensions.
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_hcp
  use particle_module
  use hydro_module
  use filename_module
  use scaling_module
  use type_module
  implicit none

  character(len=256) :: out_file    ! Name of output file
  integer :: i                      ! x-dimension counter
  integer :: j                      ! y-dimension counter
  integer :: k                      ! z-dimension counter
  integer :: nx                     ! Particle per dimension
  integer :: ny                     ! Particle per dimension
  integer :: p                      ! Particle counter
  real(kind=PR) :: dx               ! x-particle spacing
  real(kind=PR) :: dy               ! y-particle spacing
  real(kind=PR) :: xlength          ! Side-length of line/square/cube
  real(kind=PR) :: ylength          ! Side-length of line/square/cube
  real(kind=PR) :: randnumb         ! Random number
  real(kind=PR) :: rp(1:NDIM)       ! Position of particle p
  real(kind=PR) :: rorigin(1:NDIM)  ! Centre of polytrope

#ifndef DIMENSIONLESS
  write(6,*) "Compiler flag error : Only works with DIMENSIONLESS flag on"
  stop
#endif
#if NDIM != 2
  write(6,*) "Compiler flag error : Only working in 2 dimensions"
  stop
#endif


! Set parameters to default values
  call default_parameters

  write(6,*) "------------------------------"
  write(6,*) "            ic_hcp            "
  write(6,*) "------------------------------"
  write(6,*) 
  write(6,*) "Input data..."
  open(unit=1,file='hcpparams.dat',status='old')
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) nx
  read(1,*) ny
  read(1,*) xlength
  close(1)

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types

! Calculate scaling units
  call units

  dx = xlength / real(nx,PR)
  dy = dx*0.5_PR*sqrt(3.0_PR)
  ylength = dy*real(ny,PR)
  ptot = nx*ny
  write(6,*) "Number of particles : ",ptot
  write(6,*) "ylength :",ylength,"   dx :",dx,"   dy :",dy

! Allocate memory for particle data
  pgas = ptot
  call allocate_memory(.FALSE.)
  call types

! ----------------------------------------------------------------------------
#if NDIM==2
  p = 0
  do i=1,nx
     do j=1,ny
        p = p + 1
        rp(1) = xlength*(real(i - 1,PR) - 0.5_PR*real(j,PR) - 0.25_PR) &
             & / real(nx,PR)
        rp(2) = ylength*(real(j - 1,PR) - 0.5_PR) / real(ny,PR)
        if (rp(1) < 0.0_PR)  rp(1) = rp(1) + xlength
        if (rp(1) > xlength) rp(1) = rp(1) - xlength
        if (rp(2) < 0.0_PR)  rp(2) = rp(2) + ylength
        if (rp(2) > ylength) rp(2) = rp(2) - ylength
        sph(p)%r(1:NDIM) = rp(1:NDIM)
     end do
  end do
#endif
! ----------------------------------------------------------------------------

! Set all other quantities
   do p=1,ptot
      sph(p)%m = 1.0_PR / real(ptot,PR)
      sph(p)%h = 1.0_PR
      sph(p)%rho = 1.0_PR
      sph(p)%v(1:NDIM) = 0.0_PR
#if defined(HYDRO)
      sph(p)%temp = 1.0_PR
      sph(p)%sound = 1.0_PR
#endif
   end do

! Write data to file
  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  rorigin(1:NDIM) = 0.0_PR
  call write_data_debug("ICHCP.debug.dat",rorigin(1:NDIM))
#endif

! Clean up memory
  call clean_up

  stop
END PROGRAM ic_hcp
