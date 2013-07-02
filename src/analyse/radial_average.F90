! RADIAL_AVERAGE.F90
! D. A. Hubber - 25/08/2010
! ..
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM radial_average
  use particle_module
  use filename_module
  use scaling_module
  use time_module
  use type_module
  use neighbour_module
  use constant_module
  use sink_module
  use hydro_module
  implicit none

  character(len=256) :: file1               ! ..
  character(len=256) :: file1_form          ! ..
  character(len=256) :: out_file            ! Name of output file
  character(len=256) :: params_file         ! ..
  integer :: i                              ! ..
  integer :: ndiv                           ! ..
  integer :: p                              ! ..
  integer, allocatable :: nin(:)            ! ..
  real(kind=PR) :: dr(1:NDIM)               ! Relative position
  real(kind=PR) :: drad                     ! ..
  real(kind=PR) :: drmag                    ! Distance
  real(kind=PR) :: drsqd                    ! Distance squared
  real(kind=PR) :: enttot                   ! ..
  real(kind=PR) :: radmax                   ! ..
  real(kind=PR) :: rorigin(1:3)             ! ..
  real(kind=PR), allocatable :: rdata(:,:)  ! ..

! Set default parameters
  call default_parameters

! Read polytrope set-up file  
  open(unit=1,file='radial.dat',status='old')
  read(1,*) file1
  read(1,*) file1_form
  read(1,*) params_file
  read(1,*) out_file
  read(1,*) radmax
  read(1,*) ndiv
  read(1,*) rorigin(1:3)

! Read-in parameters file
  call read_parameters(params_file)
  restart = .false.

! Initialise some variables using parameters
  call initialize_seren_variables_1

! Checking compiler flags and parameter values
  call sanitycheck

! Setting up scaling units for simulation
  call units

! Reading in formatted data file
  call read_data(file1,file1_form)

! Perform all necessary SPH calculations
  call sph_setup

  drad = radmax/real(ndiv,PR)
  allocate(nin(0:ndiv))
  allocate(rdata(1:6,0:ndiv))
  do i=0,ndiv
     nin(i) = 0
     rdata(1:6,i) = 0.0_PR
     drmag = radmax*real(i,PR)/real(ndiv,PR)
     rdata(1,i) = drmag
  end do
  enttot = 0.0_PR

  do p=1,ptot
     call distance2(rorigin(1:NDIM),p,dr(1:NDIM),drsqd)
     drmag = sqrt(drsqd)
     i = int((drmag/drad) + 0.5_PR)
     if (i < 0 .or. i > ndiv) cycle
     nin(i) = nin(i) + 1
     rdata(2,i) = rdata(2,i) + rho(p)
     rdata(3,i) = rdata(3,i) + press(p)
     rdata(4,i) = rdata(4,i) + press(p)/(rho(p)**(gamma))
     rdata(5,i) = rdata(5,i) + temp(p)
     if (i <= ndiv) enttot = enttot + parray(MASS,p)*press(p)/(rho(p)**(gamma))
  end do

  do i=0,ndiv
     if (nin(i) > 0) then
        rdata(2,i) = rdata(2,i) / real(nin(i),PR)
        rdata(3,i) = rdata(3,i) / real(nin(i),PR)
        rdata(4,i) = rdata(4,i) / real(nin(i),PR)
        rdata(5,i) = rdata(5,i) / real(nin(i),PR)
     end if
  end do
  rdata(6,:) = enttot

  open(unit=1, file=out_file, status='unknown')
  do i=0,ndiv
     write(1,'(6E20.12)') rdata(1:6,i)
  end do
  close(1)

  deallocate(rdata)
  deallocate(nin)

  STOP
END PROGRAM radial_average
