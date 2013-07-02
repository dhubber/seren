! IC_DEFLECTION.F90
! R. Smith - 22/10/2011
! Place a star particle on a straight line trajectory through the centre of a 
! glass gas sphere. Gas sphere is read in from file (and mass scaled up from 
! unity to 1000.0, and all thermal properties zeroed (so run with isothermal
! equation of state) - as gas particle positions are fixed in this test). 
! Star is positioned along x-axis, with velocity (in -ve x-direction), and 
! chosen mass. These parameters are read in from params file
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_deflection
  use particle_module
  use hydro_module
  use sink_module
  use filename_module
  use periodic_module
  use type_module
  use scaling_module
  use constant_module
  implicit none

  character(len=256) :: file1_form       ! File 1 format
  character(len=256) :: file1            ! File with RANDOM_CUBE
  character(len=256) :: out_file         ! Output file name
  character(len=256) :: out_file__form   ! File 1 format
  integer :: p                           ! Particle counter
  integer :: p1                          ! No. of gas particles
  real(kind=PR) :: m1                    ! Mass of star
  real(kind=PR) :: mcloud                ! Mass of gas cloud
  real(kind=PR) :: vi                    ! Velocity (in -ve x direction)
  real(kind=PR) :: xi                    ! Initial x position on x-axis
  real(kind=PR) :: s_radius              ! Radius of 'sinks'
  real(kind=PR) :: yrot                  ! rotation angle of gas about y-axis
  real(kind=PR), allocatable :: rp(:,:)  ! Pos. of SPH part.

! Read parameters file
  call default_parameters

! Read in IC_DEFLECTION params file  
  write(6,*) "Opening ic_deflection params file"
  open(unit=1,file="deflection.dat",status='old')
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*) file1                        ! Gas input
  read(1,*) file1_form	                 ! Gas input file format
  read(1,*) m1                           ! Mass of star
  read(1,*) mcloud                       ! Mass of gas cloud
  read(1,*) xi                           ! Position on x-axis
  read(1,*) vi                           ! Velocity in -ve x direction
  read(1,*) s_radius	                 ! softening length of star
  read(1,*) yrot                         ! angle of gas rotation about y-axis
  
! Initialise some variables using parameters
  call initialize_seren_variables_1
  call units

! Read in SPH particles from file (makes this program more general...)
  call read_data(file1,file1_form)
  p1 = ptot

! make a dummy matrix of ga particle x and z positions
! in order to make a rotation
  allocate(rp(1:3,1:p1))
  do p=1,p1
    rp(1,p)=sph(p)%r(1)
    rp(3,p)=sph(p)%r(3)
  enddo

! rotate gas particles about y-axis, 
! scale up mass, and zero thermal properties)
  do p = 1,p1
     sph(p)%r(1)  = cos(yrot)*rp(1,p) + sin(yrot)*rp(3,p)
     sph(p)%r(3)  = -sin(yrot)*rp(1,p) + cos(yrot)*rp(3,p)
     sph(p)%m     = sph(p)%m*1000.0_DP
     sph(p)%temp  = 0.0_DP
     sph(p)%u     = 0.0_DP
  end do
  deallocate(rp)

! Set particle types and allocate memory
  pboundary = 0
  picm      = 0
  pgas      = p1 
  pcdm      = 0
  pdust     = 0
  pion      = 0
  stot      = 1
  ptot      = p1 
  call types

! Set star properties in Seren arrays
  sink(1)%r(2:3)    = 0.0_DP
  sink(1)%r(1)      = xi
  sink(1)%v(2:3)    = 0.0_DP
  sink(1)%v(1)      = -vi
  sink(1)%m         = m1
  sink(1)%radius    = s_radius
  sink(1)%h         = INVKERNRANGE*sink(1)%radius

! Initialise SEREN stuff...!
! ----------------------------------------------------------
  restart = .false.
  call initialize_seren_variables_2

! Write to file
  out_file = trim(adjustl(out_file))
  call write_data(out_file,out_file_form)
!#if defined(DEBUG_PLOT_DATA)
!  call write_data_debug("DEFLECTION.debug.dat",rcentre(1:NDIM))
!#endif

! Clean-up all memory
  call clean_up
  
  stop

END PROGRAM ic_deflection
