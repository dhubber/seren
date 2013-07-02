! IC_DEFLECTION.F90
! R. Smith - 22/10/2011
! Place a star particle on a straight line trajectory through the centre of a 
! glass gas cube. Gas cube is read in from file (and mass scaled up from 
! unity to 1000.0, and all thermal properties zeroed (so run with isothermal
! equation of state). Gas cube is shifted so centred on origin initially. - gas particle positions are fixed in this test. We shift gas cube in y direction to get new encounters along trajectory 
! Star is positioned along/on x-axis, with velocity (in -ve x-direction), and 
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
  real(kind=PR) :: m1                    ! Mass of star
  real(kind=PR) :: mcloud                ! Mass of gas cloud
  real(kind=PR) :: vi                    ! Velocity (in -ve x direction)
  real(kind=PR) :: xi                    ! Initial x position on x-axis
  real(kind=PR) :: s_radius              ! Radius of 'sinks'
  real(kind=PR) :: yshift                ! cube shift up y-axis
  real(kind=PR), allocatable :: rp(:,:)  ! Pos. of SPH part.

! Read parameters file
  call default_parameters

! Read in IC_DEFLECTION params file  
  write(6,*) "Opening ic_deflection params file"
  open(unit=1,file="deflection_cube.dat",status='old')
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*) file1                        ! Gas input
  read(1,*) file1_form	                 ! Gas input file format
  read(1,*) m1                           ! Mass of star
  read(1,*) mcloud                       ! Mass of gas cube
  read(1,*) xi                           ! Position on x-axis
  read(1,*) vi                           ! Velocity in -ve x direction
  read(1,*) s_radius	                 ! softening length of star
  read(1,*) yshift                       ! shift of cube up y-axis
  
! Initialise some variables using parameters
  call set_default_particle_types
  call initialize_seren_variables_1
  call units

! Read in SPH particles from file (makes this program more general...)
  call read_data(file1,file1_form)

! set periodic variables
  periodic_min(1:3)=-0.5_DP
  periodic_max(1:3)=0.5_DP

! first centre cube on origin, then shift vertically by yshift
! scale up mass, and zero thermal properties
  do p=1,ptot
     sph(p)%r(1:3) = sph(p)%r(1:3)-0.5_DP
     sph(p)%r(3)   = sph(p)%r(3) + yshift
     sph(p)%m      = sph(p)%m*mcloud
     sph(p)%temp   = 0.0_DP
#if defined(INTERNAL_ENERGY)
     sph(p)%u      = 0.0_DP
#endif
    call check_boundary_conditions(sph(p)%r(1:3),sph(p)%v(1:3))
  end do

! Set particle types and allocate memory
  pboundary = 0
  picm      = 0
  pgas      = ptot
  pcdm      = 0
  pdust     = 0
  pion      = 0
  stot      = 1
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
