! IC_SPITZER.F90
! D. A. Hubber - 30/09/2008
! ...
! in_file       : Filename of relaxed unit-sphere
! in_file_form  : Format of input sphere file 
! out_file      : Output filename
! out_file_form : Format of output file
! runit         : Length unit
! munit         : Mass unit
! rcloud        : Radius of cloud
! mcloud        : Mass of cloud
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_spitzer
  use interface_module, only : read_data,write_data
  use particle_module
  use filename_module
  use scaling_module
  use time_module
  use type_module
  use neighbour_module
  use constant_module
  use hydro_module
  implicit none
  
  character(len=256) :: out_file    ! Name of output file
  integer :: p                      ! Particle counter
  real(kind=PR) :: mcloud           ! Mass of cloud
  real(kind=PR) :: mp               ! Mass of particle p
  real(kind=PR) :: rcentre(1:NDIM)  ! Centre of polytrope
  real(kind=PR) :: rcloud           ! Radius of cloud
  real(kind=PR) :: Tgas             ! Temperature of the gas
  
  debug1("Create uniform cloud for Spitzer test [ic_spitzer.F90]")
  
! Read parameters file
  call default_parameters
  
! Read polytrope set-up file  
  open(unit=1,file='spitzer.dat',status='old')
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) in_file      
  read(1,*) in_file_form 
  read(1,*) out_file     
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) runit
  read(1,*) munit
  read(1,*) rcloud
  read(1,*) mcloud
  read(1,*) Tgas
  close(1)

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types

  in_file       = trim(adjustl(in_file))
  in_file_form  = trim(adjustl(in_file_form))
  out_file      = trim(adjustl(out_file))
  out_file_form = trim(adjustl(out_file_form))
  restart       = .false.
  run_id        = ''
  rcentre(1:3)  = 0.0_PR

! Calculate scaling variables
  call units

! Scale input parameters
  rcloud = rcloud / real(rscale,PR)
  mcloud = mcloud / real(rscale,PR)

  write(6,*) "Radius of cloud   : ",rcloud*rscale,runit
  write(6,*) "Mass of cloud     : ",mcloud*mscale,munit

! Load in uniform density sphere file
  write(6,*) "in_file : ",trim(in_file)
  call read_data(in_file,in_file_form)

  pgas = ptot
  mp = mcloud / real(pgas,PR)

! Setting up for different particle types
  call set_default_particle_types
  call types

! Initialising kernel tables
#if defined(KERNEL_TABLES)
  call tabulate_kernel_functions
#endif

  do p=1,ptot
     sph(p)%r(1:NDIM) = sph(p)%r(1:NDIM)*rcloud
     sph(p)%m         = mp
     sph(p)%v(1:NDIM) = 0.0_PR
     sph(p)%temp      = Tgas
#if defined(INTERNAL_ENERGY)
     sph(p)%u         = Pconst*sph(p)%temp/(gamma - 1.0_PR)
#endif
  end do

! Initialise all other particles arrays
  call initialize_seren_variables_2

! Initialize certain variables before first force calculation
  call initialize_sph_variables_1

! Build and stock tree
  call tree_update(nbuild,nstock)

! Make initial guesses of h either using tree or a global estimate
#if defined(BH_TREE) && !defined(CONSTANT_H)
  if (.NOT. restart) call BHhydro_hguess
#elif !defined(CONSTANT_H)
  if (.NOT. restart) call h_guess
#endif

! Calculating initial SPH quantities
  call sph_update

! Calculate thermal properties
  call initialize_thermal_properties

! Write particle data to file
  out_file = trim(adjustl(out_file))
  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  call write_data_debug("ICSPITZER.debug.dat",rcentre(1:NDIM))
#endif

! Clean-up all memory
  call clean_up
  
  stop
END PROGRAM ic_spitzer
