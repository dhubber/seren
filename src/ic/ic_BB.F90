! IC_BB.F90
! D. A. Hubber - 08/02/2010
! Prepares initial conditions for the Boss-Bodenheimer (1979) test. 
! First, read in file containing uniform density sphere of unit radius. 
! Next, add density perturbation in spherical polar coordinates and scale 
! to required units.  Finally, write to output file.  
! in_file        : Input filename
! in_file_form   : Input file format
! out_file       : Output filename
! out_file_form  : Output file format
! mass           : Mass of cloud
! munit          : Mass unit
! rcloud         : Radius of cloud
! runit          : Length unit
! temp_cloud     : Temperature of cloud
! angvel         : Angular velocity of cloud
! angvelunit     : Angular velocity unit
! mpert          : Order of azimuthal perturbation
! amp            : Amplitude of density perturbation
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_BB
  use interface_module, only : comperror,paramerror,read_data,write_data
  use particle_module
  use hydro_module
  use filename_module
  use scaling_module
  implicit none

  character(len=256) :: out_file   ! Name of output file
  character(len=50) :: pertmode    ! Perturbation mode
  character(len=50) :: velmode     ! Velocity mode
  integer :: p                     ! Particle counter
  real(kind=PR) :: amp             ! Amplitude of sinusoidal perturbation
  real(kind=PR) :: angvel          ! Initial angular velocity of cloud
  real(kind=PR) :: mpert           ! Order of perturbation
  real(kind=PR) :: mass            ! Mass of cloud
  real(kind=PR) :: mp              ! Mass of particle p
  real(kind=PR) :: r0(1:NDIM)      ! ..
  real(kind=PR) :: rcloud          ! Radius of cloud
  real(kind=PR) :: temp_cloud      ! (Isothermal) temperature of cloud

! Set default values for all parameters
  call default_parameters

  write(6,*) "------------------------------"
  write(6,*) "            ic_BB             "
  write(6,*) "------------------------------"
  write(6,*) 
  write(6,*) "Input data..."

  open(1, file="BBparams.dat", status="old", form="formatted")
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) in_file
  read(1,*) in_file_form
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) mass
  read(1,*) munit
  read(1,*) rcloud
  read(1,*) runit
  read(1,*) temp_cloud
  read(1,*) angvel
  read(1,*) angvelunit
  read(1,*) mpert
  read(1,*) amp
  close(1)

  if (mass <= 0.0_PR) call paramerror("Invalid cloud mass")
  if (rcloud <= 0.0_PR) call paramerror("Invalid cloud radius")
  if (temp_cloud <= 0.0_PR) call paramerror("Invalid cloud temperature")
  if (mpert <= 0) call paramerror("Zero/negative azimuthal perturbation order")
#if NDIM != 3
  call comperror("NDIM = 3 not selected")
#endif

! Set scaling variables that aren't specified in parameters file 
! (mainly redundant) 
  rscale     = 1.0_PR
  mscale     = 1.0_PR
  tunit      = "myr"
  vunit      = "km_s"
  aunit      = "km_s2"
  rhounit    = "m_sun_pc3"
  Punit      = "g_cms2"
  Eunit      = "GJ"
  momunit    = "m_sunkm_s"
  mu_bar     = 2.35_PR
  r0(1:NDIM) = 0.0_PR
  pertmode   = "azimuthal"
  velmode    = "angvel"
  com_frame  = .true.
  restart    = .false.

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types

! Calculate all over scaling variables
  call units

! Scale all input parameters
  angvel = angvel / angvelscale
  mass   = mass   / mscale
  rcloud = rcloud / rscale

  write(6,*) "rcloud : ",rcloud*rscale
  write(6,*) "mass   : ",mass*mscale
  write(6,*) "angvel : ",angvel*angvelscale
  write(6,*) "rho0   : ",mass*3.0_PR/(rcloud**3)/4.0_PR/PI
  write(6,*) "rhobary : ",rhobary

! Read in relaxed, uniform density sphere
  call read_data(in_file,in_file_form)

! Setting up for different particle types
  call types

! Calculate mass of individual particles (equal masses)
  mp = mass / real(ptot,PR)

! Set positions, masses and temperatures of unperturbed cloud
  do p=1,ptot
     sph(p)%r(1:NDIM) = sph(p)%r(1:NDIM)*rcloud
     sph(p)%m         = mp
     sph(p)%temp      = temp_cloud
     sph(p)%v(1:VDIM) = 0.0_PR
  end do

! Calculate COM of system, and if required, change to COM
  call COM

! Perturb positions of particles in cloud
  call add_azimuthal_density_perturbation(pertmode,2,&
       &1,ptot,amp,0.0_PR,r0(1:NDIM))

! Add solid body rotation to cloud
  call add_rotational_velocity_field(velmode,0.0_PR,&
       &0.0_PR,0.0_PR,angvel,1,ptot,r0(1:NDIM))

! Write particle data to file
  call write_data(out_file,out_file_form)

! Deallocate memory
  call clean_up

  stop
END PROGRAM ic_BB
