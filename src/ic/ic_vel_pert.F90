! IC_VEL_PERT.F90
! D. A. Hubber - 10/3/2008
! Loads file with spherically symmetric density distribution and adds 
! i) radial velocity perturbation, and/or
! ii) divergence-free turbulent velocity field
! Reads in a parameters file 'velpert.dat' with the following data : 
! in_file       : Input unit-sphere file name
! in_file_form  : Format of input sphere file 
!                 (e.g. dragon_form, dragon_unform, seren_form, seren_unform)
! out_file      : Output filename
! out_file_form : Format of output file
!                 (e.g. dragon_form, dragon_unform, seren_form, seren_unform)
! densmode      : Mode of density perturbation (not used yet)
! amp           : Amplitude of azimuthal perturbation (not used yet)
! mpert         : Azimuthal perturbation mode (not used yet)
! fenhance      : Density enhancement factor
! rstretch      : Triaxial stretching factor
! vpower        : Turbulent velocity power spectrum
! eturb         : Ratio of turb to grav energy
! ngrid         : No. of grid points per dimension for turb. vel. field
! rseed1        : Random No. seed 1 for generating turb. vel. field
! rseed2        : Random No. seed 2 for generating turb. vel. field
! velradmode    : Radial velocity mode;
!                 'energy' : Set value of radial kinetic energy / grav p.e.
!                 'dvdr'   : Set velocity gradient
! dvdr          : Radial velocity gradient
! erad          : Ratio of radial to grav energy
! velrotmode    : Rotational velocity field mode
!                 'energy' : Set value of rot. kinetic energy / grav p.e.
!                 'angvel' : Set value of angular velocity
! angmom        : Total angular momentum added to core
! angmomunit    : Angular momentum unit
! angpower      : Angular velocity power law (0 for solid-body rotation)
! angvel        : Angular velocity
! angvelunit    : Angular velocity unit
! erot          : Ratio of rot. to grav. energy
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_vel_pert
  use particle_module
  use filename_module
  use scaling_module
  use time_module
  use type_module
  implicit none

  character(len=50) :: densmode     ! Mode of density perturbation
  character(len=256) :: out_file    ! Name of output file
  character(len=50) :: velradmode   ! Radial velocity mode
  character(len=50) :: velrotmode   ! Rotational velocity mode
  integer :: iseed1                 ! Random number seed 1
  integer :: iseed2                 !   "      "     "   2
  integer :: i                      ! x-grid coordinate
  integer :: j                      ! y-grid coordinate
  integer :: k                      ! z-grid coordinate
  integer :: ngrid                  ! No. of grid points in each dimension
  integer :: p                      ! Particle counter
  real(kind=PR) :: amp              ! Amplitude of perturbation
  real(kind=PR) :: angmom           ! Angular momentum
  real(kind=PR) :: angpower         ! Angular momentum power law
  real(kind=PR) :: angvel           ! Angular velocity
  real(kind=PR) :: drmag            ! Distance
  real(kind=PR) :: drsqd            ! Distance squared
  real(kind=PR) :: dr_unit(1:NDIM)  ! Unit vector displacement
  real(kind=PR) :: dvdr             ! Dot product of dv and dr
  real(kind=PR) :: erad             ! Total radial kinetic energy
  real(kind=PR) :: erot             ! Total rotational kinetic energy
  real(kind=PR) :: eturb            ! Total turbulent kinetic energy
  real(kind=PR) :: fenhance         ! Density enhancement
  real(kind=DP) :: gpetot           ! Total gravitational potential energy
  real(kind=DP) :: keturb           ! Total kinetic energy
  real(kind=PR) :: mpert            ! Mode of perturbation
  real(kind=SP) :: mtot_temp             ! Total mass
  real(kind=PR) :: randnumb         ! Random number
  real(kind=PR) :: rcom_temp(1:NDIM)     ! Centre of mass position
  real(kind=PR) :: rp(1:NDIM)       ! Position of particle p
  real(kind=PR) :: rstretch(1:3)    ! Position scale factor for triaxial cores
  real(kind=DP) :: utot             ! Total internal energy
  real(kind=PR) :: vcom_temp(1:NDIM)     ! Velocity of centre of mass
  real(kind=PR) :: vmag             ! Speed
  real(kind=PR) :: vpower           ! Turbulent velocity power law
  real(kind=PR) :: vsqd             ! Speed squared

! Reading parameter file
  call default_parameters

! Read in data file for creating ICs.
  open(unit=1,file="velpert.dat",status="old")
  write(6,*) "------------------------------"
  write(6,*) "          ic_vel_pert         "
  write(6,*) "------------------------------"
  write(6,*) 
  write(6,*) "Reading in input parameters"
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) in_file
  read(1,*) in_file_form
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) densmode
  read(1,*) amp
  read(1,*) mpert
  read(1,*) fenhance
  read(1,*) rstretch(1),rstretch(2),rstretch(3)
  read(1,*); read(1,*); read(1,*)
  read(1,*) vpower
  read(1,*) eturb
  read(1,*) ngrid
  read(1,*) iseed1
  read(1,*) iseed2
  read(1,*); read(1,*); read(1,*)
  read(1,*) velradmode
  read(1,*) dvdr
  read(1,*) erad
  read(1,*); read(1,*); read(1,*)
  read(1,*) velrotmode
  read(1,*) angmom
  read(1,*) angmomunit
  read(1,*) angpower
  read(1,*) angvel
  read(1,*) angvelunit
  read(1,*) erot
  close(1)

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types

! Reading in formatted data file (DRAGON format)
  call read_data(in_file,in_file_form)

! Setting up for different particle types
  call types

! Initialising kernel tables
#if defined(KERNEL_TABLES)
  call tabulate_kernel_functions
#endif

! Setting up scaling units for simulation
  call units

! Read in opacity tables for radiation transport
#if defined(HYDRO) && defined(GRAVITY) && defined(RAD_WS)
  call read_eos
#endif

! Converting to dimensionless code units
  call convert_to_code_units_1(.FALSE.)
  call convert_to_code_units_2

! Density perturbations
  do p=1,ptot
     sph(p)%m = fenhance*sph(p)%m
     sph(p)%r(1:NDIM) = rstretch(1:NDIM)*sph(p)%r(1:NDIM)
  end do

! Initialize certain variables before first tree build and force calculation
  call initialize_seren_variables_2
  call initialize_sph_variables_1

! Build and stock trees for first time
  call tree_update(nbuild,nstock)

! Make initial guesses of h either using tree or a global estimate
#if defined(BH_TREE) && !defined(CONSTANT_H)
  if (.NOT. restart) call BHhydro_hguess
#elif !defined(CONSTANT_H)
  if (.NOT. restart) call h_guess
#endif

! Calculating initial SPH quantities (h, neibs, rho, etc.)
  call sph_update

! Initialize all thermal properties depending on options used
  call initialize_thermal_properties

! Calculate hydro forces on all SPH particles
#if defined(HYDRO)
  call sph_hydro_forces
#endif

! Calculate gravitational forces on all SPH particles
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  call sph_grav_forces
#endif

! Calculate gravitational forces on all sink particles
#if defined(SINKS) && defined(GRAVITY)
  call sph_sink_forces
#endif

! Calculate diagnostic properties
  call diagnostics

! First, calculate position and velocity of COM
  mtot_temp = 0.0_DP
  rcom_temp(1:NDIM) = 0.0_DP
  vcom_temp(1:NDIM) = 0.0_DP
  gpetot = 0.0_DP
  do p=1,ptot
     mtot_temp = mtot_temp + real(sph(p)%m,DP)
     rcom_temp(1:NDIM) = rcom_temp(1:NDIM) + real(sph(p)%m*sph(p)%r(1:NDIM),DP)
     vcom_temp(1:NDIM) = vcom_temp(1:NDIM) + real(sph(p)%m*sph(p)%v(1:NDIM),DP)
     gpetot = gpetot + real(sph(p)%m*sph(p)%gpot,DP)
  end do
  gpetot = 0.5_DP*gpetot
  rcom_temp(1:NDIM) = rcom_temp(1:NDIM) / mtot_temp
  vcom_temp(1:NDIM) = vcom_temp(1:NDIM) / mtot_temp

! Add turbulent velocity field to
  if (eturb > SMALL_NUMBER) then
     call add_turbulent_velocity_field(real(eturb*gpetot,PR),vpower,&
          &pgravitystart,ptot,ngrid,iseed1,iseed2)
  end if

! Add radial velocity field
  call add_radial_velocity_field(velradmode,real(erad*gpetot,PR),&
       &dvdr,1,ptot,real(rcom_temp(1:NDIM),PR))

! Add rotational velocity field
  call add_rotational_velocity_field(velrotmode,real(erot*gpetot,PR),&
       &angmom,angpower,angvel,1,ptot,real(rcom_temp(1:NDIM),PR))

! Calculate initial diagnostic quantities
#if defined(DEBUG_DIAGNOSTICS)
  call diagnostics
#endif

! Initialize other key variables (after initial force calculation)
  call initialize_sph_variables_2

! Write data to file
  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  rzero(1:NDIM) = 0.0_PR
  call write_data_debug("velpert.debug.dat",rzero(1:NDIM))
#endif

! Clean up memory
  call clean_up

  stop
END PROGRAM ic_vel_pert
