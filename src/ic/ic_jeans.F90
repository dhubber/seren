! IC_JEANS.F90
! D. A. Hubber - 25/11/2007
! Generates initial conditions for the Jeans test, i.e. applies a 
! small sinusoidal density perturbation to a relaxed cube of unit volume.
! in_file         : Input file name
! in_file_form    : Input file format
! out_file        : Output file name
! out_file_form   : Output file format
! npert           : No. of wavelengths in period box
! amp             : Amplitude of density perturbation
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_jeans
  use particle_module
  use filename_module
  use type_module
  use tree_module
  use neighbour_module
  use periodic_module
  use time_module
  use scaling_module
  implicit none

  character(len=256) :: outfile        ! IC filename
  logical :: flag                      ! Flag that solution has been found
  integer :: gridtot                   ! Total number of grid points
  integer :: i                         ! Loop counter
  integer :: j                         ! Loop counter
  integer :: npert                     ! No. of pert. wavelengths in box
  integer :: p                         ! Particle counter
  real(kind=PR) :: amp                 ! Amplitude of sinusoidal perturbation
  real(kind=PR) :: dxgrid              ! Interpolation grid spacing
  real(kind=PR) :: L                   ! Periodic box length
  real(kind=PR) :: lambda              ! Wavelength of perturbation
  real(kind=PR) :: nL                  ! waves per unit length
  real(kind=PR) :: resol               ! The resolution, i.e. lambda/h
  real(kind=PR) :: x2                  ! New position of pertubed particle
  real(kind=PR), allocatable :: xp(:)  ! Auxilary particle positions

! Read in params file
  call default_parameters

! Initial conditions parameters
  write(6,*) "Reading parameters..."
  read(5,*) in_file
  read(5,*) in_file_form
  read(5,*) out_file
  read(5,*) out_file_form
  read(5,*) npert
  read(5,*) amp

  in_file       = trim(adjustl(in_file))
  in_file_form  = trim(adjustl(in_file_form))
  out_file      = trim(adjustl(out_file))
  out_file_form = trim(adjustl(out_file_form))

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types
   
! Check flags and parameters are okay
  call sanitycheck

! Reading in formatted data file (DRAGON format)
  call read_data(in_file,in_file_form)

! Setting up for different particle types
  pgas = ptot
  call types

! Calculate scaling units
  call units

! Calculating kernel tables 
#if defined(KERNEL_TABLES)
  call tabulate_kernel_functions
#endif

! Convert to code units
  call convert_to_code_units_1(.FALSE.)
  call convert_to_code_units_2

! Create Ewald correction table
#if defined(GRAVITY) && defined(EWALD)
  call ewald_init
#endif

! Compute important quantities of perturbation
  lambda = x_per / real(npert,PR)
  resol = npert*(6.0_PR*real(pp_gather,PR)/(PI*real(ptot,PR)))**(ONETHIRD)
  write (6,*) "Wavelength of perturbation : ",lambda
  write (6,*) "Resolution - ",resol
  write (6,*) "Amplitude of perturbation : ", amp
  
  gridtot = 2000
  allocate(xp(1:gridtot))
  L = x_per
  nL = real(npert,PR)/L
  dxgrid = L / real(gridtot - 1,PR)
  
! Search for solutions to the equation x = x' + amp*sin(??)
! First sets up table of values, then finds elements that contain the 
! solution, and finally linearly interpolates the correct value.
! ----------------------------------------------------------------------------
  do i=1,gridtot
     x2 = real(i - 1,PR) / real(gridtot - 1,PR)
     xp(i) = x2 + amp*sin(TWOPI*nL*x2)/nL/TWOPI
  end do
  
! Loop over all particles
! ----------------------------------------------------------------------------
  do p=1,ptot
     sph(p)%m = 1.0_PR / real(ptot,PR)
     flag = .false.
      
     ! Now search grid for solution
     ! -----------------------------------------------------------------------
     do i=1,gridtot-1
        if (xp(i+1) > sph(p)%r(1) .AND. xp(i) < sph(p)%r(1)) then
           x2 = real(i - 1,PR)*dxgrid + &
              &dxgrid*(sph(p)%r(1) - xp(i))/(xp(i+1) - xp(i))
           sph(p)%r(1) = x2
           flag = .true.
        end if
        if (flag) exit
     end do
     ! -----------------------------------------------------------------------
      
  end do
! ----------------------------------------------------------------------------

  deallocate(xp)

! Initialize certain variables before first force calculation
  call initialize_seren_variables_2
  call initialize_sph_variables_1

! Build and stock trees for first time
  call tree_update(nbuild,nstock)

! Make initial guesses of h either using tree or a global estimate
#if defined(BH_TREE) && !defined(CONSTANT_H)
  if (.not. restart) call BHhydro_hguess
#elif !defined(CONSTANT_H)
  if (.not. restart) call h_guess
#endif

! Calculating initial SPH quantities (h, neibs, rho, press etc.)
  call sph_update

! Calculate thermal properties
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

! Calculate initial diagnostic quantities
  call diagnostics

! Initialize other key variables (after initial force calculation)
  call initialize_variables_2

! Write to files
  call write_data(out_file,out_file_form)

! Clean-up memory
  call clean_up

  stop
END PROGRAM ic_jeans
