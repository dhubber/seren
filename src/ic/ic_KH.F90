! IC_KH.F90
! D. A. Hubber - 15/3/2008
! Creates initial conditions for the Kelvin-Helmholtz instability test.
! Assumes we either load in unity sheets/cubes.
! out_file        : Output file name
! out_file_form   : Output file format
! in_file1        : 1st input file name
! in_file1_form   : 1st input file format
! in_file2        : 2nd input file name
! in_file2_form   : 2nd input file format
! p1, p2          : No. of particles in file 1, 2
! n1, n2          : No. of replicas for layers 1, 2
! vx1, vx2        : x-velocity of layer 1, 2
! rho1, rho2      : Density of layer 1, 2
! Press1, Press2  : Pressure of layer 1, 2
! x1, x2          : x-size of input file 1, 2
! y1, y2          : y-size of input file 1, 2
! z1, z2          : z-size of input file 1, 2 (redundant?)
! h_fac           : 'grad-h' h_fac
! pertmode        : Mode of perturbation
! amp             : Amplitude of perturbation
! lambda          : Wavelength of perturbation
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_KH
  use interface_module, only : check_boundary_conditions,write_data
  use particle_module
  use hydro_module
  use filename_module
  use time_module
  use periodic_module
  use type_module
  use neighbour_module, only : h_fac
  implicit none

  character(len=256) :: out_file         ! KH IC file
  character(len=50) :: pertmode          ! Perturbation type
  integer :: i                           ! Aux. counter
  integer :: k                           ! Dimension counter
  integer :: p                           ! Particle counter
  integer :: p1                          ! No. of particles in ??
  integer :: p2                          ! No. of particles in ??
  integer :: niterations                 ! No. of iterations for OSPH
  integer :: nx1                         ! No. of replicated cubes in ??
  integer :: nx2                         ! No. of replicated cubes in ??
  integer :: ny1                         ! No. of replicated cubes in ??
  integer :: ny2                         ! No. of replicated cubes in ??
  real(kind=PR) :: amp                   ! Amplitude of velocity perturbation
  real(kind=PR) :: kpert
  real(kind=PR) :: lambda                ! Wavelength of perturbation
  real(kind=PR) :: offset                ! ..
  real(kind=PR) :: mp1                   ! Mass of 'bottom' particles
  real(kind=PR) :: mp2                   ! Mass of 'top' particles
  real(kind=PR) :: Press1                ! 'Bottom' pressure
  real(kind=PR) :: Press2                ! 'Top' pressure
  real(kind=PR) :: rho1                  ! 'Bottom' density
  real(kind=PR) :: rho2                  ! 'Top' density
  real(kind=PR) :: sigmapert             ! ..
  real(kind=PR) :: T1                    ! 'Bottom' temperature
  real(kind=PR) :: T2                    ! 'Top' temperature
  real(kind=PR) :: volume                ! Volume of ..
  real(kind=PR) :: vx1                   ! 'Bottom' x-velocity 
  real(kind=PR) :: vx2                   ! 'Top' x-velocity
  real(kind=PR), allocatable :: r1(:,:)  ! Positions of 'bottom' particles
  real(kind=PR), allocatable :: r2(:,:)  ! Positions of 'top' particles
#if defined(DEBUG_PLOT_DATA)
  real(kind=PR) :: rcentre(1:NDIM)       ! Origin
#endif

! Set parameters to default values
  call default_parameters

! Read in parameter file 
  write(6,*) "Opening KHparams file"
  open(unit=1,file="KHparams.dat",status='old')
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*) nx1, nx2
  read(1,*) ny1, ny2
  read(1,*) vx1, vx2
  read(1,*) rho1, rho2
  read(1,*) Press1, Press2
  read(1,*) periodic_min(1), periodic_max(1)
  read(1,*) periodic_min(2), periodic_max(2)
  read(1,*) h_fac
  read(1,*) gamma
  read(1,*) pertmode
  read(1,*) amp
  read(1,*) lambda
  close(1)

  out_file      = trim(adjustl(out_file))
  out_file_form = trim(adjustl(out_file_form))

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types

! Checking compiler flags and parameter values
  call sanitycheck

! Setting up scaling units for simulation
  call units

  p1 = nx1*ny1
  p2 = nx2*ny2
  allocate(r1(1:NDIM,1:p1))
  allocate(r2(1:NDIM,1:p2))

  pgas = p1 + p2
  ptot = pgas


! Calculate temperature of gas
  T1 = Press1 / rho1
  T2 = Press2 / rho2
  write(6,*) "T1 :",T1
  write(6,*) "T2 :",T2

  volume = 0.5_PR*periodic_size(1)*periodic_size(2)
  mp1    = rho1*volume/real(p1,PR)
  mp2    = rho2*volume/real(p2,PR)
  kpert  = 2.0_PR*PI/lambda

! Output parameters to screen for verification
  write(6,*) "outfile :",trim(out_file),"   ",trim(out_file_form)
  write(6,*) "Number of particles :",p1,p2
  write(6,*) "vx1 :",vx1,"    vx2 ",vx2
  write(6,*) "rho1 :",rho1,"    rho2 :",rho2
  write(6,*) "Press1 :",Press1,"    Press2 :",Press2
  write(6,*) "amp :",amp,"  lambda :",lambda
  write(6,*) "Closing KHparams.dat file"
  write(6,*) "mp1 : ",mp1,"    mp2 : ",mp2

  call allocate_memory(.FALSE.)

  call add_uniform_2d_lattice(periodic_min(1),periodic_max(1),&
       &periodic_min(2),periodic_min(2) + periodic_half(2),r1,p1,nx1,ny1)
  call add_uniform_2d_lattice(periodic_min(1),periodic_max(1),&
       &periodic_min(2) + periodic_half(2),periodic_max(2),r2,p2,nx2,ny2)


! First write LHS of shock tube info 
! ----------------------------------------------------------------------------
  do p=1,p1
     sph(p)%r(1:NDIM) = r1(1:NDIM,p)
     sph(p)%r(2)      = sph(p)%r(2) + 0.5_PR*periodic_half(2)
     sph(p)%m         = mp1
     sph(p)%h         = 1.0_PR
     sph(p)%v(1:NDIM) = 0.0_PR
     sph(p)%v(1)      = vx1
     sph(p)%temp      = T1
#if defined(INTERNAL_ENERGY)
     sph(p)%u         = Pconst*T1/(gamma - 1.0_PR)
#endif
#if defined(ENTROPIC_FUNCTION)
     sph(p)%Aent      = Press1/rho1**gamma
#endif
     call check_boundary_conditions(sph(p)%r(1:NDIM),sph(p)%v(1:VDIM))
  end do


! Now write RHS of shock tube info 
! ----------------------------------------------------------------------------
  do i=1,p2
     p = p1 + i
     sph(p)%r(1:NDIM) = r2(1:NDIM,i)
     sph(p)%r(2)      = sph(p)%r(2) + 0.5_PR*periodic_half(2)
     sph(p)%m         = mp2
     sph(p)%h         = 1.0_PR
     sph(p)%v(1:NDIM) = 0.0_PR
     sph(p)%v(1)      = vx2
    sph(p)%temp       = T2
#if defined(INTERNAL_ENERGY)
     sph(p)%u         = Pconst*T2/(gamma - 1.0_PR)
#endif
#if defined(ENTROPIC_FUNCTION)
     sph(p)%Aent      = Press1/rho2**gamma
#endif
     call check_boundary_conditions(sph(p)%r(1:NDIM),sph(p)%v(1:VDIM))
  end do


! Add velocity perturbation here
! ----------------------------------------------------------------------------
  if (pertmode == 'price2008') then
     do p=1,ptot
        if (sph(p)%r(2) > periodic_max(2)) &
             &sph(p)%r(2) = sph(p)%r(2) - periodic_size(2)
        if (abs(sph(p)%r(2) - 0.25_PR) < 0.025_PR) then
           sph(p)%v(2) = amp*sin(-2.0_PR*PI*(sph(p)%r(1) + 0.5_PR)/lambda)
        else if (abs(sph(p)%r(2) + 0.25_PR) < 0.025_PR) then
           sph(p)%v(2) = amp*sin(2.0_PR*PI*(sph(p)%r(1) + 0.5_PR)/lambda)
        end if
     end do
  else if (pertmode == 'arepo') then
     sigmapert = 0.05_PR/sqrt(2.0_PR)
     do p=1,ptot
        sph(p)%v(2) = amp*sin(2.0_PR*PI*sph(p)%r(1)/lambda)*&
             &(exp(-(sph(p)%r(2) + 0.25_PR)**2/2.0_PR/sigmapert**2) +& 
             &exp(-(sph(p)%r(2) - 0.25_PR)**2/2.0_PR/sigmapert**2))
     end do
  else
     stop 'No valid velocity perturbation selected'
  end if

! Setting up for different particle types
  call types
  write(6,*) "pgas :",pgas
  write(6,*) "ptot :",ptot

! Initialising kernel tables
#if defined(KERNEL_TABLES)
  call tabulate_kernel_functions
#endif

! Initialize variables
  call initialize_seren_variables_2
  call initialize_sph_variables_1

! Build and stock tree
  call tree_update(nbuild,nstock)

! Make initial guesses of h either using tree or a global estimate
#if defined(BH_TREE) && !defined(CONSTANT_H)
  if (.not. restart) call BHhydro_hguess
#elif !defined(CONSTANT_H)
  if (.not. restart) call h_guess
#endif

#if defined(OSPH)
  niterations = 10
#else
  niterations = 1
#endif

! Iterate density-pressure relation for OSPH
! ----------------------------------------------------------------------------
  do i=1,niterations

     ! Calculating initial SPH quantities
     call sph_update
     
     ! Calculate temperature required to ensure constant pressure.
     do p=1,ptot
        sph(p)%temp  = Press1 / sph(p)%rho / Pconst
        sph(p)%u     = Pconst*sph(p)%temp/(gamma - 1.0_PR)
        sph(p)%press = (gamma - 1.0_PR)*sph(p)%rho*sph(p)%u
        sph(p)%sound = sqrt(sph(p)%press/sph(p)%rho)
#if defined(ENTROPIC_FUNCTION)
        sph(p)%Aent  = sph(p)%press/sph(p)%rho**gamma
#endif
     end do

  end do
! ----------------------------------------------------------------------------

! Calculating initial SPH quantities
  call sph_update
  call update_thermal_properties

! Calculate hydro forces on all SPH particles
  call sph_hydro_forces

! Initialize other variables
  call initialize_sph_variables_2

! Write everything to file 
  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  rcentre(1:NDIM) = 0.0_PR
  call write_data_debug("ICKH.debug.dat",rcentre(1:NDIM))
#endif

  deallocate(r1)
  deallocate(r2)

  call clean_up

  stop
END PROGRAM ic_KH
