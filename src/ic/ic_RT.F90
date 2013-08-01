! IC_RT.F90
! D. A. Hubber - 26/05/2009
! Generate initial conditions to model Rayleigh-Taylor instability.  
! Creates two gas lattices of different densities in pressure balance 
! at the boundary.  If the dense gas is on the top, then the configuration 
! is unstable to the Rayleigh-Taylor instability.  We seed a perturbation 
! at the interface following one of two possible methods
! 1) Velocity perturbation
! 2) Boundary displacement perturbation
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_RT
  use interface_module, only : read_data,write_data
  use particle_module
  use hydro_module
  use filename_module
  use time_module
  use periodic_module
  use type_module
  use scaling_module
  use neighbour_module
  implicit none

  character(len=256) :: out_file      ! Name of output file
  character(len=256) :: in_file1      ! Name of output file
  character(len=256) :: in_file1_form ! Name of output file
  integer :: i                        ! x-dimension counter
  integer :: j                        ! y-dimension counter
  integer :: k                        ! z-dimension counter
  integer :: nlayers1                 ! No. of layers of medium 1 in y
  integer :: nlayers2                 ! No. of layers of medium 2 in y
  integer :: nwall1                   ! ..
  integer :: nwall2                   ! ..
  integer :: p                        ! Particle counter
  integer :: p1                       ! No. of particles in medium 1
  integer :: p2                       ! No. of particles in medium 2
  integer :: pertmode                 ! Mode of adding perturbation
  integer :: ppd1                     ! Particle per dimension
  integer :: ppd2                     ! Particle per dimension
  real(kind=PR) :: acc_grav           ! External grav. accel. in y direction
  real(kind=PR) :: amp                ! Amplitude of velocity perturbation
  real(kind=PR) :: lambda             ! Wavelength of perturbation
  real(kind=PR) :: mp1                ! Mass of 'bottom' particles
  real(kind=PR) :: mp2                ! Mass of 'top' particles
  real(kind=PR) :: Press1             ! Constant pressure
  real(kind=PR) :: rho1               ! 'Bottom' density
  real(kind=PR) :: rho2               ! 'Top' density
  real(kind=PR) :: x1                 ! x-size of region 1
  real(kind=PR) :: x2                 ! x-size of region 2
  real(kind=PR) :: xsize              ! Inputted size of domain
  real(kind=PR) :: y1                 ! y-size of region 1
  real(kind=PR) :: y2                 ! y-size of region 2
  real(kind=PR) :: ywall1             ! ..
  real(kind=PR) :: ywall2             ! ..

#if NDIM==1 || NDIM==3
  write(6,*) "Compiler flag error : Only works in 2D"
  stop
#endif

! Reading parameter file
  call default_parameters

  write(6,*) "------------------------------"
  write(6,*) "            ic_RT             "
  write(6,*) "------------------------------"
  open(unit=1,file="RTparams.dat",status='old')
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) pertmode
  read(1,*) ppd1,ppd2
  read(1,*) nlayers1,nlayers2
  read(1,*) nwall1,nwall2
  read(1,*) rho1,rho2
  read(1,*) Press1
  read(1,*) acc_grav
  read(1,*) gamma
  read(1,*) xsize
  read(1,*) amp
  read(1,*) lambda
  read(1,*) pp_gather
  read(1,*) hmin
  read(1,*) h_fac

  p1   = ppd1*(nlayers1 + nwall1)
  p2   = ppd2*(nlayers2 + nwall2)
  ptot = p1 + p2
  write(6,*) "Number of particles : ",ptot
  write(6,*) "p1 :",p1,"   p2 :",p2,p1+p2

  x1 = xsize
  x2 = xsize
  y1 = xsize*real(nlayers1,PR)/real(ppd1,PR)
  y2 = xsize*real(nlayers2,PR)/real(ppd2,PR)
  ywall1 = xsize*real(nwall1,PR)/real(ppd1,PR)
  ywall2 = xsize*real(nwall2,PR)/real(ppd2,PR)

! Set periodic variables
  periodic_min(1)    = 0.0_PR
  periodic_max(1)    = xsize
  periodic_min(2)    = 0.0_PR - ywall1
  periodic_max(2)    = y1 + y2 + ywall2
  periodic_min(3)    = 0.0_PR
  periodic_max(3)    = 0.0_PR
  periodic_size(1:3) = periodic_max(1:3) - periodic_min(1:3)
  periodic_half(1:3) = 0.5_PR*periodic_size(1:3)

  write(6,*) periodic_min(1:3)
  write(6,*) periodic_max(1:3)
  write(6,*) periodic_size(1:3)
  write(6,*) x1,x2,y1,y2

! Initialise some variables using parameters
  call initialize_seren_variables_1

! Checking compiler flags and parameter values
  call sanitycheck

  rscale = 1.0_DP
  mscale = 1.0_DP

! Set masses of SPH particles in both regions
  mp1 = rho1*x1*(y1 + ywall1)/real(p1,PR)
  mp2 = rho2*x2*(y2 + ywall2)/real(p2,PR)

  sph(1:ptot)%ptype = GASID
  pgas = ptot


! Perturbation mode 1 (Velocity perturbation)
! ============================================================================
  if (pertmode == 1) then

     ! Allocate memory for particle data
     call allocate_memory(.FALSE.)

     ! Bottom layer of particles
     ! -----------------------------------------------------------------------
     p = 0
     do j=1,nlayers1 + nwall1
        do i=1,ppd1
           p = p + 1
           sph(p)%r(1) = xsize*(real(i,PR)-0.5_PR)/real(ppd1) + periodic_min(1)
           sph(p)%r(2) = xsize*(real(j)-0.5_PR)/real(ppd1) + periodic_min(2)
           sph(p)%m = mp1
           sph(p)%h = 1.0_PR
           sph(p)%v(1:VDIM)    = 0.0_PR
           if (j <= nwall1) sph(p)%ptype = BOUNDARYID
        end do
     end do
     
     ! Top layer of particles
     ! -----------------------------------------------------------------------
     do j=1,nlayers2 + nwall2
        do i=1,ppd2
           p = p + 1
           sph(p)%r(1) = xsize*(real(i,PR)-0.5_PR)/real(ppd2) + periodic_min(1)
           sph(p)%r(2) = xsize*(real(j,PR)-0.5_PR)/real(ppd2) &
                & + periodic_min(2) + y1 + ywall1
           sph(p)%m = mp2
           sph(p)%h = 1.0_PR
           sph(p)%v(1:VDIM)    = 0.0_PR
           if (j > nlayers2) sph(p)%ptype = BOUNDARYID
        end do
     end do
     
     ! Set velocity perturbation
     ! -----------------------------------------------------------------------
     do p=1,ptot
        
        ! Abel one
        !if (sph(p)%r(2) > 0.3_PR .and. sph(p)%r(2) < 0.7_PR) then
        !   v(2,p) = amp*(1.0_PR + cos(8.0_PR*PI*(sph(p)%r(1) + 0.25)))*&
        !        &(1.0_PR + cos(2.0_PR*PI/0.4_PR*(sph(p)%r(2) - 0.5)))/4.0_PR
        !end if

        ! Springel one
        sph(p)%v(2) = amp*(1.0_PR - cos(2.0_PR*PI*(sph(p)%r(1))/lambda))*&
             &(1.0_PR - cos(2.0_PR*PI*(sph(p)%r(2))/(y1 + y2)))

        ! Simple one
        ! if (abs(sph(p)%r(2) - y1) < 0.025) then
        !    v(2,p) = -amp*cos(2.*PI*(sph(p)%r(1))/lambda)
        ! end if
     end do



! Perturbation mode 2 (Boundary perturbation)
! ============================================================================
  else if (pertmode == 2) then

     ! Allocate memory for particle data
     call allocate_memory(.FALSE.)

     ! Set properties both sides of density boundary
     ! -----------------------------------------------------------------------
     p = 0
     do j=1,(nlayers1 + nlayers2)
        do i=1,ppd1
           p = p + 1
           sph(p)%r(1) = xsize*(real(i,PR) - 0.5_PR)/real(ppd1) - 0.5*xsize
           sph(p)%r(2) = xsize*(real(j,PR) - 0.5_PR)/real(ppd1)
           sph(p)%h = 1.0_PR
           sph(p)%v(1:VDIM) = 0.0_PR
           if (sph(p)%r(2) - y1 > -amp*cos(2.0_PR*PI*(sph(p)%r(1))/lambda)) then
              sph(p)%m = mp2
           else
              sph(p)%m = mp1
           end if
        end do
     end do 


! Load from files
! ============================================================================
  else if (pertmode == 3) then

     ! Reading in formatted data file (DRAGON format)
     call read_data(in_file1,in_file1_form)
     
     ! Set velocity perturbation
     ! -----------------------------------------------------------------------
     do p=1,ptot
        sph(p)%v(2) = amp*(1.0_PR - cos(2.0_PR*PI*(sph(p)%r(1))/lambda))*&
             &(1.0_PR - cos(2.0_PR*PI*(sph(p)%r(2))/(y1 + y2)))
     end do
 
  end if
! ============================================================================

! Sort particles in correct type order in memory
  call sort_particle_types(sph(1:ptot)%ptype)
  call types

! Setting up scaling units for simulation
  call units

! Converting to dimensionless code units
  call convert_to_code_units_1(.FALSE.)
  call convert_to_code_units_2

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

! Calculating initial SPH quantities
  call sph_update

! Calculate temperature required to ensure constant pressure
  do p=1,ptot
     sph(p)%press = Press1 + acc_grav*(sph(p)%r(2) - y1)*sph(p)%rho
     sph(p)%temp  = sph(p)%press / sph(p)%rho / Pconst
     sph(p)%u     = Pconst*sph(p)%temp/(gamma - 1.0_PR)
  end do

! Calculating initial forces/accelerations 
  call sph_hydro_forces

! Initialize other variables
  call initialize_sph_variables_2

! Write everything to file 
  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  rzero(1:NDIM) = 0.0_PR
  call write_data_debug("ICRT.debug.dat",rzero)
#endif

! Clean up memory
  call clean_up

  stop
END PROGRAM ic_RT
