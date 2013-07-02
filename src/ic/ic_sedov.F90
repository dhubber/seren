! IC_SEDOV.F90
! D. A. Hubber & M. Kaplan - 16/06/2008
! Prepare initial conditions for Sedov blastwave test.
! in_file         : Input file name
! in_file_form    : Input file format
! out_file        : Output file name
! out_file_form   : Output file format
! rho0            : Cloud density
! radius          : Cloud radius
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_sedov
  use interface_module, only : read_data,write_data,w0
  use neighbour_module, only : pp_gather
  use particle_module
  use hydro_module
  use constant_module
  use filename_module
  use scaling_module
  use type_module
  use time_module
  use kernel_module
  implicit none

  character(len=256) :: out_file          ! Name of output file
  integer :: k                            ! Dimension counter
  integer :: kern                         ! Kernel table variable
  integer :: p                            ! Particle counter
  integer :: phot                         ! No. of particles in the center
  integer :: pcold                        ! rest of the particles
  integer :: pp_out                       ! rest of the particles II
  integer, allocatable :: hotlist(:)      ! potential neighbour list
  real(kind=PR) :: center(1:NDIM)         ! Center
  real(kind=PR) :: drmag                  ! Distance
  real(kind=PR) :: kefrac                 ! ..
  real(kind=PR) :: mp                     ! Mass of particle p
  real(kind=PR) :: mtot                   ! Total mass
  real(kind=PR) :: radius                 ! ..
  real(kind=PR) :: rho0                   ! Density of initial sphere
  real(kind=PR) :: r_hot                  ! Radius of hot particle region
  real(kind=PR) :: r_neigh                ! Radius for nearest neighbours in 
                                          ! the center (pp_center*rmax/ptot)
  real(kind=PR) :: rp(1:NDIM)             ! Position of particle p
  real(kind=PR) :: ufrac                  ! ..
  real(kind=PR) :: umax                   ! ..
  real(kind=PR) :: utot_temp              ! ..
  real(kind=PR), allocatable :: utemp(:)  ! internal energy

! Reading parameter file 
  call default_parameters

  write(6,*) "------------------------------"
  write(6,*) "         ic_sedovtest         "
  write(6,*) "------------------------------"
  open(unit=1,file="sedovparams.dat",status="old")
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) in_file
  read(1,*) in_file_form
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) rho0
  read(1,*) radius
  read(1,*) kefrac
  close(1)

! Check input parameters
  if (kefrac < 0.0_PR .or. kefrac > 1.0_PR) call paramerror("kefrac out of valid range")

  in_file       = trim(adjustl(in_file))
  in_file_form  = trim(adjustl(in_file_form))
  out_file      = trim(adjustl(out_file))
  out_file_form = trim(adjustl(out_file_form))

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types

! Setting up scaling units for simulation 
  call units

! Reading in formatted data file.  File should contain unit mass sphere of 
! radius 1, so no need to scale to dimensionless units. 
  call read_data(in_file,in_file_form)

! Allocate memory for list of 'hot' particles
  allocate(hotlist(1:ptot))

! Calculating kernel tables 
#if defined(KERNEL_TABLES)
  call tabulate_kernel_functions
#endif

! Calculate average value of smoothing length and gas cloud radius to know 
! where to select boundary particles.
  r_hot = radius*(pp_gather/real(ptot,PR))**(ONETHIRD)

! Calculate correct masses by scaling input parameters
  rho0  = rho0 / rhoscale
  mtot  = (4.0_PR*PI*rho0*radius**3)/3.0_PR
  mp    = mtot / real(ptot,PR)
  pcold = 0
  phot  = 0
  ufrac = max(1.0_PR - kefrac,0.0_PR)
  write(6,*) "ufrac : ",ufrac,"   kefrac : ",kefrac

! Set all particle type variables
  pgas      = ptot
  picm      = 0
  pboundary = 0
  pcdm      = 0
  call types
  write(6,*) "rho0 : ",rho0,"   mtot : ",mtot,"   mp : ",mp

! Set some basic particle properties
  do p=1,ptot
     sph(p)%v(1:VDIM) = 0.0_PR
     sph(p)%m         = mp
     sph(p)%r(1:NDIM) = sph(p)%r(1:NDIM)*radius
  end do

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


! Determine which particles are in the core
! ----------------------------------------------------------------------------
  umax = 0.0_PR
  utot_temp = 0.0_PR
  do p=1,ptot
     rp(1:NDIM) = sph(p)%r(1:NDIM)
     drmag = sqrt(rp(1)*rp(1) + rp(2)*rp(2) + rp(3)*rp(3))
     if (drmag < r_hot) then
        hotlist(p) = 1 
        phot       = phot + 1
        sph(p)%u   = w0(KERNRANGE*drmag/r_hot)*sph(p)%m
        utot_temp  = utot_temp + sph(p)%u
        umax       = max(umax,sph(p)%u)
     else
        hotlist(p) = 0
        pcold      = pcold + 1
     end if
  end do
  umax = umax / utot_temp


! Correction for the exact values of particle numbers
! ----------------------------------------------------------------------------
  do p=1,ptot
     if (hotlist(p) == 1) then
        rp(1:NDIM) = sph(p)%r(1:NDIM)
        drmag = sqrt(rp(1)*rp(1) + rp(2)*rp(2) + rp(3)*rp(3)) + SMALL_NUMBER
        rp(1:NDIM) = rp(1:NDIM) / drmag
        sph(p)%u    = max(sph(p)%u/utot_temp/sph(p)%m,1.E-6_PR*umax/sph(p)%m)
        sph(p)%temp = (gamma - 1.0_PR)*sph(p)%u/Pconst
        sph(p)%v(1:NDIM) = sqrt(2.0_PR*kefrac*sph(p)%u)*rp(1:NDIM)
        sph(p)%u    = max(ufrac*sph(p)%u,1.E-6_PR*umax/sph(p)%m)
     else
        sph(p)%u    = 1.E-6_PR*umax/sph(p)%m
        sph(p)%temp = (gamma - 1.0_PR)*sph(p)%u/Pconst
     end if
  end do

! Calculating initial forces/accelerations 
  call sph_hydro_forces

! Initialize other variables
  call initialize_sph_variables_2
  
! Write data to file
  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  rzero(1:3) = 0.0_PR
  call write_data_debug("ICSEDOV.debug.dat",rzero)
#endif

  write(6,*) "ptot        : ",ptot
  write(6,*) "pp_gather   : ",pp_gather
  write(6,*) "mtot        : ",mtot,mp
  write(6,*) "r_hot       : ",r_hot
  write(6,*) "phot, pcold : ",phot, pcold

  call diagnostics

! Clean up memory
  deallocate(hotlist)
  call clean_up  

  stop
END PROGRAM ic_sedov
