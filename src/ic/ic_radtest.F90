! IC_RADTEST.F90
! D. A. Hubber - 8/1/2008
! Creates initial conditions for radiative transfer collapse test. 
! Sets up uniform density sphere of specified density and radius 
! surrounded by a layer of boundary particles.
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_radtest
  use particle_module
  use hydro_module
  use constant_module
  use filename_module
  use scaling_module
  use type_module
  use time_module
  use neighbour_module, only : pp_gather
  use Tprof_module, only : temp_inf
  implicit none

  character(len=256) :: out_file    ! Name of output file
  integer :: k                      ! Dimension counter
  integer :: p                      ! Particle counter
  real(kind=PR) :: drmag            ! Distance
  real(kind=PR) :: hguess           ! Estimate of h assuming uniform density
  real(kind=PR) :: mp               ! Mass of particle p
  real(kind=PR) :: mcloud             ! Total mass
  real(kind=PR) :: radiusgas        ! Radius of cloud with gas particles
  real(kind=PR) :: rmax             ! Radius of gas cloud (just gas particles)
  real(kind=PR) :: rp(1:NDIM)       ! Position of particle p
  real(kind=PR) :: rsqd             ! Distance squared
  real(kind=PR), allocatable :: rboundary(:,:)  ! Positions of boundary part.
  real(kind=PR), allocatable :: rgas(:,:)       ! Positions of gas particles

! Set default parameters
  call default_parameters

  write(6,*) "------------------------------"
  write(6,*) "         ic_radtest           "
  write(6,*) "------------------------------"
  open(unit=1,file="radparams.dat",status="old")
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) in_file
  read(1,*) in_file_form
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) mcloud
  read(1,*) munit
  read(1,*) rmax
  read(1,*) runit
  read(1,*) gamma
  close(1)

! Setting up scaling units for simulation 
  call units

! Reading in formatted data file (DRAGON snapshot).  File should contain 
! unit mass sphere of radius 1, so no need to scale to dimensionless units. 
  call read_data(in_file,in_file_form)

! Calculating kernel tables 
#if defined(KERNEL_TABLES)
  call tabulate_kernel_functions
#endif

! Calculate average value of smoothing length and gas cloud radius to know 
! where to select boundary particles.
  hguess = (pp_gather/8.0_PR/real(ptot,PR))**(ONETHIRD)
  radiusgas = 1.0_PR - KERNRANGE*hguess

  pboundary = 0
  pgas      = 0
  picm      = 0
  allocate(rboundary(1:NDIM,1:ptot))
  allocate(rgas(1:NDIM,1:ptot))

! Now sort particles into gas particles and boundary particles
! ----------------------------------------------------------------------------
  do p=1,ptot
     rp(1:NDIM) = sph(p)%r(1:NDIM)
     drmag = sqrt(rp(1)*rp(1) + rp(2)*rp(2) + rp(3)*rp(3))
     if (drmag > radiusgas) then
        pboundary = pboundary + 1
        rboundary(1:NDIM,pboundary) = rp(1:NDIM)
     else
        pgas = pgas + 1
        rgas(1:NDIM,pgas) = rp(1:NDIM)
     end if
  end do
! ----------------------------------------------------------------------------

  write(6,*) "pgas :",pgas
  write(6,*) "pboundary :",pboundary
  write(6,*) "pgas + pboundary :",pgas+pboundary
  write(6,*) "ptot :",ptot

! Now reorder particle positions for boundary and gas particles
  do p=1,pboundary
     sph(p)%r(1:NDIM) = rboundary(1:NDIM,p)
  end do
  do p=1,pgas
     parray(1:NDIM,p+pboundary) = rgas(1:NDIM,p)
  end do
  call types

! Calculate correct masses by scaling input parameters
  mcloud = mcloud / mscale
  rmax = rmax / rscale
  mp   = mcloud / real(pgas,PR)

  write(6,*) "mp :",mp
  write(6,*) "radiusgas :",radiusgas
  write(6,*) "rmax :",rmax

! Set attributes for all particles
  do p=1,ptot
     sph(p)%r(1:NDIM) = sph(p)%r(1:NDIM)*rmax/radiusgas
     sph(p)%m   = mp
     sph(p)%h   = 0.0_PR
     sph(p)%v(1:NDIM)      = 0.0_PR
     sph(p)%rho           = 3.0_PR*mcloud / (4.0_PR*PI*rmax**3)
    sph(p)%temp          = temp_inf
#if defined(INTERNAL_ENERGY)
     sph(p)%u             = Pconst*temp(p)/(gamma - 1.0_PR)
#endif
  end do

! Zero time variables
  nsteps = 0
  n = 0

! Write data to file
  call write_data(out_file,out_file_form)

! Clean up memory
  call clean_up

  deallocate(rgas)
  deallocate(rboundary)

  STOP
END PROGRAM ic_radtest

