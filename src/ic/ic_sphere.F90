! IC_SPHERE.F90
! D. A. Hubber - 12/9/2007
! Creates uniform density sphere from a relaxed 3-D cube of particles. 
! Assumes the sides of the cubes are at 0->1 in all dimensions.
! If the file containing the cube has too few particles (such that the 
! sphere would be culled), the program stops and does not produce any files.
! in_file        : Input file name
! in_file_form   : Input file format
! rcloud         : Radius of cloud (After re-scaling)
! nwant          : Required no. of particles
! out_file       : Output file name
! out_file_form  : Output file format
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_sphere
  use interface_module, only : distance3,read_data,write_data
  use particle_module
  use hydro_module
  use filename_module
  use scaling_module
  use type_module
  use time_module
  implicit none

  character(len=256) :: out_file              ! Name of outfile
  integer :: i                                ! Auxilary counter
  integer :: itmax                            ! ..
  integer :: old_ptot                         ! Original number of particles
  integer :: niterations                      ! Required number of particles
  integer :: nwant                            ! ..
  integer :: p                                ! Particle counter
  integer, allocatable :: list(:)             ! List of particles in sphere
  real(kind=PR) :: dr(1:NDIM)                 ! Relative position vector
  real(kind=PR) :: drsqd                      ! Distance squared
  real(kind=PR) :: r0(1:NDIM)                 ! Sphere centre
  real(kind=PR) :: rad                        ! Cloud scaling radius
  real(kind=PR) :: rcloud                     ! (Dimensionless) cloud radius
  real(kind=PR) :: rlow                       ! Min rad for bisection method
  real(kind=PR) :: rhigh                      ! Max rad for bisection method
  real(kind=PR), allocatable :: position(:,:) ! Temp array for part. pos.

! Set parameters to default values
  call default_parameters

  write(6,*) "------------------------------"
  write(6,*) "          ic_sphere           "
  write(6,*) "------------------------------"
  write(6,*) 
  write(6,*) "Input data..."
  open(unit=1,file='sphereparams.dat',status='old')
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) in_file
  read(1,*) in_file_form
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) rcloud
  read(1,*) nwant
  close(1)

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types

! Calculate scaling units
  call units

! Read relaxed cube file and calculate units
  call read_data(in_file,in_file_form)

! Record original particle positions in temporary array.  Also transform 
! to cube of volume 2**3 and place centre of cube at origin.
  allocate (position(1:NDIM,1:ptot))
  allocate (list(1:ptot))
  do p=1,ptot
     position(1:NDIM,p) = 2.0_PR*sph(p)%r(1:NDIM) - 1.0_PR
  end do

! Now, reallocate arrays for new replicas
  old_ptot = ptot
  call clean_up
  
! Set-up bisection method limits
  rhigh       = 2.0_PR
  rlow        = 0.0_PR
  r0(1:NDIM)  = 0.0_PR
  niterations = 0
  itmax       = 200
   
  write(6,*) "Preparing bisection method"

! Find radius which contains nwant particles using a bisection method
! ----------------------------------------------------------------------------
  do
     rad = 0.5_PR*(rhigh + rlow)
     ptot = 0
     do p=1,old_ptot
        call distance3(r0(1:NDIM),position(1:NDIM,p),dr(1:NDIM),drsqd)
        if (drsqd <= rad*rad) then 
           ptot = ptot + 1
           list(ptot) = p
        end if
     end do
     write (6,*) "rad :",rad,"   ptot :",ptot
     niterations = niterations + 1
     if (niterations == itmax) exit
     if (ptot == nwant) exit
     if (ptot > nwant) rhigh = rad
     if (ptot < nwant) rlow = rad
  end do
! ----------------------------------------------------------------------------

! Check radius is not too big so sphere is 'culled' at edges.
  if (rad > 1.0_PR) stop 'Not enough particles in cube.  Use more particles.'

! Choose required number of particles
  pgas      = ptot
  pboundary = 0
  picm      = 0
  n         = 0
  nsteps    = 0
  time      = 0.0_DP
  call allocate_memory(.FALSE.)
  call types

! Put positions in r array and normalise so rhigh = 1.0
  do p=1,ptot
     i = list(p)
     sph(p)%r(1:NDIM) = rcloud * position(1:NDIM,i) / rad
     sph(p)%m         = 1.0_PR / real(ptot,PR)
     sph(p)%h         = 0.0_PR
     sph(p)%rho       = 0.0_PR
     sph(p)%temp      = 1.0_PR
     sph(p)%v(1:NDIM) = 0.0_PR
  end do
  
! write to files
  call write_data(out_file,out_file_form)

! Clean up memory
  call clean_up

  stop
END PROGRAM ic_sphere
