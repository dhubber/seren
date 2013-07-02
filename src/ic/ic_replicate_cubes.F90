! IC_REPLICATE_CUBES.F90
! D. A. Hubber - 12/9/2007
! Reads in relaxed sheet/cube (of unit length and one corner on the origin) 
! and replicates in all directions to produce a larger cube.  If the cube 
! is relaxed, then the larger cube will also be relaxed.
! in_file         : Input file name
! in_file_form    : Input file format
! nrepeat         : No. of replicas (in each dimension)
! out_file        : Output file name
! out_file_form   : Output file format
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_replicate_cubes
  use particle_module
  use hydro_module
  use filename_module
  use scaling_module
  use type_module
  use time_module
  implicit none

  character(len=256) :: out_file              ! Output file name
  integer :: i                                ! x-direction counter 
  integer :: j                                ! y-direction counter
  integer :: k                                ! z-direction counter
  integer :: iaux                             ! Auxilary particle counter
  integer :: ncubes                           ! Total number of replicas
  integer :: nrepeat                          ! Replicas in each dimension
  integer :: old_ptot                         ! Original number of particles
  integer :: p                                ! Particle counter
  real(kind=PR), allocatable :: position(:,:) ! Position of old particles  
  
#ifndef DIMENSIONLESS
  write(6,*) "Compiler flag error : Only works with DIMENSIONLESS flag on"
  stop
#endif

! Set parameters to default values
  call default_parameters

  write(6,*) "------------------------------"
  write(6,*) "      ic_replicate_cubes      "
  write(6,*) "------------------------------"
  write(6,*) 
  write(6,*) "Input data..."
  write(6,*) "Input file : "
  read(5,*) in_file
  write(6,*) "Input file format : "
  read(5,*) in_file_form
  write(6,*) "No. of replicants in each dimension : "
  read(5,*) nrepeat
  write(6,*) "Output file : "
  read(5,*) out_file
  write(6,*) "Output file format : "
  read(5,*) out_file_form

  if (nrepeat <= 0) call paramerror("nrepeat <= 0")

  in_file       = trim(adjustl(in_file))
  in_file_form  = trim(adjustl(in_file_form))
  out_file      = trim(adjustl(out_file))
  out_file_form = trim(adjustl(out_file_form))

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types
  call units

! Read files
  call read_data(in_file,in_file_form)

! Record original particle positions in temporary array
  allocate (position(1:NDIM,1:ptot))
  do p=1,ptot
     position(1:NDIM,p) = sph(p)%r(1:NDIM)
  end do

! Now, reallocate arrays for new replicas
  old_ptot  = ptot
  call clean_up

  ncubes    = nrepeat**(NDIM)
  ptot      = old_ptot*ncubes
  pgas      = ptot
  pboundary = 0
  picm      = 0
  n         = 0
  nsteps    = 0
  iaux      = 0
  time      = 0.0_DP
  call allocate_memory(.FALSE.)
  call types

  write(6,*) "old_ptot : ",old_ptot
  write(6,*)
 
! Create nrepeat**(NDIM) cubes and resize to fit from 0 to 1.
#if NDIM==2
  do p=1,old_ptot
     do i=1,nrepeat
        do j=1,nrepeat
           iaux = iaux + 1
           sph(iaux)%r(1) = (position(1,p) + real(i,PR) - 1.0_PR) &
                & / real(nrepeat,PR)
           sph(iaux)%r(2) = (position(2,p) + real(j,PR) - 1.0_PR) &
                & / real(nrepeat,PR)
        end do
     end do
  end do
#elif NDIM==3
  do p=1,old_ptot
     do i=1,nrepeat
        do j=1,nrepeat
           do k=1,nrepeat             
              iaux = iaux + 1
              sph(iaux)%r(1) = (position(1,p) + real(i,PR) - 1.0_PR) &
                   & /real(nrepeat,PR)
              sph(iaux)%r(2) = (position(2,p) + real(j,PR) - 1.0_PR) &
                   & /real(nrepeat,PR)
              sph(iaux)%r(3) = (position(3,p) + real(k,PR) - 1.0_PR) &
                   & /real(nrepeat,PR)
           end do
        end do
     end do
  end do
#endif
  
! Set all other quantities
  do p=1,ptot
     sph(p)%m = 1.0_PR / real(ptot,PR)
     sph(p)%h = 0.0_PR
     sph(p)%v(1:NDIM)    = 0.0_PR
     sph(p)%rho         = 0.0_PR
     sph(p)%temp        = 1.0_PR
  end do
  
! Write to snapshot files
  call write_data(out_file,out_file_form)

! Clean up memory
  call clean_up

  stop
END PROGRAM ic_replicate_cubes
