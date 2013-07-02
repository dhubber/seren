! IC_COLLIDE_PLUM.F90
! R. J. Allison - 17/05/2011
! Creates two plummer spheres which collide.
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_collide_plum

  use particle_module
  use hydro_module
  use sink_module
  use filename_module
  use periodic_module
  use type_module
  use scaling_module
  use constant_module
  implicit none

!Lots of SEREN stuff...
  character(len=256) :: file1_form       ! File 1 format
  character(len=256) :: file2_form       ! File 2 format
  character(len=256) :: file1            ! File with RANDOM_CUBE
  character(len=256) :: file2            ! File with PLUMMER
  character(len=256) :: out_file         ! Output file name
  real(kind=PR) :: r_shift               ! Shift from centre for plum.
  real(kind=PR) :: v_shift1, v_shift2    ! Shift in vel for plum.s
  real(kind=PR) :: s_radius              ! Radius of 'sinks'
  integer :: p1,p2
  integer :: s1,s2
  real(kind=PR), allocatable :: r1s(:,:),r2s(:,:) !Pos of 'sinks' 
  real(kind=PR), allocatable :: v1s(:,:),v2s(:,:) !Vel of 'sinks'
  real(kind=PR), allocatable :: r1p(:,:),r2p(:,:) !Pos of SPH part.
  real(kind=PR), allocatable :: v1p(:,:),v2p(:,:) !Vel of SPH part.
  real(kind=PR), allocatable :: m1s(:),m2s(:)     !Mass of 'sinks'
  real(kind=PR), allocatable :: m1p(:),m2p(:)     !Mass of SPH part.
  real(kind=PR), allocatable :: temp1p(:),temp2p(:)!Temp of SPH part.
  real(kind=PR), allocatable :: u1p(:),u2p(:)     !Int. En. of SPH part.
  real(kind=PR) :: m_plum1,m_plum2                !Mass of plum. sph.
  integer :: p                                    !Particle counter
  integer :: s                                    !N-body counter
  real(kind=PR) :: rcentre(NDIM)

  rcentre = 0.0

! Read parameters file
! ----------------------------------------------------------------------------

! Read parameters file
  call default_parameters

! Read in IC_COLLIDE_PLUM params file  
  write(6,*) "Opening ic_collide_plum params file"
  open(unit=1,file="coll_plum.dat",status='old')
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*) file1
  read(1,*) file2
  read(1,*) file1_form
  read(1,*) file2_form
  read(1,*) r_shift     
  read(1,*) v_shift1, v_shift2
  read(1,*) s_radius
  close(1)
  
  runit = 'pc'
  vunit = 'km_s'
  munit = 'm_sun'
  rhounit = 'm_sun_pc3'

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call units

! Read in first plummer sphere
! ----------------------------------------------------------------------------
  call read_data(file1,file1_form)
  p1 = ptot
  s1 = stot
  allocate(r1s(1:NDIM,1:s1))
  allocate(v1s(1:NDIM,1:s1))
  allocate(r1p(1:NDIM,1:p1))
  allocate(v1p(1:NDIM,1:p1))
  allocate(m1s(1:s1))
  allocate(m1p(1:p1))
  allocate(temp1p(1:p1))
  allocate(u1p(1:p1))

! Fill arrays, and shift position and velocities in one dimension
! First plummer moves -ve in dimension 1 by some factor r_plum1
! and velocity shifts +ve in dimension 1 by some factor v_shift1


  m_plum1 = 0.
! Fill and shift sink array
  do s = 1,s1
     r1s(1:3,s) = sink(s)%r(1:3)/real(rscale,PR)
     v1s(1:3,s) = sink(s)%v(1:3)/real(vscale,PR)
     m1s(s)     = sink(s)%m/real(mscale,PR)
     m_plum1    = m_plum1 + m1s(s)
     
     ! Shift sinks
     r1s(1,s) = r1s(1,s) - r_shift
     v1s(1,s) = v1s(1,s) + v_shift1
  end do

! Fill and shift particle array
  do p = 1,p1
     r1p(1:3,p) = sph(p)%r(1:3)/real(rscale,PR)
     v1p(1:3,p) = sph(p)%v(1:3)/real(vscale,PR)
     m1p(p)     = sph(p)%m/real(mscale,PR)
     m_plum1    = m_plum1 + m1p(p)
     temp1p(p)  = sph(p)%temp
     u1p(p)     = sph(p)%u

     !Shift particles
     r1p(1,p) = r1p(1,p) - r_shift
     v1p(1,p) = v1p(1,p) + v_shift1
  end do

  call clean_up


! Read in second plummer sphere
! ----------------------------------------------------------------------------
  call read_data(file2,file2_form)
  p2 = ptot
  s2 = stot
  allocate(r2s(1:NDIM,1:s2))
  allocate(v2s(1:NDIM,1:s2))
  allocate(r2p(1:NDIM,1:p2))
  allocate(v2p(1:NDIM,1:p2))
  allocate(m2s(1:s2))
  allocate(m2p(1:p2))
  allocate(temp2p(1:p2))
  allocate(u2p(1:p2))

! Fill arrays, and shift position and velocities in one dimension
! Second plummer moves +ve in dimension 1 by some factor r_plum2
! and velocity shifts -ve in dimension 1 by some factor v_shift2

  m_plum2 = 0.
! Fill and shift sink array
  do s = 1,s2
     r2s(1:3,s) = sink(s)%r(1:3)/real(rscale,PR)
     v2s(1:3,s) = sink(s)%v(1:3)/real(vscale,PR)
     m2s(s)     = sink(s)%m/real(mscale,PR)
     m_plum2 = m_plum2 + m2s(s)

     ! Shift sinks
     r2s(1,s) = r2s(1,s) + r_shift
     v2s(1,s) = v2s(1,s) - v_shift2
  end do

! Fill and shift particle array
  do p = 1,p2
     r2p(1:3,p) = sph(p)%r(1:3)/real(rscale,PR)
     v2p(1:3,p) = sph(p)%v(1:3)/real(vscale,PR)
     m2p(p)     = sph(p)%m/real(mscale,PR)
     m_plum2 = m_plum2 + m2p(p)
     temp2p(p)  = sph(p)%temp
     u2p(p)     = sph(p)%u

     !Shift particles
     r2p(1,p) = r2p(1,p) + r_shift
     v2p(1,p) = v2p(1,p) - v_shift2
  end do

  call clean_up

! Combine into main arrays
! ----------------------------------------------------------------------------

! Set particle types and allocate memory
  pboundary = 0
  picm      = 0
  pgas      = p1 + p2
  pcdm      = 0
  pdust     = 0
  pion      = 0
  stot      = s1 + s2
  ptot      = p1 + p2

  call allocate_memory
  call types

! Store first plummer in memory
! Sinks
  do s = 1,s1
     sink(s)%r(1:3)    = r1s(1:3,s)
     sink(s)%v(1:3)    = v1s(1:3,s)
     sink(s)%m         = m1s(s)
     sink(s)%radius    = s_radius
     sink(s)%h         = INVKERNRANGE*sink(s)%radius
  end do

! SPH particles
  do p = 1,p1
     sph(p)%r(1:3)    = r1p(1:3,p)
     sph(p)%v(1:3)    = v1p(1:3,p)
     sph(p)%m         = m1p(p)
     sph(p)%temp      = temp1p(p)
     sph(p)%u         = u1p(p)
  end do

! Store second plummer in memory
! Sinks
  do s = 1,s2
     sink(s1+s)%r(1:3)    = r2s(1:3,s)
     sink(s1+s)%v(1:3)    = v2s(1:3,s)
     sink(s1+s)%m         = m2s(s)
     sink(s1+s)%radius    = s_radius
     sink(s1+s)%h         = INVKERNRANGE*sink(s)%radius
  end do

! SPH particles
  do p = 1,p2
     sph(p1+p)%r(1:3)    = r2p(1:3,p)
     sph(p1+p)%v(1:3)    = v2p(1:3,p)
     sph(p1+p)%m         = m2p(p)
     sph(p1+p)%temp      = temp2p(p)
     sph(p1+p)%u         = u2p(p)
  end do   

! Initialise SEREN stuff...!
! ----------------------------------------------------------

  restart = .false.
  call initialize_seren_variables_2
  out_file = trim(adjustl(out_file))

  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  call write_data_debug("ICPLUM_RCUBE.debug.dat",rcentre(1:NDIM))
#endif

! Clean-up all memory
  call clean_up
  
  stop

END PROGRAM ic_collide_plum
