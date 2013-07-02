! IC_GAS_BINARY.F90
! R. J. Allison - 19/05/2011
! Creates a single binary, and puts it in a background potential
! The system is setup with the stars at periastron, that makes it easier to
! setup the initial state of the system. 
! The users needs to input the mass of the stars, the semi-major axis, and
! the eccentricity.
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_gas_binary

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
  character(len=256) :: file1            ! File with RANDOM_CUBE
  character(len=256) :: out_file         ! Output file name
  character(len=256) :: out_file__form   ! File 1 format
  real(kind=PR) :: m1, m2                ! Mass of primary, secondary
  real(kind=PR) :: sma                   ! Semi-major axis
  real(kind=PR) :: ecc                   ! Eccentricity
  real(kind=PR) :: s_radius              ! Radius of 'sinks'
  real(kind=PR) :: magr                  ! Magnitude of binary star seperation
  real(kind=PR) :: mu                    ! Standard gravitational parameter
  real(kind=PR) :: magv                  ! Magnitude of orbital velocity
  real(kind=PR) :: r1(3),r2(3)           ! Pos. of primary, secondary
  real(kind=PR) :: v1(3),v2(3)           ! Vel. of primary, secondary
  real(kind=PR) :: COM(3)                
  real(kind=PR) :: COV(3)
  integer :: p,p1
  real(kind=PR), allocatable :: r1p(:,:)   ! Pos. of SPH part.
  real(kind=PR), allocatable :: v1p(:,:)   ! Vel. of SPH part.
  real(kind=PR), allocatable :: m1p(:)     ! Mass of SPH part.
  real(kind=PR), allocatable :: temp1p(:)  ! Temp of SPH part.
  real(kind=PR), allocatable :: u1p(:)     ! Int. En. of SPH part.
  real(kind=PR) :: rcentre(NDIM)

  rcentre = 0.0



! Read parameters file
! ----------------------------------------------------------

! Read parameters file
  call default_parameters

! Read in IC_GAS_BINARY params file  
  write(6,*) "Opening ic_gas_binary params file"
  open(unit=1,file="gas_binary.dat",status='old')
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*) file1            !Background potential input
  read(1,*) file1_form
  read(1,*) m1, m2           !Mass of primary, secondary
  read(1,*) sma              !Semi-major axis
  read(1,*) ecc              !Eccentricity
  
! Initialise some variables using parameters
  call initialize_seren_variables_1
  call units

! Setup binary system
! ----------------------------------------------------------

! We know the system is a periastron, so to find out how far 
! apart the stars are...

  magr = (sma*(1. - ecc**2)) / (1. + ecc)

! Now we know the seperation of the stars, we can find out the 
! magnitude of their velocites...

  mu = m1 + m2      !Dimensionless so G = 1
  magv = mu*((2./magr) - (1./sma))
  magv = sqrt(magv)

! Because the stars are at periastron, there is only one velocity
! component

! Set star 1 to be at zero in position and velocity

  r1 = 0.
  v1 = 0.

! Setup star 2

  r2(1) = magr     !r2(1) = magr*COS(theta) ; theta = 0
  r2(2) = 0.       !r2(2) = magr*SIN(theta) 
  r2(3) = 0.

  v2(1) = 0.       !v2(1) = magv*COS(phi) ; periastron => phi = pi/2
  v2(2) = magv     !v2(2) = magv*SIN(phi)
  v2(3) = 0.

! Shift binary system to COM and COV
! ----------------------------------------------------------

   COM(1:3)=(r1(1:3)*m1 + r2(1:3)*m2)/(m1 + m2)
   COV(1:3)=(v1(1:3)*m1 + v2(1:3)*m2)/(m1 + m2) 

   r1 = r1 - COM
   r2 = r2 - COM
   v1 = v1 - COV
   v2 = v2 - COV

! Binary system is setup!


! Setup background potential
! ----------------------------------------------------------

! Read in SPH particles from file (makes this program more
! general...

  call read_data(file1,file1_form)
  p1 = ptot
  allocate(r1p(1:NDIM,1:p1))
  allocate(v1p(1:NDIM,1:p1))
  allocate(m1p(1:p1))
  allocate(temp1p(1:p1))
  allocate(u1p(1:p1))

! Read in SPH particle info
  do p = 1,p1
     r1p(1:3,p) = sph(p)%r(1:3)/real(rscale,PR)
     v1p(1:3,p) = sph(p)%v(1:3)/real(vscale,PR)
     m1p(p)     = sph(p)%m/real(mscale,PR)
     temp1p(p)  = sph(p)%temp
     u1p(p)     = sph(p)%u
  end do

  call clean_up

! Put into Seren arrays
! ----------------------------------------------------------

! Set particle types and allocate memory
  pboundary = 0
  picm      = 0
  pgas      = p1 
  pcdm      = 0
  pdust     = 0
  pion      = 0
  stot      = 2
  ptot      = p1 

  call allocate_memory
  call types

! Put binaries into Seren arrays
! Primary
  sink(1)%r(1:3)    = r1(1:3)
  sink(1)%v(1:3)    = v1(1:3)
  sink(1)%m         = m1
  sink(1)%radius    = s_radius
  sink(1)%h         = INVKERNRANGE*sink(1)%radius

! Secondary
  sink(2)%r(1:3)    = r2(1:3)
  sink(2)%v(1:3)    = v2(1:3)
  sink(2)%m         = m2
  sink(2)%radius    = s_radius
  sink(2)%h         = INVKERNRANGE*sink(2)%radius


! SPH particles

  do p = 1,p1
     sph(p)%r(1:3)    = r1p(1:3,p)
     sph(p)%v(1:3)    = v1p(1:3,p)
     sph(p)%m         = m1p(p)
     sph(p)%temp      = temp1p(p)
     sph(p)%u         = u1p(p)
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

END PROGRAM ic_gas_binary
