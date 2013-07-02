s! IC_BINARY.F90
! D. A. Hubber - 12/05/2010
! Create a binary system from two stars loaded from separate files.  
! Often, the binaries will be polytropes created/relaxed previously.
! file1          : Input filename
! file2          : Input filename
! file1_form     : Input file format
! file2_form     : Input file format
! out_file       : Output filename
! out_file_form  : Output file format
! abin           : ..
! eccent         : ..
! corotate       : ..
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_binary
  use interface_module, only : comperror,paramerror,read_data,write_data
  use particle_module
  use hydro_module
  use filename_module
  use scaling_module
  use type_module
  implicit none

  character(len=256) :: file1_form            ! File 1 format
  character(len=256) :: file2_form            ! File 2 format
  character(len=256) :: file1                 ! File with LHS box
  character(len=256) :: file2                 ! File with RHS box
  character(len=256) :: out_file              ! Name of output file
  character(len=50) :: pertmode               ! Perturbation mode
  character(len=50) :: velmode                ! Velocity mode
  logical :: corotate                         ! Stars rotate with orbit?
  integer :: i                                ! ..
  integer :: p                                ! Particle counter
  integer :: p1                               ! Particle counter
  integer :: p2                               ! Particle counter
  real(kind=PR) :: abin                       ! Binary separation
  real(kind=PR) :: angspeed                   ! ..
  real(kind=PR) :: eccent                     ! Binary eccentricity
  real(kind=PR) :: m1                         ! Mass 1
  real(kind=PR) :: m2                         ! ..
  real(kind=PR) :: period                     ! ..
  real(kind=PR) :: rcom1(1:NDIM)              ! ..
  real(kind=PR) :: rcom2(1:NDIM)              ! ..
  real(kind=PR) :: vrot                       ! ..
  real(kind=PR), allocatable :: parray1(:,:)  ! ..
  real(kind=PR), allocatable :: parray2(:,:)  ! ..

  real(kind=PR) :: binen                      ! ..
  real(kind=PR) :: L(1:3)                  ! Angular momentum vector
  real(kind=PR) :: Lsqd                    ! Angular momentum squared
  real(kind=PR) :: q                       ! Mass-ratio
  real(kind=PR) :: reduced_mass            ! Reduced mass
  real(kind=PR) :: sma                     ! Semi-major axis

! Set default values for all parameters
  call default_parameters

  write(6,*) "------------------------------"
  write(6,*) "          ic_binary           "
  write(6,*) "------------------------------"
  write(6,*) 
  write(6,*) "Input data..."

  open(1, file="binaryparams.dat", status="old", form="formatted")
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) file1
  read(1,*) file2
  read(1,*) file1_form
  read(1,*) file2_form
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) abin
  read(1,*) eccent
  read(1,*) corotate
  close(1)

#if NDIM != 3
  call comperror("NDIM = 3 not selected")
#endif
  file1         = trim(adjustl(file1))
  file1_form    = trim(adjustl(file1_form))
  file2         = trim(adjustl(file2))
  file2_form    = trim(adjustl(file2_form))
  out_file      = trim(adjustl(out_file))
  out_file_form = trim(adjustl(out_file_form))

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types

! Calculate all over scaling variables
  call units

! Read in relaxed, uniform density sphere
  call read_data(file1,file1_form)
  p1 = ptot
  allocate(parray1(1:SMOO,1:p1))
  parray1(1:SMOO,1:p1) = parray(1:SMOO,1:p1)
  call clean_up

! Read in relaxed, uniform density sphere
  call read_data(file2,file2_form)
  p2 = ptot
  allocate(parray2(1:SMOO,1:p2))
  parray2(1:SMOO,1:p2) = parray(1:SMOO,1:p2)
  call clean_up

  write(6,*) "Now calculating centre of mass of both stars"

! Calculate position of centre of mass, and total mass, of star 1
  m1 = 0.0_DP
  rcom1(1:NDIM) = 0.0_DP
  do p=1,p1
     m1 = m1 + parray1(MASS,p)
     rcom1(1:NDIM) = rcom1(1:NDIM) + parray1(MASS,p)*parray1(1:NDIM,p)
  end do
  rcom1(1:NDIM) = rcom1(1:NDIM)/m1

! Calculate position of centre of mass, and total mass, of star 2
  m2 = 0.0_DP
  rcom2(1:NDIM) = 0.0_DP
  do p=1,p2
     m2 = m2 + parray2(MASS,p)
     rcom2(1:NDIM) = rcom2(1:NDIM) + parray2(MASS,p)*parray2(1:NDIM,p)
  end do
  rcom2(1:NDIM) = rcom2(1:NDIM)/m2

! Calculate binary properties
  vrot = sqrt((m1 + m2)/abin)
  period = TWOPI*sqrt(abin**3/(m1 + m2))

  if (corotate) then
     angspeed = vrot/abin
  else
     angspeed = 0.0_PR
  end if

  write(6,*) "Binary mass  : ",m1,m2,m1+m2
  write(6,*) "Separation   : ",abin
  write(6,*) "Period       : ",period
  write(6,*) "eccentricity : ",eccent

! Calculate angular momentum squared of binary
  reduced_mass = m1*m2 / (m1 + m2)
  L(3) = vrot*abin*reduced_mass
  Lsqd = L(3)*L(3)
        
! Now calculate all binary parameters
  binen = 0.5*reduced_mass*vrot*vrot - m1*m2/abin
  sma    = -0.5_DP*m1*m2/binen
  eccent = 1.0_DP - Lsqd/(m1 + m2)/sma/(reduced_mass**2)
  eccent = max(0.0_DP,eccent)
  eccent = sqrt(eccent)
!  write(6,*) "ECCENT :",Lsqd,m1,m2,sma,reduced_mass,eccent
  period = 2.0_DP*PI*sqrt(sma**3/(m1 + m2))
  if (m1 > m2) then
     q = m2/m1
  else
     q = m1/m2
  end if

! Now write out binary properties again
  write(6,*) "Recomputed!!"
  write(6,*) "Separation   : ",sma
  write(6,*) "Period       : ",period
  write(6,*) "eccentricity : ",eccent
  write(6,*) "vrot         : ",vrot,L(3),Lsqd

  ptot = p1 + p2
  call allocate_memory(.FALSE.)


  do i=1,p1
     p = i
     sph(p)%r(1:NDIM) = parray1(1:NDIM,i) - rcom1(1:NDIM)
     v(1,p) = m1*vrot/(m1 + m2) - sph(p)%r(2)*angspeed
     v(2,p) = sph(p)%r(1)*angspeed
     v(3,p) = 0.0_PR
     sph(p)%r(2) = sph(p)%r(2) - m1*abin/(m1 + m2)
     sph(p)%m = parray1(MASS,i)
     sph(p)%h = parray1(SMOO,i)
!     v(1,p) = -sph(p)%r(2)*vrot/abin
!     v(2,p) = sph(p)%r(1)*vrot/abin
!     v(3,p) = 0.0_PR
  end do

  do i=1,p2
     p = p1 + i
     sph(p)%r(1:NDIM) = parray2(1:NDIM,i) - rcom2(1:NDIM)
     v(1,p) = -m2*vrot/(m1 + m2) - sph(p)%r(2)*angspeed
     v(2,p) = sph(p)%r(1)*angspeed
     v(3,p) = 0.0_PR
     sph(p)%r(2) = sph(p)%r(2) + m2*abin/(m1 + m2)
     sph(p)%m = parray2(MASS,i)
     sph(p)%h = parray2(SMOO,i)
!     v(1,p) = -sph(p)%r(2)*vrot/abin
!     v(2,p) = sph(p)%r(1)*vrot/abin
!     v(3,p) = 0.0_PR
  end do


  pgas = ptot
  picm = 0

! Setting up for different particle types
  call types

! Write particle data to file
  call write_data(out_file,out_file_form)

! Deallocate memory
  call clean_up

  deallocate(parray2)
  deallocate(parray1)

  stop
END PROGRAM ic_binary
