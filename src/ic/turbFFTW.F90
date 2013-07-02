! turbsub.F90
! Andrew McLeod 03/01/2008
! Almost certainly doesn't do what it is supposed to
! which is generate a turbulent velocity field from a
! customizable power spectrum
! ----------------------------------------------------------------------------

! Input variables:
! field_type - must be 0 (natural), 1 (curl free) or 2 (divergence free)
! rmswanted - Desired rms of velocity field
! alpha - Exponent of power spectrum e.g. -4 typical for power mainly
!    on large scales - this must normally be negative
! kextent - How big a Fourier field to calculate from -kextent to kextent
!    This is over 3D i.e. (-kextent:kextent,-kextent:kextent,-kextent:kextent)
!    If the power spectrum is steep e.g. alpha of -6, you can calculate
!    for less components (components are only calculated if the magnitude of
!    the k-vector is less than kextent e.g. component (kextent,kextent,kextent)
!    will have 0 power. Also the (0,0,0) component has 0 power.
! gridsize - Size of output velocity field grid (see below)
! vfield - Output velocity field grid from 1 to gridsize (over 3D)

#include "macros.h"

subroutine turbsub(field_type,rmswanted,alpha,kextent,gridsize,vfield)
  implicit none

#include "fftw3.f"

!  integer, parameter   :: NDIM=3          ! Number of dimensions (must be 3)
  integer, parameter   :: PR=selected_real_kind(p=15) ! double precision
  integer, parameter       :: ILP = selected_int_kind(r=15)  ! Long integer
!  real(kind=PR), parameter :: PI = 3.1415926536   ! Pi
!  real(kind=PR), parameter :: TWOPI = 6.283185307   ! 2 * Pi
  integer              :: n, clock
  integer, allocatable :: seed(:)

  integer, intent(in)  :: field_type      ! 0=natural, 1=curl free, 2=div free
  real(kind=PR), intent(in)  :: rmswanted ! desired r.m.s. velocity
  real(kind=PR), intent(in)  :: alpha     ! Power spectrum index (=0 for flat)
  integer, intent(in)        :: kextent   ! Limit of k components
  integer, intent(in)        :: gridsize  ! size of velocity grid
  real(kind=PR), intent(out) :: vfield(1:NDIM,1:gridsize,1:gridsize,1:gridsize)
                                          ! Output velocity field
  complex(kind=PR), allocatable :: complexfield(:,:,:)

  integer              :: kmin       ! min in k space
  integer              :: kmax       ! max in k space
  integer (kind=ILP)   :: plan       ! FFTW plan pointer
  integer              :: shift      ! FFTW shift to rearrange Fourier format
  integer              :: i,j,k
  integer              :: n1,n2,n3,k1,k2,k3 ! Real space and fourier space loop counters
  integer              :: d          ! Dimension counter
  real(kind=PR)        :: p0         ! Power coefficient (overridden by normalisation)
  real(kind=PR)        :: F(1:NDIM)  ! Fourier component vector
  real(kind=PR)        :: unitk(1:3) ! Unit vector either parallel or perpendicular to k
  real(kind=PR), allocatable :: power(:,:,:,:), phase(:,:,:,:) ! Fourier components
  real(kind=PR)        :: Rnd(1:3),w ! Random numbers, variable in Gaussian calculation
!   real(kind=PR)        :: sumreal    ! Running sum of cosine terms
!                                      ! (real part of a fourier transform)
!   real(kind=PR)        :: theta      ! Coefficient when calculating fourier transform
  real(kind=PR)        :: rms        ! r.m.s. velocity
  logical              :: divfree, curlfree ! Selecting div-free or curl-free turbulence

  kmin=-(gridsize/2 - 1)
  kmax=gridsize/2

  p0 = 0.05    ! Power coefficient (I don't think this matters)
  divfree = .FALSE.
  curlfree = .FALSE.
  if (field_type == 1) curlfree = .TRUE.
  if (field_type == 2) divfree = .TRUE.
  if (NDIM /= 3) stop "turbsub: only NDIM==3 supported!"

  call RANDOM_SEED(size = n)
  allocate(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)
! ----------------------------------------------------------------------------------------------

  allocate(power(1:NDIM,kmin:kmax,kmin:kmax,kmin:kmax))
  allocate(phase(1:NDIM,kmin:kmax,kmin:kmax,kmin:kmax))
  allocate(complexfield(1:gridsize,1:gridsize,1:gridsize))

  vfield = 0.0

  write (6,*) "Calculating power and phase spectra"

  power(1:NDIM,:,:,:) = 0.0; phase(1:NDIM,:,:,:) = 0.0

  ! Define wave vectors in Fourier space
  ! Each wave vector has coordinates in Fourier space, random phases in
  ! three dimensions, and power in three dimensions, giving a power
  ! and a phase vector field in Fourier space
  ! With a 'natural' field type there is no coordination between x,y,z components
  ! With a div-free or curl-free field, there is a relation between
  ! coordinates in Fourier space and power vector. For a div-free field,
  ! power vectors are perpendicular to the Fourier coordinate vector, and
  ! for a curl-free field power vectors are (anti-)parallel to the Fourier
  ! coordinate vector
  ! (i,j,k) is the Fourier coordinate vector
  ! power(1:3,i,j,k) is the power vector at Fourier coordinates i,j,k
  ! phase(1:3,i,j,k) is the phase vector at Fourier coordinates i,j,k
  do k=kmin,kmax
    do j=kmin,kmax
      do i=kmin,kmax
        if (i==0.AND.j==0.AND.k==0) cycle ! Central power = 0
        !if (.NOT.(i==-j.AND.j==k)) cycle ! Along diagonal only
        !if (.NOT.(i+j==kmax(1).AND.k==1)) cycle ! Along diagonal (positive quadrant) only
        !if (.NOT.((i==0.AND.j==0).OR.(i==0.AND.k==0).OR.(j==0.AND.k==0))) cycle ! Axes only
        if (sum((/i,j,k/)**2) > kextent**2) cycle ! Sort of sphere
        ! Power value, to be multipled by random power selected from a Gaussian
        ! This is what gives the slope of the eventual power spectrum
        F = p0 * sqrt(real(i,PR)**2 + real(j,PR)**2 + real(k,PR)**2)**alpha
        do d=1,NDIM
          call random_number(Rnd(1))
          ! Random phase between -pi and pi
          !Rnd(1) = 0.0
          phase(d,i,j,k) = ((Rnd(1) * 2.0) - 1.0 ) * PI
          ! Create Gaussian distributed random numbers
          ! Marsaglia version of the Box-Muller algorithm
          ! This creates a Gaussian with a spread of 1
          ! and centred on 0 (hopefully)
          do
            call random_number(Rnd(2:3))
            Rnd(2) = 2.0 * Rnd(2) - 1.0
            Rnd(3) = 2.0 * Rnd(3) - 1.0
            w = Rnd(2)*Rnd(2) + Rnd(3)*Rnd(3)
            if (w<1.0) exit
          end do
          w = sqrt( (-2.0 * log( w ) ) / w )
          Rnd(2) = Rnd(2) * w
          Rnd(3) = Rnd(3) * w ! Don't actually need this

          ! Random power multiplied by power spectrum
          F(d) = Rnd(2) * F(d)
        end do
        ! Helmholtz decomposition!
        unitk = (/real(i,PR),real(j,PR),real(k,PR)/)
        unitk = unitk / sqrt(sum(unitk**2))
        if (curlfree) then
          ! For curl free turbulence, vector F should be parallel/anti-parallel to vector k
          power(1:NDIM,i,j,k) = unitk * dot_product(F,unitk)
        else if (divfree) then
          ! For divergence free turbulence, vector F should be perpendicular to vector k
          power(1:NDIM,i,j,k) = F(1:NDIM) - unitk * dot_product(F,unitk)
        else
          ! Natural mix of turbulence
          power(1:NDIM,i,j,k) = F(1:NDIM)
        end if
      end do
    end do
  end do

!  ! Calculate inverse discrete cosine transform at point i,j,k in real space
!   do d=1,NDIM
! !$OMP PARALLEL DO PRIVATE(n1,n2,n3,k1,k2,k3,theta,sumreal)
!     do n1=0,gridsize-1
!       !write (6,'(I0,A)',ADVANCE="NO") n1, " "
!       do n2=0,gridsize-1
!         do n3=0,gridsize-1
! 
!           sumreal = 0.0
!           ! Loop over fourier spectrum component things
!           do k1=-kextent,kextent
!             do k2=-kextent,kextent
!               do k3=-kextent,kextent
!                 theta = (real(n1*k1,PR)/real(gridsize,PR)) + &
!                       & (real(n2*k2,PR)/real(gridsize,PR)) + &
!                       & (real(n3*k3,PR)/real(gridsize,PR))
!                 sumreal = sumreal + (  power(d,k1,k2,k3) * &
!                       & cos( (TWOPI*theta) + phase(d,k1,k2,k3) )  )
!               end do
!             end do
!           end do
! 
!           vfield(d,n1+1,n2+1,n3+1) = sumreal
! 
!         end do
!       end do
!     end do
! !$OMP END PARALLEL DO
! 
!   end do

  ! FFTW time
  call dfftw_plan_dft_3d(plan, gridsize, gridsize, gridsize, complexfield, &
    & complexfield, FFTW_BACKWARD, FFTW_ESTIMATE)

  ! Remap power and phases
  shift = -kmin
  do d=1,NDIM
     power(d,:,:,:) = cshift(power(d,:,:,:),shift,1)
     power(d,:,:,:) = cshift(power(d,:,:,:),shift,2)
     power(d,:,:,:) = cshift(power(d,:,:,:),shift,3)
     phase(d,:,:,:) = cshift(phase(d,:,:,:),shift,1)
     phase(d,:,:,:) = cshift(phase(d,:,:,:),shift,2)
     phase(d,:,:,:) = cshift(phase(d,:,:,:),shift,3)
  
    complexfield = power_phase_to_complex(power(d,:,:,:),phase(d,:,:,:))
    call dfftw_execute_dft(plan, complexfield, complexfield)
    vfield(d,:,:,:) = real(complexfield)

  end do
  call dfftw_destroy_plan(plan)

  ! Normalise to some desired rms velocity
  rms = sum(sqrt(sum(vfield(1:NDIM,:,:,:)**2,1)))/real(gridsize**NDIM,PR)
  vfield = vfield * (rmswanted/rms)

  return

  contains

   complex(kind=PR) elemental function power_phase_to_complex(power,phase)
      real (kind=PR), intent(in)           :: power, phase
      real (kind=PR)                       :: re, im

      re = power * cos(phase)
      im = power * sin(phase)

      power_phase_to_complex = cmplx(re, im)

   end function power_phase_to_complex

end subroutine turbsub
