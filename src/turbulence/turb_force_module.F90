! TURB_FORCE_MODULE.F90
! A.McLeod - 24/11/2012
! Utility routines for turbulent forcing
! ============================================================================

#include "macros.h"

! ============================================================================
module turb_force_module
   use definitions
   use turbulence_module
   implicit none
   
#include "fftw3.f"
#define HERMITIAN_FIELD 1
#define POWER_LAW_TURBULENCE 1

   contains

   subroutine find_conj_pair(i,j,k,ii,jj,kk)
      integer, intent(in)     :: i,j,k         ! Grid coordinates
      integer, intent(out)    :: ii,jj,kk      ! Hermitian pair
      
      ii = TURB_GS - i
      if (ii >= TURB_GS) ii = ii - TURB_GS
      jj = TURB_GS - j
      if (jj >= TURB_GS) jj = jj - TURB_GS
      kk = TURB_GS - k
      if (kk >= TURB_GS) kk = kk - TURB_GS
      
      return
   end subroutine find_conj_pair
   
   function find_unitk(i,j,k,limit)
      real (kind=PR)           :: find_unitk(1:NDIM)
      integer, intent(in)      :: i,j,k,limit
      integer                  :: ii,jj,kk
      
      ii = i; jj = j; kk = k
      if (i > limit) ii = TURB_GS - i
      if (j > limit) jj = TURB_GS - j
      if (k > limit) kk = TURB_GS - k
      
      find_unitk = (/real(ii,PR),real(jj,PR),real(kk,PR)/)
      find_unitk = find_unitk / sqrt(sum(find_unitk**2))
      
      return
   end function find_unitk
   
   function power_spectrum(k)
      integer, intent(in)        :: k(1:NDIM)      ! Wavevector
      real(kind=PR)              :: power_spectrum ! Power value
   
      real(kind=PR)              :: k_mag          ! Wavevector magnitude
   
      ! alpha^-2 power spectrum
#if defined(POWER_LAW_TURBULENCE)
      if (all(k==0)) then
         power_spectrum = 0
         return
      end if
   
      k_mag = sqrt(real(sum(k**2),DP))

      if (k_mag > (TURB_GS/2)) then
         power_spectrum = 0
         return
      end if
   
      power_spectrum = k_mag**(-2)
#endif
      
      ! 'parabola' large-scale modes power spectrum
#if defined(PARABOLA_TURBULENCE)
      power_spectrum = 0._PR
      k_mag = sqrt(real(sum(k**2),DP))
      if ((k_mag > 1.0_PR) .AND. (k_mag < 3.0_PR)) then
          power_spectrum = 1.0 - (k_mag-2.0)**2
      end if
#endif
      
      return

   end function power_spectrum
   
   function gaussian_cmplx()
      ! Generates NDIM complex random variates, where each variate has a
      ! complex value drawn from a Gaussian distribution EITHER
      ! by assigning a Gaussian magnitude and uniformly distributed argument OR
      ! by assigning Gaussian real and imaginary components
      complex(kind=PR)     :: gaussian_cmplx(1:NDIM)
                                          ! Random gaussian values

      integer              :: d           ! Dimension counter
      real(kind=PR)        :: Rnd(1:NDIM) ! Random numbers
      real(kind=PR)        :: w           ! Aux. variable for random Gaussian
      real(kind=PR)        :: mag(1:NDIM), arg(1:NDIM) ! Magnitude and argument
      ! Create Gaussian distributed random numbers
      ! Marsaglia version of the Box-Muller algorithm
      ! This creates a Gaussian with a spread of 1
      ! and centred on 0 (hopefully)
      ! You get two independent variates each time
      
      ! Method 1: Gaussian magnitude and uniformly-random argument
      do d=1,NDIM
         do
            call random_number(Rnd(1:2))
            Rnd(1) = 2.0 * Rnd(1) - 1.0
            Rnd(2) = 2.0 * Rnd(2) - 1.0
            w = Rnd(1)**2 + Rnd(2)**2
            if (w<1.0) exit
         end do
         w = sqrt( (-2.0 * log( w ) ) / w )
         
         mag(d) = Rnd(1) * w
         ! Throwing away second random variate (something of a waste)
      end do
      
      call random_number(Rnd(1:NDIM))
      arg = (Rnd(1:NDIM) * 2.0 * PI) - 1.0
      
      gaussian_cmplx(1:NDIM) = cmplx(mag * cos(arg), mag * sin(arg))
      
      ! Method 2: Gaussian real and imaginary components
!       do d=1,NDIM
!          do
!             call random_number(Rnd(1:2))
!             Rnd(1) = 2.0 * Rnd(1) - 1.0
!             Rnd(2) = 2.0 * Rnd(2) - 1.0
!             w = Rnd(1)**2 + Rnd(2)**2
!             if (w<1.0) exit
!          end do
!          w = sqrt( (-2.0 * log( w ) ) / w )
!          
!          gaussian_cmplx(d) = cmplx(Rnd(1) * w, Rnd(2) * w)
!       end do
      
      return
   end function gaussian_cmplx
   
   subroutine add_turbulence(turb_field, dt)
      ! Take a complex Fourier field of turbulence, and add a random component
      ! from a Wiener process (in this case Gaussian perturbations) multiplied
      ! by an initial distribution of power and then projected by Helmholtz
      ! decomposition.
      ! Optionally, this field can be Hermitian i.e. purely real after the
      ! transform, but I'm not sure if this produces a purely even field...
      complex(kind=PR), intent(inout) :: turb_field(1:NDIM, 0:TURB_GS-1,&
                                                   &0:TURB_GS-1, 0:TURB_GS-1)
                                              ! Complex field to add to
      real(kind=PR), intent(in)       :: dt   ! Width of Gaussian which drives
                                              ! Wiener process

      integer              :: i, j, k         ! Loop variables
      integer              :: ii,jj,kk        ! Hermitian pair
      integer              :: halfway         ! Half of grid size
   
#if defined(HERMITIAN_FIELD)
      logical              :: hermitian_pair  ! Do we need to take conjugate?
      logical              :: own_conjg       ! Are we are own conjugate pair?
#endif
      
      complex(kind=PR)     :: complex_vec(1:NDIM) ! Complex Fourier vector
      complex(kind=PR)     :: unitk_cmplx(1:3)    ! Unit vector parallel to k
      complex(kind=PR)     :: comp_cmplx(1:NDIM)  ! Compressive component
      complex(kind=PR)     :: sol_cmplx(1:NDIM)   ! Solenoidal components
   
      halfway = TURB_GS / 2 ! integer division
   
      do k=0,TURB_GS-1
         do j=0,TURB_GS-1
            do i=0,TURB_GS-1
#if defined(HERMITIAN_FIELD)
               hermitian_pair = .FALSE.
#endif
               ! Check there is any power in this mode, else cycle
               if (power_spec(i,j,k) == 0.0_PR) cycle
               
#if defined(HERMITIAN_FIELD)
               ! Test if we are own conjugate or need to use another conjugate
               own_conjg = .FALSE.
               call find_conj_pair(i,j,k,ii,jj,kk)
               if (i==ii .AND. j==jj .AND. k==kk) own_conjg = .TRUE.
               if (k > halfway) then
                  hermitian_pair = .TRUE.
               else if (k == 0 .OR. 2*k == TURB_GS) then
                  if (j > halfway) then
                     hermitian_pair = .TRUE.
                  else if (j == 0 .OR. 2*j == TURB_GS) then
                     if (i > halfway) hermitian_pair = .TRUE.
                  end if
               end if
               if (hermitian_pair) then
                  turb_field(1:NDIM,i,j,k) = &
                     & conjg(turb_field(1:NDIM,ii,jj,kk))
                  cycle
               end if
#endif
               
               ! Ornstein-Uhlenbeck process
               ! dF(k,t) = F_0(k) P(k) dW - F(k,t) dt / T
               ! where F(k,t) is the vector fourier amplitude,
               ! F_0 is the power spectrum distribution of amplitudes,
               ! P(k) is the projection (Helmholtz decomposition),
               ! dt is the timestep of integration,
               ! and dW is the Wiener process, described by a random variate
               ! from a normal distribution with a mean of zero and a
               ! variance of dt i.e. a standard deviation of sqrt(dt)
               
               ! In this we only calculate the first term - the random
               ! growth term.

               ! Random power/phase (complex Gaussian) multiplied by power
               ! spectrum and multiplied by width dt
               complex_vec(1:NDIM) = gaussian_cmplx() * sqrt(dt) * &
                               & power_spec(i,j,k)

#if defined(HERMITIAN_FIELD)
               if (own_conjg) then
                  ! To be own conjugate, phase must be zero
                  ! (in order to be real)
                  complex_vec(1:NDIM) = abs(complex_vec(1:NDIM))
               end if
#endif
               
               ! Helmholtz decomposition!
               unitk_cmplx = cmplx(find_unitk(i,j,k,halfway))
               comp_cmplx(1:NDIM) = unitk_cmplx * &
                                    & dot_product(unitk_cmplx,complex_vec)
               sol_cmplx(1:NDIM) = complex_vec - comp_cmplx
               complex_vec = comp_cmplx*comp_frac + sol_cmplx*sol_frac
               ! note that there are two degrees of freedom for 
               ! solenoidal/transverse modes, and only one for
               ! longitudinal/compressive modes,
               ! and so for comp_frac == sol_frac == 0.5, you get
               ! F(solenoidal) == 2 * F(compressive)
               ! as per Federrath et al 2010
               turb_field(1:NDIM,i,j,k) = turb_field(1:NDIM,i,j,k) + &
                  & complex_vec(1:NDIM)

            end do
         end do
      end do
      
      return
   end subroutine add_turbulence
   
   subroutine decay_turbulence(old_turb_field, new_turb_field)
      ! Take a old complex Fourier field of turbulence, and find
      ! decay. Remove this decay from new turbulent field.
      complex(kind=PR), intent(in)    :: old_turb_field(1:NDIM, 0:TURB_GS-1,&
                                                   &0:TURB_GS-1, 0:TURB_GS-1)
                                              ! Old field for reference
      complex(kind=PR), intent(inout) :: new_turb_field(1:NDIM, 0:TURB_GS-1,&
                                                   &0:TURB_GS-1, 0:TURB_GS-1)
                                              ! Complex field to subtract from
   
     ! Ornstein-Uhlenbeck process
     ! dF(k,t) = F_0(k) P(k) dW - F(k,t) dt / T
     ! where F(k,t) is the vector fourier amplitude,
     ! F_0 is the power spectrum distribution of amplitudes,
     ! P(k) is the projection (Helmholtz decomposition),
     ! dt is the timestep of integration,
     ! and dW is the Wiener process, described by a random variate
     ! from a normal distribution with a mean of zero and a
     ! variance of dt i.e. a standard deviation of sqrt(dt)

     ! In this we only calculate the second term - the exponential
     ! decay term.
   
      new_turb_field = new_turb_field - (old_turb_field * turb_decay_frac)
      
   end subroutine decay_turbulence
   
   subroutine FFT_3D(complex_field, real_field)
      ! Transform complex field into purely real field for 3D vector field

      integer              :: d               ! Dimension counter
      integer (kind=ILP)   :: plan            ! FFTW plan pointer
      complex(kind=PR), intent(in)  :: complex_field(1:NDIM, 0:TURB_GS-1,&
                                                    &0:TURB_GS-1, 0:TURB_GS-1)
                                              ! Complex field to transform
      real(kind=PR), intent(out)    :: real_field(1:NDIM, 0:TURB_GS-1,&
                                                 &0:TURB_GS-1, 0:TURB_GS-1)
                                              ! Result of transforms
   
      complex(kind=PR), allocatable :: fftfield(:,:,:) ! Memory for FFT
      
      ! Allocate storage for performing FFTs
      allocate(fftfield(0:TURB_GS-1,0:TURB_GS-1,0:TURB_GS-1))

#if defined(QUADRUPLE_PRECISION)
      call fftwq_plan_dft_3d(plan, TURB_GS, TURB_GS, TURB_GS, fftfield, &
         & fftfield, FFTW_BACKWARD, FFTW_ESTIMATE)
#elif defined(DOUBLE_PRECISION)
      call dfftw_plan_dft_3d(plan, TURB_GS, TURB_GS, TURB_GS, fftfield, &
         & fftfield, FFTW_BACKWARD, FFTW_ESTIMATE)
#else
      call sfftw_plan_dft_3d(plan, TURB_GS, TURB_GS, TURB_GS, fftfield, &
         & fftfield, FFTW_BACKWARD, FFTW_ESTIMATE)
#endif
   
      !write (6,*) "Performing FFT..."
      do d=1,NDIM

         fftfield = complex_field(d,:,:,:)
      
#if defined(QUADRUPLE_PRECISION)
         call fftwq_execute_dft(plan, fftfield, fftfield)
#elif defined(DOUBLE_PRECISION)
         call dfftw_execute_dft(plan, fftfield, fftfield)
#else
         call sfftw_execute_dft(plan, fftfield, fftfield)
#endif
         real_field(d,:,:,:) = real(fftfield) / (turb_gs_real**NDIM)
      
         !write (6,*) "largest real component: ", maxval(abs(real(fftfield)))
         !write (6,*) "largest imaginary component: ", maxval(abs(aimag(fftfield)))
      
      end do
#if defined(QUADRUPLE_PRECISION)
      call fftwq_destroy_plan(plan)
#elif defined(DOUBLE_PRECISION)
      call dfftw_destroy_plan(plan)
#else
      call sfftw_destroy_plan(plan)
#endif
      
      deallocate(fftfield)
   
      return
   end subroutine FFT_3D
   
   function proj_rms_norm(sol_frac)
      ! Calculate the (empirically measured and fitted) normalization for
      ! projection of random vectors with solenoidal fraction 'sol_frac'
      real(kind=PR)                  :: proj_rms_norm ! Normalization constant
      real(kind=PR), intent(in)      :: sol_frac ! Solenoidal fraction
      
      proj_rms_norm = (0.797 * sol_frac**2) - (0.529 * sol_frac) + 0.568
      
      ! for reference, to maintain magnitude of vectors
      !proj_mag_norm = (0.563 * sol_frac**2) - (0.258 * sol_frac) + 0.487
      
      return
   end function proj_rms_norm
   
   function power_rms_norm(power_in)
      ! Calculate the normalizations for FFT of initial power spectrum
      real(kind=PR)             :: power_rms_norm  ! Normalization constant
      real(kind=PR), intent(in) :: power_in(0:TURB_GS-1,&
                                           &0:TURB_GS-1,0:TURB_GS-1)
      complex(kind=PR)          :: complex_field(1:NDIM, 0:TURB_GS-1,&
                                                &0:TURB_GS-1, 0:TURB_GS-1)
                                                   ! Complex field to transform
      real(kind=PR)             :: real_field(1:NDIM, 0:TURB_GS-1,&
                                             &0:TURB_GS-1, 0:TURB_GS-1)
                                                   ! Result of transforms
      integer                   :: d               ! Dimension counter
      
      do d=1,NDIM
         complex_field(d,:,:,:) = cmplx(power_in)
      end do
      
      call FFT_3D(complex_field, real_field)
      
      power_rms_norm = sqrt(sum(real_field**2)/size(real_field))
      
      return
   end function power_rms_norm

end module turb_force_module
