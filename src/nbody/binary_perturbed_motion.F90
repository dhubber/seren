! BINARY_PERTURBED_MOTION.F90
! D. A. Hubber - 9/11/2011
! Integrate motion of perturbed binary with 4th order Hermite scheme.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE binary_perturbed_motion(b)
  use sink_module
  use Nbody_module
  use time_module
  implicit none

  integer, intent(in) :: b                 ! Binary star id

  integer :: i                             ! ..
  integer :: s1                            ! ..
  integer :: s2                            ! ..
  real(kind=DP) :: adot(1:NDIM)            ! ..
  real(kind=DP) :: adot0(1:NDIM)           ! ..
  real(kind=DP) :: a2dot(1:NDIM)           ! ..
  real(kind=DP) :: a3dot(1:NDIM)           ! ..
  real(kind=DP) :: apert(1:NDIM)           ! ..
  real(kind=DP) :: arel(1:NDIM)            ! ..
  real(kind=DP) :: arel0(1:NDIM)           ! ..
  real(kind=DP) :: asqd                    ! Acceleration squared
  real(kind=DP) :: a1sqd                   ! jerk squared
  real(kind=DP) :: a2sqd                   ! 2nd derivative squared
  real(kind=DP) :: a3sqd                   ! 3rd derivative squared
  real(kind=DP) :: dt2                     ! dt*dt
  real(kind=DP) :: dt                      ! ..
  real(kind=DP) :: dt2                     ! ..
  real(kind=DP) :: dt3                     ! ..
  real(kind=DP) :: mbin                    ! ..
  real(kind=DP) :: rrel(1:NDIM)            ! ..
  real(kind=DP) :: vrel(1:NDIM)            ! ..
  real(kind=DP) :: rrel0(1:NDIM)           ! ..
  real(kind=DP) :: vrel0(1:NDIM)           ! ..
  real(kind=DP) :: taux                    ! ..
  real(kind=DP) :: tmax                    ! ..

  debug3("[binary_perturbed_motion.F90]")


! Record binary info
  s1 = binary(b)%s1
  s2 = binary(b)%s2
  rrel(1:NDIM) = binary(b)%dr(1:NDIM)
  vrel(1:NDIM) = binary(b)%dv(1:NDIM)
  mbin = binary(b)%m
  rbin(1:NDIM) = binary(b)%r(1:NDIM)
  apert(1:NDIM) = binary(b)%apert(1:NDIM)

! Calculate initial acceleration of each component
  drsqd = dot_product(rrel(1:NDIM),rrel(1:NDIM))
  drmag = sqrt(drsqd) + SMALL_NUMBER_DP
  invdrmag = 1.0_DP / drmag
  drdt = dot_product(vrel(1:NDIM),rrel(1:NDIM))*invdrmag
  arel(1:NDIM) = mbin*dr(1:NDIM)*invdrmag**3
  adot(1:NDIM) = mbin*dv(1:NDIM)*invdrmag**3 &
       & - 3.0_DP*mpp*drdt*dr(1:NDIM)*invdrmag**4

! ..
  rrel0(1:NDIM) = rrel(1:NDIM)
  vrel0(1:NDIM) = vrel(1:NDIM)
  arel0(1:NDIM) = arel(1:NDIM)
  adot0(1:NDIM) = adot(1:NDIM)

! Calculate initial timestep
  call nbody_timestep_size(s1,dt2)
  call nbody_timestep_size(s2,dt3)
  dt = min(dt2,dt3)
#if defined(DEBUG_BINARY_PERTURBED_MOTION)
  write(6,*) "Initial binary motion :",b
  write(6,*) "dt :",dt,"    Period :",binary(b)%period
#endif


! Main loop over internal integrated motion
! ============================================================================
  do

     ! Prediction step of motion
     rrel(1:NDIM) = rrel(1:NDIM) + vrel(1:NDIM)*dt + 0.5_DP*arel(1:NDIM)*dt2
     vrel(1:VDIM) = vrel(1:VDIM) + arel(1:VDIM)*dt

     ! P(EC)^n loop
     ! -----------------------------------------------------------------------
     do i=1,3

        ! Calculate acceleration of each component, including the perturbation
        drsqd = dot_product(rrel(1:NDIM),rrel(1:NDIM))
        drmag = sqrt(drsqd) + SMALL_NUMBER_DP
        invdrmag = 1.0_DP / drmag
        drdt = dot_product(vrel(1:NDIM),rrel(1:NDIM))*invdrmag
        arel(1:NDIM) = mbin*dr(1:NDIM)*invdrmag**3 + apert(1:NDIM)
        adot(1:NDIM) = mbin*dv(1:NDIM)*invdrmag**3 &
             & - 3.0_DP*mpp*drdt*dr(1:NDIM)*invdrmag**4
        
        ! Higher-order derivatives
        a2dot(1:NDIM) = (-6.0_DP*(arel0(1:NDIM) - arel(1:NDIM))&
             & - dt*(4.0_DP*adot0(1:NDIM) + 2.0_DP*adot(1:NDIM))) / dt2
        a3dot(1:NDIM) = (12.0_DP*(arel0(1:NDIM) - arel(1:NDIM))&
             & + 6.0_DP*dt*(adot0(1:NDIM) + adot(1:NDIM))) / dt3

        ! Correction term for motion
#if defined(NBODY_HERMITE4_TS)

#else        
        rrel(1:NDIM) = rrel(1:NDIM) + a2dot(1:NDIM)*dt2*dt2/24.0_DP &
             & + a3dot(1:NDIM)*dt3*dt2/120.0_DP
        vrel(1:NDIM) = vrel(1:NDIM) + a2dot0(1:NDIM)*dt3/6.0_DP &
             & + a3dot(1:NDIM)*dt2*dt2/24.0_DP
#endif

     end do
     ! -----------------------------------------------------------------------
     
     ! ..
     rrel0(1:NDIM) = rrel(1:NDIM)
     vrel0(1:NDIM) = vrel(1:NDIM)
     arel0(1:NDIM) = arel(1:NDIM)
     adot0(1:NDIM) = adot(1:NDIM)



     ! Add time and find out if we have reached end of complete orbit.
     ! If so, exit loop


     ! Calculate timestep
     asqd  = real(dot_product(arel(1:NDIM),arel(1:NDIM)),DP)
     a1sqd = real(dot_product(adot(1:NDIM),adot(1:NDIM)),DP)
     a2sqd = real(dot_product(a2dot(1:NDIM),a2dot(1:NDIM)),DP)
     a3sqd = real(dot_product(a3dot(1:NDIM),a3dot(1:NDIM)),DP)
     dt    = (sqrt(asqd*a2sqd) + a1sqd) / (sqrt(a1sqd*a3sqd) + a2sqd)
     dt    = nbody_timemult*sqrt(dt2)
     dt2   = dt*dt
     dt3   = dt2*dt

  end do
! ============================================================================


! Record binary properties
  binary(b)%dr(1:NDIM) = rrel(1:NDIM)
  binary(b)%dv(1:NDIM) = vrel(1:NDIM)

! Record star properties
  star(s1)%r(1:NDIM) = rbin(1:NDIM) + star(s1)%m*dr(1:NDIM)/mbin
  star(s2)%r(1:NDIM) = rbin(1:NDIM) - star(s2)%m*dr(1:NDIM)/mbin


  return
END SUBROUTINE binary_perturbed_motion
