! NBODY_BINARY_TIMESTEP_SIZE.F90
! D. A. Hubber - 9/11/2011
! Calculate appropriate N-body timestep for binary b.
! Currently uses Aarseth timestep for 4th-order Hermite scheme.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_binary_timestep_size(s,dt)
  use Nbody_module
  implicit none

  integer, intent(in) :: b             ! Binary id
  real(kind=DP), intent(out) :: dt     ! Physical timestep

  integer :: s1                        ! ..
  integer :: s2                        ! ..
  real(kind=DP) :: asqd                ! Acceleration squared
  real(kind=DP) :: a1sqd               ! jerk squared
  real(kind=DP) :: a2sqd               ! 2nd derivative squared
  real(kind=DP) :: a3sqd               ! 3rd derivative squared
  real(kind=DP) :: dt2                 ! dt*dt
  real(kind=DP) :: acom(1:NDIM)        ! ..
  real(kind=DP) :: adotcom(1:NDIM)     ! ..
  real(kind=DP) :: a2dotcom(1:NDIM)    ! ..
  real(kind=DP) :: a3dotcom(1:NDIM)    ! ..

  s1 = binary(b)%s1
  s2 = binary(b)%s2

! Timestep for 4-th order Hermite scheme
! ----------------------------------------------------------------------------
#if defined(NBODY_HERMITE4)

  ! Calculate various acceleration derivatives for COM
  acom(1:NDIM) = (star(s1)%m*star(s1)%a(1:NDIM) + &
       & star(s2)%m*star(s2)%a(1:NDIM)) / binary(b)%m
  adotcom(1:NDIM) = (star(s1)%m*star(s1)%adot(1:NDIM) + &
       & star(s2)%m*star(s2)%adot(1:NDIM)) / binary(b)%m
  a2dotcom(1:NDIM) = (star(s1)%m*star(s1)%a2dot(1:NDIM) + &
       & star(s2)%m*star(s2)%a2dot(1:NDIM)) / binary(b)%m
  a3dotcom(1:NDIM) = (star(s1)%m*star(s1)%a3dot(1:NDIM) + &
       & star(s2)%m*star(s2)%a3dot(1:NDIM)) / binary(b)%m

  ! Now calculate Aarseth timestep for COM
  asqd  = real(dot_product(acom,acom),DP)
  a1sqd = real(dot_product(adotcom,adotcom),DP)
  a2sqd = real(dot_product(a2dotcom,a2dotcom),DP)
  a3sqd = real(dot_product(a3dotcom,a3dotcom),DP)
  dt2   = (sqrt(asqd*a2sqd) + a1sqd) / (sqrt(a1sqd*a3sqd) + a2sqd)
  dt    = nbody_timemult*sqrt(dt2)

  ! Ensure timestep of COM is not greater than the period of the binary
  dt = min(dt,binary(b)%period)

#if defined(DEBUG_TIMESTEP_SIZE)
  write(6,*) "Timestep for binary ",b,dt,binary(b)%period
  write(6,*) "asqd : ",asqd,a1sqd,a2sqd,a3sqd
#endif

! ----------------------------------------------------------------------------
#elif defined(NBODY_HERMITE6)

#endif
! ----------------------------------------------------------------------------

  return
END SUBROUTINE nbody_binary_timestep_size
