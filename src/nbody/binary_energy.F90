! BINARY_ENERGY.F90
! D. A. Hubber - 12/3/2008
! Calculates total energy (kinetic plus gravitational) between 
! a) two single stars, 
! b) a binary and a single star, 
! c) two binary stars
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE binary_energy(s1,s2,h1,h2,m1,m2,r1,r2,v1,v2,binen)
  use interface_module, only : distance3_dp,wpot
  use definitions
  use kernel_module
  use sink_module, only : stot
  implicit none

  integer,intent(in) :: s1                 ! id of system 1
  integer,intent(in) :: s2                 ! id of system 2
  real(kind=DP),intent(in) :: h1           ! Smoothing length of system 1
  real(kind=DP),intent(in) :: h2           ! Smoothing length of system 2
  real(kind=DP),intent(in) :: m1           ! Mass of system 1
  real(kind=DP),intent(in) :: m2           ! Mass of system 2
  real(kind=DP),intent(in) :: r1(1:NDIM)   ! Position of system 1
  real(kind=DP),intent(in) :: r2(1:NDIM)   ! Position of system 2
  real(kind=DP),intent(in) :: v1(1:NDIM)   ! Velocity of system 1
  real(kind=DP),intent(in) :: v2(1:NDIM)   ! Velocity of system 2
  real(kind=DP),intent(out) :: binen       ! 2-body energy
  
  real(kind=DP) :: dr(1:NDIM)              ! Relative displacement vector  
  real(kind=DP) :: drmag                   ! Distance
  real(kind=DP) :: drsqd                   ! Distance squared
  real(kind=DP) :: dv(1:NDIM)              ! Relative velocity vector
  real(kind=DP) :: dvsqd                   ! Relative velocity squared
  real(kind=DP) :: reduced_mass            ! Reduced mass

  call distance3_dp(r1(1:NDIM),r2(1:NDIM),dr(1:NDIM),drsqd)
  drmag = sqrt(drsqd) + SMALL_NUMBER_DP
  dv(1:NDIM) = v2(1:NDIM) - v1(1:NDIM)
  dvsqd = dot_product(dv(1:NDIM),dv(1:NDIM))
  reduced_mass = m1*m2 / (m1 + m2)

! Calculate total energy depending on whether this is standard gravity or 
! kernel-softened gravity. 
! ----------------------------------------------------------------------------
#ifndef N_BODY
! Only calculated kernel-softened gravity if both components are 
! stars (as opposed to other binary systems). 
  if (s1 <= stot .and. s2 <= stot .and. &
       & (drmag <= KERNRANGE*h1 .or. drmag <= KERNRANGE*h2)) then
     binen = 0.5_DP*(reduced_mass*dvsqd - &
          &m1*m2*(wpot(real(drmag/h1,PR))/h1 + wpot(real(drmag/h2,PR))/h2))
  else
     binen  = 0.5_DP*reduced_mass*dvsqd - m1*m2/drmag
  end if

#else
  binen = 0.5_DP*reduced_mass*dvsqd - m1*m2/drmag
#endif
! ----------------------------------------------------------------------------
  
  return
END SUBROUTINE binary_energy

