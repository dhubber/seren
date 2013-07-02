! GRAVITY_HERMITE4.F90
! D. A. Hubber - 28/6/2008
! Calculate gravitational acceleration and 1st time derivative of the 
! acceleration (the jerk) of particle p due to particle pp.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE gravity_hermite4(invhp,hpp,mpp,rp,rpp,vp,vpp,&
                          &atemp,adottemp,dpotp)
  use interface_module, only : distance3_dp,ewald_force,w0,wgrav,wpot
  use definitions
  use kernel_module
  implicit none

  real(kind=DP), intent(in)  :: hpp               ! Smoothing length of pp
  real(kind=DP), intent(in)  :: invhp             ! Smoothing length of p
  real(kind=DP), intent(in)  :: mpp               ! Mass of pp
  real(kind=DP), intent(in)  :: rp(1:NDIM)        ! Position of p
  real(kind=DP), intent(in)  :: rpp(1:NDIM)       ! Position of pp
  real(kind=DP), intent(in)  :: vp(1:NDIM)        ! Velocity of p
  real(kind=DP), intent(in)  :: vpp(1:NDIM)       ! Velocity of pp
  real(kind=DP), intent(out) :: atemp(1:NDIM)     ! Acceleration vector
  real(kind=DP), intent(out) :: adottemp(1:NDIM)  ! Accel derivative
  real(kind=DP), intent(out) :: dpotp             ! Gravitational potential

  real(kind=DP) :: dr(1:NDIM)         ! Relative position vector
  real(kind=DP) :: drdt               ! Rate of change of distance
  real(kind=DP) :: drmag              ! Distance
  real(kind=DP) :: drsqd              ! Distance squared
  real(kind=DP) :: dv(1:NDIM)         ! Relative velocity vector
  real(kind=DP) :: grav               ! aux variable
  real(kind=DP) :: invdrmag           ! 1.0 / drmag
  real(kind=DP) :: invhpp             ! 1.0 / hpp
  real(kind=DP) :: wmean              ! wmean
#if defined(EWALD)
  real(kind=DP) :: eaccel(1:NDIM)     ! Ewald grav. acceleration
#endif

! Calculate relative displacement vector between bodies
#if defined(EWALD)
  call distance3_dp(rp(1:NDIM),rpp(1:NDIM),dr(1:NDIM),drsqd)
#else
  dr(1:NDIM) = rpp(1:NDIM) - rp(1:NDIM)
  drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif
  drmag = sqrt(drsqd) + SMALL_NUMBER_DP
  invdrmag = 1.0_DP / drmag
  dv(1:NDIM) = vpp(1:NDIM) - vp(1:NDIM)
  drdt = dot_product(dv(1:NDIM),dr(1:NDIM))*invdrmag

! Calculate contribution of gravity from particle p
  if (invdrmag > real(INVKERNRANGE,DP)*invhp) then
     grav = invdrmag*invhp*invhp*real(wgrav(drmag*invhp),DP)
     dpotp = invhp*real(wpot(drmag*invhp),DP)
     wmean = real(w0(drmag*invhp)*invhp**(NDIM),DP)
  else
     grav = invdrmag**3
     dpotp = invdrmag
     wmean = 0.0_DP
  end if

! Calculate contribution of gravity from particle pp
  if (drmag < real(KERNRANGE,DP)*hpp) then
     invhpp = 1.0_DP / hpp
     grav = grav + invdrmag*invhpp*invhpp*real(wgrav(drmag*invhpp),DP)
     dpotp = dpotp + invhpp*real(wpot(drmag*invhpp),DP)
     wmean = wmean + real(w0(drmag*invhpp)*invhpp**(NDIM),DP)
  else
     grav = grav + invdrmag**3
     dpotp = dpotp + invdrmag
  end if

! Record acceleration in output vector      
  atemp(1:NDIM) = 0.5_DP*mpp*grav*dr(1:NDIM)
  adottemp(1:NDIM) = 0.5_DP*mpp*grav*dv(1:NDIM) - &
       & 1.5_DP*mpp*drdt*grav*invdrmag*dr(1:NDIM) + &
       & TWOPI_DP*mpp*drdt*wmean*invdrmag*dr(1:NDIM)
  dpotp = 0.5_DP*mpp*dpotp
  
! Add Ewald periodic correction if required
#if defined(EWALD)
  call ewald_force(dr,mpp,eaccel)
  atemp(1:NDIM) = atemp(1:NDIM) + eaccel(1:NDIM)
#endif

#if defined(DEBUG_HERMITE4)
  write(6,*) "grav_hermite : ",grav,dpotp,drdt,drdt
  write(6,*) "rp           : ",rp(1:NDIM),vp(1:NDIM)
  write(6,*) "rpp          : ",rpp(1:NDIM),vpp(1:NDIM),hpp,mpp
  write(6,*) "dr           : ",dr(1:NDIM),drmag,invdrmag,invhp
#endif

  return
END SUBROUTINE gravity_hermite4
