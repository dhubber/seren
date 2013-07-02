! GRAVITY_HERMITE4_GRADH_MEANH.F90
! D. A. Hubber - 28/6/2008
! Calculate gravitational acceleration and 1st time derivative of the 
! acceleration (the jerk) of particle p due to particle pp.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE gravity_hermite4_gradh_meanh(hp,hpp,mpp,rp,rpp,vp,vpp,&
                          &zo_p,zo_pp,atemp,adottemp,dpotp)
  use interface_module, only : distance3_dp,ewald_force,w0,wgrav,wpot,wzeta
  use definitions
  use kernel_module
  implicit none

  real(kind=PR), intent(in)  :: hp                ! Smoothing length of p
  real(kind=PR), intent(in)  :: hpp               ! Smoothing length of pp
  real(kind=DP), intent(in)  :: mpp               ! Mass of pp
  real(kind=DP), intent(in)  :: rp(1:NDIM)        ! Position of p
  real(kind=DP), intent(in)  :: rpp(1:NDIM)       ! Position of pp
  real(kind=DP), intent(in)  :: vp(1:NDIM)        ! Velocity of p
  real(kind=DP), intent(in)  :: vpp(1:NDIM)       ! Velocity of pp
  real(kind=PR), intent(in)  :: zo_p              ! Zeta/Omega for p
  real(kind=PR), intent(in)  :: zo_pp             ! Zeta/Omega for 
  real(kind=DP), intent(out) :: atemp(1:NDIM)     ! Acceleration vector
  real(kind=DP), intent(out) :: adottemp(1:NDIM)  ! Accel derivative
  real(kind=DP), intent(out) :: dpotp             ! Gravitational potential

  real(kind=DP) :: dr(1:NDIM)         ! Relative position vector
  real(kind=DP) :: drdt               ! Rate of change of distance
  real(kind=DP) :: drmag              ! Distance
  real(kind=DP) :: drsqd              ! Distance squared
  real(kind=DP) :: dv(1:NDIM)         ! Relative velocity vector
  real(kind=DP) :: grav               ! aux variable
  real(kind=DP) :: hmean              ! Mean smoothing length
  real(kind=DP) :: invdrmag           ! 1.0 / drmag
  real(kind=DP) :: invh               ! 1.0 / h
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
  hmean = 0.5_DP*(hp + hpp)

! Calculate contribution of gravity from particle
  if (drmag < real(KERNRANGE,DP)*hmean) then
     invh = 1.0_DP / hmean
     grav = invdrmag*invh*invh*real(wgrav(drmag*invh),DP)
     dpotp = invh*real(wpot(drmag*invh),DP)
     wmean = real(w0(drmag*invh),DP)*real(invh**(NDIMPR),DP)
  else
     grav = invdrmag**3
     dpotp = invdrmag
     wmean = 0.0_DP
  end if

! Calculate meanh grad-h correction term for grav. force for p
  if (drmag < KERNRANGE*hp) then
     invh = 1.0_PR / hp
     grav = grav + 0.5_PR*invdrmag*zo_p*wzeta(drmag*invh)*invh**4
  end if

! Calculate meanh grad-h correction term for grav. force for pp
  if (drmag < KERNRANGE*hpp) then
     invh = 1.0_PR / hpp
     grav = grav + 0.5_PR*invdrmag*zo_pp*wzeta(drmag*invh)*invh**4
  end if

! Record acceleration and jerk vectors
  atemp(1:NDIM) = mpp*grav*dr(1:NDIM)
  adottemp(1:NDIM) = mpp*grav*dv(1:NDIM) - &
       & 3.0_DP*mpp*drdt*grav*invdrmag*dr(1:NDIM) + &
       & 2.0_DP*TWOPI_DP*mpp*drdt*wmean*invdrmag*dr(1:NDIM)
  dpotp = mpp*dpotp
  
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
END SUBROUTINE gravity_hermite4_gradh_meanh
