! NBODY_HERMITE4_EXTRA_TERMS.F90
! D. A. Hubber - 23/6/2008
! Calculate the 3rd and 4th order acceleration time-derivatives explicitly 
! for the first step only in the Hermite integration scheme for all stars.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_hermite4_extra_terms
  use particle_module
  use sink_module
  use kernel_module
  use time_module
  use Nbody_module
  use type_module
  use scaling_module, only : tscale
  implicit none

  integer :: ndead                    ! Number of newly accreted particles
  integer :: p                        ! Particle counter
  integer :: s                        ! Sink counter
  integer :: ss                       ! Secondary sink counter
  integer :: sinkid_p(1:SMAX)         ! Auxilary sorting arrays
  integer, allocatable :: deadlist(:) ! List of dead (accreted) particles
  integer, allocatable :: sinkid(:)   ! Sink that particle p is most bound to
  real(kind=DP) :: adottemp(1:NDIM)   ! Aux. variable for calculating jerk
  real(kind=DP) :: a2dot(1:NDIM)      ! 2nd derivative of acceleration
  real(kind=DP) :: a2dottemp(1:NDIM)  ! Aux var for calculating 2nd derivative
  real(kind=DP) :: a3dot(1:NDIM)      ! 3rd derivative of acceleration
  real(kind=DP) :: atemp(1:NDIM)      ! Auxilary accel variable 
  real(kind=PR) :: atemp_pr(1:NDIM)   ! ..
  real(kind=DP) :: afac               ! Aux var. for 2nd/3rd derivative calc
  real(kind=DP) :: bfac               !  "   "
  real(kind=DP) :: cfac               !  "   "
  real(kind=PR) :: dpotp              ! Aux. sink pot variable
  real(kind=DP) :: dr(1:NDIM)         ! Relative position vector (DP)
  real(kind=DP) :: drsqd              ! Distance squared (DP)
  real(kind=PR) :: dv(1:NDIM)         ! Relative velocity vector
  real(kind=DP) :: dv2(1:NDIM)        ! Relative velocity vector (DP)
  real(kind=PR) :: energy_p(1:SMAX)   ! Energy of p-s system
  real(kind=PR) :: invhs              ! Smoothing length of sink s
  real(kind=PR) :: invhp              ! 1 / hp
  real(kind=DP) :: invdrmag           ! 1 / drmag
  real(kind=DP) :: invdrsqd           ! 1 / drsqd
  real(kind=PR) :: mp                 ! Mass of particle p
  real(kind=PR) :: ms                 ! Mass of sink s
  real(kind=DP) :: potp               ! Potential of sink s
  real(kind=PR) :: reduced_mass       ! Reduced mass of sink-particle system
  real(kind=PR) :: rp(1:NDIM)         ! Position of particle p
  real(kind=PR) :: rs(1:NDIM)         ! Position of sink s
  real(kind=DP) :: rs2(1:NDIM)        ! Position of sink s (DP)
  real(kind=DP) :: rss(1:NDIM)        ! Position of sink ss
  real(kind=PR) :: vp(1:NDIM)         ! Velocity of particle p
  real(kind=PR) :: vs(1:NDIM)         ! Velocity of sink s
  real(kind=DP) :: vsqd               ! Velocity squared (DP)

  debug2("[nbody_hermite4_extra_terms.F90]")

! Loop over all stars
! ----------------------------------------------------------------------------
  do s=1,stot

     ! Zero arrays
     a2dot(1:NDIM) = 0.0_DP
     a3dot(1:NDIM) = 0.0_DP
     potp = 0.0_DP

     ! Store local copies of sink position, mass and smoothing length 
     rs2(1:NDIM) = star(s)%r(1:NDIM)
     invhs = 1.0_DP / star(s)%h

     ! Loop over all other sinks and calculate net gravitational acceleration
     ! -----------------------------------------------------------------------
     do ss=1,stot
        if (s == ss) cycle
        rss(1:NDIM) = star(ss)%r(1:NDIM)
#if defined(MEANH_GRAVITY)
        call gravity_hermite4_meanh(0.5_DP*(star(s)%h + star(ss)%h),&
             &star(ss)%m,star(s)%r(1:NDIM),star(ss)%r(1:NDIM),&
             &star(s)%v(1:NDIM),star(ss)%v(1:NDIM),atemp(1:NDIM),&
             &adottemp(1:NDIM),dpotp)
#else
        call gravity_hermite4(invhs,star(ss)%h,star(ss)%m,&
             &rs(1:NDIM),rss(1:NDIM),vs(1:NDIM),&
             &star(ss)%v(1:NDIM),atemp(1:NDIM),adottemp(1:NDIM),dpotp)
#endif
        call distance3_dp(rs2(1:NDIM),rss(1:NDIM),dr(1:NDIM),drsqd)
        invdrsqd = 1.0_DP/drsqd
        invdrmag = sqrt(invdrsqd)
        dv2(1:NDIM) = star(ss)%v(1:NDIM) - star(s)%v(1:NDIM)
        vsqd = dot_product(dv2(1:NDIM),dv2(1:NDIM))
        afac = dot_product(dv2(1:NDIM),dr(1:NDIM))*invdrsqd
        bfac = vsqd*invdrsqd + afac*afac + invdrsqd*&
             & dot_product(dr(1:NDIM),star(ss)%a(1:NDIM) - star(s)%a(1:NDIM))
        cfac = 3.0_DP*invdrsqd*dot_product(dv2(1:NDIM),&
             & star(ss)%a(1:NDIM) - star(s)%a(1:NDIM)) + &
             & invdrsqd*dot_product(dr(1:NDIM),&
             & star(ss)%adot(1:NDIM) - star(s)%adot(1:NDIM)) + &
             & afac*(3.0_DP*bfac - 4.0_DP*afac*afac)
        a2dottemp(1:NDIM) = -star(ss)%m*(star(s)%a(1:NDIM) - &
             & star(ss)%a(1:NDIM))*invdrmag*invdrsqd - &
             & 6.0_DP*afac*adottemp(1:NDIM) - 3.0_DP*bfac*atemp(1:NDIM)
        a2dot(1:NDIM) = a2dot(1:NDIM) + a2dottemp(1:NDIM)
        a3dot(1:NDIM) = a3dot(1:NDIM) - star(ss)%m*(star(s)%adot(1:NDIM) - &
             & star(ss)%adot(1:NDIM))*invdrmag*invdrsqd - &
             & 9.0_DP*afac*a2dottemp(1:NDIM) - 9.0_DP*bfac*adottemp(1:NDIM) &
             & - 3.0_DP*cfac*atemp(1:NDIM)
     end do
     ! -----------------------------------------------------------------------

     ! Record acceleration and potential energy in main array
     star(s)%a0(1:NDIM)     = star(s)%a(1:NDIM)
     star(s)%adot0(1:NDIM)  = star(s)%adot(1:NDIM)
     star(s)%a2dot(1:NDIM)  = a2dot(1:NDIM)
     star(s)%a2dot0(1:NDIM) = a2dot(1:NDIM)
     star(s)%a3dot(1:NDIM)  = a3dot(1:NDIM)
#if defined(FORCE_SPLITTING)
     star(s)%a2dotstar(1:NDIM) = a2dot(1:NDIM)
     star(s)%a3dotstar(1:NDIM) = a3dot(1:NDIM)
#endif

#if defined(DEBUG_HERMITE4)
     write(6,*) "Higher order terms for star :",s
     write(6,*) "a     : ",star(s)%a(1:NDIM)
     write(6,*) "adot  : ",star(s)%adot(1:NDIM)
     write(6,*) "a2dot : ",star(s)%a2dot(1:NDIM)
     write(6,*) "a3dot : ",star(s)%a3dot(1:NDIM)
#endif

  end do
! ----------------------------------------------------------------------------


  return
END SUBROUTINE nbody_hermite4_extra_terms
