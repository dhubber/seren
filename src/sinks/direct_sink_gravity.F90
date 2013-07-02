! SINK_GRAVITY.F90
! D. A. Hubber - 23/02/2010
! Calculates gravitational accelerations for particle p (or sink -p) 
! due to sinks only by direct sum.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE direct_sink_gravity(p,hp,rp,agravp,potp)
  use interface_module, only : gravity_gradh_meanh,gravity_nbody,gravity_meanh
  use definitions
  use particle_module
  use sink_module
  implicit none

  integer, intent(in) :: p                      ! Id of current particle
  real(kind=PR), intent(in) :: hp               ! Smoothing length of p
  real(kind=PR), intent(in) :: rp(1:NDIM)       ! Position of particle p
  real(kind=DP), intent(out) :: agravp(1:NDIM)  ! Gravitational accelertation
  real(kind=DP), intent(out) :: potp            ! Gravitational potential

  integer       :: s                 ! sink particle identifier
  integer       :: ss                ! sink particle identifier
  real(kind=PR) :: atemp(1:NDIM)     ! temp. grav acceleration var
  real(kind=PR) :: dpotp             ! grav potential of p due to pp
#if defined(GRAD_H_SPH) && defined(SELF_GRAVITY)
  real(kind=PR) :: zeta_p            ! grad-h zeta correction for particle p
  real(kind=PR) :: zeta_pp           ! grad-h zeta correction for particle pp
#endif

  debug3("Calculating gravitational force [sink_gravity.F90] for particle",p)

! Initialise gravitational acceleration and potential to zero
  agravp(1:NDIM) = 0.0_DP
  potp = 0.0_DP

! Set sink id if required
  s = -1
  if (p < 0) s = -p
#if defined(GRAD_H_SPH) && defined(SELF_GRAVITY)
  if (p < 0) then
     zeta_p = 0.0_PR
  else if (p <= ptot) then
     zeta_p = sph(p)%zo
  end if
  zeta_pp = 0.0_PR
#endif

! Add contributions due to all sinks
! ----------------------------------------------------------------------------
  do ss=1,stot
     if (s == ss) cycle
#if defined(N_BODY)
     call gravity_nbody(real(sink(ss)%m,PR),rp(1:NDIM),&
          &sink(s)%r(1:NDIM),atemp,dpotp)
#elif defined(GRAD_H_SPH) && defined(MEANH_GRAVITY)
     call gravity_meanh(0.5_PR*(hp + sink(ss)%h),real(sink(ss)%m,PR),&
          &rp(1:NDIM),sink(ss)%r(1:NDIM),atemp,dpotp)
#elif defined(GRAD_H_SPH)
     call gravity_sph(1.0_PR/hp,sink(ss)%h,real(sink(ss)%m,PR),&
          &rp(1:NDIM),sink(ss)%r(1:NDIM),atemp(1:NDIM),dpotp)
#elif defined(MEANH_GRAVITY)
     call gravity_meanh(0.5_PR*(hp + sink(ss)%h),real(sink(ss)%m,PR),&
          &rp(1:NDIM),sink(ss)%r(1:NDIM),atemp,dpotp)
#else
     call gravity_sph(1.0_PR/hp,sink(ss)%h,real(sink(ss)%m,PR),&
          &rp(1:NDIM),sink(ss)%r(1:NDIM),atemp(1:NDIM),dpotp)
#endif
     agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),DP)
     potp = potp + real(dpotp,DP)
  end do
! ----------------------------------------------------------------------------

  return
END SUBROUTINE direct_sink_gravity
