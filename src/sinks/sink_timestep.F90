! SINK_TIMESTEP.F90
! D. A. Hubber - 23/2/2007
! Calculates optimum timestep for sink particle s
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sink_timestep(s,dt_ideal)
  use interface_module, only : distance3_dp
  use sink_module
  use time_module, only : sink_mult
  implicit none

  integer, intent(in) :: s                ! Sink identifier
  real(kind=DP), intent(out) :: dt_ideal  ! Min. of computed timesteps

  integer :: ss                           ! Secondary sink loop counter
  real(kind=DP) :: acc_mag                ! Magnitude of acceleration
  real(kind=DP) :: as(1:NDIM)             ! acceleration of sink s
  real(kind=DP) :: dist                   ! Distance
  real(kind=DP) :: dr(1:NDIM)             ! Displacement
  real(kind=DP) :: drsqd                  ! Distance squared
  real(kind=DP) :: dt1                    ! rads/acc_mag timestep condition
  real(kind=DP) :: dt2                    ! rads/vel_mag timestep condition
  real(kind=DP) :: dv(1:VDIM)             ! Relative velocity
  real(kind=DP) :: rads                   ! Radius of sink s
  real(kind=DP) :: rs(1:NDIM)             ! Position of sink s
  real(kind=DP) :: rss(1:NDIM)            ! Position of sink ss
  real(kind=DP) :: vel_mag                ! Magnitude of velocity
  real(kind=DP) :: vs(1:NDIM)             ! velocity of sink s

  rads       = real(sink(s)%radius,DP)
  rs(1:NDIM) = real(sink(s)%r(1:NDIM),DP)
  as(1:NDIM) = real(sink(s)%a(1:NDIM),DP)
  vs(1:NDIM) = real(sink(s)%v(1:NDIM),DP)
  acc_mag    = sqrt(dot_product(as(1:NDIM),as(1:NDIM)))
  vel_mag    = sqrt(dot_product(vs(1:NDIM),vs(1:NDIM)))

  dt1 = sink_mult*sqrt(rads/(acc_mag + SMALL_NUMBER_DP))
  dt2 = sink_mult*rads/(vel_mag + SMALL_NUMBER_DP)
  dt_ideal = min(dt1,dt2)

  if (stot > 1) then
     dt2 = BIG_NUMBER_DP
     do ss=1,stot
        if (s == ss) cycle
        rss(1:NDIM) = real(sink(ss)%r(1:NDIM),DP)
        call distance3_dp(rs(1:NDIM),rss(1:NDIM),dr(1:NDIM),drsqd)
        dist = max(0.5_DP*(rads + real(sink(ss)%radius,DP)),sqrt(drsqd))
        dv(1:NDIM) = vs(1:VDIM) - real(sink(ss)%v(1:VDIM),DP)
        vel_mag = sqrt(dot_product(dv(1:VDIM),dv(1:VDIM)))
        dt2 = sink_mult*dist/(vel_mag + SMALL_NUMBER)
        dt_ideal = min(dt_ideal,dt2)
     end do
  end if

#if defined(DEBUG_SINK_TIMESTEP)
  write(6,*) "sink timestep :",s,rads,acc_mag,vel_mag,dt1,dt2,dt_ideal
#endif  

  return
END SUBROUTINE sink_timestep
