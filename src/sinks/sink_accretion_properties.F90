! SINK_ACCRETION_PROPERTIES.F90
! D. A. Hubber - 8/12/2008
! Calculates an approximation of the accretion rate.  Next calculates 
! the properties of the unresolved protostar in the sink using simple 
! model equations (e.g. Stamatellos, Whitworth, Boyd & Goodwin 2005).
! N.B. UNITS AND SCALING NOT CORRECT YET
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sink_accretion_properties(s,maccreted)
  use Eos_module,only : rad_const
  use sink_module
  use time_module
  use scaling_module
  use constant_module
  implicit none

  integer, intent(in) :: s                 ! sink id
  real(kind=DP), intent(in) :: maccreted   ! mass accreted in current step

  integer :: i                          ! Aux. counter
  real(kind=DP) :: mtemp(1:DMDT_RANGE)  ! Mass aux. variable
  real(kind=DP) :: summ                 ! m summation (for least-squares)
  real(kind=DP) :: sumt                 ! t summation (for least-squares)
  real(kind=DP) :: sumtm                ! m*t summation (for least-squares)
  real(kind=DP) :: sumtsqd              ! t*t summation (for least-squares)
  real(kind=DP) :: ttemp(1:DMDT_RANGE)  ! Time aux. variables

! If the sink has accreted any mass this timestep, shuffle all values down 
! one array element (therefore removing the oldest)
!  if (maccreted > SMALL_NUMBER_DP) then
     do i=DMDT_RANGE-1,1,-1
        sink(s)%macc(i+1) = sink(s)%macc(i)
        sink(s)%tacc(i+1) = sink(s)%tacc(i)
     end do
     if (sink(s)%ncreate /= nsteps) sink(s)%macc(1) = maccreted
!  end if
  sink(s)%tacc(1) = time


! Calculate least-squares fit to accretion rate
  mtemp(1:DMDT_RANGE) = 0.0_DP
  ttemp(1:DMDT_RANGE) = 0.0_DP
  sumt    = 0.0_DP
  sumtsqd = 0.0_DP
  sumtm   = 0.0_DP
  summ    = 0.0_DP
  do i=DMDT_RANGE-1,1,-1
     mtemp(i) = mtemp(i+1) + sink(s)%macc(i)
     ttemp(i) = 0.5*(sink(s)%tacc(i+1) + sink(s)%tacc(i)) - &
          &sink(s)%tacc(DMDT_RANGE)
     sumt     = sumt + ttemp(i)
     sumtsqd  = sumtsqd + ttemp(i)**2
     sumtm    = sumtm + ttemp(i)*mtemp(i)
     summ     = summ + mtemp(i)
  end do
  sink(s)%dmdt = (real(DMDT_RANGE-1,DP)*sumtm - summ*sumt) / &
       &(real(DMDT_RANGE-1,DP)*sumtsqd - sumt*sumt)
  if (sink(s)%dmdt < 0.0_DP) sink(s)%dmdt = 0.0_DP


! Calculate properties of unresolved protostar
  sink(s)%star_radius = 4.0_DP*r_sun / (rscale*r_SI)  !star_radius

!! luminosity of protostar 
!#if defined(EPISODIC_ACCRETION)
!  call episodic_accretion_model(s)
!#else
!  sink(s)%Mstar = sink(s)%m
!  sink(s)%dmdt_star = sink(s)%dmdt
!#endif

! Calculate properties of unresolved protostar 
  if (time - sink(s)%tcreate > feedback_tdelay .and. &
       &real(sink(s)%m,DP)>feedback_minmass) then
     !sink(s)%luminosity_old = sink(s)%luminosity
     sink(s)%luminosity = &
          & (real(sink(s)%m,DP)*mscale*m_SI/m_sun)**3/Lscale + &
          & f_accretion*(1.0_PR - sink(s)%star_radius/(2.0_PR*sink(s)%radius))&
          & *sink(s)%dmdt*real(sink(s)%m,DP)/sink(s)%star_radius
     sink(s)%temperature = (sink(s)%luminosity / (4.0_DP*rad_const*PI*&
          & sink(s)%star_radius*sink(s)%star_radius))**(0.25_DP) 
!#if defined(STAR_SIMPLE_HEATING)
!     if (sink(s)%luminosity > 0.0_PR .and. abs(sink(s)%luminosity - &
!          &sink(s)%luminosity_old)/sink(s)%luminosity > 0.01_PR) then
!        sync_flag = .true.
!	sync_steps = 0
!     end if
!#endif
  end if

  return
END SUBROUTINE sink_accretion_properties
