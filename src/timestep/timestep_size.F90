! TIMESTEP_SIZE.F90
! C. P. Batty & D. A. Hubber - 29/3/2007
! Calculates ideal timestep for particle p based on multiple possible 
! criterion.  Calculates minimum of:
! 1. Acceleration timestep, dt = sqrt(h / accel)
! 2. Courant condition,     dt = h / (h*div_v + sound)
!    or (with viscosity)    dt = h / (sound + h*div_v_p + 
!                                     TVISC_FAC*(alpha*sound + beta*h*div_v_p)
! 3. Energy condition,      dt = u / (dudt + SMALL_NUMBER)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE timestep_size(p,dt)
  use definitions
  use type_module
  use hydro_module
  use sink_module
  use neighbour_module, only : hmin
  use particle_module, only : sph
  use time_module, only : accel_mult,courant_mult
  use periodic_module, only : periodic_half_minval
#if defined(TURBULENT_FORCING)
  use turbulence_module, only : turb_dt
#endif
  implicit none

  integer, intent(in)        :: p    ! Particle counter
  real(kind=DP), intent(out) :: dt   ! Step size for particle p

  real(kind=DP) :: amag              ! Magnitude of acceleration
  real(kind=DP) :: ap(1:NDIM)        ! Local copy of acceleration
  real(kind=DP) :: div_v_p           ! Local copy of velocity divergence
  real(kind=DP) :: hp                ! Local copy of smoothing length
  real(kind=DP) :: tacc              ! Acceleration timestep
  real(kind=DP) :: tcour             ! Courant time
#if defined(MHD)
  real(kind=DP) :: tmag              ! Magnetic courant time
#endif
  real(kind=DP) :: vsignal           ! Signal speed of particle p
#if defined(PERIODIC)
  real(kind=DP) :: vp(1:VDIM)        ! Local copy of velocity
  real(kind=DP) :: vmag              ! Magnitude of velocity
  real(kind=DP) :: tperiodic         ! Time to cross smallest 1/2 of box
#endif
#if defined(ARTIFICIAL_VISCOSITY)
  real(kind=DP) :: alpha_p           ! Local copy of alpha visc. value 
  real(kind=DP) :: beta_p            ! Local copy of beta visc. value
#endif
#if defined(ENERGY_EQN) && defined(EXPLICIT_COOLING_HEATING) && !defined(RAD_WS)
  real(kind=DP) :: tenergy           ! Energy time
#endif
#if defined(SINKS)
  integer :: s                       ! Sink counter
  real(kind=DP) :: drsqd             ! Distance squared
  real(kind=DP) :: dr(1:NDIM)        ! Relative position vector
#endif

  debug3("Calculating timestep [timestep_size.F90] for particle ", p)

! Make local copies of important particle properties
  ap(1:NDIM) = real(sph(p)%a(1:NDIM),DP)
  hp = real(sph(p)%h,DP)
#if defined(PERIODIC)
  vp(1:VDIM) = real(sph(p)%v(1:VDIM),DP)
  vmag = sqrt(dot_product(vp(1:VDIM),vp(1:VDIM)))
#endif
  div_v_p = real(abs(sph(p)%div_v),DP)
#if defined(HYDRO)
  vsignal = real(sph(p)%sound,DP)
#if defined(ARTIFICIAL_VISCOSITY) && defined(VISC_TD) && defined(VISC_BALSARA)
  alpha_p = real(sph(p)%talpha,DP)*real(sph(p)%balsara,DP)
  beta_p  = 2.0_DP*alpha_p
#elif defined(ARTIFICIAL_VISCOSITY) && defined(VISC_TD)
  alpha_p = real(sph(p)%talpha,DP)
  beta_p  = 2.0_DP*alpha_p
#elif defined(ARTIFICIAL_VISCOSITY) && defined(VISC_BALSARA)
  alpha_p = real(alpha*sph(p)%balsara,DP)
  beta_p  = real(beta*sph(p)%balsara,DP)
#elif defined(ARTIFICIAL_VISCOSITY) && defined(VISC_PATTERN_REC)
  alpha_p = real(alpha*sph(p)%pattrec,DP)
  beta_p = real(beta*sph(p)%pattrec,DP)
#elif defined(ARTIFICIAL_VISCOSITY) 
  alpha_p = real(alpha,DP)
  beta_p  = real(beta,DP)
#endif
#endif
  amag = sqrt(dot_product(ap(1:NDIM),ap(1:NDIM)))
  if (p > phydroend) vsignal = SMALL_NUMBER_DP

! Acceleration condition on timestep (Always calculated)
  tacc = accel_mult * sqrt(hp / (amag + SMALL_NUMBER_DP))
  dt = tacc

! Courant condition (with or without artificial viscosity), or without 
! gravity (no sound speed, but velocity divergence).
#if defined(HYDRO) && defined(ARTIFICIAL_VISCOSITY)
  tcour = courant_mult * hp / (vsignal + hp*div_v_p + &
       & TVISC_FAC*(alpha_p*vsignal + beta_p*hp*div_v_p))
  dt = min(dt,tcour)
#elif defined(HYDRO) && !defined(ARTIFICIAL_VISCOSITY)
  tcour = courant_mult * hp / (vsignal + hp*div_v_p)
  dt = min(dt,tcour)
#elif !defined(HYDRO)
  tcour = courant_mult / (div_v_p + SMALL_NUMBER)
  dt = min(dt,tcour)
#endif

! Magnetic courant time
#if defined(MHD)
  tmag = courant_mult * sph(p)%B_t_signal
  dt = min(dt,tmag)
#endif

! Periodic boundary timestep - prevents particles crossing more than
! one half of a box length in a single timestep
#if defined(PERIODIC)
  tperiodic = periodic_half_minval / (vmag + SMALL_NUMBER)
  dt = min(dt,tperiodic)
#endif

! Internal energy timestep
#if defined(HYDRO) && defined(EXPLICIT_COOLING_HEATING) && !defined(RAD_WS)
  if (p < phydroend) then
     tenergy  = min(accel_mult,courant_mult) * real(sph(p)%u,DP) / &
          & (real(abs(sph(p)%dudt),DP) + SMALL_NUMBER_DP)
     dt = min(dt,tenergy)
  end if
#endif

! If particle is inside a sink, set timestep equal to sink timestep
#if defined(SINKS) && defined(SMOOTH_ACCRETION)
  do s=1,stot
     call distance2(sink(s)%r(1:NDIM),p,dr(1:NDIM),drsqd)
     if (drsqd <= sink(s)%radius*sink(s)%radius) then
        dt = 0.4_DP*sink(s)%h / (vsignal + sink(s)%h*div_v_p + &
             & TVISC_FAC*(alpha_p*vsignal + beta_p*sink(s)%h*div_v_p))
     end if
  end do
#endif

! For now, particles must have timesteps shorter than turbulent update time
#if defined(TURBULENT_FORCING)
  dt = min(dt, turb_dt)
#endif

#if defined(DEBUG_TIMESTEP_SIZE)
  if (dt < 0.0_PR) then
     write(6,*) "Minimum timestep,      dt        : ",dt,p
     write(6,*) "Acceleration timestep, tacc      : ",tacc,amag
#if defined(PERIODIC)
     write(6,*) "Periodic timestep,     tperiodic : ", tperiodic
#endif
#if defined(HYDRO)
     write(6,*) "Courant timestep,      tcour     : ",tcour
#if defined(HYDRO) && defined(EXPLICIT_COOLING_HEATING) && !defined(RAD_WS)
     write(6,*) "Energy timestep,       tenergy   : ",tenergy,up,dudt_p
#endif
#endif
     stop 'Negative timestep!!'
  end if
#endif

  return
END SUBROUTINE timestep_size
