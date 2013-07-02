! RECORD_PARTICLE_DATA.F90
! D. A. Hubber - 06/03/2010
! Records and returns all float variables for particle p in a single array 
! (alldata) relative to some specified origin (rorigin).  Also returns no. 
! of filled elements in array (nelements)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE record_particle_data(p,nelements,alldata,rorigin)
  use particle_module
  use hydro_module
  use neighbour_module
  use scaling_module
  use filename_module
  use time_module
  use type_module
  use HP_module
#if defined(RAD_WS)
  use Eos_module
#endif
  implicit none

  integer, intent(in) :: p                      ! I.d. of particle
  integer, intent(out) :: nelements             ! No. of filled elements
  real(kind=PR), intent(out) :: alldata(1:100)  ! Particle data array
  real(kind=PR), intent(in) :: rorigin(1:NDIM)  ! Origin for output

  integer :: k                      ! Dimension counter
  real(kind=PR) :: a_rad            ! magnitude of acceleration
  real(kind=PR) :: dr(1:NDIM)       ! relative displacement vector
  real(kind=PR) :: drmag            ! Magnitude of distance
  real(kind=PR) :: drsqd            ! Distance squared
  real(kind=PR) :: dr_unit(1:NDIM)  ! Unit vector
  real(kind=DP) :: dt               ! timestep
  real(kind=PR) :: dv_unit(1:NDIM)  ! Unit vector in direction of velocity
  real(kind=PR) :: rp(1:NDIM)       ! Particle position
  real(kind=PR) :: vmag             ! Magnitude of velocity
  real(kind=PR) :: v_rad            ! Radial component of velocity
  real(kind=PR) :: vtemp(1:NDIM)    ! temp vector for scalar product call
#if defined(DEBUG_FORCES)
  real(kind=PR) :: ag_rad           ! radial gravitational acceleration
  real(kind=PR) :: ah_rad           ! radial hydro acceleration
  real(kind=PR) :: av_rad           ! radial viscous acceleration
  real(kind=PR) :: av_vel           ! viscous accel in direction of velocity
  real(kind=PR) :: aw_rad           ! radial wind acceleration
#endif

  nelements      = 0
  alldata(1:100) = 0.0_PR
  rp = sph(p)%r
#if defined(USE_MPI) && defined(PERIODIC)
  call unwrap_particle_position(rp)
#endif
! Calculate unit vector here
  call distance3(rorigin(1:NDIM),rp,dr(1:NDIM),drsqd)
  drmag = sqrt(drsqd) + SMALL_NUMBER
  if (drmag < SMALL_NUMBER) then
     dr_unit(1:NDIM) = 0.0_PR
  else
     dr_unit(1:NDIM) = dr(1:NDIM) / drmag
  end if
  vmag = sqrt(dot_product(sph(p)%v(1:VDIM),sph(p)%v(1:VDIM))) + SMALL_NUMBER
  if (vmag < SMALL_NUMBER) then
     dv_unit(1:NDIM) = 0.0_PR
  else
     dv_unit(1:NDIM) = sph(p)%v(1:NDIM) / vmag
  end if
  
! Positions (always written)
  do k=1,NDIM
     nelements = nelements + 1
     alldata(nelements) = rp(k)*real(rscale,PR)
  end do
  
! Velocities (always written)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = sph(p)%v(k)*real(vscale,PR)
     vtemp(k) = alldata(nelements)
  end do
  v_rad = dot_product(vtemp,dr_unit)

! Smoothed velocities
#if defined(SMOOTHED_VELOCITY)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = v_smooth(k,p)*real(vscale,PR)
     vtemp(k) = alldata(nelements)
  end do
#endif
  
! Accelerations (always written)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = sph(p)%a(k)*real(ascale,PR)
     vtemp(k) = alldata(nelements)
  end do
  a_rad = dot_product(vtemp,dr_unit)
  
! Hydro accelerations
#if defined(DEBUG_FORCES)
#if defined(HYDRO)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = sph(p)%a_hydro(k)*real(ascale,PR)
     vtemp(k) = alldata(nelements)
  end do
  ah_rad = dot_product(vtemp,dr_unit)
#endif
  
! Gravitational acceleration
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = sph(p)%a_grav(k)*real(ascale,PR)
     vtemp(k) = alldata(nelements)
  end do
  ag_rad = dot_product(vtemp,dr_unit)
#endif

! Viscous acceleration
#if defined(ARTIFICIAL_VISCOSITY)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = sph(p)%a_visc(k)*real(ascale,PR)
     vtemp(k) = alldata(nelements)
  end do
  av_rad = dot_product(vtemp(1:NDIM),dr_unit(1:NDIM))
  av_vel = dot_product(vtemp(1:NDIM) - &
       &av_rad*dr_unit(1:NDIM),dv_unit(1:NDIM))
#endif

! Wind acceleration
#if defined(STELLAR_WIND)
  do k=1,VDIM
     nelements = nelements + 1
     alldata(nelements) = sph(p)%a_wind(k)*real(ascale,PR)
     vtemp(k) = alldata(nelements)
  end do
  aw_rad = dot_product(vtemp(1:NDIM),dr_unit(1:NDIM))
#endif
#endif
  
! Smoothing length (always written)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%h*real(rscale,PR)
  
! Density (always written)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%rho*real(rhoscale*rhocgs,PR)
  
! div_v (always written)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%div_v/real(tscale,PR)
  
! timestep (always written)
  call timestep_size(p,dt)
  nelements = nelements + 1
  alldata(nelements) = real(dt,PR)*real(tscale,PR)
  
#if defined(HYDRO)
! Temperature
  nelements = nelements + 1
  alldata(nelements) =sph(p)%temp
  
! Sound speed
  nelements = nelements + 1
  alldata(nelements) = sph(p)%sound*real(vscale,PR)
  
! Pressure
  nelements = nelements + 1
  alldata(nelements) = sph(p)%press*real(Pscale,PR)
  
! Specific internal energy
  nelements = nelements + 1
  alldata(nelements) = sph(p)%u*real(uscale,PR)

! Compressional heating rate
#if defined(ENERGY_EQN)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%dudt*real(uscale/tscale,PR)
#endif

! Entropic function
#if defined(ENTROPIC_FUNCTION)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%Aent
  nelements = nelements + 1
  alldata(nelements) = sph(p)%dAdt
#endif
#endif
  
! Radial distance (always written)
  nelements = nelements + 1
  alldata(nelements) = drmag*real(rscale,PR)
  
! Other radial values
  nelements = nelements + 1
  alldata(nelements) = v_rad
  nelements = nelements + 1
  alldata(nelements) = a_rad
#if defined(DEBUG_FORCES)
#if defined(HYDRO)
  nelements = nelements + 1
  alldata(nelements) = ah_rad
#endif
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  nelements = nelements + 1
  alldata(nelements) = ag_rad
  nelements = nelements + 1
  alldata(nelements) = sph(p)%gpot*real(Escale/mscale,PR)
#endif
#if defined(ARTIFICIAL_VISCOSITY)
  nelements = nelements + 1
  alldata(nelements) = av_rad
  nelements = nelements + 1
  alldata(nelements) = av_vel
#endif
#if defined(STELLAR_WIND)
  nelements = nelements + 1
  alldata(nelements) = aw_rad
#endif
#endif
  
#if defined(GRAD_H_SPH)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%omega
  nelements = nelements + 1
  alldata(nelements) = h_fac*(sph(p)%m/sph(p)%rho)**(INVNDIM) / &
       &sph(p)%h
#if defined(SELF_GRAVITY)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%zo
#endif
#endif
  
! Smoothing length (always written)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%m*real(mscale,PR)
  
! Particle type
  nelements = nelements + 1
  if (pboundary > 0 .and. p <= pboundary) then
     alldata(nelements) = BOUNDARYID
  else if (picm > 0 .and. p <= pboundary + picm) then
     alldata(nelements) = ICMID
  else if (pgas > 0 .and. p <= pboundary + picm + pgas) then
     alldata(nelements) = GASID
  end if
  
! Radiative transfer routines
#if defined(RAD_WS) && defined(DEBUG_RAD) 
  nelements = nelements + 1
  alldata(nelements) = rad_info(1,p)*real(kappascale,PR)
  nelements = nelements + 1
  alldata(nelements) = rad_info(2,p)*real(kappascale,PR)
  nelements = nelements + 1
  alldata(nelements) = rad_info(3,p)
  nelements = nelements + 1
  alldata(nelements) = rad_info(4,p)
  nelements = nelements + 1
  alldata(nelements) = &
       &rad_info(5,p)*real(mscale*mcgs/rscale/rcgs/rscale/rcgs,PR)
  nelements = nelements + 1
  alldata(nelements) = &
       &rad_info(6,p)/sph(p)%m*real(Escale*Ecgs/mscale/mcgs,PR)
  nelements = nelements + 1
  alldata(nelements) = &
       &rad_info(7,p)*real(Escale*Ecgs/tscale/tcgs/mscale/mcgs,PR)
  nelements = nelements + 1
  alldata(nelements) = &
       &rad_info(8,p)*real(Escale*Ecgs/tscale/tcgs/mscale/mcgs,PR)
  nelements = nelements + 1
  alldata(nelements) = &
       &rad_info(9,p)*real(Escale*Ecgs/tscale/tcgs/mscale/mcgs,PR)
  nelements = nelements + 1
  alldata(nelements) = &
       &rad_info(10,p)*real(Escale*Ecgs/tscale/tcgs/mscale/mcgs,PR)
#endif
  
! Flux-limited diffusion info
#if defined(DIFFUSION)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%ueq*real(uscale,PR)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%dt_therm*tscale
  nelements = nelements + 1
  alldata(nelements) = sph(p)%du_dt_diff*real(uscale/tscale,PR)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%k_cond*Escale*Ecgs/rscale/rcgs/tscale
#endif
  
#if defined(DEBUG_SINK_BOUNDARY_PROPERTIES)
  nelements = nelements + 1
  alldata(nelements) = a_rad
  nelements = nelements + 1
  alldata(nelements) = &
       &dot_product(sph(p)%a(1:NDIM) - a_rad*dr_unit(1:NDIM),dr_unit(1:NDIM))
#endif
  
#if defined(VISC_TD)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%talpha
#endif
  
#if defined(VISC_BALSARA)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%balsara
#endif

#if defined(VISC_PATTERN_REC)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%pattrec
#endif
  
#if defined(DIV_A)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%div_a
#endif

#if defined(SIGNAL_VELOCITY)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%vsigmax*real(vscale,PR)
#endif

#if defined(SM2012_SPH)
  nelements = nelements + 1
  alldata(nelements) = sph(p)%q*real(uscale,PR)*real(rhoscale,PR)
#endif

#if defined(DEBUG_HP_WALK_ALL_RAYS)
  nelements = nelements + 1
  alldata(nelements) = whichHPlevel(p)
#endif

! Time
  nelements = nelements + 1
  alldata(nelements) = real(time*tscale,PR)
  

  return
END SUBROUTINE record_particle_data
