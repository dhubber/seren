! RAD_WS_UPDATE.F90
! D. A. Hubber - 14/7/2008
! Update radiation transport quantities.
! (Ref : Stamatellos et al. 2007; Forgan et al. 2009)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE rad_ws_update
  use interface_module
  use particle_module
  use hydro_module
  use Eos_module
  use type_module
  use time_module
  implicit none

  integer :: acctot                   ! No. of particles on acc. step
  integer :: i                        ! Aux. particle counter
  integer(kind=ILP) :: dn             ! Int. timestep since start of timestep
  integer :: p                        ! particle counter
  integer, allocatable :: acclist(:)  ! List of particles on acc. step
  real(kind=PR) :: dt                 ! Timestep
  real(kind=PR) :: dt_new             ! Latest timestep
  real(kind=PR) :: dt_old             ! Previous timestep
  real(kind=PR) :: mu_bar_p           ! Mean gas particle mass for p
#if defined(OPENMP)
  integer :: chunksize                ! Data packet size for dynamic OpenMP
#endif

  debug2("Updating radiation transport quantities [rad_ws_update.F90]")
  debug_timing("RAD_WS")

! For block timesteps, first make a list of all hydro SPH particles on 
! an acceleration step, and then parallelize over that list.
  acctot = 0
  allocate(acclist(1:ptot))
!  do p=pgasstart,pgasend
  do p=1,pgasend
     if (sph(p)%accdo) then
        acctot = acctot + 1
        acclist(acctot) = p
     end if
  end do
#if defined(OPENMP)
  chunksize = int(CHUNKFRAC*real(acctot,PR)) + 1
#endif


! Calculate flux-limited diffusion terms
! ----------------------------------------------------------------------------
#if defined(DIFFUSION)
 if (acctot > 0) then
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
        call find_idens(sph(p)%rho,sph(p)%idens)
        call find_itemp(sph(p)%temp,sph(p)%itemp)
        call conductivity(p)
     end do
     !$OMP END PARALLEL DO
 end if
#endif
! ----------------------------------------------------------------------------


! Calculate itemp, idens and column density to infinity of all particles
! ----------------------------------------------------------------------------
#if defined(DIFFUSION)
  if (acctot > 0) then
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) &
    !$OMP PRIVATE(p,dt,dt_old,dt_new)
     do i=1,acctot
        p = acclist(i)

        dt_old = real(sph(p)%laststep,PR)
        dt_new = real(2**(level_step - sph(p)%nlevel),PR)*real(timestep,PR)
#if defined(EULER) || defined(LEAPFROG_KDK) 
        dt = dt_new
#elif defined(RUNGE_KUTTA2)
        dt = 0.5_PR*dt_new
#elif defined(LEAPFROG_DKD)
        dt = 0.5_PR*(dt_old + dt_new)
#endif
        call diffusion(p,dt)
     end do
     !$OMP END PARALLEL DO
  end if
#endif
! ----------------------------------------------------------------------------


! Perform implicit integration of internal energy
! ----------------------------------------------------------------------------
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) &
!$OMP PRIVATE(dn,dt,mu_bar_p) 
  do p=pgasstart,pgasend

     if (sph(p)%accdo) call find_equilibrium_temp_ws(p)

     ! Calculate time since last force computation.
     dn = n - sph(p)%nlast
     dt = real(timestep,PR)*real(dn,PR)
     
#if defined(IONIZING_UV_RADIATION)
     if (sph(p)%temp < sph(p)%tempmin) then
        sph(p)%temp = sph(p)%tempmin
        call find_idens(sph(p)%rho,sph(p)%idens)
        call find_itemp(sph(p)%temp,sph(p)%itemp)
        sph(p)%u = eosenergy(sph(p)%rho,sph(p)%temp,sph(p)%idens,sph(p)%itemp)
        sph(p)%ueq = sph(p)%u
     end if
#endif
     
     ! Perform implicit integration depending on timestep.
     if (sph(p)%dt_therm <= SMALL_NUMBER) then
        sph(p)%u = sph(p)%u_old
     else if (dt < 40.0_PR*sph(p)%dt_therm) then
        sph(p)%u = sph(p)%u_old*exp(-dt/sph(p)%dt_therm) &
             & + sph(p)%ueq*(1.0_PR - exp(-dt/sph(p)%dt_therm))
     else if (dt >= 40.0_PR*sph(p)%dt_therm) then
        sph(p)%u = sph(p)%ueq
     end if
     
     ! Now update all other thermal properties
     call find_idens(sph(p)%rho,sph(p)%idens)
     call find_temp_from_energy(sph(p)%idens,sph(p)%u,sph(p)%itemp,sph(p)%temp)
     mu_bar_p = eosmu(sph(p)%rho,sph(p)%temp,sph(p)%idens,sph(p)%itemp)
     sph(p)%press = Pconst2*sph(p)%temp*sph(p)%rho/mu_bar_p
     sph(p)%sound = sqrt(sph(p)%press/sph(p)%rho)

  end do
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------------

  deallocate(acclist)

  return
END SUBROUTINE rad_ws_update
