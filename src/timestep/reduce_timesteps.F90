! REDUCE_TIMESTEPS.F90
! D. A. Hubber - 18/6/2011
! Reduces the timestep of particles if 
! i)  the 'synchronise_all' flag is set in order to synchronise and re-set all 
!     particle timesteps in the event (e.g. feedback or supernovae explosion) 
!     which can significantly reduce the timesteps of some fraction of the 
!     SPH particles in the simulation.
! ii) Reduce the timesteps of any particles which have neighbours with 
!     relatively low timesteps.  Only reduces timesteps of particles when the 
!     new timestep is correctly synchronized, and the particle has not just 
!     started a new timestep (which is handled in sph_timestep.F90 and other 
!     such routines).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE reduce_timesteps
  use interface_module, only : reduce_particle_timestep
  use particle_module, only : ptot,sph
  use seren_sim_module
  use neighbour_module
  use sink_module
  use Nbody_module
  use time_module
  implicit none

  integer :: p                   ! Particle counter
#if defined(SINKS)
  integer :: s                   ! Sink counter
#endif

  debug2("Reduce timesteps of all flagged particles [reduce_timesteps.F90]")

  ! Synchronise ALL particles (SPH and stars) if 'syncronise_all' flag 
  ! has been activated.
  ! ==========================================================================
  if (synchronise_all) then

     ! Reduce timestep of all SPH particles (if applicable)
     ! -----------------------------------------------------------------------
#if defined(SPH_SIMULATION) || defined(NBODY_SPH_SIMULATION)
     if (sph_sim .or. nbody_sph_sim) then
        !$OMP PARALLEL DO
        do p=1,ptot
           call reduce_particle_timestep(p,nlevels-1)
        end do
        !$OMP END PARALLEL DO
     end if
#endif     

     ! If using sinks in SPH simulation, reduce sink timesteps here.
     ! -----------------------------------------------------------------------
#if defined(SPH_SIMULATION) && defined(SINKS)
     if (sph_sim) then
        if (n /= nlast_sinks) &
             &laststep_sinks = timestep*real(n - nlast_sinks,DP)
        do s=1,stot
           sink(s)%rold(1:NDIM) = sink(s)%r(1:NDIM)
           sink(s)%vold(1:NDIM) = sink(s)%v(1:NDIM)
#if defined(LEAPFROG_KDK)
           if (n /= nlast_sinks) sink(s)%r(1:NDIM) = sink(s)%rold(1:NDIM) &
                & + sink(s)%vold(1:NDIM)*real(laststep_sinks,PR) + &
                & 0.5_PR*sink(s)%aold(1:NDIM)*real(laststep_sinks,PR)**2
#endif
        end do
        nlast_sinks = n
     end if
#endif

     ! If using sinks in SPH simulation, reduce sink timesteps here.
     ! -----------------------------------------------------------------------
#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
     if (nbody_sph_sim .or. nbody_sim) then
        do s=1,stot
           star(s)%nlast = n
        end do
     end if
#endif

     ! Calculate new timesteps for all particles in next timestep level update
     nresync = n


  ! If activated, reduce particle timesteps to within TIMESTEP_LEVEL_DIFF_MAX
  ! levels of its neighbours minimum timestep.
  ! ==========================================================================
#if defined(CHECK_NEIGHBOUR_TIMESTEPS) && defined(IMMEDIATE_TIMESTEP_REDUCTION)
  else

     ! Loop over all particles and reduce the timestep of any SPH particles  
     ! if the timesteps are synchronised with the new timestep level.
     ! -----------------------------------------------------------------------
#if defined(SPH_SIMULATION) || defined(NBODY_SPH_SIMULATION)
     if (sph_sim .or. nbody_sph_sim) then
        !$OMP PARALLEL DO
        do p=1,ptot
           if (sph(p)%nminneib > level_max) then
              write(6,*) "Problem with nminneib > level_max : ",&
                   &sph(p)%nminneib,level_max
              stop
           end if
           if (sph(p)%nminneib - sph(p)%nlevel > TIMESTEP_LEVEL_DIFF_MAX &
                & .and. mod(n,2**(level_step - sph(p)%nminneib + &
                & TIMESTEP_LEVEL_DIFF_MAX)) == 0_ILP) then
              call reduce_particle_timestep(p,&
                   &sph(p)%nminneib - TIMESTEP_LEVEL_DIFF_MAX)
           end if
        end do
        !$OMP END PARALLEL DO
     end if
#endif
     ! -----------------------------------------------------------------------
#endif

  end if
  ! ==========================================================================

  synchronise_all = .false.

  return
END SUBROUTINE reduce_timesteps
