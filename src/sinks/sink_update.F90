! SINK_UPDATE.F90
! D. A. Hubber - 22/04/2010
! Control subroutine for searching for new sinks and accreting SPH particles 
! to existing sinks.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sink_update
  use time_module, only : n,nsteps,nsearchnext,nsearchstep,time
  use sink_module, only : nlast_sinks,stot
  implicit none
  integer            :: psink       ! Candidate sink particle number

  debug2("Searching for new sinks and accreting to old sinks [sink_update.F90]")

! Create sink particles once density reaches required threshold.
! Only check at end of sink (i.e. minimum) timestep.
#if defined(SELF_GRAVITY)
  if (nsteps == nsearchnext) then
     ! Search step, there may be a new sink
     call sink_search(psink)
#if defined(USE_MPI)
     call sink_transfer_particles(psink)
     call create_sink(psink)
#else
     if (psink > 0) call create_sink(psink)
#endif
  else
#if defined(USE_MPI)
     call sink_transfer_particles(-2)
#endif
  end if
  if (nsteps == nsearchnext) then
     nsearchnext = nsearchnext + nsearchstep
  end if
#endif

! Remove particles that are now bound to sink particles.
! Ensure accretion only occurs at the end of a full sink timestep.
  if (n == nlast_sinks .and. stot > 0) then
#if !defined(SELF_GRAVITY) && defined(USE_MPI)
     call sink_transfer_particles(-1)
#endif
     call accrete_particles
  end if

  return
END SUBROUTINE sink_update
