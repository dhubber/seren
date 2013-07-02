! ADVANCE_SINK_EULER.F90
! D. A. Hubber - 19/3/2007
! Advance positions and velocities of all sink particles using 1st 
! order Euler integration scheme. 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE advance_sink_euler
  use particle_module
  use time_module
  use sink_module
#if defined(USE_MPI)
  use mpi_communication_module
#endif
  implicit none

  integer :: dn                    ! Integer timestep since start of timestep
  integer :: nfull                 ! Full integer timestep
  integer :: s                     ! Sink particle counter
  real(kind=PR) :: dt              ! Physical time since beginning of timestep
  real(kind=PR) :: rs(1:NDIM)      ! Local copy of position 
  real(kind=PR) :: vs(1:NDIM)      ! Local copy of velocity

  debug2("Advancing sink particles [advance_sink_euler.F90]")

! Calculate integer and physical time intervals since beginning of step
  dn    = n - nlast_sinks
!  nfull = nstep_sinks
  nfull = 2**(level_step -
  debug2("Advancing sink particles [advance_sink_euler.F90]")

! Calculate integer and physical time intervals since beginning of step
  dn    = n - nlast_sinks
!  nfull = nstep_sinks
  nfull = 2**(level_step - nlevel_sinks)
  dt    = real(timestep,PR)*real(dn,PR)
 nlevel_sinks)
  dt    = real(timestep,PR)*real(dn,PR)

! Loop over all sink particles 
! ----------------------------------------------------------------------------
  do s=1,stot
#if defined(USE_MPI)
     if (sink(s)%domain /= rank) cycle
#endif

     ! Advance sink particle s
     ! --------------------------------------------------------------------

     ! Full Euler integration step
     rs(1:NDIM) = sink(s)%rold(1:NDIM) + sink(s)%vold(1:NDIM)*dt
     vs(1:NDIM) = sink(s)%vold(1:NDIM) + sink(s)%a(1:NDIM)*dt
        
#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(rs(1:NDIM),vs(1:VDIM))
#endif

     ! If end of timestep, record as 'old' values
     if (dn == nfull) then
        sink(s)%rold(1:NDIM) = rs(1:NDIM)
        sink(s)%vold(1:NDIM) = vs(1:NDIM)
     end if

     ! Record positions and velocities in arrays for all cases 
     sink(s)%r(1:NDIM) = rs(1:NDIM)
     sink(s)%v(1:NDIM) = vs(1:NDIM)

  end do
! ----------------------------------------------------------------------------

! Update quantities for next timestep/force calculation
  if (dn == nfull) then
     accdo_sinks    = .true.
     nlast_sinks    = n
     laststep_sinks = timestep*real(nfull,DP)
  end if

  return
END SUBROUTINE advance_sink_euler
