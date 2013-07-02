! ADVANCE_SINK_LEAPFROG_KDK.F90
! D. A. Hubber - 19/3/2007
! Advance positions and velocities of all sink particles using 2nd order 
! Leapfrog kick-drift-kick integration scheme.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE advance_sink_leapfrog_kdk
  use particle_module
  use time_module
  use sink_module
#if defined(USE_MPI)
  use mpi_communication_module
#endif
  implicit none

  integer(kind=ILP) :: dn         ! Integer timestep since start of timestep
  integer(kind=ILP) :: nfull      ! Full integer timestep
  integer(kind=ILP) :: nhalf      ! Half integer timestep
  integer :: s                    ! Sink particle counter
  real(kind=PR) :: dt             ! Physical time since beginning of timestep
  real(kind=PR) :: dt_half        ! Time since halfstep
  real(kind=PR) :: rs(1:NDIM)     ! Local copy of position 
  real(kind=PR) :: vs(1:NDIM)     ! Local copy of velocity

  debug2("Advancing sink particles [advance_sink_LFKDK.F90]")


! Calculate integer and physical time intervals since beginning of step
  accdo_sinks = .false.
  dn          = n - nlast_sinks
  nfull       = 2**(level_step - nlevel_sinks)
  nhalf       = nfull / 2
  dt          = real(timestep,PR)*real(dn,PR)
  dt_half     = real(timestep,PR)*real(dn - nhalf,PR)


! Loop over all sink particles 
! ----------------------------------------------------------------------------
  do s=1,stot
#if defined(USE_MPI)
     if (sink(s)%domain /= rank) cycle
#endif
     if (.not. sink(s)%static) then
        
        ! Advance sink particle s
        ! --------------------------------------------------------------------
        rs(1:NDIM) = sink(s)%rold(1:NDIM) + &
             &sink(s)%vold(1:NDIM)*dt + 0.5_PR*sink(s)%a(1:NDIM)*dt*dt
        vs(1:NDIM) = sink(s)%vold(1:NDIM) + sink(s)%a(1:NDIM)*dt
        
#if defined(BOUNDARY_CONDITIONS)
        call check_boundary_conditions(rs(1:NDIM),vs(1:VDIM))
#endif
        
        ! If end of timestep, record as 'old' values
        if (dn == nfull) then
           !sink(s)%aold(1:NDIM) = sink(s)%a(1:NDIM)
           sink(s)%rold(1:NDIM) = rs(1:NDIM)
           sink(s)%vold(1:NDIM) = vs(1:NDIM)
        end if
        
        ! Record positions and velocities in arrays for all cases 
        sink(s)%r(1:NDIM) = rs(1:NDIM)
        sink(s)%v(1:NDIM) = vs(1:NDIM)
     end if
  end do
! ----------------------------------------------------------------------------
  
! Update quantities for next timestep/force calculation
  if (dn == nfull) then
     accdo_sinks    = .true.
     nlast_sinks    = n
     laststep_sinks = timestep*real(nfull,DP)
  end if

  return
END SUBROUTINE advance_sink_leapfrog_kdk
