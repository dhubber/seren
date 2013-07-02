! NBODY_TIMESTEPS.F90
! D. A. Hubber - 6/9/2010
! Computes ideal timesteps for all N-body particles, then quantises to lower 
! block timestep level.  Allows the timestep of a particle, and the global 
! minimum, to increase by one level if synchronised correctly.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_timesteps
  use time_module
  use Nbody_module
  use scaling_module
  use sink_module
  implicit none

  integer :: s                           ! Sink particle counter
  integer :: s_update                    ! No. of timesteps updated
  integer(kind=ILP) :: level             ! Timestep level
  integer(kind=ILP) :: level_old         ! Old timestep level
  integer(kind=ILP) :: nlastlevel        ! Timestep level
  real(kind=DP) :: dt_min                ! Minimum step size
  real(kind=DP) :: dt                    ! Step size of particle p
  real(kind=DP), allocatable :: step(:)  ! Step size storage
#if defined(DEBUG_BLOCK_TIMESTEPS)
  integer :: pmin                               ! Particle with dt_min
  integer(kind=ILP), allocatable :: ninlevel(:) ! Number of particles in level
#endif

  debug2("Calculating new/updated timesteps [nbody_timesteps.F90]")
  debug_timing("NBODY_TIMESTEPS")
  s_update = 0
  allocate(step(1:SMAX+nbin))

! Check if we are at a resynchronisation step
! ============================================================================
  if (n == nresync) then

      debug2("Resynchronising block timesteps [nbody_timesteps.F90]")

      ! Wrap n back around to 0 for resync step
      n = 0_ILP
      dt_min = BIG_NUMBER_DP
      
      ! If integrating binary COMs, calculate timesteps
      do s=1,stot
         star(s)%nlast = 0_ILP
#if defined(BINARY_COM_MOTION)
         b = star(s)%binud
         if (b > 0) then
            call nbody_binary_timestep_size(b,dt)
         else
            call nbody_timestep_size(s,dt)
         end if
#else
         call nbody_timestep_size(s,dt)
#endif
         step(s) = dt
         dt_min = min(dt,dt_min)
      end do
      
#if defined(DEBUG_BLOCK_TIMESTEPS)
      write(6,*) "dt_min : ",dt_min*tscale
#endif

      ! Compute min/max timestep values for resync 
      level_max = nlevels - 1_ILP
#if defined(RESTRICTED_TIMESTEP_LEVELS)
      level     = int(INVLOGETWO*log(dt_fixed / dt_min),ILP) + 1_ILP
      dt_min    = dt_fixed / 2.0_DP**(level)
#elif defined(FIXED_TIMESTEP_LEVELS)
      dt_min    = dt_fixed / 2.0_DP**(level_max)
#endif
      dt_max    = dt_min * 2.0_DP**(level_max)

#if defined(DEBUG_BLOCK_TIMESTEPS)
      write(6,*) "dt_min : ",dt_min*tscale,"    dt_max : ",dt_max*tscale
      write(6,*) "level_max : ",level_max
#endif
      
      ! Assign new step levels to all stars
      ! ----------------------------------------------------------------------
      do s=1,stot
         dt    = step(s)
         level = min(int(INVLOGETWO*log(dt_max/dt),ILP) + 1_ILP, level_max)
         level = max(level,0_ILP)
         star(s)%nlevel = level
#if defined(DEBUG_BLOCK_TIMESTEPS)
         if (level > level_max) then
            write(6,*) "Error : level > level_max : ",&
                 &level,level_max,s,dt,dt_min
            stop
         end if
#endif
      end do

      level_step = level_max
      nresync    = 2**(level_step)
      timestep   = dt_max / real(nresync,DP)
      nstepsize  = 1_ILP
      nmaxsteps  = nmaxsteps + 1_ILP


! ============================================================================
! If not resync time, check if any timesteps need to be recomputed
! ============================================================================
  else

     debug2("Recalculating particle stepsizes [nbody_timesteps.F90]")
     level_old = level_max

     ! Loop stars and find which ones need their timesteps to be updated.
     ! -----------------------------------------------------------------------
     do s=1,stot

        ! Skip if timestep not finished
        if (n == star(s)%nlast) then
           s_update = s_update + 1
           nlastlevel = star(s)%nlevel

           ! If integrating binary COMs, calculate timesteps
#if defined(BINARY_COM_MOTION)
           b = star(s)%binud
           if (b > 0) then
              call nbody_binary_timestep_size(b,dt)
           else
              call nbody_timestep_size(s,dt)
           end if
#else
           call nbody_timestep_size(s,dt)
#endif
           level = max(int(INVLOGETWO*log(dt_max/dt),ILP) + 1_ILP, 0_ILP)
           
           ! Allow timestep to go up one level if synchronised correctly.
           ! Else, allow it to fall to any lower timetep.
           if (level < nlastlevel .and. nlastlevel > 1_ILP .and. &
                mod(n,2**(level_step - (nlastlevel - 1_ILP))) == 0_ILP) then
              star(s)%nlevel = nlastlevel - 1_ILP
           else if (level > nlastlevel) then
              star(s)%nlevel = level
           else
              star(s)%nlevel = nlastlevel
           end if
        end if
        
     end do
     ! -----------------------------------------------------------------------


     ! If no stars are updated, do not update timesteps
     ! ----------------------------------------------------------------------
     if (s_update > 0) then
         
        ! Calculate maximum level
        level_max = 0
        do s=1,stot
           level_max = max(level_max,star(s)%nlevel)
        end do
         
        ! If minimum stepsize is synchronized, then increase if possible
        if (level_max <= level_old - 1_ILP .and. level_old > 1_ILP .and. &
             & mod(n,2_ILP**(level_step - (level_old - 1_ILP))) == 0_ILP) then
#if defined(DEBUG_BLOCK_TIMESTEPS)
           write(6,*) "Timestep level removed : ",&
                &level_old,level_old-1,level_max
#endif
           level_max = level_old - 1_ILP
        else if (level_max < level_old) then
           level_max = level_old
        end if
        
#if defined(DEBUG_BLOCK_TIMESTEPS)
        write(6,*) "level_max : ",level_max,"     level_old : ",level_old
#endif
        
        ! Rescale all integer step variables if level has changed
        ! -------------------------------------------------------------------
        if (level_max > level_old) then
           do s=1,stot
              star(s)%nlast = star(s)%nlast*(2**(level_max - level_old))
           end do
#if defined(DEBUG_BLOCK_TIMESTEPS)
           write(6,*) "Rescaled n : ",&
                &level_max,level_old,n,n*(2_ILP**(level_max - level_old))
#endif
           n = n * (2_ILP**(level_max - level_old))
           nresync = nresync*(2_ILP**(level_max - level_old))
        else if (level_max < level_old) then
           do s=1,stot
              star(s)%nlast = star(s)%nlast / (2**(level_old - level_max))
           end do
#if defined(DEBUG_BLOCK_TIMESTEPS)
           write(6,*) "Rescaled n : ",&
                &level_max,level_old,n,n/(2_ILP**(level_old - level_max))
#endif
           n = n / (2_ILP**(level_old - level_max))
           nresync = nresync / (2_ILP**(level_old - level_max))
        end if
        
        level_step = level_max 
        timestep = dt_max / real(nresync,DP)
        
     end if
     ! ----------------------------------------------------------------------
     
  end if
! ===========================================================================


#if defined(DEBUG_BLOCK_TIMESTEPS)
  write(6,*) "n      : ",n,"    nresync : ",nresync
  write(6,*) "nsteps : ",nsteps,"   time : ",time*tscale
  if (n == 0 .or. s_update > 0) then
     allocate(ninlevel(0:level_max))
     ninlevel(0:level_max) = 0
     do s=1,stot
        level = star(s)%nlevel
        ninlevel(level) = ninlevel(level) + 1_ILP
     end do
     do level=0,level_max
        write(6,*) "No. of particles in level ",level," : ",ninlevel(level)
     end do
     if (ninlevel(level_max) == 0 .and. stot == 0) then
        write(6,*) "No particles in highest level : ",ninlevel(level_max)
        do s=1,stot
           if (star(s)%nlevel == level_max) &
                &write(6,*) "But there's one here ... ",s,star(s)%nlevel
        end do
        stop
     end if
     deallocate(ninlevel)
  end if
#endif

  if (allocated(step)) deallocate(step)

  return
END SUBROUTINE nbody_timesteps
