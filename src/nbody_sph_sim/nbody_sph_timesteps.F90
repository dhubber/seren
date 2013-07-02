! SPH_SIM_TIMESTEPS.F90
! D. A. Hubber - 6/9/2010
! Computes ideal timesteps for all particles, then quantises to lower 
! block timestep level.  Allows the timestep of a particle, and the global 
! minimum, to increase by one level if synchronised correctly.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_sph_timesteps
  use time_module
  use particle_module
  use hydro_module
  use scaling_module
  use sink_module
  use Nbody_module
  implicit none

  integer :: p                                  ! Particle counter
  integer :: p_update                           ! No. of timesteps updated
  integer :: s                                  ! Sink particle counter
  integer(kind=ILP) :: level                    ! Timestep level
  integer(kind=ILP) :: level_old                ! Old timestep level
  integer(kind=ILP) :: nlastlevel               ! Timestep level
  real(kind=DP) :: dt_min                       ! Minimum step size
  real(kind=DP) :: dt                           ! Step size of particle p
  real(kind=DP), allocatable :: step(:)         ! Step size storage
#if defined(FIXED_SPH_TIMESTEP_LEVEL)
  real(kind=DP) :: dt_sph                       ! SPH particle timestep
#endif
#if defined(DEBUG_BLOCK_TIMESTEPS)
  integer(kind=ILP), allocatable :: ninlevel(:) ! Number of particles in level
  integer :: pmin                               ! Particle with dt_min
  integer :: smin                               ! ..
  real(kind=DP) :: asqd                         ! Acceleration squared
  real(kind=DP) :: a1sqd                        ! jerk squared
  real(kind=DP) :: a2sqd                        ! 2nd derivative squared
  real(kind=DP) :: a3sqd                        ! 3rd derivative squared
  real(kind=DP) :: dt2                          ! dt*dt
#endif

  debug_timing("NBODY_SPH_TIMESTEPS")
  p_update = 0
#if defined(FORCE_SPLITTING)
  sphminstep = .false.
#endif

! Check if we are at a resynchronisation step
! ============================================================================
  if (n == nresync) then

      debug2("Resynchronising block timesteps [nbody_sph_timesteps.F90]")

      ! Wrap n back around to 0 for resync step
      n = 0_ILP
      allocate(step(1:ptot+stot))
      dt_min = BIG_NUMBER_DP
      
      ! Loop over all hydro and gas particles and find ideal time steps,
      ! and also the minimum of all particle timesteps
      ! ----------------------------------------------------------------------
      !$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(MIN:dt_min) PRIVATE(dt)
      do p=1,ptot
         sph(p)%nlast = 0_ILP
         call timestep_size(p,dt)
         step(p) = dt
         dt_min = min(dt,dt_min)
      end do
      !$OMP END PARALLEL DO
#if defined(FIXED_SPH_TIMESTEP_LEVEL)
      dt_sph = dt_min
#endif

      do s=1,stot
         star(s)%nlast = 0_ILP
         call nbody_sph_timestep_size(s,dt)
         step(s + ptot) = dt
         dt_min = min(dt,dt_min)
      end do
      
#if defined(DEBUG_BLOCK_TIMESTEPS)
      write(6,*) "dt_min : ",dt_min*tscale
#if defined(FIXED_SPH_TIMESTEP_LEVEL)
      write(6,*) "dt_sph : ",dt_sph*tscale
#endif
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
      write(6,*) "level_max : ",&
           &int(INVLOGETWO*log(dt_max / minval(step)),ILP) + 1_ILP
#endif
      
      ! Loop over hydro and gas particles and assign step levels
      ! ----------------------------------------------------------------------
#if defined(FIXED_SPH_TIMESTEP_LEVEL)
      level = min(int(INVLOGETWO*log(dt_max/dt_sph),ILP) + 1_ILP, level_max)
      level = max(level,0_ILP)
      do p=1,ptot
         sph(p)%nlevel = level
#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
         sph(p)%nminneib = level
#endif
      end do
      ! ----------------------------------------------------------------------
#else
      do p=1,ptot
         dt        = step(p)
         level     = min(int(INVLOGETWO*log(dt_max/dt),ILP) + 1_ILP, level_max)
         level     = max(level,0_ILP)
         sph(p)%nlevel = level
#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
         sph(p)%nminneib = level
#endif
#if defined(DEBUG_BLOCK_TIMESTEPS)
         if (level > level_max) then
            write(6,*) "Error : level > level_max : ",&
                 &level,level_max,p,dt,dt_min
            stop
         end if
#endif
      end do
#endif
      ! ----------------------------------------------------------------------


      ! Loop over hydro and gas particles and assign step levels
      ! ----------------------------------------------------------------------
      do s=1,stot
         dt    = step(s + ptot)
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


      ! 2 integer steps required per real time step for Runge-Kutta,
      ! leapfrog or predictor-corrector, else 1 integer timestep for Euler.
#if defined(EULER) || defined(LEAPFROG_KDK)
      level_step = level_max
#else
      level_step = level_max + 1_ILP
#endif
      nresync = 2**(level_step)
      timestep = dt_max / real(nresync,DP)
      nstepsize = 1
      nmaxsteps = nmaxsteps + 1_ILP
      synchronise_all = .false.
#if defined(FORCE_SPLITTING)
      sphminstep = .true.
#endif

      ! Free up memory for temp arrays
      deallocate(step)


! ============================================================================
! If not resync time, check if any timesteps need to be recomputed
! ============================================================================
  else

     debug2("Recalculating particle stepsizes [sph_timesteps.F90]")
     level_old = level_max

     ! Loop over all hydro and gas particles and find which particles 
     ! need their timesteps to be updated.
     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dt,level,nlastlevel) &
     !$OMP REDUCTION(+ : p_update)
     do p=1,ptot

        ! Skip if timestep not finished
        if (n == sph(p)%nlast) then
           p_update = p_update + 1
           nlastlevel = sph(p)%nlevel
           
           ! Find ideal timestep for particle p
           call timestep_size(p,dt)
           level = max(int(INVLOGETWO*log(dt_max/dt),ILP) + 1_ILP, 0_ILP)
           
           ! Allow timestep to go up one level if synchronised correctly.
           ! Else, allow it to fall to any lower timetep.
           if (level < nlastlevel .and. nlastlevel > 1_ILP .and. &
                mod(n,2**(level_step - (nlastlevel - 1_ILP))) == 0_ILP) then
              sph(p)%nlevel = nlastlevel - 1_ILP
           else if (level > nlastlevel) then
              sph(p)%nlevel = level
           else
              sph(p)%nlevel = nlastlevel
           end if
        end if
        
     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------

#if defined(FORCE_SPLITTING)
     if (p_update > 0) sphminstep = .true.
#endif

     ! Force all SPH particles onto the same level if selected
#if defined(FIXED_SPH_TIMESTEP_LEVEL)
     if (p_update > 0) then
        dt_sph = BIG_NUMBER_DP
        !$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(MIN:dt_sph) PRIVATE(dt)
        do p=1,ptot
           call timestep_size(p,dt)
           dt_sph = min(dt,dt_sph)
        end do
        !$OMP END PARALLEL DO
        level = min(int(INVLOGETWO*log(dt_max/dt_sph),ILP) + 1_ILP, level_max)
        level = max(level,0_ILP)
        do p=1,ptot
           sph(p)%nlevel = level
#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
           sph(p)%nminneib = level
#endif
        end do
     end if
#endif


     ! Loop over all hydro and gas particles and find which particles 
     ! need their timesteps to be updated.
     ! -----------------------------------------------------------------------
     do s=1,stot
        
        ! Skip if timestep not finished
        if (n == star(s)%nlast) then
           p_update = p_update + 1
           nlastlevel = star(s)%nlevel
           
           ! Find ideal timestep for particle p
           call nbody_sph_timestep_size(s,dt)
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


     ! If no particles or sinks are updated, do not update timesteps
     ! -----------------------------------------------------------------------
     if (p_update > 0) then
        
        ! Calculate maximum level
        level_max = 0        
        do p=1,ptot
           level_max = max(level_max,sph(p)%nlevel)
        end do
        do s=1,stot
           level_max = max(level_max,star(s)%nlevel)
        end do
        
        ! If minimum stepsize is synchronized, then increase if possible
        if (level_max <= level_old - 1_ILP .and. level_old > 1_ILP .and. &
             & mod(n,2**(level_step - (level_old - 1_ILP))) == 0_ILP) then
#if defined(DEBUG_BLOCK_TIMESTEPS)
           write(6,*) "Timestep level removed : ",&
                &level_old,level_old-1,level_max
#endif
           level_max = level_old - 1_ILP
        else if (level_max < level_old) then
           level_max = level_old
        end if

        ! If neighbour's timesteps are more than TIMESTEP_DIFF_MAX times 
        ! times smaller, immediately reduce the timestep if possible. 
        ! Also, skip update if particle is neighbouring a sink.
        ! -------------------------------------------------------------------
#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
        do p=1,ptot
           sph(p)%nminneib = min(sph(p)%nminneib,level_max)
           if (sph(p)%nminneib - sph(p)%nlevel > TIMESTEP_LEVEL_DIFF_MAX .and. &
                & sph(p)%nlast == n .and. &
                & mod(n,2**(level_step - sph(p)%nminneib + &
                & TIMESTEP_LEVEL_DIFF_MAX)) == 0_ILP) then
              sph(p)%nlevel = min(sph(p)%nminneib - TIMESTEP_LEVEL_DIFF_MAX,level_max)
              sph(p)%nminneib = sph(p)%nlevel
           else if (sph(p)%nlast == n) then
              sph(p)%nminneib = sph(p)%nlevel
           end if
        end do
#endif
        
        ! Rescale all integer step variables if level has changed
        ! -------------------------------------------------------------------
        if (level_max > level_old) then
           sph(1:ptot)%nlast = sph(1:ptot)%nlast*(2**(level_max - level_old))
           star(1:stot)%nlast = star(1:stot)%nlast*(2**(level_max - level_old))
           n = n*(2_ILP**(level_max - level_old))
           nresync = nresync*(2_ILP**(level_max - level_old))
#if defined(SINKS)
           nlast_sinks = nlast_sinks*(2**(level_max - level_old))
#endif
        else if (level_max < level_old) then
           sph(1:ptot)%nlast = sph(1:ptot)%nlast / (2**(level_old - level_max))
           star(1:stot)%nlast = star(1:stot)%nlast &
                & / (2_ILP**(level_old - level_max))
           n = n / (2_ILP**(level_old - level_max))
           nresync = nresync / (2_ILP**(level_old - level_max))
#if defined(SINKS)
           nlast_sinks = nlast_sinks / (2_ILP**(level_old - level_max))
#endif
        end if

        ! 2 integer steps required per real time step for Runge-Kutta,
        ! leapfrog or predictor-corrector, else 1 integer timestep for Euler.
#if defined(EULER) || defined(LEAPFROG_KDK)
        level_step = level_max
#else
        level_step = level_max + 1_ILP
#endif
        timestep = dt_max / real(nresync,DP)

      end if
      ! ----------------------------------------------------------------------

  end if
! ============================================================================


#if defined(DEBUG_BLOCK_TIMESTEPS)
  write(6,*) "n      : ",n,"    nresync : ",nresync
  write(6,*) "nsteps : ",nsteps,"   time : ",time*tscale
  write(6,*) "dt_max : ",dt_max*tscale,"   timestep : ",timestep*tscale
  if (n == 0 .or. p_update > 0) then
     dt_min = BIG_NUMBER_DP
     smin = 1
     do p=1,ptot
        call timestep_size(p,dt)
        dt_min = min(dt,dt_min)
     end do
     write(6,*) "SPH dt_min  : ",dt_min*tscale
     dt_min = BIG_NUMBER_DP
     do s=1,stot
        call nbody_timestep_size(s,dt)
!        write(6,*) "STUFF : ",s,smin,dt,dt_min
        if (dt < dt_min) smin = s
        dt_min = min(dt,dt_min)
     end do
     write(6,*) "Star dt_min : ",dt_min*tscale
!     write(6,*) "STAR HERE : ",star(smin)
     asqd  = real(dot_product(star(smin)%a,star(smin)%a),DP)
     a1sqd = real(dot_product(star(smin)%adot,star(smin)%adot),DP)
     a2sqd = real(dot_product(star(smin)%a2dot,star(smin)%a2dot),DP)
     a3sqd = real(dot_product(star(smin)%a3dot,star(smin)%a3dot),DP)
     dt2   = (sqrt(asqd*a2sqd) + a1sqd) / (sqrt(a1sqd*a3sqd) + a2sqd)
     write(6,*) "Timestep for star ",smin,dt2,nbody_timemult*sqrt(dt2)*tscale
     write(6,*) "asqd : ",asqd,a1sqd,a2sqd,a3sqd
     allocate(ninlevel(0:level_max))
     ninlevel(0:level_max) = 0
     do p=1,ptot
        level = sph(p)%nlevel
        ninlevel(level) = ninlevel(level) + 1_ILP
     end do
     do level=0,level_max
        write(6,*) "No. of SPH particles in level  ",level," : ",ninlevel(level)
     end do
     ninlevel(0:level_max) = 0
     do s=1,stot
        level = star(s)%nlevel
        ninlevel(level) = ninlevel(level) + 1_ILP
     end do
     do level=0,level_max
        write(6,*) "No. of star particles in level ",level," : ",ninlevel(level)
     end do
     deallocate(ninlevel)
  end if
!read(5,*) level
#endif


! Advance time and integer times
  n      = n + 1
  nsteps = nsteps + 1
  time   = time + timestep

! Check that something hasn't gone wrong with times
  if (time < 0.0_DP .or. timestep < 0.0_DP) then
     write(6,*) "Something wrong with timesteps : ",time*tscale,timestep*tscale
     stop
  end if

!  if (timestep*tscale < 1.0D-08) stop 'dt going wrong!!'

  return
END SUBROUTINE nbody_sph_timesteps
