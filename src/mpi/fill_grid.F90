! FILL_GRID.F90
! A. McLeod - 23/10/2008
! Fills a grid with approximate particle density information
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine fill_grid(dx, minr, minimal)
    use definitions
    use mpi_communication_module
    use mpi_decomposition_module
    use particle_module
    use type_module, only : phydrostart
    use time_module
    use neighbour_module, only : pp_gather
    use periodic_module
#ifdef SINKS
    use sink_module
#endif
    implicit none

    real(kind=PR), intent(in) :: minr(1:NDIM)  ! Min particle position
    real(kind=PR), intent(in) :: dx(1:NDIM)    ! Grid spacing
    logical, intent(in)       :: minimal       ! Are we using minimal particles?
    
    integer :: gridslot(1:NDIM)        ! Position of the ptcl in the grid
    integer :: p                       ! Loop counter
    real(kind=PR) :: particleweight    ! 'Weight' given to particle
    real(kind=PR) :: cost              ! 'Cost' given to particle
    real(kind=PR) :: acc_count         ! Np. of expected acceleration steps
    real(kind=PR) :: maxadd            ! Maximum per particle
#ifdef SINKS
    integer :: s                       ! Sink counter
#endif
    
    
    !maxadd = max(1._PR,real(totalptot/pp_gather,PR)) ! Integer division
    maxadd = BIG_NUMBER
    MPIpdensitygrid = 0.0_PR
    
    ! Loop over all SPH particles that require work
    ! (i.e. assume boundary particles are no work).
    ! ------------------------------------------------------------------------
    do p=phydrostart,ptot
       if (minimal) then
          gridslot = int((minimal_sph(p)%r - minr) / dx)
       else
          gridslot = int((sph(p)%r - minr) / dx)
       end if
#if defined(LOW_MEM)
       cost = 1.0_PR
#else
       cost = 0.0_PR
#endif
#ifdef DEBUG2
       if (any(gridslot < 0)) then
          write (6,*) "gridslot = ", gridslot
          write (6,*) "minr = ", minr
          write (6,*) "sph(p)%r = ", sph(p)%r
          write (6,*) "dx = ", dx
          if (minimal) then
             write (6,*) "minimal_sph(p)%r - minr = ", &
                         &minimal_sph(p)%r - minr
             write (6,*) "[minimal_sph(p)%r - minr] / dx= ", &
                         &(minimal_sph(p)%r - minr)/dx
          else
             write (6,*) "sph(p)%r - minr = ", sph(p)%r - minr
             write (6,*) "[sph(p)%r - minr] / dx= ", (sph(p)%r - minr)/dx
          end if
       end if
#endif
       ! Estimate number of times particle will have accdo=.TRUE.
       ! in the next loadbalance_nsteps
       if (minimal) then
          particleweight = pp_gather*COST_HYDRO
       else
#if defined(EULER) || defined(LEAPFROG_KDK)
          acc_count = real(loadbalance_nsteps,PR) * &
               & (COST_BASE + real(2**(level_step - sph(p)%nlevel),PR) )
#else
          acc_count = real(loadbalance_nsteps,PR) * &
               & (COST_BASE + real(2**(level_step - sph(p)%nlevel - 1_ILP),PR) )
#endif
#if !defined(LOW_MEM)
#if HYDRO
          cost = cost + sph(p)%hydro_calcs * COST_HYDRO
#endif
#if defined(SELF_GRAVITY)
          cost = cost + sph(p)%gravity_calcs * COST_GRAVITY
#endif
#endif
          particleweight = cost * acc_count
       end if
       
       ! Add 1 + 1/(int. timestep size) for each particle in each grid cell
       ! This is to approximate the effect of multiple particle timesteps
#if NDIM==3
       MPIpdensitygrid(gridslot(1),gridslot(2),gridslot(3)) = &
            & MPIpdensitygrid(gridslot(1),gridslot(2),gridslot(3))&
            & + min(maxadd,particleweight)
#elif NDIM==2
       MPIpdensitygrid(gridslot(1),gridslot(2)) = &
            & MPIpdensitygrid(gridslot(1),gridslot(2))&
            & + min(maxadd,particleweight)
#else
       MPIpdensitygrid(gridslot(1)) = &
            &MPIpdensitygrid(gridslot(1)) &
            & + min(maxadd,particleweight)
#endif
    end do
    ! ------------------------------------------------------------------------
    
    
    ! If using sink particles, estimate number of times sink will have 
    ! accdo = .true. in the next loadbalance_nsteps
    ! ------------------------------------------------------------------------
#ifdef SINKS
    if (minimal) then
       particleweight = COST_SINK
    else
#ifndef EULER
       acc_count = real(loadbalance_nsteps,PR) * &
            & (COST_BASE + real(2**(level_step - nlevel_sinks - 1_ILP),PR) )
#else
       acc_count = real(loadbalance_nsteps,PR) * &
            & (COST_BASE + real(2**(level_step - nlevel_sinks),PR) )
#endif
       cost = COST_SINK
       particleweight = min(maxadd,cost*acc_count)
    end if
    
    do s=1,stot
       if (sink(s)%domain /= rank) cycle
       gridslot = int((sink(s)%r - minr) / dx)
       
#if NDIM==3
       MPIpdensitygrid(gridslot(1),gridslot(2),gridslot(3)) = &
            & MPIpdensitygrid(gridslot(1),gridslot(2),gridslot(3))&
            & + particleweight
#elif NDIM==2
       MPIpdensitygrid(gridslot(1),gridslot(2)) = &
            & MPIpdensitygrid(gridslot(1),gridslot(2))&
            & + particleweight
#else
       MPIpdensitygrid(gridslot(1)) = &
            &MPIpdensitygrid(gridslot(1)) &
            & + particleweight
#endif
    end do
#endif
    ! ------------------------------------------------------------------------
    
    return
end subroutine fill_grid
