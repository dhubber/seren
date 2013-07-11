! CHECK_NEIGHBOUR_TIMESTEPS.F90
! D. A. Hubber - 2/6/2009
! Check through the neighbour lists of all active particles and calculate 
! the minimum neighbour timestep of all their neighbours.  Also flag any
! particles which are 'neighbouring' sink particles so they can be set to
! the minimum timestep.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE check_neighbour_timesteps
  use interface_module, only : distance2,get_neib_on_fly
  use particle_module, only : ptot,sph
  use neighbour_module
  use seren_sim_module
#if defined(SINKS)
  use sink_module
  use Nbody_module
#endif
  implicit none

  integer :: i                             ! Aux. particle counter
  integer :: ilist                         ! ..
  integer :: p                             ! Particle id
  integer :: pp                            ! Neighbour id
  integer :: pp_numb                       ! No. of neighbours
  integer, allocatable :: pp_templist(:)   ! List of neighbours
  integer(kind=ILP) :: levelp              ! Timestep level of particle p
  real(kind=PR) :: dr(1:NDIM)              ! Relative displacement vector
  real(kind=PR) :: drsqd                   ! Distance squared
  real(kind=PR) :: hrangesqd               ! Kernel extent squared
  real(kind=PR) :: rp(1:NDIM)              ! Position of particle p
#if defined(SINKS)
  integer :: s                             ! Sink id
#endif

  debug2("Calculate min. neighbour timesteps [check_neighbour_timesteps.F90]")
  debug_timing("CHECK_NEIBSTEPS")

  allocate(pp_templist(1:ptot))

! Loop over all active particles
! ----------------------------------------------------------------------------
  do p=1,ptot
     if (.not. sph(p)%accdo) cycle
     levelp = sph(p)%nlevel
     rp(1:NDIM) = sph(p)%r(1:NDIM)
     hrangesqd = KERNRANGESQD*sph(p)%h*sph(p)%h

#if defined(NEIGHBOUR_LISTS)
     pp_numb = pptot(p)
     if (pp_numb <= pp_limit) then
        pp_templist(1:pp_numb) = pplist(1:pp_numb,p)
     else
        if (allocated(pp_templist)) deallocate(pp_templist)
        call get_neib_on_fly(p,pp_numb,ptot,pp_templist,&
             &sph(p)%r(1:NDIM),sph(p)%h)
     end if
#else
     if (allocated(pp_templist)) deallocate(pp_templist)
     call get_neib_on_fly(p,pp_numb,ptot,pp_templist,&
          &sph(p)%r(1:NDIM),sph(p)%h)
#endif

     ! Loop over all neighbours and check if timestep of particle p is 
     ! the maximum of its neighbours timestep levels.
     do i=1,pp_numb
        pp = pp_templist(i)
        if (p == pp) cycle
        call distance2(rp(1:NDIM),p,dr(1:NDIM),drsqd)
        if (drsqd >= hrangesqd .and. &
             & drsqd >= KERNRANGESQD*sph(pp)%h*sph(pp)%h) cycle
        if (levelp > sph(pp)%nminneib) sph(pp)%nminneib = levelp
     end do

  end do
! ----------------------------------------------------------------------------


! Loop over all active star particles
! ----------------------------------------------------------------------------
#if defined(NBODY_SPH_SIMULATION)
  if (nbody_sph_sim) then
     do s=1,stot
        if (.not. star(s)%accdo) cycle
        levelp = star(s)%nlevel
        
        ! Copy neighbour lists from array or walk the tree
        call get_neib_on_fly(-s,pp_numb,ptot,pp_templist,&
             &star(s)%r(1:NDIM),0.0_PR)
        
        ! Loop over all neighbours and check if star timestep level is 
        ! the maximum of its neighbours timestep levels.
        do i=1,pp_numb
           pp = pp_templist(i)
           if (levelp > sph(pp)%nminneib .and. sph(pp)%nminneib /= -1_ILP) &
                & sph(pp)%nminneib = levelp
        end do

     end do
  end if
#endif
! ----------------------------------------------------------------------------


! Also check sinks if required.  Mark nminneib as -1 to signify that
! minimum timestep (i.e. same as sinks) must be imposed.
! ----------------------------------------------------------------------------
#if defined(SPH_SIMULATION) && defined(SINKS)
  if (sph_sim) then
     do p=1,ptot
        do s=1,stot
           call distance2(sink(s)%r(1:NDIM),p,dr(1:NDIM),drsqd)
           if (drsqd <= REXTENTSQD*sink(s)%radius*sink(s)%radius) &
                & sph(p)%nminneib = nlevel_sinks
                !& sph(p)%nminneib = min(sph(p)%nminneib,nlevel_sinks)
        end do
     end do
  end if
#endif
! ----------------------------------------------------------------------------

  if (allocated(pp_templist)) deallocate(pp_templist)

  return
END SUBROUTINE check_neighbour_timesteps
