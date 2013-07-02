! SEARCH_GHOST_PARTICLES.F90
! D. A. Hubber - 09/05/2011
! Search through all particles in the periodic box to find which ones must be 
! replicated as ghost particles.  Once flagged, the new ghost particles are 
! copied from the originals with basic properties (e.g. mass, velocity, etc..).
! ============================================================================

#include "macros.h"
#if defined(USE_MPI)
#define GHOSTRANGE 1.0_PR
#else
#define GHOSTRANGE 1.4_PR
#endif

! ============================================================================
SUBROUTINE search_ghost_particles
  use interface_module, only : active_particle_list
  use particle_module
  use periodic_module
  use type_module
  implicit none

  integer :: p                        ! Particle counter
#if !defined(USE_MPI)
  integer :: acctot                   ! Total no. of active particles
  integer, allocatable :: acclist(:)  ! List of particles on acc. step
#endif
  
  debug2("Search for new ghosts near boundaries [search_ghost_particles.F90]")

#if !defined(USE_MPI)
  ! Prepare list of particles which are on an acceleration step
  allocate(acclist(1:ptot))
  call active_particle_list(acctot,acclist,sphmask)
  deallocate(acclist)
  
  ! If there are no active particles, return without creating any ghost ptcls.
  if (acctot == 0) then
     return
  end if
#endif

  ! Repeat ghost search in case we need to reallocate memory
  ! ==========================================================================
  do
     pghost = 0
     pghostmax = pmax - ptot

     ! First create ghost in x-direction.
     ! -----------------------------------------------------------------------
#if defined(PERIODIC_X)
     do p=1,ptot+pghost
        if (sph(p)%r(1) < ghost_search_min(1) + GHOSTRANGE*KERNRANGE*sph(p)%h) &
             &call create_ghost_particle(p,1,&
             &sph(p)%r(1) + periodic_size(1),sph(p)%v(1),sph(p)%m)
        if (sph(p)%r(1) > ghost_search_max(1) - GHOSTRANGE*KERNRANGE*sph(p)%h) &
             &call create_ghost_particle(p,1,&
             &sph(p)%r(1) - periodic_size(1),sph(p)%v(1),sph(p)%m)
     end do
#endif
     
     ! Now create ghosts in y-direction (including any ghosts from x-ghosts).
     ! -----------------------------------------------------------------------
#if NDIM > 1 && defined(PERIODIC_Y)
     if (pghost <= pghostmax) then
        do p=1,ptot+pghost
            if (sph(p)%r(2) < ghost_search_min(2) + GHOSTRANGE*KERNRANGE*sph(p)%h) &
                &call create_ghost_particle(p,2,&
                &sph(p)%r(2) + periodic_size(2),sph(p)%v(2),sph(p)%m)
            if (sph(p)%r(2) > ghost_search_max(2) - GHOSTRANGE*KERNRANGE*sph(p)%h) &
                &call create_ghost_particle(p,2,&
                &sph(p)%r(2) - periodic_size(2),sph(p)%v(2),sph(p)%m)
        end do
     end if
#endif
     
     ! Finally create z-ghosts (including any from x-ghosts and y-ghosts).
     ! -----------------------------------------------------------------------
#if NDIM==3 && defined(PERIODIC_Z)
     if (pghost <= pghostmax) then
        do p=1,ptot+pghost
            if (sph(p)%r(3) < ghost_search_min(3) + GHOSTRANGE*KERNRANGE*sph(p)%h) &
                &call create_ghost_particle(p,3,&
                &sph(p)%r(3) + periodic_size(3),sph(p)%v(3),sph(p)%m)
            if (sph(p)%r(3) > ghost_search_max(3) - GHOSTRANGE*KERNRANGE*sph(p)%h) &
                &call create_ghost_particle(p,3,&
                &sph(p)%r(3) - periodic_size(3),sph(p)%v(3),sph(p)%m)
        end do
     end if
#endif

     ! If there's already enough memory, continue with the next part of code.
     if (pghost <= pghostmax) exit

     ! If there's not enough memory, reallocate the main arrays
     write(6,*) "Not enough memory for ghosts : ",pghost,pghostmax
     pghostmax = ptot + int(1.5_DP*real(min(pghost,pghostmax),DP))
     pghost = 0 ! So that expand does not attempt to copy ghosts
     call expand(ptot + pghostmax)
#if defined(BH_TREE)
     call expand_ghost_tree(pghostmax)
#endif
     !stop

  end do
  ! ==========================================================================

  ! Update counter of purely periodic (i.e. non-MPI) ghosts
  pperiodic = pghost

#if defined(DEBUG_GHOST_PARTICLES)
  write(6,*) "pghost : ",pghost
#endif

  return
END SUBROUTINE search_ghost_particles
