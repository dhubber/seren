! SEARCH_EXTERNAL_GHOST_PARTICLES.F90
! D. A. Hubber - 09/05/2011
! Search through all particles in the periodic box to find which ones must be 
! replicated as ghost particles on other tasks. Only store the location and
! original IDs to save memory.
! ============================================================================

#include "macros.h"
#if defined(USE_MPI)
#define GHOSTRANGE 1.0_PR
#else
#define GHOSTRANGE 1.4_PR
#endif
#define GKR (GHOSTRANGE*KERNRANGE)

! ============================================================================
SUBROUTINE search_external_ghost_particles
  use interface_module, only : active_particle_list
  use particle_module
  use periodic_module
  use type_module
  use mpi_communication_module
  implicit none

  integer :: p                        ! Particle counter
  integer :: d                        ! Domain counter
  integer :: porig                    ! Original particle
  integer :: pghost_external          ! External ghosts for this task
  integer :: pghost_external_prev     ! External ghosts (previous dimensions)
  integer :: pghost_external_max      ! Array size of temp_ghosts

  type(ghost_sph_reference_type), allocatable :: ghosts(:), tmp_ghosts(:)
  
  debug2("Search for new ghosts near boundaries &
         &[search_external_ghost_particles.F90]")

  allocate(ghosts(1:ptot)) ! This really should be enough most of the time
  pghost_external_max = ptot
  
  if (allocated(ghost_sph_reference)) deallocate(ghost_sph_reference)
  allocate(ghost_sph_reference(0:lastrank))

  ! For each domain, find ghost particles
  ! ==========================================================================
  do d=0, lastrank
     if (d==rank) cycle
     
#if defined(PERIODIC_X)
     ghost_search_min(1) = bbmax(1,d) - periodic_size(1)
     ghost_search_max(1) = bbmin(1,d) + periodic_size(1)
#endif
#if defined(PERIODIC_Y)
     ghost_search_min(2) = bbmax(2,d) - periodic_size(2)
     ghost_search_max(2) = bbmin(2,d) + periodic_size(2)
#endif
#if defined(PERIODIC_Z)
     ghost_search_min(3) = bbmax(3,d) - periodic_size(3)
     ghost_search_max(3) = bbmin(3,d) + periodic_size(3)
#endif

     pghost_external = 0
     pghost_external_prev = 0

     ! First create ghost in x-direction.
     ! -----------------------------------------------------------------------
#if defined(PERIODIC_X)
     do p=1,ptot
        if (sph(p)%r(1) < ghost_search_min(1) + GKR*sph(p)%h) then
           if (pghost_external == pghost_external_max) call expand_temp_ghosts
           pghost_external = pghost_external + 1
           call create_external_ghost_particle(p,sph(p)%r,&
              & ghosts(pghost_external))
           ghosts(pghost_external)%r(1) = &
              & ghosts(pghost_external)%r(1) + periodic_size(1)
        end if
        if (sph(p)%r(1) > ghost_search_max(1) - GKR*sph(p)%h) then
           if (pghost_external == pghost_external_max) call expand_temp_ghosts
           pghost_external = pghost_external + 1
           call create_external_ghost_particle(p,sph(p)%r,&
              & ghosts(pghost_external))
           ghosts(pghost_external)%r(1) = &
              & ghosts(pghost_external)%r(1) - periodic_size(1)
        end if
     end do
     pghost_external_prev = pghost_external
#endif
     
     ! Now create ghosts in y-direction (including any ghosts from x-ghosts).
     ! -----------------------------------------------------------------------
#if NDIM > 1 && defined(PERIODIC_Y)
     do p=1,ptot
        if (sph(p)%r(2) < ghost_search_min(2) + GKR*sph(p)%h) then
           if (pghost_external == pghost_external_max) call expand_temp_ghosts
           pghost_external = pghost_external + 1
           call create_external_ghost_particle(p,sph(p)%r,&
              & ghosts(pghost_external))
           ghosts(pghost_external)%r(2) = &
              & ghosts(pghost_external)%r(2) + periodic_size(2)
        end if
        if (sph(p)%r(2) > ghost_search_max(2) - GKR*sph(p)%h) then
           if (pghost_external == pghost_external_max) call expand_temp_ghosts
           pghost_external = pghost_external + 1
           call create_external_ghost_particle(p,sph(p)%r,&
              & ghosts(pghost_external))
           ghosts(pghost_external)%r(2) = &
              & ghosts(pghost_external)%r(2) - periodic_size(2)
        end if
     end do
     do p=1,pghost_external_prev
        porig = ghosts(p)%porig
        if (ghosts(p)%r(2) < ghost_search_min(2) + GKR*sph(porig)%h) then
           if (pghost_external == pghost_external_max) call expand_temp_ghosts
           pghost_external = pghost_external + 1
           call create_external_ghost_particle(porig,ghosts(p)%r,&
              & ghosts(pghost_external))
           ghosts(pghost_external)%r(2) = &
              & ghosts(pghost_external)%r(2) + periodic_size(2)
        end if
        if (ghosts(p)%r(2) > ghost_search_max(2) - GKR*sph(porig)%h) then
           if (pghost_external == pghost_external_max) call expand_temp_ghosts
           pghost_external = pghost_external + 1
           call create_external_ghost_particle(porig,ghosts(p)%r,&
              & ghosts(pghost_external))
           ghosts(pghost_external)%r(2) = &
              & ghosts(pghost_external)%r(2) - periodic_size(2)
        end if
     end do
     pghost_external_prev = pghost_external
#endif
     
     ! Finally create z-ghosts (including any from x-ghosts and y-ghosts).
     ! -----------------------------------------------------------------------
#if NDIM==3 && defined(PERIODIC_Z)
     do p=1,ptot
        if (sph(p)%r(3) < ghost_search_min(3) + GKR*sph(p)%h) then
           if (pghost_external == pghost_external_max) call expand_temp_ghosts
           pghost_external = pghost_external + 1
           call create_external_ghost_particle(p,sph(p)%r,&
              & ghosts(pghost_external))
           ghosts(pghost_external)%r(3) = &
              & ghosts(pghost_external)%r(3) + periodic_size(3)
        end if
        if (sph(p)%r(3) > ghost_search_max(3) - GKR*sph(p)%h) then
           if (pghost_external == pghost_external_max) call expand_temp_ghosts
           pghost_external = pghost_external + 1
           call create_external_ghost_particle(p,sph(p)%r,&
              & ghosts(pghost_external))
           ghosts(pghost_external)%r(3) = &
              & ghosts(pghost_external)%r(3) - periodic_size(3)
        end if
     end do
     do p=1,pghost_external_prev
        porig = ghosts(p)%porig
        if (ghosts(p)%r(3) < ghost_search_min(3) + GKR*sph(porig)%h) then
           if (pghost_external == pghost_external_max) call expand_temp_ghosts
           pghost_external = pghost_external + 1
           call create_external_ghost_particle(porig,ghosts(p)%r,&
              & ghosts(pghost_external))
           ghosts(pghost_external)%r(3) = &
              & ghosts(pghost_external)%r(3) + periodic_size(3)
        end if
        if (ghosts(p)%r(3) > ghost_search_max(3) - GKR*sph(porig)%h) then
           if (pghost_external == pghost_external_max) call expand_temp_ghosts
           pghost_external = pghost_external + 1
           call create_external_ghost_particle(porig,ghosts(p)%r,&
              & ghosts(pghost_external))
           ghosts(pghost_external)%r(3) = &
              & ghosts(pghost_external)%r(3) - periodic_size(3)
        end if
     end do
#endif

#if defined(DEBUG_GHOST_PARTICLES)
     write(6,*) "pghost_external (domain ", d, ") : ",pghost_external
#endif

     allocate(ghost_sph_reference(d)%data(1:pghost_external))
     ghost_sph_reference(d)%data(1:pghost_external) = ghosts(1:pghost_external)

  end do
  ! ==========================================================================

  return

  contains
  
  subroutine expand_temp_ghosts
  
     allocate(tmp_ghosts(1:pghost_external_max*10))
     tmp_ghosts(1:pghost_external) = ghosts(1:pghost_external)
     deallocate(ghosts)
     pghost_external_max = pghost_external_max * 10
     allocate(ghosts(1:pghost_external_max))
     ghosts(1:pghost_external) = tmp_ghosts(1:pghost_external)
     deallocate(tmp_ghosts)
     return
  
  end subroutine expand_temp_ghosts
END SUBROUTINE search_external_ghost_particles
