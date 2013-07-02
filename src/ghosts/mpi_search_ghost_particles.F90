! MPI_SEARCH_GHOST_PARTICLES.F90
! A. McLeod - 09/05/2011
! Search through all particles in the periodic box to find which ones must be 
! replicated as periodic ghost particles. Then find particles needed as
! periodic ghost particles for other domains.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE mpi_search_ghost_particles(preemptive_pass)
   use particle_module
   use periodic_module
   use type_module
   use mpi_communication_module
   implicit none

   logical, intent(in) :: preemptive_pass         ! Is this a preemptive pass?
!    real(kind=PR) :: all_bbmin(1:NDIM)             ! min(bbmin) across tasks
!    real(kind=PR) :: all_bbmax(1:NDIM)             ! max(bbmax) across tasks
  
   debug2("Search for new ghosts [mpi_search_ghost_particles.F90]")
   debug_timing("GHOSTS PERIODIC")
   
!    all_bbmin = minval(bbmin,dim=2)
!    all_bbmax = maxval(bbmax,dim=2)

   ! Find externally required ghost particles
   call search_external_ghost_particles
   
   ! Set up to find locally needed periodic ghost particles
#if defined(PERIODIC_X)
   ghost_search_min(1) = bbmax(1,rank) - periodic_size(1)
   ghost_search_max(1) = bbmin(1,rank) + periodic_size(1)
#endif
#if defined(PERIODIC_Y)
   ghost_search_min(2) = bbmax(2,rank) - periodic_size(2)
   ghost_search_max(2) = bbmin(2,rank) + periodic_size(2)
#endif
#if defined(PERIODIC_Z)
   ghost_search_min(3) = bbmax(3,rank) - periodic_size(3)
   ghost_search_max(3) = bbmin(3,rank) + periodic_size(3)
#endif

   ! Find ghosts
   call search_ghost_particles

#if defined(BH_TREE)
   if (pghost > 0) then
      call BHghost_build
      call BHhydro_stock(cmax_ghost,ctot_ghost,&
           &ltot_ghost,first_cell_ghost,last_cell_ghost,BHghost)
   end if
   
   if (.NOT. preemptive_pass) then
      ! Save copy of the periodic ghost tree, we can restore it later
      ctot_periodic = ctot_ghost
      ltot_periodic = ltot_ghost
      first_cell_periodic = first_cell_ghost
      last_cell_periodic = last_cell_ghost
      allocate(BHperiodic(0:ctot_periodic))
      BHperiodic(0:ctot_periodic) = BHghost(0:ctot_periodic)
   end if
#endif

   return
END SUBROUTINE mpi_search_ghost_particles
