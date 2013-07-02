! MPI_SCATTER_GHOSTS.F90
! A. McLeod - 18/11/2011
! Calculate size of active particle box and active particle smoothing
! kernel bounding box, and use this to send gather and scatter ghosts
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE mpi_scatter_ghosts
   use interface_module, only : active_particle_list,get_neib
   use particle_module
!    use neighbour_module
   use type_module
!    use timing_module
! #if defined(OPENMP)
!    use omp_lib
! #endif
   use mpi_communication_module
   use domain_comparisons
   use tree_module
!    use mpi
   implicit none
! 
   integer :: acctot                   ! No. of particles on acc. step
   integer :: acchydrotot              ! No. of active hydro particles
   integer :: i                        ! Aux. particle counter
   integer :: p                        ! particle counter
   integer, allocatable :: acclist(:)  ! List of particles on acc. step
#if defined(OPENMP)
   integer :: chunksize                ! Loop chunksize for OMP
#endif

   debug2("Exporting scatter ghosts [mpi_sph_update.F90]")

   debug_timing("MPI_SCATTER_GHOSTS")
   
! Prepare list of particles which are on an acceleration step
! ----------------------------------------------------------------------------
   allocate(acclist(1:ptot))
   call active_particle_list(acctot,acclist,sphmask)
#if defined(OPENMP)
   chunksize = int(CHUNKFRAC*real(acctot,PR)/real(omp_get_max_threads(),PR)) + 1
#endif

  ! Find a new set of ghost particles, including scatter neighbours,
  ! for use with hydro and other routines
  ! -----------------------------------------------------------------------

! Update (periodic) ghost particle properties
! -----------------------------------------------------------------------
   pghost = pperiodic

#if defined(BH_TREE) && defined(PERIODIC)
   ! Restore neighbour tree of ONLY periodic ghost particles
   ! otherwise we get references to old MPI ghosts in get_ghosts' search
!       call BHghost_build
!       call BHhydro_stock(cmax_ghost,ctot_ghost,&
!          &ltot_ghost,first_cell_ghost,last_cell_ghost,BHghost)
   BHghost(0:ctot_periodic) = BHperiodic
   deallocate(BHperiodic)
   ctot_ghost = ctot_periodic
   ltot_ghost = ltot_periodic
   first_cell_ghost = first_cell_periodic
   last_cell_ghost = last_cell_periodic
   if (pghost > 0) then
      call BHhydro_update_hmax(cmax_ghost,ctot_ghost,&
         &ltot_ghost,first_cell_ghost,last_cell_ghost,BHghost)
   end if
#endif
   
   ! Find final size of bounding boxes containing only active particles
   call maxmin_hydro_filter(bbmax(1:NDIM,rank),bbmin(1:NDIM,rank),&
         & acclist(1:acctot),acctot,.TRUE.)
   call maxmin_hydro_filter(activemax(1:NDIM,rank),activemin(1:NDIM,rank),&
         & acclist(1:acctot),acctot,.FALSE.)
   debug_timing("MPI_SCATTER_GHOSTS")
   
   ! Share the active bounding boxes with other tasks
   call broadcastboundingboxes(bbmin,bbmax)
   call broadcastboundingboxes(activemin,activemax)
   
#if defined(PERIODIC)
   ! Find externally required periodic ghost particles
   call search_external_ghost_particles
#endif
   
   call get_ghosts(.TRUE.)
   last_pghost = max(pghost,last_pghost)
   ! Save the current number of ghosts (local periodic and MPI ghosts for
   ! scatter), if larger than the number of ghosts used for h_gather, to
   ! estimate memory requirements in transfer_particles

#if defined(BH_TREE)
   debug_timing("GHOST_TREE 2")
   ! Rebuild neighbour tree of ghost particles
   if (pghost > pperiodic) then
      call BHghost_build
      call BHhydro_stock(cmax_ghost,ctot_ghost,&
         &ltot_ghost,first_cell_ghost,last_cell_ghost,BHghost)
      call BHhydro_update_hmax(cmax_ghost,ctot_ghost,&
         &ltot_ghost,first_cell_ghost,last_cell_ghost,BHghost)
   end if
#endif

   ! Rebuild neighbour lists for all particles
   ! -----------------------------------------------------------------------
#if defined(NEIGHBOUR_LISTS)
   debug2("Calculating neighbour lists for all particles")
   debug_timing("GET_NEIB")

   ! Reconstruct active particle list for hydro particles only
   call active_particle_list(acchydrotot,acclist,hydromask)
#if defined(OPENMP)
   chunksize = int(CHUNKFRAC*real(acchydrotot,PR)/&
        &real(omp_get_max_threads(),PR)) + 1
#endif

   if (acchydrotot > 0) then
      !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) &
      !$OMP DEFAULT(SHARED) PRIVATE(p)
      do i=1,acchydrotot
         p = acclist(i)
         call get_neib(p)
      end do
      !$OMP END PARALLEL DO
   end if
#endif

     ! Mean-h zeta for self-gravity
     ! -----------------------------------------------------------------------
#if defined(GRAD_H_SPH) && defined(SELF_GRAVITY) && defined(MEANH_GRAVITY)
   debug2("Calculating meanh zeta for all particles")
   debug_timing("MEANH_ZETA")

   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p)
   do i=1,acctot
      p = acclist(i)
      call meanh_zeta(p)
   end do
   !$OMP END PARALLEL DO
#endif
     ! -----------------------------------------------------------------------

! ============================================================================

   if (allocated(acclist)) deallocate(acclist)

  return
END SUBROUTINE mpi_scatter_ghosts
