! CLEAN_UP.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Frees up memory after simulation finishes (i.e. deallocates main arrays). 
! All variables should be deallocated in reverse order to which they are 
! originally allocated in allocate_memory.F90 (or in other routines).   
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE clean_up
  use particle_module
  use neighbour_module
  use hydro_module
  use kernel_module
  use tree_module
  use Nbody_module
  use time_module
  use sink_module
#if defined(RAD_WS)
  use Eos_module
#endif
#if defined(HEALPIX)
  use HP_module
#endif
#if defined(USE_MPI)
  use mpi_communication_module
#endif
  implicit none

  debug2("Deallocating all global arrays [clean_up.F90]")
  debug_timing("CLEAN_UP")


! Deallocate all sink arrays only if we are not using the N-body integrator
! (If using N-body integrator, these arrays are deallocated instead in 
! nbody_clean_up.F90)
! ----------------------------------------------------------------------------
#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  if (allocated(star)) deallocate(star)
#if defined(BINARY_STATS)
  if (allocated(binary)) deallocate(binary)
#endif
#endif

! Sink arrays
! ----------------------------------------------------------------------------
#if defined(SINKS) || defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  if (allocated(sink)) deallocate(sink)
#endif

! Kernel arrays
! ----------------------------------------------------------------------------
#if defined(KERNEL_TABLES)
  if (allocated(wzetatable)) deallocate(wzetatable)
  if (allocated(womegatable)) deallocate(womegatable)
  if (allocated(wpottable)) deallocate(wpottable)
  if (allocated(wgravtable)) deallocate(wgravtable)
  if (allocated(w2table)) deallocate(w2table)
  if (allocated(w1table)) deallocate(w1table)
  if (allocated(w0table)) deallocate(w0table)
#endif


! Deallocate all particle arrays if allocated (just check main array here)
! ============================================================================
  if (allocated(sph)) then

     ! Ionization routines
     ! -----------------------------------------------------------------------
#if defined(HEALPIX)
#if defined(DEBUG_HP_WALK_ALL_RAYS)
     if (allocated(whichHPlevel)) deallocate(whichHPlevel)
#endif
     if (allocated(HPray)) deallocate(HPray)
#endif

     ! Radiative transport approximation arrays
     ! -----------------------------------------------------------------------
#if defined(RAD_WS) && defined(DEBUG_RAD)
     if (allocated(rad_info)) deallocate(rad_info)
#endif

     ! Tree arrays
     ! -----------------------------------------------------------------------
#if defined(BH_TREE)
     if (allocated(BHstock)) deallocate(BHstock)
     if (allocated(whichchild)) deallocate(whichchild)
     if (allocated(BHnextptcl)) deallocate(BHnextptcl)
     if (allocated(cellof)) deallocate(cellof)

     if (allocated(BHskeleton)) deallocate(BHskeleton)
     if (allocated(BHhydro)) deallocate(BHhydro)
#if defined(GHOST_PARTICLES)
     if (allocated(BHghost)) deallocate(BHghost)
#endif
#if defined(SELF_GRAVITY)
     if (allocated(BHgrav)) deallocate(BHgrav)
#endif

#if defined(USE_MPI) && defined(SELF_GRAVITY)
     if (allocated(BHlocal_grav)) deallocate(BHlocal_grav)
     if (allocated(BHremote_grav)) deallocate(BHremote_grav)
     if (allocated(treesend)) deallocate(treesend)
     if (allocated(treerecv)) deallocate(treerecv)
#endif

#endif

#if defined(NEIGHBOUR_LISTS)
     if (allocated(pplist)) deallocate(pplist)
     if (allocated(pptot)) deallocate(pptot)
#endif
     if (allocated(parray)) deallocate(parray)
     if (allocated(sph)) deallocate(sph)
     if (allocated(minimal_sph)) deallocate(minimal_sph)

  end if
! ============================================================================


  return
END SUBROUTINE clean_up
