! EXPAND.F90
! A. McLeod - 30/07/2008
! Subroutine to expand data arrays
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine expand(newsize)
  use definitions
  use particle_module
  use hydro_module
  use type_module
  use time_module
  use neighbour_module
#if defined(BH_TREE)
  use tree_module
#endif
#ifdef SINKS
  use sink_module
#ifdef SINKCORRECTION
  use sink_correction_module
#endif
#endif
#if defined(RAD_WS) || defined(RAD_WS_HYBRID)
  use Eos_module
#endif
  use interface_module, only : err_stop
#if defined(USE_MPI)
  use mpi_communication_module
#endif
  implicit none

  integer, intent(in)             :: newsize          ! New size to expand to
  integer                         :: oldsize          ! Old size to expand from
  type(sph_particle), allocatable :: scratch(:)       ! Array for particle data
#if defined(BH_TREE)
#if defined(SELF_GRAVITY)
  type(BHgrav_node), allocatable  :: BHgrav_temp(:)   ! Array for gravity tree
#endif
  type(BHhydro_node), allocatable :: BHhydro_temp(:)  ! BH hydro tree array
#endif
#if defined(NEIGHBOUR_LISTS)
  integer, allocatable :: pptot_temp(:)        ! number of neighbours
  integer, allocatable :: pplist_temp(:,:)     ! list of neighbours
#endif
  integer              :: ierr                 ! Error value for allocation

  debug_timing("EXPAND_ALL_ARRAYS")

#if defined(GHOST_PARTICLES)
  oldsize = ptot + pghost
#else
  oldsize = ptot
#endif
  if (newsize > pmax) then
#if defined(USE_MPI)
     write (6,'(A,I0,A,I0,A,I0)') "Expanding all arrays in rank ", rank, &
     & " from size ", oldsize, "(limit ", pmax, ") to size ", newsize
#else
     write (6,'(A,I0,A,I0,A,I0)') "Expanding all arrays from size ", &
     & oldsize, "(limit ", pmax, ") to size ", newsize
#endif
  else
#if defined(USE_MPI)
     write (6,'(A,I0,A,I0,A,I0)') "Shrinking all arrays in rank ", rank, &
     & " from size ", oldsize, "(limit ", pmax, ") to size ", newsize
#else
     write (6,'(A,I0,A,I0,A,I0)') "Shrinking all arrays from size ", &
     & oldsize, "(limit ", pmax, ") to size ", newsize
#endif
  end if

#if defined(F2003)
  allocate(scratch(1:newsize), stat=ierr)
#else
  allocate(scratch(1:oldsize), stat=ierr)
#endif
  if (ierr /= 0) call err_stop("Unable to allocate scratch in expand!")
  scratch(1:oldsize) = sph(1:oldsize)
#if defined(F2003)
  call move_alloc(scratch, sph)
#else
  deallocate(sph)
  allocate(sph(1:newsize), stat=ierr)
  if (ierr /= 0) call err_stop("Unable to allocate sph in expand!")
  sph(1:oldsize) = scratch(1:oldsize)
  deallocate(scratch)
#endif
  
#if defined(NEIGHBOUR_LISTS)
#if defined(F2003)
  allocate(pptot_temp(1:newsize), stat=ierr)
#else
  allocate(pptot_temp(1:oldsize), stat=ierr)
#endif
  if (ierr /= 0) call err_stop("Unable to allocate pptot_temp in expand!")
  pptot_temp(1:oldsize) = pptot(1:oldsize)
#if defined(F2003)
  call move_alloc(pptot_temp, pptot)
#else
  deallocate(pptot)
  allocate(pptot(1:newsize), stat=ierr)
  if (ierr /= 0) call err_stop("Unable to allocate pptot in expand!")
  pptot(1:oldsize) = pptot_temp(1:oldsize)
  deallocate(pptot_temp)
#endif
  
#if defined(F2003)
  allocate(pplist_temp(1:pp_limit,1:newsize), stat=ierr)
#else
  allocate(pplist_temp(1:pp_limit,1:oldsize), stat=ierr)
#endif
  if (ierr /= 0) call err_stop("Unable to allocate pplist_temp in expand!")
  pplist_temp(1:pp_limit,1:oldsize) = pplist(1:pp_limit,1:oldsize)
#if defined(F2003)
  call move_alloc(pplist_temp, pplist)
#else
  deallocate(pplist)
  allocate(pplist(1:pp_limit,1:newsize), stat=ierr)
  if (ierr /= 0) call err_stop("Unable to allocate pplist in expand!")
  pplist(1:pp_limit,1:oldsize) = pplist_temp(1:pp_limit,1:oldsize)
  deallocate(pplist_temp)
#endif
#endif
  
  ! Reset pmax to new size
  pmax = newsize

  ! --------------------------------------------------------------------------

#if defined(BH_TREE)
  debug2("Allocating tree arrays...")
  if (LEAFMAX > 8) then
     cmax_grav  = int((2*(pgas + pcdm)*PMAXMULT) / (LEAFMAX - 4))
     cmax_hydro = (2*pmax) / (LEAFMAX - 4)
  else
     cmax_grav  = int(2*(pgas + pcdm)*PMAXMULT)
     cmax_hydro = 2*pmax
  end if
  cmax_grav = max(ceiling(PMAXMULT*ctot_grav), cmax_grav)
  cmax_hydro = max(ceiling(PMAXMULT*ctot_hydro), cmax_hydro)
  cmax_skeleton = cmax_hydro

#if defined(SELF_GRAVITY)
#if defined(F2003)
  allocate(BHgrav_temp(0:cmax_grav), stat=ierr)
#else
  allocate(BHgrav_temp(0:ctot_grav), stat=ierr)
#endif
  if (ierr /= 0) call err_stop("Unable to allocate BHgrav_temp in expand!")
  BHgrav_temp(0:ctot_grav) = BHgrav(0:ctot_grav)
#if defined(F2003)
  call move_alloc(BHgrav_temp, BHgrav)
#else
  deallocate(BHgrav)
  allocate(BHgrav(0:cmax_grav), stat=ierr)
  if (ierr /= 0) call err_stop("Unable to allocate BHgrav in expand!")
  BHgrav(0:ctot_grav) = BHgrav_temp(0:ctot_grav)
  deallocate(BHgrav_temp)
#endif
#endif
#if defined(F2003)
  allocate(BHhydro_temp(0:cmax_hydro), stat=ierr)
#else
  allocate(BHhydro_temp(0:ctot_hydro), stat=ierr)
#endif
  if (ierr /= 0) call err_stop("Unable to allocate BHhydro_temp in expand!")
  BHhydro_temp(0:ctot_hydro) = BHhydro(0:ctot_hydro)
#if defined(F2003)
  call move_alloc(BHhydro_temp, BHhydro)
#else
  deallocate(BHhydro)
  allocate(BHhydro(0:cmax_hydro), stat=ierr)
  if (ierr /= 0) call err_stop("Unable to allocate BHhydro in expand!")
  BHhydro(0:ctot_hydro) = BHhydro_temp(0:ctot_hydro)
  deallocate(BHhydro_temp)
#endif

  deallocate(BHskeleton)
  allocate(BHskeleton(0:cmax_skeleton), stat=ierr)
  if (ierr /= 0) call err_stop("Unable to allocate BHskeleton in expand!")

  deallocate(cellof)
  deallocate(BHnextptcl)
  deallocate(whichchild)
  deallocate(BHstock)
  allocate(cellof(1:pmax), stat=ierr)
  if (ierr /= 0) call err_stop("Unable to allocate cellof in expand!")
  allocate(BHnextptcl(1:pmax), stat=ierr)
  if (ierr /= 0) call err_stop("Unable to allocate BHnextptcl in expand!")
  allocate(whichchild(1:pmax), stat=ierr)
  if (ierr /= 0) call err_stop("Unable to allocate whichchild in expand!")
  allocate(BHstock(0:cmax_hydro), stat=ierr)
  if (ierr /= 0) call err_stop("Unable to allocate BHstock in expand!")
#endif

  return
end subroutine expand


subroutine expand_ghost_tree(pghost_in)
  use definitions
  use particle_module
  use tree_module
  use interface_module, only : err_stop
  
  integer, intent(in)             :: pghost_in        ! Number of ghosts
                                                      ! (or estimate)
  type(BHhydro_node), allocatable :: BHghost_temp(:)  ! BH ghost tree array
  integer              :: ierr                 ! Error value for allocation
  
  if (LEAFMAX > 8) then
     cmax_ghost = int((2*pghost_in*PMAXMULT) / (LEAFMAX - 4))
  else
     cmax_ghost = int(2*pghost_in*PMAXMULT)
  end if
  cmax_ghost = max(cmax_ghost, ceiling(ctot_ghost*PMAXMULT))
  
#if defined(F2003)
  allocate(BHghost_temp(0:cmax_ghost), stat=ierr)
#else
  allocate(BHghost_temp(0:ctot_ghost), stat=ierr)
#endif
  if (ierr /= 0) then
     call err_stop("Unable to allocate BHghost_temp in expand_ghost!")
  end if
  BHghost_temp(0:ctot_ghost) = BHghost(0:ctot_ghost)
#if defined(F2003)
  call move_alloc(BHghost_temp, BHghost)
#else
  deallocate(BHghost)
  allocate(BHghost(0:cmax_ghost), stat=ierr)
  if (ierr /= 0) call err_stop("Unable to allocate BHghost in expand_ghost!")
  BHghost(0:ctot_ghost) = BHghost_temp(0:ctot_ghost)
  deallocate(BHghost_temp)
#endif

end subroutine expand_ghost_tree
