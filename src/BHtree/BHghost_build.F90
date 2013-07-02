! BHGHOST_BUILD.F90
! D. A. Hubber - 18/03/2011
! Builds Barnes-Hut tree for a given distribution of particles. Only builds 
! tree for ghost particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHghost_build
  use definitions
  use tree_module
  use particle_module
  use type_module
  !use mpi_communication_module
  implicit none

  integer :: c                             ! cell counter
  integer :: i                             ! aux. counter
  integer :: nptcls                        ! no. of particles
  integer :: p                             ! particle counter
  integer, allocatable :: plist(:)         ! particle list
  real(kind=PR), allocatable :: r(:,:)     ! particle positions

  debug2("[BHghost_build.F90]")

  allocate(plist(1:pghost))
  allocate(r(1:NDIM,1:pghost))

  nptcls = pghost

  ! Create particle id list and particle positions  (Not using 
  ! create_particle_list.F90 now since we don't have a ghost mask yet).
  do i=1,pghost
     p = ptot + i
     plist(i) = p
     r(1:NDIM,i) = sph(p)%r(1:NDIM)
  end do

  ! Build skeleton of tree using reduced particle list
  call BH_build_skeleton(nptcls,plist(1:nptcls),r(1:NDIM,1:nptcls))

  ! Now copy skeleton of tree to main ghost tree
  ctot_ghost = ctot_skeleton
  ltot_ghost = ltot_skeleton

  first_cell_ghost(0:ltot_ghost) = first_cell_skeleton(0:ltot_ghost)
  last_cell_ghost(0:ltot_ghost) = last_cell_skeleton(0:ltot_ghost)

  do c=0,ctot_ghost
     BHghost(c)%leaf = BHskeleton(c)%leaf
     BHghost(c)%nextcell = BHskeleton(c)%nextcell
     BHghost(c)%ifopen = BHskeleton(c)%ifopen
     BHghost(c)%plist(1:LEAFMAX) = BHskeleton(c)%plist(1:LEAFMAX)
  end do

#if defined(DEBUG_GHOST_PARTICLES)
  write(6,*) "ltot_ghost : ",ltot_ghost,"   ctot_ghost : ",ctot_ghost
#endif

  deallocate(r)
  deallocate(plist)

  return 
END SUBROUTINE BHghost_build
