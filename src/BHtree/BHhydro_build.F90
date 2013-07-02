! BHHYDRO_BUILD.F90
! D. A. Hubber - 22/01/2008
! Builds Barnes-Hut tree for a given distribution of particles. Only builds 
! tree for hydro particles (i.e. boundary, icm and 
! self-gravitating gas particles). 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHhydro_build
  use interface_module, only : create_particle_list
  use definitions
  use tree_module
  use particle_module
  use type_module
  implicit none

  integer :: c                             ! cell counter
  integer :: i                             ! aux. counter
  integer :: nptcls                        ! no. of particles
  integer :: p                             ! particle counter
  integer, allocatable :: plist(:)         ! particle list
  real(kind=PR), allocatable :: r(:,:)     ! particle positions
#if defined(REORDER_TREE)
  integer :: caux                          ! aux. tree cell counter
#endif

  allocate(plist(1:ptot))
  allocate(r(1:NDIM,1:ptot))

  ! Create list of self-gravitating particles
  call create_particle_list(nptcls,plist,sphmask)

  ! Create list of particle positions
  do i=1,nptcls
     p = plist(i)
     r(1:NDIM,i) = sph(p)%r(1:NDIM)
  end do

  ! Build skeleton of tree using reduced particle list
  call BH_build_skeleton(nptcls,plist(1:nptcls),r(1:NDIM,1:nptcls))

  ! Now copy skeleton of tree to main hydro tree
  ctot_hydro = ctot_skeleton
  ltot_hydro = ltot_skeleton

  first_cell_hydro(0:ltot_hydro) = first_cell_skeleton(0:ltot_hydro)
  last_cell_hydro(0:ltot_hydro) = last_cell_skeleton(0:ltot_hydro)

  ! Reorder tree to walk-order
#if defined(REORDER_TREE)
  call BH_reorder_tree(hydroorder)
  do c=0,ctot_hydro
     caux = hydroorder(c)
     BHhydro(caux)%leaf = BHskeleton(c)%leaf
     BHhydro(caux)%nextcell = hydroorder(BHskeleton(c)%nextcell)
     BHhydro(caux)%ifopen = hydroorder(BHskeleton(c)%ifopen)
     BHhydro(caux)%plist(1:LEAFMAX) = BHskeleton(c)%plist(1:LEAFMAX)
  end do
#else
  do c=0,ctot_hydro
     BHhydro(c)%leaf = BHskeleton(c)%leaf
     BHhydro(c)%nextcell = BHskeleton(c)%nextcell
     BHhydro(c)%ifopen = BHskeleton(c)%ifopen
     BHhydro(c)%plist(1:LEAFMAX) = BHskeleton(c)%plist(1:LEAFMAX)
  end do
#endif

  deallocate(r)
  deallocate(plist)

  return 
END SUBROUTINE BHhydro_build
