! BH_REORDER_TREE.F90
! D. A. Hubber - 16/04/2011
! Reorder the tree cells so they are in 'walk' order in the main arrays.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BH_reorder_tree(neworder)
  use particle_module
  use hydro_module
  use time_module
  use type_module
  use tree_module
  use filename_module
#if defined(HEALPIX)
  use HP_module, only : HPtot
#endif
  implicit none
  
  integer, intent(out) :: neworder(0:cmax_skeleton)  ! New ids of cells
  integer :: c                                       ! Cell counter
  integer :: caux                                    ! Aux. cell counter
  integer, allocatable :: dummylist(:)               ! List of cells in order

  debug2("Order BH tree to walk-order [BH_reorder_tree.F90]")

  allocate(dummylist(0:cmax_skeleton))

  do c=0,cmax_skeleton
     neworder(c) = c
  end do

! Make dummy list of particles in tree order
! ============================================================================
  c = 0
  caux = 0
  do
     !write(6,*) "Reordering tree ids : ",c,caux,ctot_skeleton
     dummylist(caux) = c
     if (BHskeleton(c)%leaf > 0) then
        c = BHskeleton(c)%nextcell
     else if (BHskeleton(c)%leaf == 0) then
        c = BHskeleton(c)%ifopen
     else
        c = BHskeleton(c)%nextcell
     end if
     if (c > ctot_skeleton) exit
     caux = caux + 1
  end do
! ============================================================================

  do c=0,ctot_skeleton
     neworder(dummylist(c)) = c
  end do

! Free allocated memory
  deallocate(dummylist)

  return
END SUBROUTINE BH_reorder_tree
