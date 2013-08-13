! TREE_UPDATE.F90
! D. A. Hubber - 21/9/2007
! Control subroutine to call tree build and tree stock routines.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE tree_update(nbuild,nstock)
#if defined(BH_TREE)
  use interface_module, only : BHhydro_stock
#endif
  use particle_module, only : ptot
  use time_module, only : nsteps
  use tree_module
  use type_module
  implicit none

  integer(kind=ILP), intent(in) :: nbuild    ! Tree build step no.
  integer(kind=ILP), intent(in) :: nstock    ! Tree stock step no.


! Barnes-Hut tree subroutine calls
! ============================================================================
#if defined(BH_TREE)
  debug2("Updating BH tree [tree_update.F90]")
  debug_timing("BH_TREE")

! Build new hydro tree
! ----------------------------------------------------------------------------
  if (nbuild == nsteps .or. build_tree) then
     call BHhydro_build
#if defined(REORDER_PARTICLES)
     call BH_reorder_particles
#endif
  end if

! Stock trees on every acceleration step
  call BHhydro_stock(cmax_hydro,ctot_hydro,&
       &ltot_hydro,first_cell_hydro,last_cell_hydro,BHhydro)


! Build new gravity tree
! ----------------------------------------------------------------------------
#if defined(SELF_GRAVITY)
  if (nbuild == nsteps .or. build_tree) then

     ! If all particles are gas or cdm particles, copy the hydro tree 
     ! to the gravity tree.  Else, build the gravity tree from scratch.
     if (pgas + pcdm == ptot) then
        call copy_BHhydro_to_BHgrav
     else
        call BHgrav_build
     end if

  end if

! Stock trees on every acceleration step
  call BHgrav_stock
#endif


#endif
! ============================================================================


  return
END SUBROUTINE tree_update
