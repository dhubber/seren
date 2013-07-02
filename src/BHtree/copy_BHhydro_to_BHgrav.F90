! COPY_BHHYDRO_TO_BHGRAV.F90
! D. A. Hubber - 19/7/2008
! Copies tree structure of BHgrav tree to BHhydro tree.  Used when all 
! of the SPH particles are self-gravitating gas particles therefore 
! avoiding duplicating an identical tree build for both hydro and gravity.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE copy_BHhydro_to_BHgrav
  use tree_module
  implicit none

  integer :: c                     ! Cell counter

  debug2("Copying hydro to grav tree [copy_BHhydro_to_BHgrav.F90]")

  ctot_grav = ctot_hydro
  ltot_grav = ltot_hydro
  
  if (cmax_grav < ctot_hydro) then
     debug2("Expanding gravity tree before copy")
     cmax_grav = cmax_hydro
     deallocate(BHgrav)
     allocate(BHgrav(0:cmax_grav))
  end if

! Copy level information for stocking
  first_cell_grav(0:ltot_hydro) = first_cell_hydro(0:ltot_hydro)
  last_cell_grav(0:ltot_hydro)  = last_cell_hydro(0:ltot_hydro)

! Now copy all necessary lists and variables
!$OMP PARALLEL DO DEFAULT(SHARED)
  do c=0,ctot_hydro
     BHgrav(c)%leaf = BHhydro(c)%leaf
     BHgrav(c)%nextcell = BHhydro(c)%nextcell
     BHgrav(c)%ifopen = BHhydro(c)%ifopen
     BHgrav(c)%plist(1:LEAFMAX) = BHhydro(c)%plist(1:LEAFMAX)
  end do
!$OMP END PARALLEL DO

  return
END SUBROUTINE copy_BHhydro_to_BHgrav
