! BHGRAV_SLICE.F90
! A. McLeod - 23/07/2008
! Slices off the top layers of the BH Gravity tree
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHgrav_slice()
  use mpi_communication_module
  use particle_module
  use tree_module
  implicit none

  integer           :: c                      ! Cell identifier

  debug_timing("BH_GRAV_SLICE")
  debug2("Slicing off top layers of gravity tree [BHgrav_slice.F90]")
  
  ctot_localgrav = last_cell_grav(min(remotetreedepth,ltot_grav))
  
  BHlocal_grav(0:ctot_localgrav) = BHgrav(0:ctot_localgrav)

   do c=0,ctot_localgrav
      ! Test if nextcell to end of tree is sufficiently large (must be > cmax_remotegrav)
!       if (BHlocal_grav(c)%nextcell > ctot_localgrav) then
!          BHlocal_grav(c)%nextcell = cmax_remotegrav + 1
!       end if
      ! No longer required - nextcell would be LARGEST_INT for 'finished'
      ! Test if cell is a leaf cell
      if (BHlocal_grav(c)%leaf > 0) then
         ! Leaf cell, set ifopen = -1
         BHlocal_grav(c)%ifopen = -1
      else
         ! Test if opening cell takes us outside of top remotetreedepth levels
         if (BHlocal_grav(c)%ifopen > ctot_localgrav) then
            BHlocal_grav(c)%ifopen = -1 ! Needs exporting to open, set ifopen = -1
         end if
      end if

   end do

   return
END SUBROUTINE BHgrav_slice
