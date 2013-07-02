! TURB_CHECK_TIME.F90
! A.McLeod - 26/11/2012
! Check whether we need to update to a new turbulent field
! Repeats if necessary (not necessary if using sensible options for fixed_dt)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE turb_check_time
   use definitions
   use time_module
   use turbulence_module
   implicit none

   do
      if (time >= turb_next_time) then
         call turb_next_field
      else
         exit
      end if
   end do
   
   return
  
END SUBROUTINE turb_check_time