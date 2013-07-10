! QUIT_SEREN.F90
! D. A. Hubber - 17/09/2011
! Print error message and then quits the program.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE quit_seren(errmsg)
  implicit none

  character(len=*), intent(in) :: errmsg     ! Error message

  write(6,*) errmsg
  stop

  return
END SUBROUTINE quit_seren
