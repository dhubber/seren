! ERRORS.F90
! A. McLeod - 25/02/2013
! Allows a simple interface for errors that leave a return value in the shell
! ============================================================================

#include "macros.h"
! Next bit only needed while people still on F95 compilers!
#if !defined(F2003)
#define ERROR_UNIT 0
#endif

! ============================================================================
SUBROUTINE err_stop(message)
#if defined(F2003)
  use, intrinsic :: ISO_FORTRAN_ENV
#endif
  implicit none
  character (LEN=*), intent(in), optional    :: message      ! Error message
  
  integer                                    :: err_unit     ! Unit to write to
  
  err_unit = ERROR_UNIT
  
  if (present(message)) then
     write (err_unit, *) message
  end if
  
  stop 1
  
  return
END SUBROUTINE err_stop
