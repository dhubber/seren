! PARAMSTORE.F90
! C. P. Batty & D. A. Hubber - 1/4/2007
! Writes compiler flag and simulation parameters to file as record.  
! Gives all information required to reproduce a simulation. Also writes 
! important information about the columns of various ASCII output files.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE paramstore(store_file)
  use interface_module, only : write_debug_column_info,write_makefile_options
  use definitions
  use filename_module
  implicit none

  character(len=256), intent(in) :: store_file  ! params file name

  integer :: i                                  ! params counter

  debug1("Storing simulation information : "//trim(store_file)//" [paramstore.F90]")

  open(1, file=store_file, status="unknown", form="formatted")


! Makefile options
! ----------------------------------------------------------------------------
  call write_makefile_options(1)


! Parameter information
! ----------------------------------------------------------------------------
  write(1,*) "=========="
  write(1,*) "PARAMETERS"
  write(1,*) "=========="
  do i=1,nparams
     select case(params(i)%var_type)
     case('i'); write(1,*) trim(params(i)%var_name)," = ",params(i)%var_i
     case('j'); write(1,*) trim(params(i)%var_name)," = ",params(i)%var_j
     case('r'); write(1,*) trim(params(i)%var_name)," = ",params(i)%var_r
     case('d'); write(1,*) trim(params(i)%var_name)," = ",params(i)%var_d
     case('c'); write(1,*) trim(params(i)%var_name)," = ",trim(params(i)%var_c)
     case('u'); write(1,*) trim(params(i)%var_name)," = ",trim(params(i)%var_u)
     case('l'); write(1,*) trim(params(i)%var_name)," = ",params(i)%var_l
     case default; stop
     end select
  end do


! Debug plot data column index 
! ----------------------------------------------------------------------------
  call write_debug_column_info(1)


! Close file
  close(1)


  return
END SUBROUTINE paramstore
