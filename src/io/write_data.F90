! WRITE_DATA.F90
! D. A. Hubber - 10/8/2007
! Calls various subroutines to write snapshot files
! ============================================================================ 

#include "macros.h"

! ============================================================================ 
SUBROUTINE write_data(filename,file_form)
  implicit none

  character(len=*), intent(in) :: filename    ! snapshot name
  character(len=*), intent(in) :: file_form   ! snapshot format

  debug2("Calling write data subroutines [write_data.F90]")

  select case(file_form)

! DRAGON data format
! ----------------------------------------------------------------------------
#if defined(DRAGON_OUTPUT)
  case ("dragon_form","df")
     call write_data_dragon_form(filename) 
  case ("dragon_unform","du")
     call write_data_dragon_unform(filename)
#endif

! SEREN data format
! ----------------------------------------------------------------------------
#if defined(SEREN_OUTPUT)
  case ("seren_form","sf")
     call write_data_seren_form(filename)
  case ("seren_unform","su")
     call write_data_seren_unform(filename)
#endif

! Simple column data ASCII format
! ----------------------------------------------------------------------------
#if defined(ASCII_OUTPUT)
  case ("ascii")
     call write_data_ascii(filename)
#endif

! If no recognised file format
! ----------------------------------------------------------------------------
  case default
     write(6,*) "Output file format unrecognised or not selected in Makefile"&
          &,trim(file_form)
  end select

  return
END SUBROUTINE write_data
