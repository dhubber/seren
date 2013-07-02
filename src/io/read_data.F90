! READ_DATA.F90
! C. P. Batty & D. A. Hubber - 11/9/2007
! Calls various subroutines to read data from snapshot files.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_data(filename,file_form,decomp_read)
  implicit none

  character(len=*), intent(in) :: filename    ! snapshot name
  character(len=*), intent(in) :: file_form   ! snapshot format
  logical, intent(in), optional :: decomp_read ! MPI decomposition read
  
  logical                  :: decomp_read_pass

  debug2("Calling read data subroutines [read_data.F90]")

#if defined(USE_MPI)
  decomp_read_pass = .FALSE.
  if (present(decomp_read)) then
     if (decomp_read) decomp_read_pass = .TRUE.
  end if
#else
  decomp_read_pass = .FALSE.
#endif

  select case(file_form)

! DRAGON data format
! ----------------------------------------------------------------------------
#if defined(DRAGON_INPUT)
  case ("dragon_form","df")
     call read_data_dragon_form(filename,decomp_read_pass)
  case ("dragon_unform","du")
     call read_data_dragon_unform(filename,decomp_read_pass)
#endif

! SEREN data format
! ----------------------------------------------------------------------------
#if defined(SEREN_INPUT)
  case ("seren_form","sf")
     call read_data_seren_form(filename,decomp_read_pass)
  case ("seren_unform","su")
     call read_data_seren_unform(filename,decomp_read_pass)
#endif

! Simple column data ASCII format
! ----------------------------------------------------------------------------
#if defined(ASCII_INPUT)
  case ("ascii")
     call read_data_ascii(filename,decomp_read_pass)
#endif

! If no recognised file format
! ----------------------------------------------------------------------------
  case default
     write(6,*) "Input file format unrecognised or not selected in Makefile :"&
          &,trim(file_form)
     stop
  end select


  return
END SUBROUTINE read_data
