! READ_TURB_FIELDS.F90
! A.McLeod - 24/11/2012
! Read turbulent fields from binary files for restarting
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_turb_fields
   use definitions
   use filename_module
   use turbulence_module
   use scaling_module
   implicit none
   
   character(len=16)    :: precision_str    ! String containing precision
   
   debug2("Reading turbulence fields [read_turb_fields.F90]")
   
#if defined(DOUBLE_PRECISION)
   precision_str = 'DOUBLE_PRECISION'
#else
   precision_str = 'SINGLE_PRECISION'
#endif
   
   turb_file_header = trim(run_dir)//'turb_fields.dat'
   
   ! Read header file of important information
   open(1,file=turb_file_header,status="old",form="formatted")
   read(1,*) nturbtemp
   read(1,*) precision_str
   read(1,*) turb_file_last
   read(1,*) turb_last_time
   read(1,*) turb_file_next
   read(1,*) turb_next_time
   close(1)
   ! We don't care about the rest of the stuff in the file
   ! as it should be set in the params file
   
   turb_last_time = turb_last_time / tscale
   turb_next_time = turb_next_time / tscale
   
   ! Read turbulent fields
   open(1,file=trim(run_dir)//trim(turb_file_last),status="unknown",&
        &access="stream",form="unformatted")
   read(1) turb_last
   close(1)
   
   open(1,file=trim(run_dir)//trim(turb_file_next),status="unknown",&
        access="stream",form="unformatted")
   read(1) turb_next
   close(1)
   
   ! No need to change files yet
   turb_changed = .FALSE.

   return
  
END SUBROUTINE read_turb_fields