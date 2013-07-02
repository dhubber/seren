! WRITE_TURB_FIELDS.F90
! A.McLeod - 24/11/2012
! Write turbulent fields to binary files for restarting
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_turb_fields
   use definitions
   use filename_module
   use turbulence_module
   use scaling_module
   implicit none
   
   character(len=8)     :: file_ext       ! String version of nturbtemp
   character(len=16)    :: precision_str  ! String containing precision
   character(len=50000) :: file_buffer    ! Buffer for header file
   character(len=1)     :: c              ! Mostly pointless variable
   
   ! Don't bother writing out files if nothing has changed
   if (.NOT. turb_changed) return
   
   debug2("Writing turbulence fields [write_turb_fields.F90]")
   
#if defined(DOUBLE_PRECISION)
   precision_str = 'DOUBLE_PRECISION'
#else
   precision_str = 'SINGLE_PRECISION'
#endif
   
   write(file_ext,"(I0)") nturbtemp
   
   turb_file_last = 'turb_last.'//trim(file_ext)//'.dat'
   turb_file_next = 'turb_next.'//trim(file_ext)//'.dat'
   turb_file_header = trim(run_dir)//'turb_fields.dat'

   ! Write turbulent fields
   open(1,file=trim(run_dir)//trim(turb_file_last),status="unknown",&
        &access="stream",form="unformatted")
   write(1) turb_last
   close(1)
   
   open(1,file=trim(run_dir)//trim(turb_file_next),status="unknown",&
        access="stream",form="unformatted")
   write(1) turb_next
   close(1)
   
   ! Write header file of important information
   write(file_buffer,*) nturbtemp, NEW_LINE(c),&
                      & precision_str, NEW_LINE(c),&
                      & trim(turb_file_last), NEW_LINE(c),&
                      & turb_last_time*tscale, NEW_LINE(c),&
                      & trim(turb_file_next), NEW_LINE(c),&
                      & turb_next_time*tscale, NEW_LINE(c),&
                      & TURB_GS, TURB_GS, TURB_GS, NEW_LINE(c),&
                      & turb_min*real(rscale,PR), NEW_LINE(c),&
                      & turb_max*real(rscale,PR), NEW_LINE(c),&
                      & turb_T*tscale, turb_dt*tscale, NEW_LINE(c),&
                      & comp_frac, sol_frac, NEW_LINE(c)
   
   ! Now write atomically (safest)
   open(1,file=turb_file_header,status="unknown",access='stream',&
        &form="formatted")
   write(1,*) trim(file_buffer)
   close(1)
   
   ! Update record-keeping on what files have been output, and that
   ! we have now made an output
   if (nturbtemp==1) then
      nturbtemp = 2
   else
      nturbtemp = 1
   end if
   turb_changed = .FALSE.

   return
  
END SUBROUTINE write_turb_fields