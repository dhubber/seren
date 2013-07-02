! TURB_NEXT_FIELD.F90
! A.McLeod - 24/11/2012
! Generate next turbulent field
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE turb_next_field
   use definitions
   use time_module
   use turbulence_module
   use turb_force_module
   implicit none
   
   debug_timing("TURB_NEXT_FIELD")
   debug2("Calculating next turbulent field [turb_next_field.F90]")

MPI_ROOT
   ! Set turbulent field times
   turb_last_time = turb_next_time
   turb_next_time = turb_last_time + turb_dt

   ! Copy current field, and add more turbulence
   turb_last = turb_next
   afield_last = afield_next
   call add_turbulence(turb_next, turb_dt)

   ! Subtract decay of turbulence
   call decay_turbulence(turb_last, turb_next)

   ! Fourier transform
   call FFT_3D(turb_next, afield_next)
   afield_next = afield_next * turb_norm * turb_rms

   ! Set turb_changed to trigger output at next snapshot
   turb_changed = .TRUE.
MPI_END

#if defined(USE_MPI)
   call mpi_share_turb_fields
#endif

   return
  
END SUBROUTINE turb_next_field