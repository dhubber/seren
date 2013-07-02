! UNWRAP_PARTICLE_POSITION.F90
! A. McLeod - 23/08/12
! With some MPI periodic options, particles may not always be within periodic
! boxes. This subroutine unwraps particle positions for output purposes etc.
! ============================================================================

#include "macros.h"

! ============================================================================

subroutine unwrap_particle_position(rp)
   use particle_module
   use periodic_module
   implicit none
   real(kind=PR), intent(inout) :: rp(1:NDIM)  ! Particle position
   
#if defined(USE_MPI) && defined(PERIODIC)

#if defined(PERIODIC_X)
   if (rp(1) < periodic_min(1)) then
      rp(1) = rp(1) + periodic_size(1)
   else if (rp(1) > periodic_max(1)) then
      rp(1) = rp(1) - periodic_size(1)
   end if
#endif

#if defined(PERIODIC_Y) && (NDIM == 2 || NDIM == 3)
   if (rp(2) < periodic_min(2)) then
      rp(2) = rp(2) + periodic_size(2)
   else if (rp(2) > periodic_max(2)) then
      rp(2) = rp(2) - periodic_size(2)
   end if
#endif

#if defined(PERIODIC_Z) && (NDIM == 3)
   if (rp(3) < periodic_min(3)) then
      rp(3) = rp(3) + periodic_size(3)
   else if (rp(3) > periodic_max(3)) then
      rp(3) = rp(3) - periodic_size(3)
   end if
#endif

#else
   continue ! Do nothing
#endif
   
   return
end subroutine unwrap_particle_position