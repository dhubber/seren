! NBODY_END_STEP.F90
! D. A. Hubber - 24/11/2011
! Set all quantities at the end of a star's timestep after all P(EC)^n 
! iterative steps have been completed.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_end_step
  use Nbody_module
  use sink_module
  use time_module
  implicit none

  integer :: s                    ! Star particle counter
  real(kind=PR) :: dt             ! Physical time since beginning of timestep

! Loop over all stars in simulation
! ----------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED)
  do s=1,stot
     star(s)%accdo = .false.

     ! If at end of the step, record star properties
     ! -----------------------------------------------------------------------
     if (n == star(s)%nlast) then

        ! Record 'old' values of star properties
        star(s)%rold(1:NDIM)  = star(s)%r(1:NDIM)
        star(s)%vold(1:NDIM)  = star(s)%v(1:NDIM)
        star(s)%a0(1:NDIM)    = star(s)%a(1:NDIM)
        star(s)%adot0(1:NDIM) = star(s)%adot(1:NDIM)
#if defined(FORCE_SPLITTING)
        if (sphminstep) star(s)%accsph = .true.
        star(s)%a0star(1:NDIM)    = star(s)%astar(1:NDIM)
        star(s)%adot0star(1:NDIM) = star(s)%adotstar(1:NDIM)
        star(s)%a0sph(1:NDIM)    = star(s)%asph(1:NDIM)
        star(s)%adot0sph(1:NDIM) = star(s)%adotsph(1:NDIM)
#endif
     end if
     ! -----------------------------------------------------------------------

  end do
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------------

  return
END SUBROUTINE nbody_end_step
