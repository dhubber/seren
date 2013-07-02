! HP_INVERSE_POSITIONS
! T. Bisbas - 15/9/2008
! Converts vector from rotated frame to original HEALPix frame.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE HP_inverse_positions(rvec,isource)
  use definitions
  use HP_module, only : HPsource
  implicit none

  integer, intent(in) :: isource                  ! id of ionization source
  real(kind=DP), intent(inout) :: rvec(1:NDIM)    ! Vector in original frame

  real(kind=DP)                :: rtemp(1:NDIM)   ! Aux. vector

  rtemp(1) = HPsource(isource)%arot(1,1)*rvec(1) + &
       & HPsource(isource)%arot(2,1)*rvec(2) + &
       & HPsource(isource)%arot(3,1)*rvec(3)
  rtemp(2) = HPsource(isource)%arot(1,2)*rvec(1) + &
       & HPsource(isource)%arot(2,2)*rvec(2) + &
       & HPsource(isource)%arot(3,2)*rvec(3)
  rtemp(3) = HPsource(isource)%arot(1,3)*rvec(1) + &
       & HPsource(isource)%arot(2,3)*rvec(2) + &
       & HPsource(isource)%arot(3,3)*rvec(3)
  rvec(1:3) = rtemp(1:3)
  
  return
END SUBROUTINE HP_inverse_positions
