! EIGENVALUE_MAC.F90
! D. A. Hubber - 18/10/2008
! Calculates maximum eigenvalue of quadrupole moment tensor for 
! eigenvalue MAC.  Since the Q matrix is traceless, the characteristic 
! equation is a depressed cubic of the form y^3 + py + q = 0.  Using a 
! solution developed by Vieta (e.g. Geometric Construction, Martin 1998), 
! the maximum magnitude of the three solutions is given by y = sqrt(-4p/3).
! We use this directly to calculate the maximum eigenvalue for the MAC.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE eigenvalue_mac(qc,mac)
  use definitions
  use tree_module, only : abserror
  implicit none

  real(kind=DP), intent(in) :: qc(1:NQUAD)  ! Quadrupole moment tensor of c
  real(kind=PR), intent(out) :: mac         ! MAC value

  real(kind=DP) :: lambda                   ! max eigenvalue
  real(kind=DP) :: p                        ! aux. quadrupole variable

! Compute q-matrix variable depending on dimensionality
#if NDIM==2
  p = qc(1)*qc(3) - (qc(1) + qc(3))*(qc(1) + qc(3)) - qc(2)*qc(2)
#elif NDIM==3
  p = qc(1)*qc(3) - (qc(1) + qc(3))*(qc(1) + qc(3)) &
       & - qc(2)*qc(2) - qc(4)*qc(4) - qc(5)*qc(5)
#endif

! If p is zero (or more due to rounding errors), set MAC to zero.  
! Otherwise, calculate maximum eigenvalue and return
  if (p >= 0.0_DP) then
     mac = 0.0_PR
  else
     lambda = 2.0_DP*sqrt(-p/3.0_DP)
     mac = real((0.5_DP*lambda/abserror)**(TWOTHIRDS_DP),PR)
  end if

  return
END SUBROUTINE eigenvalue_mac
