! BHGRAV_NODE_ACCEL.F90
! A. McLeod & D. A. Hubber - 23/01/2008
! Computes gravitational force exerted on at position rp due to the tree 
! cell 'node'.  Also computes higher order multipole terms if activated.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHgrav_node_accel(node,rp,atemp,dpot)
  use interface_module, only : distance3,ewald_force,&
       &gravity_gradh,gravity_nbody,gravity_sph,gravity_meanh,wpot
  use definitions
  use tree_module
  use particle_module
  implicit none

  type(BHgrav_node) :: node                     ! Grav. node
  real(kind=PR), intent(in) :: rp(1:NDIM)       ! Position of particle p
  real(kind=PR), intent(out) :: atemp(1:NDIM)   ! Grav. accelertation
  real(kind=PR), intent(out) :: dpot            ! Grav. potential

  real(kind=PR) :: dr(1:NDIM)       ! Relative displacement vector
  real(kind=PR) :: drsqd            ! Distance squared
  real(kind=PR) :: invdrmag         ! ( 1 / drmag )
  real(kind=PR) :: invdrsqd         ! ( 1 / drsqd )
#if defined(QUADRUPOLE) || defined(OCTUPOLE)
  real(kind=PR) :: invdr5           ! ( 1 / drmag^5 )
#endif
#if defined(QUADRUPOLE)
  real(kind=PR) :: qscalar          ! Inner product of quad tensor
#endif
#if defined(OCTUPOLE)
  real(kind=PR) :: sscalar          ! Scalar component of octupole terms
#endif
#if defined(EWALD)
  real(kind=PR) :: eaccel(1:NDIM)   ! Ewald grav. acceleration
#endif

  ! Calculate distance depending on whether we use Ewald gravity
#if defined(EWALD)
  call distance3(node%r(1:NDIM),rp(1:NDIM),dr(1:NDIM),drsqd)
#else
  dr(1:NDIM) = rp(1:NDIM) - node%r(1:NDIM)
  drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif
  
  invdrsqd = 1.0_PR / (drsqd + SMALL_NUMBER)
  invdrmag = sqrt(invdrsqd)
  atemp(1:NDIM) = -node%m*invdrsqd*invdrmag*dr(1:NDIM)
  dpot = node%m*invdrmag

  ! Add quadrupole moment correction terms
#if defined(QUADRUPOLE)
  invdr5 = invdrsqd*invdrsqd*invdrmag
#if NDIM==2
  qscalar = node%q(1)*dr(1)*dr(1) + node%q(3)*dr(2)*dr(2) &
       & + 2.0_PR*node%q(2)*dr(1)*dr(2)
  atemp(1) = atemp(1) + (node%q(1)*dr(1) + node%q(2)*dr(2))&
       &*invdr5 - 2.5_PR*qscalar*dr(1)*invdr5*invdrsqd
  atemp(2) = atemp(2) + (node%q(2)*dr(1) + node%q(3)*dr(2))&
       &*invdr5 - 2.5_PR*qscalar*dr(2)*invdr5*invdrsqd
#elif NDIM==3
  qscalar = node%q(1)*dr(1)*dr(1) + node%q(3)*dr(2)*dr(2) &
       & - (node%q(1) + node%q(3))*dr(3)*dr(3) &
       & + 2.0_PR*(node%q(2)*dr(1)*dr(2) &
       & + node%q(4)*dr(1)*dr(3) + node%q(5)*dr(2)*dr(3))
  atemp(1) = atemp(1) + (node%q(1)*dr(1) + node%q(2)*dr(2)&
       & + node%q(4)*dr(3))*invdr5 &
       & - 2.5_PR*qscalar*dr(1)*invdr5*invdrsqd
  atemp(2) = atemp(2) + (node%q(2)*dr(1)+node%q(3)*dr(2)&
       & + node%q(5)*dr(3))*invdr5 &
       & - 2.5_PR*qscalar*dr(2)*invdr5*invdrsqd
  atemp(3) = atemp(3) + (node%q(4)*dr(1) + node%q(5)*dr(2)&
       & - (node%q(1)+node%q(3))*dr(3))*invdr5 &
       & - 2.5_PR*qscalar*dr(3)*invdr5*invdrsqd
#endif
  dpot = dpot + 0.5_PR*qscalar*invdr5 
#endif
  
  ! Add octupole moment correction terms
#if defined(OCTUPOLE)
  sscalar = node%s(1)*dr(1)*dr(1)*dr(1) + &
       & node%s(2)*dr(1)*dr(1)*dr(2) + &
       & node%s(3)*dr(2)*dr(2)*dr(1) + &
       & node%s(4)*dr(2)*dr(2)*dr(2) + &
       & node%s(5)*dr(3)*dr(3)*dr(1) + &
       & node%s(6)*dr(3)*dr(3)*dr(2) + &
       & node%s(7)*dr(3)*dr(3)*dr(3) + &
       & node%s(8)*dr(1)*dr(1)*dr(3) + &
       & node%s(9)*dr(2)*dr(2)*dr(3) + &
       & node%s(10)*dr(1)*dr(2)*dr(3)
  atemp(1) = atemp(1) + 0.5_PR*(3.0_PR*node%s(1)*dr(1)*dr(1) + &
       & 2.0_PR*node%s(2)*dr(1)*dr(2) + &
       & 2.0_PR*node%s(8)*dr(1)*dr(3) + &
       & node%s(3)*dr(2)*dr(2) + node%s(5)*dr(3)*dr(3) + &
       & node%s(10)*dr(2)*dr(3))*invdr5*invdrsqd - &
       & 3.5_PR*sscalar*dr(1)*invdr5*invdrsqd*invdrsqd
  atemp(2) = atemp(2) + 0.5_PR*(3.*node%s(4)*dr(2)*dr(2) + &
       & 2.0_PR*node%s(3)*dr(1)*dr(2) + &
       & 2.0_PR*node%s(9)*dr(2)*dr(3) + &
       & node%s(2)*dr(1)*dr(1) + node%s(6)*dr(3)*dr(3) + &
       & node%s(10)*dr(1)*dr(3))*invdr5*invdrsqd - &
       & 3.5_PR*sscalar*dr(2)*invdr5*invdrsqd*invdrsqd
  atemp(3) = atemp(3) + 0.5_PR*(3.*node%s(7)*dr(3)*dr(3) + &
       & 2.0_PR*node%s(5)*dr(1)*dr(3) + &
       & 2.0_PR*node%s(6)*dr(2)*dr(3) + &
       & node%s(8)*dr(1)*dr(1) + node%s(9)*dr(2)*dr(2) + &
       & node%s(10)*dr(1)*dr(2))*invdr5*invdrsqd - &
       & 3.5_PR*sscalar*dr(3)*invdr5*invdrsqd*invdrsqd
  dpot = dpot + 0.5_PR*sscalar*invdr5*invdrsqd
#endif
  
  ! Add Ewald correction force if required
#if defined(EWALD)
  call ewald_force(dr,node%m,eaccel)
  atemp(1:NDIM) = atemp(1:NDIM) + eaccel(1:NDIM)
#endif
  
  return
END SUBROUTINE BHgrav_node_accel
