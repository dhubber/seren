! MEMORY_POOL.F90
! D. A. Hubber - ??
! ..
! ============================================================================

#include "macros.h"

! ============================================================================
MODULE memory_pool
  use definitions

  integer :: l1_size
  integer :: l1_free
  integer :: l2_size
  integer :: l2_free
  integer :: i1_size
  integer :: i1_free
  integer :: i2_size
  integer :: i2_free
  integer :: r1_size
  integer :: r1_free
  integer :: r2_size
  integer :: r2_free
  integer :: d1_size
  integer :: d1_free
  integer :: d2_size
  integer :: d2_free


  ! Generic arrays
  ! --------------------------------------------------------------------------
  logical, allocatable :: logical1(:)
  logical, allocatable :: logical2(:,:)
  integer, allocatable :: integer1(:)
  integer, allocatable :: integer2(:,:)
  real(kind=PR), allocatable :: real1(:)
  real(kind=PR), allocatable :: real2(:,:)
  real(kind=DP), allocatable :: double1(:)
  real(kind=DP), allocatable :: double2(:,:)


  ! SEREN-specific arrays
  ! --------------------------------------------------------------------------
  !type(BHhydro_node), allocatable :: BHhydro1(:)
  !type(BHgrav_node), allocatable :: BHgrav1(:)
  !type(sink_node), allocatable :: sink(:)
  !type(star_node), allocatable :: star(:)


contains


END MODULE memory_pool
