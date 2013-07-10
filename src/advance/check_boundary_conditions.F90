! CHECK_BOUNDARY_CONDITIONS.F90
! C. P. Batty & D. A. Hubber - 4/2/2009
! Check if particles have crossed any periodic, mirror or wall boundaries.  
! If so, correct the particle's position and velocity.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE check_boundary_conditions(rp,vp)
  use definitions
  use periodic_module
#if defined(USE_MPI)
  use mpi_communication_module
#endif
  implicit none

  real(kind=PR), intent(inout) :: rp(1:NDIM)  ! Position of particle/sink
  real(kind=PR), intent(inout) :: vp(1:VDIM)  ! Velocity of particle/sink

#if defined(USE_MPI)
  real(kind=PR) :: wrap_min(1:3)              ! Point at which to wrap particles
  real(kind=PR) :: wrap_max(1:3)              ! Point at which to wrap particles
#endif


! Reposition particle for periodic boundary conditions
! ----------------------------------------------------------------------------
#if defined(USE_MPI)
  wrap_min = periodic_mpi_min
  wrap_max = periodic_mpi_max
#if defined(PERIODIC_X)
  wrap_min(1) = max(domain_bbmax(1,rank) - periodic_size(1), wrap_min(1))
  wrap_max(1) = min(domain_bbmin(1,rank) + periodic_size(1), wrap_max(1))
  wrap_min(1) = min(periodic_min(1), wrap_min(1))
  wrap_max(1) = max(periodic_max(1), wrap_max(1))
  if (rp(1) < wrap_min(1)) rp(1) = rp(1) + periodic_size(1)
  if (rp(1) > wrap_max(1)) rp(1) = rp(1) - periodic_size(1)
#endif
#if defined(PERIODIC_Y) && (NDIM == 2 || NDIM == 3)
  wrap_min(2) = max(domain_bbmax(2,rank) - periodic_size(2), wrap_min(2))
  wrap_max(2) = min(domain_bbmin(2,rank) + periodic_size(2), wrap_max(2))
  wrap_min(2) = min(periodic_min(2), wrap_min(2))
  wrap_max(2) = max(periodic_max(2), wrap_max(2))
  if (rp(2) < wrap_min(2)) rp(2) = rp(2) + periodic_size(2)
  if (rp(2) > wrap_max(2)) rp(2) = rp(2) - periodic_size(2)
#endif
#if defined(PERIODIC_Z) && (NDIM == 3)
  wrap_min(3) = max(domain_bbmax(3,rank) - periodic_size(3), wrap_min(3))
  wrap_max(3) = min(domain_bbmin(3,rank) + periodic_size(3), wrap_max(3))
  wrap_min(3) = min(periodic_min(3), wrap_min(3))
  wrap_max(3) = max(periodic_max(3), wrap_max(3))
  if (rp(3) < wrap_min(3)) rp(3) = rp(3) + periodic_size(3)
  if (rp(3) > wrap_max(3)) rp(3) = rp(3) - periodic_size(3)
#endif

#else

#if defined(PERIODIC_X)
  if (rp(1) < periodic_min(1)) rp(1) = rp(1) + periodic_size(1)
  if (rp(1) > periodic_max(1)) rp(1) = rp(1) - periodic_size(1)
#endif
#if defined(PERIODIC_Y)
#if NDIM == 2 || NDIM == 3
  if (rp(2) < periodic_min(2)) rp(2) = rp(2) + periodic_size(2)
  if (rp(2) > periodic_max(2)) rp(2) = rp(2) - periodic_size(2)
#endif
#endif
#if defined(PERIODIC_Z)
#if NDIM == 3
  if (rp(3) < periodic_min(3)) rp(3) = rp(3) + periodic_size(3)
  if (rp(3) > periodic_max(3)) rp(3) = rp(3) - periodic_size(3)
#endif
#endif

#endif
! ----------------------------------------------------------------------------


! Reposition particle for mirror or wall boundary conditions
! ----------------------------------------------------------------------------
#if defined(MIRROR_X_LHS) || defined(WALL_X_LHS)
  if (rp(1) < periodic_min(1)) then
     rp(1) = periodic_min(1) + (periodic_min(1) - rp(1))
     vp(1) = -vp(1)
  end if
#endif
#if defined(MIRROR_X_RHS) || defined(WALL_X_RHS)
  if (rp(1) > periodic_max(1)) then
     rp(1) = periodic_max(1) - (rp(1) - periodic_max(1));
     vp(1) = -vp(1)
  end if
#endif
#if NDIM==2 || NDIM==3
#if defined(MIRROR_Y_LHS) || defined(WALL_Y_LHS)
  if (rp(2) < periodic_min(2)) then
     rp(2) = periodic_min(2) + (periodic_min(2) - rp(2))
     vp(2) = -vp(2)
  end if
#endif
#if defined(MIRROR_Y_RHS) || defined(WALL_Y_RHS)
  if (rp(2) > periodic_max(2)) then
     rp(2) = periodic_max(2) - (rp(2) - periodic_max(2))
     vp(2) = -vp(2)
  end if
#endif
#endif
#if NDIM==3
#if defined(MIRROR_Z_LHS) || defined(WALL_Z_LHS)
  if (rp(3) < periodic_min(3)) then
     rp(3) = periodic_min(3) + (periodic_min(3) - rp(3))
     vp(3) = -vp(3)
  end if
#endif
#if defined(MIRROR_Z_RHS) || defined(WALL_Z_RHS)
  if (rp(3) > periodic_max(3)) then
     rp(3) = periodic_max(3) - (rp(3) - periodic_max(3))
     vp(3) = -vp(3)
  end if
#endif
#endif
! ----------------------------------------------------------------------------


  return
END SUBROUTINE check_boundary_conditions
