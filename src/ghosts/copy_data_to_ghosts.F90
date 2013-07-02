! COPY_DATA_TO_GHOSTS.F90
! D. A. Hubber - 10/05/2011
! Copy all extra data (except the position and velocity) to all ghost 
! particles from their original 'parent' particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE copy_data_to_ghosts
  use particle_module
  use periodic_module
  use type_module
  implicit none

  integer :: i                       ! Aux. particle counter
  integer :: p                       ! Ghost particle id
  integer :: porig                   ! Id of original particle
  real(kind=PR) :: rp(1:NDIM)        ! Position of ghost particle
  real(kind=PR) :: vp(1:NDIM)        ! Velocity of ghost particle

  debug2("Copy all particle data to ghosts from originals [copy_data_to_ghosts.F90]")

  ! Loop over all new (periodic) ghost particles
  ! --------------------------------------------------------------------------
  do i=1,pperiodic
     p = ptot + i
     porig = sph(p)%porig

     ! Create temporary copies of position and velocity of ghost
     rp(1:NDIM) = sph(p)%r(1:NDIM)
     vp(1:VDIM) = sph(p)%v(1:NDIM)

     ! Now copy all data from original particle to ghost
     sph(p) = sph(porig)

     ! Copy back position and velocity info to ghost
     sph(p)%r(1:NDIM) = rp(1:NDIM)
     sph(p)%v(1:VDIM) = vp(1:NDIM)
     sph(p)%porig = porig
     sph(p)%accdo = .false.

  end do
  ! --------------------------------------------------------------------------

  return
END SUBROUTINE copy_data_to_ghosts
