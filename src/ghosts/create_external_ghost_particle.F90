! CREATE_GHOST_PARTICLE.F90
! D. A. Hubber - 09/05/2011
! Create a new ghost particle, pghost (with id ptot + pghost in main arrays), 
! from particle p.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE create_external_ghost_particle(p,rp,ghost)
  use particle_module
  use type_module
  use mpi_communication_module
  implicit none

  integer, intent(in)       :: p               ! particle id
  real(kind=PR), intent(in) :: rp(1:NDIM)      ! ghost position
  type(ghost_sph_reference_type), intent(out) :: ghost ! Ghost reference

  ! Set up ghost particle refence
  ghost%porig       = p
  ghost%r           = rp

  return
END SUBROUTINE create_external_ghost_particle
