! CREATE_GHOST_PARTICLE.F90
! D. A. Hubber - 09/05/2011
! Create a new ghost particle, pghost (with id ptot + pghost in main arrays), 
! from particle p.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE create_ghost_particle(p,k,rk,vk,mp)
  use particle_module
  use type_module
  implicit none

  integer, intent(in) :: p               ! particle id
  integer, intent(in) :: k               ! ghost dimension (??)
  real(kind=PR), intent(in) :: rk        ! k-position of ghost
  real(kind=PR), intent(in) :: vk        ! k-velocity of ghost
  real(kind=PR), intent(in) :: mp        ! ghost particle mass

  ! Increment ghost particle counter.  If counter has exceeded maximum 
  ! allowed number, stop adding ghosts and flag up.
  pghost = pghost + 1
  if (pghost > pghostmax) return

  ! Create new ghost and copy all information from original particle, p, 
  ! except the position and velocity which are different.  Also, set accdo 
  ! to false since we will copy, not calculate, ghost properties.
  sph(ptot + pghost)       = sph(p)
  sph(ptot + pghost)%r(k)  = rk
  sph(ptot + pghost)%v(k)  = vk
  sph(ptot + pghost)%m     = mp
  sph(ptot + pghost)%accdo = .false.

  ! Record id of original particle in porig array for later copying.
  ! (N.B. using current id, not actual original id, so a bit mis-leading, 
  ! but easier to work with).
  if (p > ptot) then
     sph(ptot + pghost)%porig = sph(p)%porig
  else
     sph(ptot + pghost)%porig = p
  end if

  return
END SUBROUTINE create_ghost_particle
