! HP_CALCULATE_BASIS_VECTOR
! T. Bisbas & D. A. Hubber - 09/11/2010
! Calculate basis vector and then transform to rotated frame
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE HP_calculate_basis_vector(level,ipix,rbasis,isource)
  use interface_module, only : pix2vec_nest
  use definitions
  use healpix_types
  use HP_module, only : pix2x,pix2y,HPsource
  implicit none

  integer, intent(in) :: isource                ! id of ionizing source
  integer, intent(in) :: level                  ! Current level
  integer(kind=I4B), intent(in) :: ipix         ! Pixel id
  real(kind=PR), intent(out) :: rbasis(1:NDIM)  ! Unit vector of pixel ipix

  integer(kind=I4B) :: nside                    ! ????
  real(kind=DP) :: vector(1:3)                  ! Unit vector

! Calculate basis vector in original HEALPix frame
  nside = 2**(int(level,I4B))
  call pix2vec_nest(nside,ipix,pix2x,pix2y,vector)

! Rotate basis vector to new frame of reference using rotation matrix
  rbasis(1) = real(HPsource(isource)%arot(1,1)*vector(1) + &
       & HPsource(isource)%arot(1,2)*vector(2) + &
       & HPsource(isource)%arot(1,3)*vector(3),PR)
  rbasis(2) = real(HPsource(isource)%arot(2,1)*vector(1) + &
       & HPsource(isource)%arot(2,2)*vector(2) + &
       & HPsource(isource)%arot(2,3)*vector(3),PR)
  rbasis(3) = real(HPsource(isource)%arot(3,1)*vector(1) + &
       & HPsource(isource)%arot(3,2)*vector(2) + &
       & HPsource(isource)%arot(3,3)*vector(3),PR)
  
  return
END SUBROUTINE HP_calculate_basis_vector


