! TURB_FORCE_APPLY.F90
! A.McLeod - 26/11/2012
! Apply turbulent forcing to a particle. Find particle position in grid, and
! use trilinear interpolation in space and linear interpolation in time
! to determine forcing 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE turb_force_apply(p, aturb)
   use definitions
   use time_module
   use type_module
   use particle_module
   use turbulence_module
   implicit none

   integer, intent(in)        :: p               ! id of current particle
   real(kind=PR), intent(out) :: aturb(1:NDIM)   ! Turbulent forcing
   
   real(kind=PR)              :: r(1:NDIM)       ! Position in grid
   real(kind=PR)              :: dr1(1:NDIM)     ! Position in cell
   real(kind=PR)              :: dr2(1:NDIM)     ! 1 - position in cell
   integer                    :: bmin(1:NDIM)    ! 'top-left' corner of box
   integer                    :: bmax(1:NDIM)    ! 'bottom-right' corner of box
   real(kind=PR)              :: cube_vals(1:NDIM,1:8) 
                                                 ! Values of turbulent forcing
                                                 ! from top left to bottom right
   real(kind=PR)              :: interp(1:NDIM,1:8) ! Interpolation weights
   
   ! Find position of particle in 'grid' space
   r = sph(p)%r
#if defined(USE_MPI) && defined(PERIODIC)
   call unwrap_particle_position(r)
#endif
   r = (r - turb_min) / turb_space
   
   ! If particles are outside turb_box, no forcing occurs.
   if (any(r < 0.0_PR)) then
      aturb = 0.0_PR
      return
   else if (any(r > turb_gs_real)) then
      aturb = 0.0_PR
      return
   end if
   
   ! Find ids of top left and bottom right corner of encompassing cube
   ! These can be the same if a particle is right on the boundary, but won't
   ! affect the result.
   bmin = floor(r)
   bmax = ceiling(r)
   
   ! Right-hand edge is equal to left-hand edge (periodic),
   ! which gives TURB_GS + 1 interpolation points
   where (bmin==TURB_GS) bmin=0 ! this only happens for r==TURB_GS
   where (bmax==TURB_GS) bmax=0
   
   dr1 = r - real(floor(r),PR)
   dr2 = 1.0_PR - dr1
   
   ! Find cube values
   cube_vals(:,1) = afield_now(:, bmin(1), bmin(2), bmin(3))
   cube_vals(:,2) = afield_now(:, bmax(1), bmin(2), bmin(3))
   cube_vals(:,3) = afield_now(:, bmin(1), bmax(2), bmin(3))
   cube_vals(:,4) = afield_now(:, bmax(1), bmax(2), bmin(3))
   cube_vals(:,5) = afield_now(:, bmin(1), bmin(2), bmax(3))
   cube_vals(:,6) = afield_now(:, bmax(1), bmin(2), bmax(3))
   cube_vals(:,7) = afield_now(:, bmin(1), bmax(2), bmax(3))
   cube_vals(:,8) = afield_now(:, bmax(1), bmax(2), bmax(3))

   ! Find interpolation values
   interp(:,1) = dr2(1) * dr2(2) * dr2(3)
   interp(:,2) = dr1(1) * dr2(2) * dr2(3)
   interp(:,3) = dr2(1) * dr1(2) * dr2(3)
   interp(:,4) = dr1(1) * dr1(2) * dr2(3)
   interp(:,5) = dr2(1) * dr2(2) * dr1(3)
   interp(:,6) = dr1(1) * dr2(2) * dr1(3)
   interp(:,7) = dr2(1) * dr1(2) * dr1(3)
   interp(:,8) = dr1(1) * dr1(2) * dr1(3)

   ! Find interpolated value
   aturb = sum(interp * cube_vals, dim=2)

   return
  
END SUBROUTINE turb_force_apply
