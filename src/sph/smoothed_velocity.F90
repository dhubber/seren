! SMOOTHED_VELOCITY.F90
! D. A. Hubber - 20/11/2010
! Calculates smoothed velocity at position of particle p
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE smoothed_velocity(p,v_smooth_p)
  use interface_module, only : distance2,gather_neib_on_fly,w0
  use neighbour_module
  use hydro_module
  use kernel_module
  use particle_module
  implicit none

  integer, intent(in) :: p                          ! id of particle
  real(kind=PR), intent(out) :: v_smooth_p(1:VDIM)  ! Smoothed velocity of p

  integer :: i                            ! auxilary neighbour counter 
  integer :: kern                         ! kernel table element for p
  integer :: pp                           ! neighbouring particles (p')
  integer :: pp_max                       ! ..
  integer :: pp_numb                      ! number of neighbours 
  integer, allocatable :: pp_templist(:)  ! temp. list of neighbours
  real(kind=PR) :: div_v_p                ! Local copy of velocity divergence
  real(kind=PR) :: dr(1:NDIM)             ! relative position vector
  real(kind=PR) :: dv(1:VDIM)             ! Relative velocity vector
  real(kind=PR) :: dvdr                   ! Scalar product of dv and dr
  real(kind=PR) :: drmag                  ! magnitude of separation
  real(kind=PR) :: drsqd                  ! separation squared
  real(kind=PR) :: hfactor                ! invhp ^ NDIM
  real(kind=PR) :: hp                     ! Smoothing length of particle p
  real(kind=PR) :: invdrmag               ! ( 1 / drmag )
  real(kind=PR) :: invhp                  ! ( 1 / hp )
  real(kind=PR) :: mpp                    ! mass of neighbour pp
  real(kind=PR) :: rp(1:NDIM)             ! position of particle p
  real(kind=PR) :: skern                  ! 0.5 * (r/h) * KERNTOT for p
  real(kind=PR) :: vp(1:VDIM)             ! Local copy of velocity of p

  debug3("Calculating smoothed velocity [smoothed_velocity.F90] for particle ",p)

! Create local copies and initialize all other quantities
  rp(1:NDIM) = sph(p)%r(1:NDIM)
  hp = sph(p)%h
  invhp = 1.0_PR / hp
  vp(1:VDIM) = sph(p)%v(1:VDIM)
  hfactor = invhp**(NDIMPR)
  v_smooth_p(1:VDIM) = sph(p)%v(1:VDIM)*sph(p)%m*hfactor*w0(0.0_PR)
#if defined(NEIGHBOUR_LISTS)
  pp_max = max(LISTSIZE,ptot)
  pp_numb = pptot(p)
  if (pp_numb <= pp_limit) then 
     allocate(pp_templist(1:pp_max))
     pp_templist(1:pp_numb) = pplist(1:pp_numb,p)
  else
     call gather_neib_on_fly(p,pp_max,pp_numb,&
          &pp_templist,rp(1:NDIM),KERNRANGE*sph(p)%h)
  end if
#else
  call gather_neib_on_fly(p,pp_max,pp_numb,&
       &pp_templist,rp(1:NDIM),KERNRANGE*sph(p)%h)
#endif


! Now loop over all neighbours to find contributions
! ----------------------------------------------------------------------------
  do i=1,pp_numb
     pp = pp_templist(i)
     dv(1:VDIM) = sph(pp)%v(1:VDIM) - vp(1:VDIM)
     mpp = sph(pp)%m
     call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
     dvdr  = dot_product(dv(1:NDIM),dr(1:NDIM))
     drmag = sqrt(drsqd) + SMALL_NUMBER
     if (drmag >= KERNRANGE*hp) cycle
     invdrmag = 1.0_PR / drmag
     v_smooth_p(1:VDIM) = v_smooth_p(1:VDIM) + &
          &sph(pp)%v(1:VDIM)*mpp*hfactor*w0(drmag*invhp)

  end do
! ----------------------------------------------------------------------------

! Normalise and store SPH quantities in main arrays
  v_smooth_p(1:VDIM) = v_smooth_p(1:VDIM) / sph(p)%rho

  if (allocated(pp_templist)) deallocate(pp_templist)

  return
END SUBROUTINE smoothed_velocity
