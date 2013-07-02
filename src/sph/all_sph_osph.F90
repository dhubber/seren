! ALL_SPH_OSPH.F90
! D. A. Hubber - 29/01/2011
! Calculates all SPH quantities (density, velocity divergence and velocity 
! curl if needed) for particle p using the gather method.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE all_sph_osph(p)
  use interface_module, only : distance2,gather_neib_on_fly
  use neighbour_module
  use hydro_module
  use kernel_module
  use particle_module
  implicit none

  integer, intent(in) :: p              ! id of current particle

  integer :: i                          ! auxilary neighbour counter 
  integer :: kern                       ! kernel table element for p
  integer :: pp                         ! neighbouring particles (p')
  integer :: pp_numb                    ! number of neighbours 
  integer :: pp_templist(1:LISTSIZE)    ! temp. list of neighbours
  real(kind=PR) :: div_v_p              ! Local copy of velocity divergence
  real(kind=PR) :: dr(1:NDIM)           ! relative position vector
  real(kind=PR) :: dv(1:VDIM)           ! Relative velocity vector
  real(kind=PR) :: dvdr                 ! Scalar product of dv and dr
  real(kind=PR) :: drmag                ! magnitude of separation
  real(kind=PR) :: drsqd                ! separation squared
  real(kind=PR) :: hfactor              ! invhp ^ NDIM
  real(kind=PR) :: hp                   ! Smoothing length of particle p
  real(kind=PR) :: invAp                ! (1 / Aent) for p
  real(kind=PR) :: invdrmag             ! (1 / drmag)
  real(kind=PR) :: invhp                ! (1 / hp)
  real(kind=PR) :: mpp                  ! mass of neighbour pp
  real(kind=PR) :: rhotemp              ! local value of density
  real(kind=PR) :: rp(1:NDIM)           ! position of particle p
  real(kind=PR) :: skern                ! 0.5 * (r/h) * KERNTOT for p
  real(kind=PR) :: vp(1:VDIM)           ! Local copy of velocity of p
#if defined(SINKS) && defined(GRAVITY)
  real(kind=PR) :: gpotmin              ! Minimum potential of neighbours 
#endif
#if defined(VISC_BALSARA)
  real(kind=PR) :: curl_v_p(1:3)        ! Curl of velocity vector (always 3-D)
  real(kind=PR) :: dvXdr(1:3)           ! Cross product of dv and dr
#endif

  debug3("Calculating SPH quantities [all_sph_osph.F90] for particle ",p)

! Calculate smoothing length here
! (Inefficient use of neighbour lists at present, but can improve in future).
#if defined(CONSTANT_H)
  sph(p)%h = hmin
#else
  call h_gather(p,sph(p)%h,sph(p)%r(1:NDIM))
#endif

! Store local copies and initialize all other quantities
  rp(1:NDIM) = sph(p)%r(1:NDIM)
  hp = sph(p)%h
  invhp = 1.0_PR / hp
  hfactor = invhp**(NDIMPR)
  rhotemp = sph(p)%m*w0(0)*hfactor
  vp(1:VDIM) = sph(p)%v(1:VDIM)
  div_v_p = 0.0_PR
  invAp = 1.0_PR/sph(p)%Aent
#if defined(SINKS)
  gpotmin = 0.0_PR
#endif
#if defined(VISC_BALSARA)
  curl_v_p(1:3) = 0.0_PR
  dvXdr(1:3) = 0.0_PR
#endif

! Copy over neighbour lists if activiated.  If not, or list is too small, 
! then obtain neighbour list on-the-fly (e.g. via treewalk).
#if defined(NEIGHBOUR_LISTS)
  pp_numb = pptot(p)
  if (pp_numb <= pp_limit) then 
     allocate(pp_templist(1:pp_numb))
     pp_templist(1:pp_numb) = pplist(1:pp_numb,p)
  else
     pp_numb = 0
  end if
#endif
  if (pp_numb == 0) call gather_neib_on_fly(p,pp_max,&
       &pp_numb,pp_templist,rp(1:NDIM),KERNRANGE*sph(p)%h)


! Now loop over all neighbours to find contributions
! ----------------------------------------------------------------------------
  do i=1,pp_numb
     pp = pp_templist(i)
     if (p == pp) cycle
     dv(1:VDIM) = sph(pp)%v(1:VDIM) - vp(1:VDIM)
     mpp = sph(pp)%m
     call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
     dvdr  = dot_product(dv(1:NDIM),dr(1:NDIM))
     drmag = sqrt(drsqd) + SMALL_NUMBER
     if (drmag >= KERNRANGE*hp) cycle
     invdrmag = 1.0_PR / drmag
     skern = HALFKERNTOT*drmag*invhp
     kern  = int(skern)
     kern  = min(kern,KERNTOT)
#if defined(ENERGY_EQN)
     rhotemp = rhotemp + mpp*hfactor*w0(kern)*sph(pp)%u/sph(p)%u
     div_v_p = div_v_p - mpp*dvdr*w4(kern)*hfactor*invhp*&
          &invdrmag*sph(pp)%u/sph(p)%u
#elif defined(ENTROPY_EQN)
     rhotemp = rhotemp + mpp*hfactor*w0(kern)*&
          &(sph(pp)%Aent*invAp)**(1.0_PR/gamma)
     div_v_p = div_v_p - mpp*dvdr*w4(kern)*hfactor*invhp*invdrmag*&
          &(sph(pp)%Aent*invAp)**(1.0_PR/gamma)
#endif
#if defined(SINKS) && defined(GRAVITY)
     gpotmin = min(-gpot(pp),gpotmin)
#endif
#if defined(VISC_BALSARA)
#if NDIM==3
     dvXdr(1) = dv(2)*dr(3) - dv(3)*dr(2)
     dvXdr(2) = dv(3)*dr(1) - dv(1)*dr(3)
#endif
#if NDIM==2 || NDIM==3
     dvXdr(3) = dv(1)*dr(2) - dv(2)*dr(1)
     curl_v_p(1:3) = curl_v_p(1:3) - mpp*dvXdr(1:3)*invdrmag* &
              & w1(kern)*hfactor*invhp
#endif
#endif
  end do
! ----------------------------------------------------------------------------


! Normalise and store SPH quantities in main arrays
  sph(p)%rho    = rhotemp
  sph(p)%invrho = 1.0_PR / rhotemp
  sph(p)%invh   = invhp
  div_v_p       = div_v_p / rhotemp
  sph(p)%div_v  = div_v_p
  sph(p)%drhodt = -rhotemp * div_v_p
#if defined(SINKS) && defined(GRAVITY)
  if (gpotmin > -sph(p)%gpot) then
     sph(p)%ispotmin = .true.
  else
     sph(p)%ispotmin = .false.
  endif
#endif
 
! Store magnitude of curl in balsara array for now, and calculate 
! balsara factor after thermal once we know sound speed.
#if defined(VISC_BALSARA)
  curl_v_p(1:3) = curl_v_p(1:3) / rhotemp
  sph(p)%balsara = sqrt(dot_product(curl_v_p(1:3),curl_v_p(1:3)))
#endif  

  return
END SUBROUTINE all_sph_osph
