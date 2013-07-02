! ALL_SPH.F90
! C. P. Batty & D. A. Hubber - 12/12/2006
! Calculates all SPH quantities (density, velocity divergence and velocity 
! curl if needed) for particle p using the gather method.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE all_sph(p,typemask)
  use interface_module, only : distance2,gather_neib_on_fly,&
       &specific_internal_energy,w0,w2,womega,wzeta
  use neighbour_module
  use hydro_module
  use kernel_module
  use particle_module
  use type_module
  implicit none

  integer, intent(in) :: p               ! id of current particle
  logical, optional, intent(in) :: typemask(1:ntypes) ! part. types to include?

  integer :: i                           ! Auxilary neighbour counter 
  integer :: pp                          ! Neighbouring particles (p')
  integer :: pp_max                      ! Max. no. of neighbours in array
  integer :: pp_numb                     ! Number of neighbours 
  integer, allocatable :: pp_templist(:) ! Temp. list of neighbours
  real(kind=PR) :: div_v_p               ! Local copy of velocity divergence
  real(kind=PR) :: dr(1:NDIM)            ! Relative position vector
  real(kind=PR) :: dv(1:VDIM)            ! Relative velocity vector
  real(kind=PR) :: dvdr                  ! Scalar product of dv and dr
  real(kind=PR) :: drmag                 ! magnitude of separation
  real(kind=PR) :: hfactor               ! invhp ^ NDIM
  real(kind=PR) :: hrangesqd             ! Kernel range
  real(kind=PR) :: invhp                 ! ( 1 / hp )
  real(kind=PR) :: mpp                   ! mass of neighbour pp
  real(kind=PR) :: rhotemp               ! local value of density
  real(kind=PR) :: rp(1:NDIM)            ! Position of particle p
  real(kind=PR) :: up                    ! Specific internal energy of p
  real(kind=PR) :: vp(1:VDIM)            ! Velocity of particle p
#if defined(SM2012_SPH)
  real(kind=PR) :: q_p                   ! Internal energy density
#endif
#if defined(GRAD_H_SPH)
  real(kind=PR) :: omega_p               ! Omega correction factor for p
#endif
#if defined(SELF_GRAVITY) && defined(GRAD_H_SPH)
  real(kind=PR) :: zeta_p                ! local value of zeta
#endif
#if defined(SINKS) && defined(SELF_GRAVITY)
  real(kind=PR) :: gpotmin               ! Minimum potential of neighbours 
#endif
#if defined(DIV_A)
  real(kind=PR) :: ap(1:NDIM)            ! accel. of particle p
  real(kind=PR) :: div_a_p               ! local value of div_a
#endif
#if defined(SMOOTHED_VELOCITY)
  real(kind=PR) :: v_smooth_p(1:VDIM)    ! smoothed velocity of particle p
#endif
#if defined(VISC_BALSARA)
  real(kind=PR) :: curl_v_p(1:3)         ! Curl of velocity vector (always 3-D)
  real(kind=PR) :: dvXdr(1:3)            ! Cross product of dv and dr
  real(kind=PR) :: mag_curl_v            ! Magnitude iof curl v
#endif
#if defined(VISC_PATTERN_REC)
  real(kind=PR) :: mv2, negcc, vneib, rneib, dotprod
  real(kind=PR) :: ratio, modr, modv, rcount
#endif

  debug3("Calculating SPH quantities [all_sph.F90] for particle ",p)

! Set constant smoothing length here
#if defined(CONSTANT_H)
  sph(p)%h = hmin
#endif

! Store local copies and initialize all other quantities
  pp_numb = 0
  rp(1:NDIM) = sph(p)%r(1:NDIM)
  invhp = 1.0_PR / sph(p)%h
  hfactor = invhp**(NDIM)
  hrangesqd = KERNRANGESQD*sph(p)%h*sph(p)%h
  rhotemp = sph(p)%m*hfactor*w0(0.0_PR)
  vp(1:VDIM) = sph(p)%v(1:VDIM)
  div_v_p = 0.0_PR
  up = sph(p)%u
#if defined(SM2012_SPH)
  q_p = sph(p)%m*up*hfactor*w0(0.0_PR)
#endif
#if defined(GRAD_H_SPH)
  omega_p = sph(p)%m*hfactor*invhp*womega(0.0_PR)
#endif
#if defined(SELF_GRAVITY) && defined(GRAD_H_SPH)
  zeta_p = sph(p)%m*invhp*invhp*wzeta(0.0_PR)
#endif
#if defined(DIV_A)
  div_a_p = 0.0_PR
#endif
#if defined(SMOOTHED_VELOCITY)
  v_smooth_p(1:VDIM) = sph(p)%v(1:VDIM)*mp*hfactor_p*w0(0.0_PR)
#endif
#if defined(VISC_BALSARA)
  curl_v_p(1:3) = 0.0_PR
#endif
#if defined(RTSPH)
#endif
#if defined(VISC_PATTERN_REC)
  modr   = sqrt(sph(p)%r(1)**2 + sph(p)%r(2)**2)
  modv   = sqrt(v(1,p)**2 + v(2,p)**2)
  mv2    = modr*modv**2
  negcc  = 0.0_PR
  rcount = 0.0_PR
#endif
#if defined(SINKS) && defined(SELF_GRAVITY)
  gpotmin = 0.0_PR
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


! Now loop over all gather neighbours
! ----------------------------------------------------------------------------
  do i=1,pp_numb
     pp = pp_templist(i)
     if (pp == p) cycle
#if defined(PERIODIC) && !defined(GHOST_PARTICLES)
     call distance2(rp(1:NDIM),pp,dr(1:NDIM),drmag)
#else
     dr(1:NDIM) = sph(pp)%r(1:NDIM) - rp(1:NDIM)
     drmag = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif
     if (drmag > hrangesqd) cycle
     drmag = sqrt(drmag) + SMALL_NUMBER
     mpp = sph(pp)%m
     dv(1:NDIM) = sph(pp)%v(1:NDIM) - vp(1:NDIM)
     dvdr = dot_product(dv(1:NDIM),dr(1:NDIM))
#if defined(SM2012_SPH)
     q_p = q_p + mpp*sph(pp)%u*hfactor*w0(drmag*invhp)
#endif
#if defined(RTSPH)
     rhotemp = rhotemp + mpp*hfactor*w0(drmag*invhp)*&
          &specific_internal_energy(pp)/up
#else
     rhotemp = rhotemp + mpp*hfactor*w0(drmag*invhp)
#endif
     div_v_p = div_v_p - mpp*dvdr*w2(drmag*invhp)*hfactor*invhp/drmag
#if defined(GRAD_H_SPH)
     omega_p = omega_p + mpp*hfactor*invhp*womega(drmag*invhp)
#endif
#if defined(DIV_A)
     div_a_p = div_a_p - mpp*w2(drmag*invhp)*hfactor*invhp*&
          & dot_product(sph(pp)%a(1:NDIM) - ap(1:NDIM),dr(1:NDIM))/drmag
#endif
#if defined(SMOOTHED_VELOCITY)
     v_smooth_p(1:VDIM) = v_smooth_p(1:VDIM) + &
          & sph(pp)%v(1:VDIM)*mpp*hfactor*w0(drmag*invhp)
#endif
#if defined(SELF_GRAVITY) && defined(GRAD_H_SPH)
     zeta_p = zeta_p + mpp*invhp*invhp*wzeta(drmag*invhp)
#endif
#if defined(VISC_BALSARA)
#if NDIM==3
     dvXdr(1) = dv(2)*dr(3) - dv(3)*dr(2)
     dvXdr(2) = dv(3)*dr(1) - dv(1)*dr(3)
#endif
#if NDIM==2 || NDIM==3
     dvXdr(3) = dv(1)*dr(2) - dv(2)*dr(1)
#endif
     curl_v_p(1:3) = curl_v_p(1:3) - &
          &(mpp * (dvXdr(1:3) / drmag) * (w2(drmag*invhp)*hfactor*invhp))
#endif
#if defined(VISC_PATTERN_REC)
     rcount = rcount + 1.0_PR
     vneib  = sqrt(v(1,pp)**2 + v(2,pp)**2)
     rneib  = sqrt(sph(pp)%r(1)**2 + sph(pp)%r(2)**2)
     
     ! Two checks. If either fail, neighbour is anomalous
     dotprod = (v(1,pp)*sph(pp)%r(1) + v(2,pp)*sph(pp)%r(2))/vneib/rneib
     ratio   = mv2/vneib/vneib/rneib    
     if (dotprod > 0.03_PR .or. dotprod < -0.03_PR .or. &
          &ratio > 1.03_PR .or. ratio < 0.97_PR) negcc = negcc + 1
#endif
#if defined(SINKS) && defined(SELF_GRAVITY)
     gpotmin = min(-sph(pp)%gpot,gpotmin)
#endif
  end do
! ----------------------------------------------------------------------------


! Store various SPH quantities in main arrays.  For the Balsara switch, 
! store magnitude of curl in balsara array for now, and calculate 
! balsara factor after thermal once we know sound speed
! ----------------------------------------------------------------------------
  sph(p)%rho    = rhotemp
#if !defined(LOW_MEM)
  sph(p)%invrho = 1.0_PR / rhotemp
  sph(p)%invh   = invhp
  sph(p)%div_v  = div_v_p*sph(p)%invrho
#else
  sph(p)%div_v  = div_v_p/rhotemp
#endif
#if defined(SM2012_SPH)
  sph(p)%q      = q_p
#endif
#if defined(GRAD_H_SPH)
  sph(p)%omega  = 1.0_PR + (sph(p)%h*omega_p) / (NDIMPR*rhotemp)
#endif
#if defined(SELF_GRAVITY) && defined(GRAD_H_SPH)
  sph(p)%zo    = -1.0_PR*(sph(p)%h*zeta_p) / (NDIMPR*sph(p)%rho*sph(p)%omega)
#endif
#if defined(DIV_A)
  sph(p)%div_a = div_a_p*sph(p)%invrho
#endif
#if defined(SMOOTHED_VELOCITY)
  sph(p)%v_smooth(1:VDIM) = v_smooth_p(1:VDIM)*sph(p)%invrho
#endif
#if defined(VISC_BALSARA)
  curl_v_p(1:3) = curl_v_p(1:3)*sph(p)%invrho
  mag_curl_v = sqrt(dot_product(curl_v_p(1:3),curl_v_p(1:3)))
  sph(p)%balsara = mag_curl_v
#endif
#if defined(VISC_PATTERN_REC)
  if (rcount > 0) sph(p)%pattrec = negcc/rcount
#endif
#if defined(SINKS) && defined(SELF_GRAVITY)
  if (gpotmin > -sph(p)%gpot) then
     sph(p)%ispotmin = .true.
  else
     sph(p)%ispotmin = .false.
  endif
#endif

  if (allocated(pp_templist)) deallocate(pp_templist)

  return
END SUBROUTINE all_sph
