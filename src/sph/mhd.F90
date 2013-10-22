! MHD.F90
! A. McLeod - 21/10/2013
! Calculates MHD stuff for particle p
! ============================================================================

#include "macros.h"
#define BETA_RESIST 1.0
#define BETA_TENSILE 1.0
#define CLEAN_FAST_WAVE_SPEED 1.0
#if NDIM==1
#define SIGMA_DAMP 'no idea what this should be'
#elif NDIM==2
#define SIGMA_DAMP 0.25
#else
#define SIGMA_DAMP 1.0
#endif

! ============================================================================
SUBROUTINE mhd(p)
  use interface_module, only : distance2,get_neib_on_fly,w1
  use particle_module
  use neighbour_module
  use hydro_module
  use kernel_module
!   use time_module
!   use type_module
  implicit none

  integer, intent(in) :: p               ! particle identifier

  integer :: i                           ! neighbour counter
  integer :: pp                          ! neighbouring particle number
  integer :: pp_numb                     ! number of neighbours
  integer, allocatable :: pp_templist(:) ! temp. list of neighbours
  real(kind=PR) :: alpha_p               ! resistivity alpha for particle p
  real(kind=PR) :: alpha_pp              ! resistivity alpha for particle pp
  real(kind=PR) :: amag_temp(1:VDIM)     ! mag acceleration of particle p
  real(kind=PR) :: B_p(1:BDIM)           ! magnetic field of particle p
  real(kind=PR) :: B_pp(1:BDIM)          ! magnetic field of particle pp
  real(kind=PR) :: B_t_signal            ! Signal time between particles
  real(kind=PR) :: Bdr                   ! B dot drunit
  real(kind=PR) :: c_h                   ! Hyperbolic cleaning speed
  real(kind=PR) :: dBdt_temp(1:BDIM)     ! dBdt for particle p
  real(kind=PR) :: dBdt_diss(1:BDIM)     ! dBdt resistivity for particle p
  real(kind=PR) :: dBdt_damp(1:BDIM)     ! dBdt hyperbolic damping, p
  real(kind=PR) :: dphi_dt_temp          ! scalar damping field for particle p
  real(kind=PR) :: dphi_dt_econs         ! phi energy cons. term for particle p
  real(kind=PR) :: drmag                 ! magnitude of separation
  real(kind=PR) :: dr_unit(1:BDIM)       ! Unit displacement vector
  real(kind=PR) :: dv(1:VDIM)            ! Relative velocity
  real(kind=PR) :: dvdr                  ! Scalar product of dv and dr
  real(kind=PR) :: hmean                 ! (h(p) + h(pp)) / 2
  real(kind=PR) :: hfactor_p             ! invhp ^ (NDIM + 1)
  real(kind=PR) :: hfactor_pp            ! invhpp ^ (NDIM + 1)
  real(kind=PR) :: hp                    ! smoothing length of particle p
  real(kind=PR) :: hpp                   ! smoothing length of neighbour pp
  real(kind=PR) :: invhp                 ! ( 1 / hp )
  real(kind=PR) :: invhpp                ! ( 1 / hpp )
  real(kind=PR) :: invdrmag              ! ( 1 / drmag )
  real(kind=PR) :: factors               ! all scalar factors (temporary)
  real(kind=PR) :: maxwell_p(1:BDIM, 1:BDIM) ! Maxwell stress tensor, p
  real(kind=PR) :: maxwell_pp(1:BDIM, 1:BDIM) ! Maxwell stress tensor, pp
  real(kind=PR) :: mpp                   ! mass of neighbour pp
  real(kind=PR) :: pfactor_p             ! 1.0 / (rho**2 * omega)
  real(kind=PR) :: pfactor_pp            ! 1.0 / (rho**2 * omega) for pp
  real(kind=PR) :: phi_p                 ! scalar damping field for p
  real(kind=PR) :: phi_pp                ! scalar damping field for pp
  real(kind=PR) :: rho_p                 ! density of particle p
  real(kind=PR) :: rho_pp                ! density of particle pp
  real(kind=PR) :: rhomean               ! mean density of particles p and pp
  real(kind=PR) :: rp(1:NDIM)            ! position of particle p
  real(kind=PR) :: sound_p               ! sound speed of particle p
  real(kind=PR) :: sound_pp              ! sound speed of particle pp
  real(kind=PR) :: tensile_temp          ! tensile instability correction, p
  real(kind=PR) :: vp(1:VDIM)            ! velocity of particle p
  real(kind=PR) :: vsig_1_p              ! first part of signal velocity, p
  real(kind=PR) :: vsig_2_p              ! second part of signal velocity, p
  real(kind=PR) :: vsig_1_pp             ! first part of signal velocity, pp
  real(kind=PR) :: vsig_2_pp             ! second part of signal velocity, pp
  real(kind=PR) :: vsig_p                ! signal velocity, p
  real(kind=PR) :: vsig_pp               ! signal velocity, pp
  real(kind=PR) :: vsig                  ! signal velocity, combined (timestep)
  real(kind=PR) :: vsig_resist           ! signal velocity, combined (resist.)
#if !defined(GRAD_H_SPH)
  real(kind=PR) :: wmean                 ! (W(p) + W(pp)) / 2
#endif

  debug3("Calculating mhd variables and forces [mhd.F90] for particle ", p)

! Create local copies of important properties of particle p
  pp_numb    = 0
  rp(1:NDIM) = sph(p)%r(1:NDIM)
  hp         = sph(p)%h
  invhp      = 1.0_PR / hp
  hfactor_p  = invhp**(NDIMPLUS1)
  vp(1:NDIM) = sph(p)%v(1:NDIM)
  sound_p    = sph(p)%sound
  rho_p      = sph(p)%rho
  B_p        = sph(p)%B
  phi_p      = sph(p)%phi
  alpha_p    = sph(p)%alpha_resist
  maxwell_p  = maxwell_stress_tensor(B_p)
#if defined(GRAD_H_SPH)
  pfactor_p  = 1.0_PR / (rho_p**2 * sph(p)%omega)
#else
  pfactor_p  = 1.0_PR / rho_p**2
#endif
  
  vsig_1_p   = sound_p**2 + sum(B_p**2)/(rho_p*mu_0)
  vsig_2_p   = 2.0_PR / sqrt(rho_p * mu_0)

! Zero arrays
  dr_unit = 0.0_PR
  amag_temp = 0.0_PR
  tensile_temp = 0.0_PR
  dBdt_temp = 0.0_PR
  dBdt_diss = 0.0_PR
  dBdt_damp = 0.0_PR
  dphi_dt_temp = 0.0_PR
  dphi_dt_econs = 0.0_PR
  B_t_signal = BIG_NUMBER

! Copy neighbour lists if recorded in arrays, or recompute the neighbour list
! by walking the tree
#if defined(NEIGHBOUR_LISTS)
  pp_numb = pptot(p)
  if (pp_numb <= pp_limit) then 
     allocate(pp_templist(1:pp_numb))
     pp_templist(1:pp_numb) = pplist(1:pp_numb,p)
  else
     pp_numb = 0
  end if
#endif
  if (pp_numb == 0) call get_neib_on_fly(p,pp_numb,&
       &ptot,pp_templist,rp(1:NDIM),KERNRANGE*hp)


! Loop over all neighbours, summing all contributions
! ============================================================================
  do i=1,pp_numb
     pp = pp_templist(i)
     if (p == pp) cycle

     ! Create local copies for neighbour pp
!      wmean    = 0.0_PR
     rho_pp      = sph(pp)%rho
     mpp         = sph(pp)%m
     hpp         = sph(pp)%h
     sound_pp    = sph(pp)%sound
     B_pp        = sph(pp)%B
     phi_pp      = sph(pp)%phi
     alpha_pp    = sph(pp)%alpha_resist
     maxwell_pp  = maxwell_stress_tensor(B_pp)
     vsig_1_pp   = sound_pp**2 + sum(B_pp**2)/(rho_pp*mu_0)
     vsig_2_pp   = 2.0_PR / sqrt(rho_pp * mu_0)
     rhomean     = 0.5 * (rho_p + rho_pp)
     
#if defined(PERIODIC) && !defined(GHOST_PARTICLES)
     call distance2(rp(1:NDIM),pp,dr_unit(1:NDIM),drmag)
#else
     dr_unit(1:NDIM) = sph(pp)%r(1:NDIM) - rp(1:NDIM)
     drmag = dot_product(dr_unit(1:NDIM),dr_unit(1:NDIM))
#endif
     if (drmag >= KERNRANGESQD*hp*hp .and. drmag >= KERNRANGESQD*hpp*hpp) cycle
     drmag = sqrt(drmag) + SMALL_NUMBER
     dv(1:NDIM) = sph(pp)%v(1:NDIM) - vp(1:NDIM)
     dvdr = dot_product(dv,dr_unit)
     invdrmag = 1.0_PR / drmag
     dr_unit(1:NDIM) = dr_unit(1:NDIM)*invdrmag
     
     ! Signal speed/time
     hmean = 0.5 * (hp + hpp)
     Bdr = dot_product(B_p, dr_unit)
     vsig_p = 0.5_PR * (sqrt(vsig_1_p + vsig_2_p*Bdr) +&
                       &sqrt(vsig_1_p - vsig_2_p*Bdr))
     Bdr = dot_product(B_pp, dr_unit)
     vsig_pp = 0.5_PR * (sqrt(vsig_1_pp + vsig_2_pp*Bdr) +&
                        &sqrt(vsig_1_pp - vsig_2_pp*Bdr))
     vsig_resist = vsig_p + vsig_pp
     vsig = vsig_resist + BETA_RESIST * abs(dvdr)
     B_t_signal = min(B_t_signal, hmean / vsig)
     
     ! Magnetic force (plus symmetric operator for induction equation damping)
#if defined(GRAD_H_SPH)
     ! Add gather-neighbour contribution
     if (drmag < KERNRANGE*hp) then
        factors = mpp * pfactor_p * hfactor_p * w1(drmag*invhp)
        amag_temp(1:VDIM) = amag_temp(1:VDIM) + matmul(maxwell_p, dr_unit) *&
             & factors
        tensile_temp = tensile_temp + factors * dot_product(B_p, dr_unit)
        dBdt_damp = dBdt_damp + factors * phi_p
     end if

     ! Add scatter-neighbour contribution
     if (drmag < KERNRANGE*hpp) then
        invhpp = 1.0_PR / hpp
        hfactor_pp = invhpp**(NDIMPLUS1)
        pfactor_pp  = 1.0_PR / (rho_pp**2 * sph(pp)%omega)
        factors = mpp * pfactor_pp * hfactor_pp * w1(drmag*invhpp)
        amag_temp(1:VDIM) = amag_temp(1:VDIM) + matmul(maxwell_pp, dr_unit) *&
             & factors
        tensile_temp = tensile_temp + factors * dot_product(B_pp, dr_unit)
        dBdt_damp = dBdt_damp + factors * phi_p
     end if
#else
     ! Add two kernel contributions
     if (drmag < KERNRANGE*hp) wmean = 0.5_PR*hfactor_p*w1(drmag*invhp)
     if (drmag < KERNRANGE*hpp) wmean = wmean + &
          &0.5_PR*invhpp**(NDIMPLUS1)*w1(drmag*invhpp)
     pfactor_pp  = 1.0 / rho_pp**2
     factors = mpp * wmean
     amag_temp(1:VDIM) = amag_temp(1:VDIM) + factors*&
          &(pfactor_p*matmul(maxwell_p, dr_unit) + &
           &pfactor_pp*matmul(maxwell_pp, dr_unit))
     tensile_temp = tensile_temp + factors*&
          &(pfactor_p*dot_product(B_p, dr_unit) + &
            pfactor_pp*dot_product(B_pp, dr_unit))
#endif

     ! Induction equation, hyperbolic damping
     factors = mpp * hfactor_p * w1(drmag*invhp)
     dBdt_temp = dBdt_temp + factors * (dv * dot_product(B_p, dr_unit) - &
                                       &B_p * dot_product(dv, dr_unit))
     
     dphi_dt_temp = dphi_dt_temp + factors * dot_product(B_pp - B_p, dr_unit)
     
     dphi_dt_econs = dphi_dt_econs + factors * dot_product(dv, dr_unit)
     
     factors = factors * 0.5 * (alpha_p + alpha_pp) / rhomean**2
     dBdt_diss = dBdt_diss + factors * vsig_resist * (B_pp - B_p)
     
  end do
! ============================================================================

  amag_temp = -amag_temp + (BETA_TENSILE * B_p * tensile_temp)

  ! Fastest cleaning velocity / shortest smoothing length
  sph(p)%B_t_signal = B_t_signal / CLEAN_FAST_WAVE_SPEED
  
! Record magnetic acceleration in main array
  sph(p)%a(1:VDIM) = sph(p)%a(1:VDIM) + amag_temp(1:VDIM)
#if defined(DEBUG_FORCES)
  sph(p)%a_mag(1:VDIM) = amag_temp(1:VDIM)
#endif

  ! Induction equation
#if defined(GRAD_H_SPH)
  dBdt_temp = dBdt_temp / (rho_p * sph(p)%omega)
#else
  dBdt_temp = dBdt_temp / rho_p
#endif
  sph(p)%dBdt = dBdt_temp + rho_p * (dBdt_diss + dBdt_damp)
  
  ! Hyperbolic cleaning
  
  ! Take the cleaning speed from smoothing length and minimum signal time
  c_H = hp / sph(p)%B_t_signal
#if defined(GRAD_H_SPH)
  dphi_dt_temp = c_H**2 * dphi_dt_temp / (sph(p)%omega * rho_p)
  dphi_dt_econs = 0.5_PR * phi_p * dphi_dt_econs / (sph(p)%omega * rho_p)
#else
  dphi_dt_temp = c_H**2 * dphi_dt_temp / rho_p
  dphi_dt_econs = 0.5_PR * phi_p * dphi_dt_econs / rho_p
#endif
  sph(p)%dphi_dt = dphi_dt_temp - (phi_p*SIGMA_DAMP*c_H/hp) + dphi_dt_econs

#if defined(DIV_B)
  sph(p)%div_B = -dphi_dt_temp / c_H**2
#endif

  if (allocated(pp_templist)) deallocate(pp_templist)
  
  return

  contains
  
  function maxwell_stress_tensor(B)
     ! Horribly unoptimized...
     real(kind=PR)             :: maxwell_stress_tensor(1:BDIM, 1:BDIM)
     real(kind=PR), intent(in) :: B(1:BDIM)
     real(kind=PR)             :: Bsqd
     integer                   :: i, j
     
     Bsqd = sum(B**2)
     do i=1,BDIM
        do j=1,BDIM
           if (i==j) then
              maxwell_stress_tensor(i,j) = B(i) * B(j) - 0.5 * Bsqd
           else
              maxwell_stress_tensor(i,j) = B(i) * B(j)
           end if
        end do
     end do
     
     maxwell_stress_tensor = maxwell_stress_tensor / mu_0
  
     return
  end function

END SUBROUTINE mhd
