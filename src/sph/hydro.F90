! HYDRO.F90
! C. P. Batty & D. A. Hubber - 27/3/2007
! Calculates hydro acceleration of particle p
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE hydro(p)
  use interface_module, only : cooling_rate,distance2,effective_gamma,&
       &get_neib_on_fly,isothermal_riemann_solver,riemann_solver,w1
  use particle_module
  use neighbour_module
  use hydro_module
  use kernel_module
  use time_module
  implicit none

  integer, intent(in) :: p               ! particle identifier

  integer :: i                           ! neighbour counter
  integer :: pp                          ! neighbouring particle number
  integer :: pp_numb                     ! number of neighbours
  integer, allocatable :: pp_templist(:) ! temp. list of neighbours
  real(kind=PR) :: ahydro_temp(1:NDIM)   ! hydro acceleration of particle p
  real(kind=PR) :: dr(1:NDIM)            ! vector displacements (p'-p)
  real(kind=PR) :: drmag                 ! magnitude of separation
  real(kind=PR) :: drsqd                 ! separation squared
  real(kind=PR) :: dr_unit(1:NDIM)       ! Unit displacement vector
  real(kind=PR) :: dv(1:NDIM)            ! Relative velocity
  real(kind=PR) :: dvdr                  ! Scalar product of dv and dr
  real(kind=PR) :: hfactor_p             ! invhp ^ NDIM
  real(kind=PR) :: hp                    ! smoothing length of particle p
  real(kind=PR) :: hpp                   ! smoothing length of neighbour pp
  real(kind=PR) :: invhp                 ! ( 1 / hp )
  real(kind=PR) :: invhpp                ! ( 1 / hpp )
  real(kind=PR) :: invdrmag              ! ( 1 / drmag )
  real(kind=PR) :: pfactor_p             ! press / rho / rho ( / omega)
  real(kind=PR) :: mpp                   ! mass of neighbour pp
  real(kind=PR) :: rho_p                 ! density of particle p
  real(kind=PR) :: rho_pp                ! density of particle pp
  real(kind=PR) :: rp(1:NDIM)            ! position of particle p
  real(kind=PR) :: sound_p               ! sound speed of particle p
  real(kind=PR) :: up                    ! Specific internal energy of p
  real(kind=PR) :: vp(1:NDIM)            ! velocity of particle p
  real(kind=PR) :: wmean                 ! (W(p) + W(pp)) / 2
#if defined(ARTIFICIAL_VISCOSITY) || defined(ARTIFICIAL_CONDUCTIVITY)
  real(kind=PR) :: rhomean               ! mean density
  real(kind=PR) :: vsignal               ! signal velocity
#if defined(ENERGY_EQN)
  real(kind=PR) :: dudt_diss             ! dissipative dudt
#endif
#endif
#if defined(ARTIFICIAL_VISCOSITY)
  real(kind=PR) :: ahydro_diss(1:NDIM)   ! dissipative hydro accel
  real(kind=PR) :: alpha_mean            ! mean vlaue of alpha
#if defined(VISC_AB)
  real(kind=PR) :: mu                    ! pseudo divergence term
  real(kind=PR) :: visc_factor           ! artificial viscosity
#endif
#if defined(VISC_TD)
  real(kind=PR) :: talpha_p              ! TD value of alpha for p
#endif
#if defined(VISC_BALSARA)
  real(kind=PR) :: balsara_p             ! Balsara factor for particle p
#endif
#if defined(VISC_PATTERN_REC)
  real(kind=PR) :: pattrec_p             ! Pattern recognition factor for p
#endif
#endif
#if defined(SM2012_SPH)
  real(kind=PR) :: q_p                   ! Energy density of p
#endif
#if defined(ENERGY_EQN)
  real(kind=PR) :: dudt_temp             ! dudt of particle p
#if defined(ARTIFICIAL_CONDUCTIVITY) 
  real(kind=PR) :: vsignal_u             ! Second signal velocity for cond
#endif
#endif

  debug3("Calculating hydro forces [hydro.F90] for particle ", p)

! Create local copies of important properties of particle p
  rp(1:NDIM) = sph(p)%r(1:NDIM)
  hp         = sph(p)%h
  invhp      = 1.0_PR / hp
  hfactor_p  = invhp**(NDIMPLUS1)
  vp(1:NDIM) = sph(p)%v(1:NDIM)
  sound_p    = sph(p)%sound
  rho_p      = sph(p)%rho
  up         = sph(p)%u
#if defined(SM2012_SPH)
  q_p        = sph(p)%q
#endif
#if defined(EXTERNAL_PRESSURE)
  pfactor_p  = (sph(p)%press - Pext)*sph(p)%invrho*sph(p)%invrho
#else
  pfactor_p  = sph(p)%press*sph(p)%invrho*sph(p)%invrho
#endif
#if defined(ARTIFICIAL_VISCOSITY) && defined(VISC_TD)
  talpha_p = sph(p)%talpha
#endif
#if defined(ARTIFICIAL_VISCOSITY) && defined(VISC_BALSARA)
  balsara_p = sph(p)%balsara
#endif
#if defined(ARTIFICIAL_VISCOSITY) && defined(VISC_PATTERN_REC)
  pattrec_p = sph(p)%pattrec
#endif

! Zero accumulation variables
  ahydro_temp(1:NDIM) = 0.0_PR
#if defined(ARTIFICIAL_VISCOSITY)
  ahydro_diss(1:VDIM) = 0.0_PR
#endif
#if defined(ENERGY_EQN)
  sph(p)%dudt = 0.0_PR
  dudt_temp = 0.0_PR
#if defined(ARTIFICIAL_VISCOSITY) || defined(ARTIFICIAL_CONDUCTIVITY)
  dudt_diss = 0.0_PR
#endif
#endif

#if defined(NEIGHBOUR_LISTS)
  pp_numb = pptot(p)
  if (pp_numb <= pp_limit) then 
     allocate(pp_templist(1:pp_numb))
     do i=1,pp_numb
        pp_templist(i) = pplist(i,p)
     end do
  else
     call get_neib_on_fly(p,pp_numb,ptot,pp_templist,rp(1:NDIM),KERNRANGE*hp)
  end if
#else
  call get_neib_on_fly(p,pp_numb,ptot,pp_templist,rp(1:NDIM),KERNRANGE*hp)
#endif


! Loop over all neighbours, summing each (p'-p) acceleration component
! ============================================================================
  do i=1,pp_numb
     pp = pp_templist(i)

     ! Create local copies for neighbour pp
     wmean = 0.0_PR
     rho_pp = sph(pp)%rho
     mpp = sph(pp)%m
     hpp = sph(pp)%h
     invhpp = 1.0_PR/hpp
#if defined(PERIODIC)
     call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
#else
     dr(1:NDIM) = sph(pp)%r(1:NDIM) - rp(1:NDIM)
     drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif
     if (drsqd >= KERNRANGESQD*hp*hp .and. drsqd >= KERNRANGESQD*hpp*hpp) cycle
     drmag = sqrt(drsqd) + SMALL_NUMBER
     invdrmag = 1.0_PR / drmag
     dr_unit(1:NDIM) = dr(1:NDIM)*invdrmag
     dv(1:NDIM) = sph(pp)%v(1:NDIM) - vp(1:NDIM)
     dvdr = dot_product(dv(1:NDIM),dr(1:NDIM))

     ! Add two kernel contributions
     if (drmag < KERNRANGE*hp) wmean = 0.5_PR*hfactor_p*w1(drmag*invhp)
     if (drmag < KERNRANGE*hpp) wmean = wmean + &
          &0.5_PR*invhpp**(NDIMPLUS1)*w1(drmag*invhpp)


     ! Compute hydro pressure forces
     ! -----------------------------------------------------------------------
#if defined(SM2012_SPH)
     ahydro_temp(1:NDIM) = ahydro_temp(1:NDIM) + 0.5_PR*(gamma - 1.0_PR)*&
            &up*sph(pp)%u*mpp*(1.0_PR/q_p + 1.0_PR/sph(pp)%q)*dr_unit(1:NDIM)*&
            &(w1(drmag*invhp)*hfactor_p + w1(drmag*invhpp)*invhpp**(NDIMPLUS1))
#elif defined(EXTERNAL_PRESSURE)
     ahydro_temp(1:NDIM) = ahydro_temp(1:NDIM) + mpp*wmean*&
          &(pfactor_p + (sph(pp)%press - Pext)/&
          &(rho_pp*rho_pp))*dr_unit(1:NDIM)
#else
     ahydro_temp(1:NDIM) = ahydro_temp(1:NDIM) + mpp*wmean*&
          &(pfactor_p + sph(pp)%press/(rho_pp*rho_pp))*dr_unit(1:NDIM)
#endif

     ! Compute rate of change of internal energy (Need to check minus sign)
#if defined(ENERGY_EQN) && defined(SM2012_SPH)
     dudt_temp = dudt_temp + (gamma - 1.0_PR)*up*sph(pp)%u*mpp*&
          &dvdr*wmean*invdrmag/q_p
#endif

     ! Add all artificial dissipation terms
     ! -----------------------------------------------------------------------
#if defined(ARTIFICIAL_VISCOSITY) || defined(ARTIFICIAL_CONDUCTIVITY)
     rhomean = 0.5_PR*(rho_p + rho_pp)

     ! Only add artificial viscosity terms if particles are approaching
     if (dvdr < 0.0_PR) then

        ! Calculate effective value of alpha for viscosity 
        ! (i.e. time-dependent viscosity term plus Balsara switch)
#if defined(ARTIFICIAL_VISCOSITY) && defined(VISC_TD) && defined(VISC_BALSARA)
        alpha_mean = 0.25_PR*(talpha_p + sph(pp)%talpha)*(balsara_p + sph(pp)%balsara)
#elif defined(ARTIFICIAL_VISCOSITY) && defined(VISC_TD)
        alpha_mean = 0.5_PR*(talpha_p + sph(pp)%talpha)
#elif defined(ARTIFICIAL_VISCOSITY) && defined(VISC_BALSARA)
        alpha_mean = 0.5_PR*alpha*(balsara_p + sph(pp)%balsara)
#elif defined(ARTIFICIAL_VISCOSITY) && defined(VISC_PATTERN_REC)
        alpha_mean = 0.5_PR*alpha*(pattrec_p + sph(pp)%pattrec)
#elif defined(ARTIFICIAL_VISCOSITY) 
        alpha_mean = alpha
#endif

        ! Viscous pressure term 
#if defined(ARTIFICIAL_VISCOSITY) && defined(VISC_AB)
        vsignal = 0.5_PR*(sound_p + sph(pp)%sound)
        mu = 0.5_PR*(hp + hpp)*dvdr / &
             &(drsqd + 0.25_PR*ETA_SQD*(hp + hpp)*(hp + hpp))
        visc_factor = (-alpha_mean*vsignal*mu + &
             &2.0_PR*alpha_mean*mu*mu)/rhomean
        ahydro_diss(1:NDIM) = ahydro_diss(1:NDIM) + &
             & mpp*visc_factor*wmean*dr_unit(1:NDIM)
#elif defined(ARTIFICIAL_VISCOSITY) && defined(VISC_MON97)
        vsignal = sound_p + sph(pp)%sound - beta*dvdr*invdrmag
        ahydro_diss(1:NDIM) = ahydro_diss(1:NDIM) - mpp*alpha_mean*vsignal* &
             & dvdr*invdrmag*dr_unit(1:NDIM)*wmean/rhomean
#endif
        
        ! Viscous heating term
#if defined(ENERGY_EQN) && defined(VISC_AB)
        dudt_diss = dudt_diss + 0.5_PR*mpp*visc_factor*wmean*dvdr*invdrmag
#elif defined(ENERGY_EQN) && defined(VISC_MON97)
        dudt_diss = dudt_diss - 0.5_PR*mpp*alpha_mean*invdrmag*invdrmag* &
             &vsignal*wmean*(dot_product(dv(1:NDIM),dr(1:NDIM))**2)/rhomean
#endif
     end if

     ! Artificial conductivity term
#if defined(ENERGY_EQN) && defined(ARTIFICIAL_CONDUCTIVITY)
#if defined(COND_PRICE2008)
     vsignal = sqrt(abs(sph(p)%press - sph(pp)%press) / rhomean)
#elif defined(COND_WADSLEY2008)
     vsignal = abs(dvdr*invdrmag)
#endif
     dudt_diss = dudt_diss + &
          &mpp*alpha_cond*vsignal*(up - sph(pp)%u)*wmean/rhomean
#endif

#endif
     ! -----------------------------------------------------------------------

  end do
! ============================================================================

#if defined(ARTIFICIAL_VISCOSITY)
  ahydro_temp(1:VDIM) = ahydro_temp(1:VDIM) + ahydro_diss(1:VDIM)
#endif
#if defined(DEBUG_FORCES) && defined(ARTIFICIAL_VISCOSITY)
  sph(p)%a_visc(1:VDIM) = ahydro_diss(1:VDIM)
#endif

! Record hydro acceleration in main array
  sph(p)%a(1:NDIM) = sph(p)%a(1:NDIM) + ahydro_temp(1:NDIM)
#if defined(DEBUG_FORCES)
  sph(p)%a_hydro(1:VDIM) = ahydro_temp(1:VDIM)
#endif

! Record artificial viscosity variables
#if defined(ARTIFICIAL_VISCOSITY) && defined(VISC_TD)
  sph(p)%dalpha_dt = (C_1*(alpha_min - sph(p)%talpha)*sound_p/hp) + &
       & max(-sph(p)%div_v,0.0_PR)*(alpha - sph(p)%talpha)
#endif


! Dissipation terms in entropy formulation
! ----------------------------------------------------------------------------
#if defined(ENTROPIC_FUNCTION) && defined(ENERGY_EQN)
#if defined(ARTIFICIAL_VISCOSITY) || defined(ARTIFICIAL_CONDUCTIVITY)
  sph(p)%dAdt = (gamma - 1.0_PR)*dudt_diss/rho_p**(gamma - 1.0_PR)
#else
  sph(p)%dAdt = 0.0_PR
#endif

! Add main contribution to change in internal energy (plus artificial 
! conductivity if applicable) and record in main array
! ----------------------------------------------------------------------------
#elif defined(ENERGY_EQN) 
#if defined(ARTIFICIAL_VISCOSITY) || defined(ARTIFICIAL_CONDUCTIVITY)
  sph(p)%dudt = sph(p)%dudt + dudt_diss
#endif
#if defined(SM2012_SPH)
  sph(p)%dudt = sph(p)%dudt + dudt_temp
#else
  sph(p)%dudt = sph(p)%dudt - sph(p)%press*sph(p)%div_v/rho_p
#endif
#if defined(COOLING_HEATING) && defined(EXPLICIT_COOLING_HEATING)
  sph(p)%dudt = sph(p)%dudt + cooling_rate(p)
#endif
#endif

  deallocate(pp_templist)
  
! MPI loadbalancing instrumentation
#if defined(USE_MPI) && !defined(LOW_MEM)
  if (p > 0) sph(p)%hydro_calcs = sph(p)%hydro_calcs + real(pp_numb,PR)
#endif

  return
END SUBROUTINE hydro
