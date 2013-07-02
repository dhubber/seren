! HYDRO_GRADH.F90
! C. P. Batty & D. A. Hubber - 27/3/2007
! Calculates hydro acceleration of particle p
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE hydro_gradh(p)
  use interface_module, only : distance2,cooling_rate,effective_gamma,&
       &get_neib_on_fly,isothermal_riemann_solver,riemann_solver,w1
  use particle_module
  use neighbour_module
  use hydro_module
  use kernel_module
  use time_module
  use type_module
  implicit none

  integer, intent(in) :: p               ! particle identifier

  integer :: i                           ! neighbour counter
  integer :: pp                          ! neighbouring particle number
  integer :: pp_numb                     ! number of neighbours
  integer, allocatable :: pp_templist(:) ! temp. list of neighbours
  real(kind=PR) :: ahydro_temp(1:NDIM)   ! hydro acceleration of particle p
  real(kind=PR) :: alpha_mean            ! mean value of alpha
  real(kind=PR) :: drmag                 ! magnitude of separation
  real(kind=PR) :: dr_unit(1:NDIM)       ! Unit displacement vector
  real(kind=PR) :: dv(1:NDIM)            ! Relative velocity
  real(kind=PR) :: dvdr                  ! Scalar product of dv and dr
  real(kind=PR) :: hfactor_p             ! invhp ^ (NDIM + 1)
  real(kind=PR) :: hfactor_pp            ! invhpp ^ (NDIM + 1)
  real(kind=PR) :: hp                    ! smoothing length of particle p
  real(kind=PR) :: hpp                   ! smoothing length of neighbour pp
  real(kind=PR) :: invhp                 ! ( 1 / hp )
  real(kind=PR) :: invhpp                ! ( 1 / hpp )
  real(kind=PR) :: invdrmag              ! ( 1 / drmag )
  real(kind=PR) :: pfactor_p             ! press / rho / rho / omega
  real(kind=PR) :: pfactor_pp            ! press / rho / rho / omega for pp
  real(kind=PR) :: mpp                   ! mass of neighbour pp
  real(kind=PR) :: rho_p                 ! density of particle p
  real(kind=PR) :: rho_pp                ! density of particle pp
  real(kind=PR) :: rp(1:NDIM)            ! position of particle p
  real(kind=PR) :: sound_p               ! sound speed of particle p
  real(kind=PR) :: up                    ! Specific internal energy of p
  real(kind=PR) :: vp(1:NDIM)            ! velocity of particle p
  real(kind=PR) :: wmean                 ! (W(p) + W(pp)) / 2
#if defined(SIGNAL_VELOCITY)
  real(kind=PR) :: vsigmax_p             ! Max. signal velocity of p
#endif
#if defined(RIEMANN_SOLVER)
  real(kind=PR) :: gamma_eff             ! Effective ratio of specific heats
  real(kind=PR) :: Pstar                 ! Mean pressure variable
  real(kind=PR) :: vstar                 ! Mean velocity variable
#endif
#if defined(ARTIFICIAL_VISCOSITY) || defined(ARTIFICIAL_CONDUCTIVITY)
  real(kind=PR) :: invrhomean            ! 1 / mean density
  real(kind=PR) :: vsignal               ! signal velocity
#if defined(ENERGY_EQN)
  real(kind=PR) :: dudt_diss             ! dissipative dudt
#endif
#endif
#if defined(ARTIFICIAL_VISCOSITY)
  real(kind=PR) :: ahydro_diss(1:NDIM)   ! dissipative hydro accel
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
  real(kind=PR) :: pattrec_p             ! Patern recognition factor for p
#endif
#endif
integer :: j
  debug3("Calculating hydro forces [hydro_gradh.F90] for particle ", p)

! Create local copies of important properties of particle p
  pp_numb    = 0
  rp(1:NDIM) = sph(p)%r(1:NDIM)
  hp         = sph(p)%h
  invhp      = 1.0_PR / hp
  hfactor_p  = invhp**(NDIMPLUS1)
  vp(1:NDIM) = sph(p)%v(1:NDIM)
  sound_p    = sph(p)%sound
  rho_p      = sph(p)%rho
  up         = sph(p)%u
#ifndef RIEMANN_SOLVER
  pfactor_p  = sph(p)%press / rho_p / rho_p / sph(p)%omega
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

! Zero arrays
  ahydro_temp(1:NDIM) = 0.0_PR
#if defined(SIGNAL_VELOCITY)
  vsigmax_p = 0.0_PR
#endif
#if defined(ARTIFICIAL_VISCOSITY)
  ahydro_diss(1:VDIM) = 0.0_PR
#endif
#if defined(ENERGY_EQN)
  sph(p)%dudt = 0.0_PR
#if defined(ARTIFICIAL_VISCOSITY) || defined(ARTIFICIAL_CONDUCTIVITY)
  dudt_diss = 0.0_PR
#endif
#endif

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


! Loop over all neighbours, summing all contributions to the hydro accel.
! and other hydro quantities
! ============================================================================
  do i=1,pp_numb
     pp = pp_templist(i)
     if (p == pp) cycle

     ! Create local copies for neighbour pp
     wmean  = 0.0_PR
     rho_pp = sph(pp)%rho
     mpp    = sph(pp)%m
     hpp    = sph(pp)%h
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
#if defined(SIGNAL_VELOCITY)
     vsigmax_p = max(vsigmax_p,sound_p,&
          &sound_p + sph(pp)%sound - beta*dvdr*invdrmag)
#endif

     ! Add gather-neighbour contribution
     if (drmag < KERNRANGE*hp) then
        wmean = 0.5_PR*hfactor_p*w1(drmag*invhp)
        ahydro_temp(1:NDIM) = ahydro_temp(1:NDIM) + &
             & mpp*pfactor_p*hfactor_p*w1(drmag*invhp)*dr_unit(1:NDIM)
     end if

     ! Add scatter-neighbour contribution
     if (drmag < KERNRANGE*hpp) then
        invhpp = 1.0_PR / hpp
        hfactor_pp = invhpp**(NDIMPLUS1)
        pfactor_pp = sph(pp)%press / rho_pp / rho_pp / sph(pp)%omega
        wmean  = wmean + 0.5_PR*hfactor_pp*w1(drmag*invhpp)
        ahydro_temp(1:NDIM) = ahydro_temp(1:NDIM) + &
             & mpp*pfactor_pp*hfactor_pp*w1(drmag*invhpp)*dr_unit(1:NDIM)
     end if


     ! Add all artificial dissipation terms
     ! -----------------------------------------------------------------------
#if defined(ARTIFICIAL_VISCOSITY) || defined(ARTIFICIAL_CONDUCTIVITY)
     invrhomean = 2.0_PR/(rho_p + rho_pp)

     ! Only add artificial viscosity if particles are approaching
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
             &(drmag*drmag + 0.25_PR*ETA_SQD*(hp + hpp)*(hp + hpp))
        visc_factor = (-alpha_mean*vsignal*mu + 2.0_PR*alpha_mean*mu*mu)*&
             &invrhomean
        ahydro_diss(1:NDIM) = ahydro_diss(1:NDIM) + &
             & mpp*visc_factor*wmean*dr_unit(1:NDIM)
#elif defined(ARTIFICIAL_VISCOSITY) && defined(VISC_MON97)
        vsignal = sound_p + sph(pp)%sound - beta*dvdr*invdrmag
        ahydro_diss(1:NDIM) = ahydro_diss(1:NDIM) - mpp*alpha_mean*vsignal* &
             & dvdr*invdrmag*dr_unit(1:NDIM)*wmean*invrhomean
#endif
        
        ! Viscous heating term
#if defined(ENERGY_EQN) && defined(VISC_AB)
        dudt_diss = dudt_diss + 0.5_PR*mpp*visc_factor*wmean*dvdr*invdrmag
#elif defined(ENERGY_EQN) && defined(VISC_MON97)
        dudt_diss = dudt_diss - 0.5_PR*mpp*alpha_mean*drmag*invdrmag* &
             & vsignal*wmean*(dot_product(dv,dr_unit)**2)*invrhomean
#endif
     end if

     ! Artificial conductivity term
#if defined(EXPLICIT_ENERGY_EQN) && defined(ARTIFICIAL_CONDUCTIVITY)
#if defined(COND_PRICE2008)
     vsignal = sqrt(abs(sph(p)%press - sph(pp)%press)*invrhomean)
#elif defined(COND_WADSLEY2008)
     vsignal = abs(dvdr*invdrmag)
#endif
     dudt_diss = dudt_diss + &
          & mpp*alpha_cond*vsignal*(up - sph(pp)%u)*wmean*invrhomean
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
#if defined(SIGNAL_VELOCITY)
  sph(p)%vsigmax = vsigmax_p
#endif

! Record hydro acceleration in main array
  sph(p)%a(1:NDIM) = sph(p)%a(1:NDIM) + ahydro_temp(1:NDIM)
#if defined(DEBUG_FORCES)
  sph(p)%a_hydro(1:NDIM) = ahydro_temp(1:NDIM)
#endif

! Record artificial viscosity variables
#if defined(ARTIFICIAL_VISCOSITY) && defined(VISC_TD)
  sph(p)%dalpha_dt = (C_1*(alpha_min - sph(p)%talpha)*sound_p*invhp) + &
       & max(-sph(p)%div_v,0.0_PR)*(alpha - sph(p)%talpha)
#endif

! Dissipation terms in entropy formulation
! ----------------------------------------------------------------------------
#if defined(ENTROPIC_FUNCTION)
#if defined(ARTIFICIAL_VISCOSITY) || defined(ARTIFICIAL_CONDUCTIVITY)
  sph(p)%dAdt = (gamma - 1.0_PR)*dudt_diss/rho_p**(gamma - 1.0_PR)
#else
  sph(p)%dAdt = 0.0_PR
#endif

! Add main contribution to change in internal energy (plus artificial 
! conductivity if applicable) and record in main array
! ----------------------------------------------------------------------------
#elif defined(ENERGY_EQN) 
  sph(p)%dudt = sph(p)%dudt - sph(p)%press*sph(p)%div_v/sph(p)%omega/rho_p
#if defined(ARTIFICIAL_VISCOSITY) || defined(ARTIFICIAL_CONDUCTIVITY)
  sph(p)%dudt = sph(p)%dudt + dudt_diss
#endif
#if defined(COOLING_HEATING) && defined(EXPLICIT_COOLING_HEATING)
  sph(p)%dudt = sph(p)%dudt + cooling_rate(p)
#endif
#endif

! MPI loadbalancing instrumentation
#if defined(USE_MPI) && !defined(LOW_MEM)
  if (p > 0) sph(p)%hydro_calcs = sph(p)%hydro_calcs + real(pp_numb,PR)
#endif

#if defined(DEBUG_VERIFY_ALL)
  do pp=1,ptot+pghost
     if (pp == p) cycle
     dr_unit(1:NDIM) = sph(p)%r(1:NDIM) - sph(pp)%r(1:NDIM)
     drmag = sqrt(dot_product(dr_unit(1:NDIM),dr_unit(1:NDIM)))
     if (drmag <= KERNRANGE*sph(p)%h .or. drmag <= KERNRANGE*sph(pp)%h) then
        j = 0
        do i=1,pp_numb
           if (pp == pp_templist(i)) j = 1
        end do
        if (j == 0) then
           write(6,*) "Error finding neighbour in tree : ",&
                &p,pp,ptot,pghost,drmag/(KERNRANGE*sph(p)%h),drmag/(KERNRANGE*sph(pp)%h),sph(p)%h/sph(pp)%h
           stop
        end if
        !write(6,*) "Particle verified : ",p,pp,ptot,pghost,rp(1:NDIM)
     end if
  end do
#endif

  if (allocated(pp_templist)) deallocate(pp_templist)
  
  return
END SUBROUTINE hydro_gradh
