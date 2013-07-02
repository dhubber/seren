! DIFFUSION.F90
! Dimitris Stamatellos - 02/09/2009
! Calculates thermal conductivity of p
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE diffusion(p,dt)
  use interface_module, only : distance2,getkappa,get_neib_on_fly,w2
  use particle_module
  use neighbour_module
  use kernel_module
  use hydro_module
  use constant_module
  use eos_module
  use scaling_module
  implicit none

  integer, intent(in) :: p               ! particle identifier
  real(kind=PR), intent(in) :: dt        ! ..
  integer :: i                           ! neighbour counter
  integer :: pp                          ! neighbouring particle number
  integer :: pp_numb                     ! number of neighbours
  integer, allocatable :: pp_templist(:) ! temp. list of neighbours
  real(kind=PR) :: dr(1:NDIM)            ! Relative displacement vector
  real(kind=PR) :: drmag                 ! magnitude of separation
  real(kind=PR) :: dr_unit(1:NDIM)       ! Unit displacement vector
  real(kind=PR) :: dt_diff_pp            ! energy diff. timescale from p to pp
  real(kind=PR) :: dudt_diff_p           ! energy diff. rate from p
  real(kind=PR) :: dudt_diff_pp          ! energy diff. rate from pp
  real(kind=PR) :: du_diff_p             ! energy diff.  from p
  real(kind=PR) :: du_diff_pp            ! energy diff. from pp
  real(kind=PR) :: du_p_pp               ! ..
  real(kind=PR) :: hfactor_p             ! invhp ^ NDIM
  real(kind=PR) :: hfactor_pp            ! invhpp ^ NDIM
  real(kind=PR) :: hp                    ! smoothing length of particle p
  real(kind=PR) :: hpp                   ! smoothing length of particle pp
  real(kind=PR) :: invhp                 ! ( 1 / hp )
  real(kind=PR) :: invhpp                ! ( 1 / hpp )
  real(kind=PR) :: invdrmag              ! ( 1 / drmag )
  real(kind=PR) :: kcond_p               ! ..
  real(kind=PR) :: rho_p                 ! density of particle p
  real(kind=PR) :: rho_pp                ! density of particle pp
  real(kind=PR) :: rp(1:NDIM)            ! position of particle p
  real(kind=PR) :: up                    ! specific internal energy of p
  real(kind=PR) :: vp(1:NDIM)            ! velocity of particle p
  real(kind=PR) :: wmean                 ! (W(p) + W(pp)) / 2
  real(kind=PR) :: kappaT                ! ..
  real(kind=PR) :: kappaT_pp             ! ..
  real(kind=PR) :: kappapT               ! ..
  real(kind=PR) :: kappapT_pp            ! ..
  real(kind=PR) :: kapparT               ! ..
  real(kind=PR) :: kapparT_pp            ! ..
  real(kind=PR) :: tau_p                 ! ..
  real(kind=PR) :: tau_pp                ! ..
  real(kind=PR) :: temp_p                ! ..
  real(kind=PR) :: drsqd                 ! separation squared

! Create local copies of important properties of particle p
  rp(1:NDIM) = sph(p)%r(1:NDIM)
  hp         = sph(p)%h
  invhp      = 1.0_PR / hp
  hfactor_p  = invhp**(NDIM)
  vp(1:NDIM) = sph(p)%v(1:NDIM)
  rho_p      = sph(p)%rho
  up         = sph(p)%u
  kcond_p    = sph(p)%k_cond
  temp_p     = sph(p)%temp

#if defined(NEIGHBOUR_LISTS)
  pp_numb = pptot(p)
  if (pp_numb <= pp_limit) then 
     allocate(pp_templist(1:pp_numb))
     do i=1,pp_numb
        pp_templist(i) = pplist(i,p)
     end do
  else
     allocate(pp_templist(1:ptot))
     call get_neib_on_fly(p,pp_numb,ptot,pp_templist,rp(1:NDIM),hp)
  end if
#else
  allocate(pp_templist(1:ptot))
  call get_neib_on_fly(p,pp_numb,ptot,pp_templist,rp(1:NDIM),hp)
#endif


! Initialise variables
  du_diff_p = 0.0_PR
  dudt_diff_p = 0.0_PR
  call getkappa(rho_p,sph(p)%temp,sph(p)%idens,kappaT,kapparT,kappapT)
  tau_p = sqrt(sph(p)%column2)*kappaT


! Loop over all neighbours and sum terms
! ----------------------------------------------------------------------------
  do i=1,pp_numb
     pp = pp_templist(i)
     if (p == pp) cycle
 
     ! Create local copies for neighbour pp     
     rho_pp = sph(pp)%rho
     hpp = sph(pp)%h
     invhpp = 1.0_PR / hpp
     hfactor_pp = invhpp**(NDIM)
     call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
     drmag = sqrt(drsqd) + SMALL_NUMBER
     if (drmag > KERNRANGE*hp) cycle
     invdrmag = 1.0_PR / drmag
     dr_unit(1:NDIM) = dr(1:NDIM)*invdrmag

     wmean = 0.5_PR*(hfactor_p*invhp*w2(drmag*invhp) + &
          &hfactor_pp*invhpp*w2(drmag*invhpp))

     call getkappa(rho_pp,sph(pp)%temp,sph(pp)%idens,&
          &kappaT_pp,kapparT_pp,kappapT_pp)
     tau_pp = sqrt(sph(pp)%column2)*kappaT_pp

     ! Skip particle if optically thin
     if ((kapparT_pp*rho_pp + kapparT*rho_p)*drmag < 2.0_PR) cycle
     
     ! transfer of energy from particle pp to particle p (sign is correct)
     dudt_diff_pp = 4.0_PR*(sph(pp)%m/(rho_pp*rho_p)) * &
             &(kcond_p*sph(pp)%k_cond/(kcond_p + sph(pp)%k_cond)) * &
             &(temp_p - sph(pp)%temp)*wmean*invdrmag

     ! calculate diffusion timescale
     dt_diff_pp = 0.5_PR*abs((up + sph(pp)%u)/dudt_diff_pp)

     ! calculate energy flow to particle p assuming a constant 
     ! flow rate during the dynamical timestep dt
     du_diff_pp = dudt_diff_pp*dt

     ! energy flow cannot be larger than the energy difference 
     ! between the two particles
     du_p_pp = 0.5_PR*(sph(pp)%u - up)
     du_diff_pp = min(abs(du_diff_pp),abs(du_p_pp))*(du_p_pp/abs(du_p_pp))
     du_diff_p = du_diff_p + du_diff_pp 
     
  end do
! ----------------------------------------------------------------------------

  sph(p)%du_dt_diff = du_diff_p/dt
 
  deallocate(pp_templist)

  return
END SUBROUTINE diffusion
