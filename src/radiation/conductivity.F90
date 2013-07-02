! CONDUCTIVITY.F90
! Dimitris Stamatellos 2/09/09
! Calculates thermal conductivity of p
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE conductivity(p)
  use interface_module, only : distance2,getkappa,get_neib_on_fly,w0,w2
  use particle_module
  use neighbour_module
  use kernel_module
  use hydro_module
  use constant_module
  use eos_module
  use scaling_module
  implicit none

  integer, intent(in) :: p               ! particle identifier

  integer :: i                           ! neighbour counter
  integer :: pp                          ! neighbouring particle number
  integer :: pp_numb                     ! number of neighbours
  integer, allocatable :: pp_templist(:) ! temp. list of neighbours
  real(kind=PR) :: dr(1:NDIM)            ! Relative displacement vector
  real(kind=PR) :: drmag                 ! magnitude of separation
  real(kind=PR) :: drsqd                 ! separation squared
  real(kind=PR) :: dr_unit(1:NDIM)       ! Unit displacement vector
  real(kind=PR) :: hfactor               ! invhp ^ NDIM
  real(kind=PR) :: hp                    ! smoothing length of particle p
  real(kind=PR) :: hpp                   ! smoothing length of neighbour pp
  real(kind=PR) :: invhp                 ! ( 1 / hp )
  real(kind=PR) :: invdrmag              ! ( 1 / drmag )
  real(kind=PR) :: lambda_diff           ! ..
  real(kind=PR) :: rho_p                 ! density of particle p
  real(kind=PR) :: rho_pp                ! density of particle pp
  real(kind=PR) :: rp(1:NDIM)            ! position of particle p
  real(kind=PR) :: wmean                 ! (W(p) + W(pp)) / 2
  real(kind=PR) :: radenergy             ! radiation energy 
  real(kind=PR) :: radenergygrad(1:NDIM) ! radiation energy gradient
  real(kind=PR) :: R_diff                ! R diffusion term 
  real(kind=PR) :: kappaT                ! Pseudo opacity
  real(kind=PR) :: kappapT               ! ..
  real(kind=PR) :: kapparT               ! ..

! Create local copies of important properties of particle p
  rp(1:NDIM) = sph(p)%r(1:NDIM)
  hp         = sph(p)%h
  invhp      = 1.0_PR / hp
  hfactor    = invhp**(NDIM)
  rho_p      = sph(p)%rho
#if defined(NEIGHBOUR_LISTS)
  pp_numb = pptot(p)
  if (pp_numb <= pp_limit) then 
     allocate(pp_templist(1:pp_numb))
     pp_templist(1:pp_numb) = pplist(1:pp_numb,p)
  else
     allocate(pp_templist(1:ptot))
     call get_neib_on_fly(p,pp_numb,ptot,pp_templist,rp(1:NDIM),hp)
  end if
#else
  allocate(pp_templist(1:ptot))
  call get_neib_on_fly(p,pp_numb,ptot,pp_templist,rp(1:NDIM),hp)
#endif

! Initialise summation variables, including self-contribution
  radenergygrad(1:NDIM) = 0.0_PR
  radenergy = (sph(p)%m/rho_p)*sph(p)%temp**4*hfactor*w0(0.0_PR)  


! Loop over all neighbours and sum terms
! ----------------------------------------------------------------------------
  do i=1,pp_numb
     pp = pp_templist(i)
     if (p == pp) cycle
     
     ! Create local copies for neighbour pp
     rho_pp  = sph(pp)%rho
     call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
     drmag = sqrt(drsqd) + SMALL_NUMBER
     if (drmag > KERNRANGE*hp) cycle
     invdrmag = 1.0_PR / drmag
     dr_unit(1:NDIM) = dr(1:NDIM)*invdrmag
     wmean = hfactor*invhp*w2(drmag*invhp)
 
     radenergygrad(1:NDIM) = radenergygrad(1:NDIM) + (sph(pp)%m/rho_p)&
          &*(sph(pp)%temp**4 - sph(p)%temp**4)*wmean*dr_unit(1:NDIM)
  
     radenergy = radenergy + &
          &(sph(pp)%m/rho_pp)*sph(pp)%temp**4*hfactor*w0(drmag*invhp)
  end do
! ----------------------------------------------------------------------------

  call getkappa(sph(p)%rho,sph(p)%temp,sph(p)%idens,kappaT,kapparT,kappapT)


  if (dot_product(radenergygrad,radenergygrad) == 0.0_PR) then
     R_diff = BIG_NUMBER
  else
     R_diff = sqrt(dot_product(radenergygrad,radenergygrad))&
          &/radenergy/rho_p/kapparT
  end if
  lambda_diff = (2.0_PR + R_diff) / (6.0_PR + 3.0_PR*R_diff + R_diff**2)
  sph(p)%k_cond = 16.0_PR*rad_const*lambda_diff*sph(p)%temp**3/rho_p/kapparT

  deallocate(pp_templist)

  return
END SUBROUTINE conductivity
