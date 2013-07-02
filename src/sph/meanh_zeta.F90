! MEANH_ZETA.F90
! D. A. Hubber - 20/06/2011
! Calculates value of the 'grad-h' gravitational force correction term, 
! zeta, for particle p.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE meanh_zeta(p)
  use interface_module, only : distance2,get_neib_on_fly,wzeta
  use neighbour_module
  use hydro_module
  use kernel_module
  use particle_module
#if defined(NBODY_SPH_SIMULATION) && defined(SINKS)
  use sink_module
  use Nbody_module
#endif
  implicit none

  integer, intent(in) :: p                ! id of particle

  integer :: i                            ! auxilary neighbour counter 
  integer :: pp                           ! neighbouring particles (p')
  integer :: pp_max                       ! ..
  integer :: pp_numb                      ! number of neighbours 
  integer, allocatable :: pp_templist(:)  ! temp. list of neighbours
  real(kind=PR) :: dr(1:NDIM)             ! relative position vector
  real(kind=PR) :: drmag                  ! magnitude of separation
  real(kind=PR) :: drsqd                  ! separation squared
  real(kind=PR) :: hp                     ! Smoothing length of particle p
  real(kind=PR) :: invhmean               ! 2 / (hp + hpp)
  real(kind=PR) :: mpp                    ! mass of neighbour pp
  real(kind=PR) :: rp(1:NDIM)             ! position of particle p
  real(kind=PR) :: zeta_p                 ! zeta for particle p
#if defined(NBODY_SPH_SIMULATION) && defined(SINKS)
  integer :: s                            ! Star counter
#endif

  debug3("Calculating zeta for meanh gravity [meanh_zeta.F90] for particle ",p)

  ! Create local copies and initialize all other quantities
  rp(1:NDIM) = sph(p)%r(1:NDIM)
  hp = sph(p)%h
  pp_numb = 0

  ! Find, or copy particle neighbour lists
  ! --------------------------------------------------------------------------
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

  invhmean = 1.0_PR/sph(p)%h
  zeta_p = sph(p)%m*invhmean*invhmean*wzeta(0.0_PR)

  ! Now loop over all neighbours to find contributions
  ! --------------------------------------------------------------------------
  do i=1,pp_numb
     pp = pp_templist(i)
     if (p == pp) cycle
     mpp = sph(pp)%m
     invhmean = 2.0_PR/(hp + sph(pp)%h)
     call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
     drmag = sqrt(drsqd) + SMALL_NUMBER
     if (drmag*invhmean > KERNRANGE) cycle
     zeta_p = zeta_p + mpp*invhmean*invhmean*wzeta(drmag*invhmean)
  end do
  ! --------------------------------------------------------------------------

  ! Now loop over all sink/star particles to add contributions
  ! --------------------------------------------------------------------------
#if defined(NBODY_SPH_SIMULATION) && defined(SINKS)
  do s=1,stot
     mpp = star(s)%m
     invhmean = 2.0_PR/(hp + star(s)%h)
     call distance3(rp(1:NDIM),star(s)%r(1:NDIM),dr(1:NDIM),drsqd)
     drmag = sqrt(drsqd) + SMALL_NUMBER
     if (drmag*invhmean > KERNRANGE) cycle
     zeta_p = zeta_p + mpp*invhmean*invhmean*wzeta(drmag*invhmean)
  end do
#endif
  ! --------------------------------------------------------------------------

  ! Normalise and store zeta in main array.
  sph(p)%zo = -(hp*zeta_p) / (NDIMPR*sph(p)%rho*sph(p)%omega)

  if (allocated(pp_templist)) deallocate(pp_templist)

  return
END SUBROUTINE meanh_zeta
