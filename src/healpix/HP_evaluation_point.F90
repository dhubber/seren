! H_DENSITY_EVALUATION_POINT.F90
! D. A. Hubber - 20/10/2010
! Calculates the properties at the evaluation point either by interpolation 
! or by performing a tree walk.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE HP_evaluation_point(rep,rmag_ep,rsqd_ep,rsource,pprev,&
     &p,plast,nextptcl,rhoep,hep,rtree,htree,rhotree,treeep)
  use definitions
  use particle_module
  use hydro_module
  use HP_module
  implicit none

  logical, intent(inout) :: treeep                ! Is e.p. from tree-walk?
  integer, intent(in) :: nextptcl(1:ptot)         ! Linked list of particles
  integer, intent(in) :: p                        ! id of particle after e.p.
  integer, intent(in) :: plast                    ! id of last particle in list
  integer, intent(in) :: pprev                    ! id of particle before e.p.
  real(kind=PR), intent(in) :: rep(1:NDIM)        ! Pos. of evaluation point
  real(kind=PR), intent(in) :: rsource(1:NDIM)    ! Pos. of ionizing source
  real(kind=PR), intent(in) :: rmag_ep            ! Dist. of e.p. from source
  real(kind=PR), intent(in) :: rsqd_ep            ! Dist. sqd from source
  real(kind=PR), intent(inout) :: hep             ! h at e.p.
  real(kind=PR), intent(inout) :: rhoep           ! Density at e.p.
  real(kind=PR), intent(inout) :: rtree(1:NDIM)   ! Position of tree walk
  real(kind=PR), intent(inout) :: htree           ! h from tree walk
  real(kind=PR), intent(inout) :: rhotree         ! rho from tree walk

  integer :: plow                       ! Aux. particle id 1
  integer :: phigh                      ! Aux. particle id 2
  real(kind=PR) :: dr1(1:NDIM)          ! Relative displacement vector 1
  real(kind=PR) :: dr2(1:NDIM)          ! Relative displacement vector 2
  real(kind=PR) :: dr3(1:NDIM)          ! ..
  real(kind=PR) :: drsqdlow             ! Dist. sqd of plow from source
  real(kind=PR) :: drsqdhigh            ! Dist. sqd of phigh from source
  real(kind=PR) :: drsqdtree            ! Dist. sqd of treewalk from source
  real(kind=PR) :: hhigh                ! ..
  real(kind=PR) :: hlow                 ! ..
  real(kind=PR) :: rhigh                ! ..
  real(kind=PR) :: rlow                 ! ..
  real(kind=PR) :: rhohigh              ! ..
  real(kind=PR) :: rholow               ! ..


  plow = pprev
  phigh = p
  drsqdlow = 0.0_PR  !BIG_NUMBER
  if (plow /= -1) call distance2(rsource(1:NDIM),plow,dr1(1:NDIM),drsqdlow)
!  if (treeep) call distance3(rsource(1:NDIM),rtree(1:NDIM),&
!       &dr3(1:NDIM),drsqdtree)
  call distance3(rsource(1:NDIM),rtree(1:NDIM),dr3(1:NDIM),drsqdtree)

! Use linked lists to find particles immediatly before and after 
! evaluation point
  do
     call distance2(rsource(1:NDIM),phigh,dr2(1:NDIM),drsqdhigh)
     if (drsqdhigh >= rsqd_ep) exit
     if (phigh == plast) exit
     drsqdlow = drsqdhigh
     plow     = phigh
     phigh    = nextptcl(phigh)
  end do

! Is tree-walk e.p. the lower bound?
! ----------------------------------------------------------------------------
  if (drsqdtree >= drsqdlow .and. drsqdtree < rsqd_ep .and. treeep) then
     rholow  = rhotree
     hlow    = htree
     rhohigh = sph(phigh)%rho
     hhigh   = sph(phigh)%h
     rlow    = sqrt(drsqdtree)
     rhigh   = sqrt(drsqdhigh)
     plow    = 0

! Is tree-walk e.p. the upper bound?
! ----------------------------------------------------------------------------
  else if (drsqdtree >= rsqd_ep .and. drsqdtree < drsqdhigh .and. &
       &plow /= -1 .and. treeep) then
     rholow  = sph(plow)%rho
     hlow    = sph(plow)%h
     rhohigh = rhotree
     hhigh   = htree
     rlow    = sqrt(drsqdlow)
     rhigh   = sqrt(drsqdtree)

! If two particles form the bounds
! ----------------------------------------------------------------------------
  else if (plow /= phigh .and. plow /= -1) then
     rholow  = sph(plow)%rho
     hlow    = sph(plow)%h
     rhohigh = sph(phigh)%rho
     hhigh   = sph(phigh)%h
     rlow    = sqrt(drsqdlow)
     rhigh   = sqrt(drsqdhigh)

! If only one particle is a bound so we perform a tree-walk
! ----------------------------------------------------------------------------
  else 
     rholow  = sph(phigh)%rho
     hlow    = sph(phigh)%h
     rhohigh = sph(phigh)%rho
     hhigh   = sph(phigh)%h
     rlow    = BIG_NUMBER
     rhigh   = sqrt(drsqdhigh)

  end if
! ----------------------------------------------------------------------------


! If evaluation point is inbetween two particles, interpolate the density at 
! the evaluation point.  
! ----------------------------------------------------------------------------
  if (rlow < rmag_ep .and. rhigh >= rmag_ep .and. &
       & plow /= phigh .and. plow /= -1) then

     ! Check distance between particles and e.p. is not too great.  
     ! If so, perform a tree walk instead to obtain h and rho
     if (abs(rhigh - rlow) < 0.5_PR*f4*KERNRANGE*(hlow + hhigh)) then
        rhoep = rholow + (rmag_ep - rlow)*(rhohigh - rholow)/(rhigh - rlow)
        hep = hlow*(rholow/rhoep)**(1.0_PR/real(NDIM,PR))
     else
        hep = 0.5_PR*(hlow + hhigh)
        rhoep = 0.5_PR*(rholow + rhohigh)
        call HP_rhoh_ep(rep(1:NDIM),hep,rhoep)
        treeep        = .true.
        htree         = hep
        rhotree       = rhoep
        rtree(1:NDIM) = rep(1:NDIM)
     end if

! Else, calculate h and rho via a tree walk.
! ----------------------------------------------------------------------------
  else 
     hep = sph(p)%h
     rhoep = sph(p)%rho
     call HP_rhoh_ep(rep(1:NDIM),hep,rhoep)
     treeep        = .true.
     htree         = hep
     rhotree       = rhoep
     rtree(1:NDIM) = rep(1:NDIM)

  end if
! ----------------------------------------------------------------------------

  return
END SUBROUTINE HP_evaluation_point
