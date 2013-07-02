! HP_WALK_RAY.F90
! T. Bisbas & D. A. Hubber - 15/9/2008
! Walk healpix ray i from beginning until the ray has 
! i) reached the splitting criterion, or 
! ii) been terminated due to satisfying some criterion.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE HP_walk_ray(i,isource,level,res,nep,rmag_ep,&
     &rsource,nextptcl,prevptcl)
  use interface_module, only : distance2,distance3,HP_rhoh_ep,&
       &HP_calculate_basis_vector,HP_evaluation_point,&
       &HP_propagate_UV_radiation
  use particle_module
  use hydro_module
  use HP_module
  implicit none
  
  integer, intent(in) :: i                      ! Ray id
  integer, intent(in) :: isource                ! Source id
  integer, intent(in) :: level                  ! current level
  integer, intent(in) :: nextptcl(1:ptot)       ! forward linked list
  integer, intent(in) :: prevptcl(1:ptot)       ! backwards linked list
  integer, intent(inout) :: nep                 ! No. of eval points
  real(kind=PR), intent(out) :: res             ! ray level resolution
  real(kind=PR), intent(out) :: rmag_ep         ! distance of e.p. from source
  real(kind=PR), intent(in) :: rsource(1:NDIM)  ! Position of HEALPix source

  logical :: treeep                    ! Is last e.p. a tree-walk?
  integer :: clevel                    ! Optimal level from ray-split criterion
  integer :: iaux                      ! Aux. integer variable
  integer :: nnotdone                  ! No. of nodes not finished
  integer :: p                         ! particle id
  integer :: pprev                     ! id of particle from previous e.p.
  real(kind=PR) :: dr(1:NDIM)          ! relative displacement vector
  real(kind=PR) :: drsqd               ! distance squared
  real(kind=PR) :: fh                  ! Next integration step-size
  real(kind=PR) :: hep                 ! Smoothing length of the e.p.
  real(kind=PR) :: hp                  ! Smoothing length of particle p
  real(kind=PR) :: htree               ! Smoothing length computed by tree 
  real(kind=PR) :: intfunc             ! integral contribution of current step
  real(kind=PR) :: oldfh               ! previous value of fh
  real(kind=PR) :: func1               ! Function integrand at start of step
  real(kind=PR) :: func2               ! Function integrand at mid. of step
  real(kind=PR) :: func3               ! Function integrand at end of step
  real(kind=PR) :: ray_integral        ! ionization integral
  real(kind=PR) :: raux                ! Aux. real variable
  real(kind=PR) :: rbasis(1:NDIM)      ! Rotated HEALPix basis vector for ray i
  real(kind=PR) :: rep(1:NDIM)         ! position of evaluation point
  real(kind=PR) :: rep_prev(1:NDIM)    ! position of previous eval. point
  real(kind=PR) :: rhoep               ! density at evaluation point
  real(kind=PR) :: rhotree             ! density computed by tree walk
  real(kind=PR) :: rtree(1:NDIM)       ! Position of last tree-walk
  real(kind=PR) :: rsqd_ep             ! distance squared of e.p. from source
#if defined(DEBUG_HP_WALK_RAY)
  integer :: j                         ! Aux. counter
#endif

! Initialise all important variables for this ray
  clevel   = level
  nep      = 0
  nnotdone = 0
  p        = HPray(i)%rayend
  pprev    = prevptcl(p)
  rep(1:NDIM) = HPray(i)%rep(1:NDIM)
  call distance2(rsource(1:NDIM),p,dr(1:NDIM),rsqd_ep)
  raux     = sqrt(rsqd_ep)
  call distance3(rsource(1:NDIM),rep(1:NDIM),dr(1:NDIM),rsqd_ep)
  rmag_ep  = sqrt(rsqd_ep)
  raux     = raux - rmag_ep
  res      = 0.0_PR
#if defined(IONIZING_UV_RADIATION)
  ray_integral = HPray(i)%integral
#endif

! Set tree walk variables
  treeep        = .false.
  htree         = 0.0_PR
  rhotree       = 0.0_PR
  rtree(1:NDIM) = 0.0_PR
  
  if (HPray(i)%rhoep < SMALL_NUMBER) then
     rhoep = sph(p)%rho
     hep = sph(p)%h
     call HP_evaluation_point(rep(1:NDIM),rmag_ep,rsqd_ep,rsource(1:NDIM),&
          &pprev,p,HPray(i)%last,nextptcl,rhoep,hep,&
          &rtree(1:NDIM),htree,rhotree,treeep)
  else
     rhoep = HPray(i)%rhoep
     hep = HPray(i)%hep
  end if
  
  raux    = min(hep,max(sph(p)%h,raux))
  fh      = f1*raux
#if defined(TRAPEZOIDAL_RULE)
  intfunc = (rhoep**2)*rsqd_ep
  oldfh   = fh
  func1   = intfunc
#elif defined(SIMPSONS_RULE)
  intfunc = (rhoep**2)*rsqd_ep
  oldfh   = fh
  func1   = intfunc
#endif
  
! Calculate required level for current evaluation point.  If more resolution  
! is needed, the current level is skipped completely for this ray.
  if (rsqd_ep > 0.0_PR) then
     clevel = int(log10(rmag_ep/(f2*raux))*INVLOG10TWO + 0.5_PR)
     if (clevel < level) clevel = level
     if (clevel > lmax_hp) clevel = lmax_hp
     res = max(res,(log10(rmag_ep/(f2*hep))*INVLOG10TWO + 0.5_PR))
  end if
  
! Calculate rotated HEALPix basis vector
  call HP_calculate_basis_vector(level,HPray(i)%ipix,rbasis(1:NDIM),isource)

#if defined(DEBUG_HP_WALK_RAY)
  write(6,*) "==========================================================="
  write(6,*) "Walking ray : ",i,level,"    isource : ",isource
  write(6,*) "p : ",p,"     pprev : ",pprev,"    next : ",nextptcl(p)
  write(6,*) "rsource : ",rsource(1:NDIM)
  write(6,*) "rbasis  : ",rbasis(1:NDIM)
  write(6,*) "rp      : ",sph(p)%r(1:NDIM) - rsource(1:NDIM)
  write(6,*) "rmag_ep : ",rmag_ep,"        dist. : ",raux
  write(6,*) "rhoep   : ",rhoep,"        hep   : ",hep
  write(6,*) "clevel  : ",clevel,"      res : ",res
  do j=1,ptot
     if (nextptcl(j) /= -1) write(6,*) "final : ",nextptcl(j),j
     if (nextptcl(j) /= -1) exit
  end do
#endif
  
  
! Follow each ray until we reach the first manditory split, 
! or are forced to split early to a high density region.
! ============================================================================
  do while (clevel <= level)
     nep = nep + 1
     rep_prev(1:NDIM) = rep(1:NDIM)
     rep(1:NDIM) = rep(1:NDIM) + fh*real(rbasis(1:NDIM),PR)
     call distance3(rsource(1:NDIM),rep(1:NDIM),dr(1:NDIM),rsqd_ep)
     rmag_ep = sqrt(rsqd_ep)

     call HP_evaluation_point(rep(1:NDIM),rmag_ep,rsqd_ep,&
          &rsource(1:NDIM),pprev,p,HPray(i)%last,nextptcl,rhoep,hep,&
          &rtree(1:NDIM),htree,rhotree,treeep)

     ! Use smoothing length of nearby particles (via the linked-list)
     ! to determine the optimal position of next eval. point
     iaux = p
     do
        hp = sph(iaux)%h
        call distance2(rsource(1:NDIM),iaux,dr(1:NDIM),drsqd)
        raux = sqrt(drsqd) - rmag_ep
        if (raux > f1*hp .or. raux > f1*hep) exit
        if (iaux == HPray(i)%last) exit
        iaux = nextptcl(iaux)
     end do
     raux = min(hep,max(hp,raux))
     fh = f1*raux

#if defined(DEBUG_HP_WALK_RAY)
     write(6,*) "-----------------------------------------------------------"
     write(6,*) "nep   : ",nep,"   p : ",p,"   pprev : ",pprev,"  iaux : ",iaux
     write(6,*) "rhoep : ",rhoep,"   hep : ",hep,"   htree : ",htree
     write(6,*) "rmag_ep : ",rmag_ep,"   raux : ",raux,hep,hp,sqrt(drsqd)-rmag_ep
#endif
     
     ! Calculate required level of next eval. point
     if (rsqd_ep > 0.0_PR) then
        clevel = int(log10(rmag_ep/(f2*raux))*INVLOG10TWO + 0.5_PR)
        if (clevel < level) clevel = level
        if (clevel > lmax_hp) clevel = lmax_hp
        res = max(res,(log10(rmag_ep/(f2*hep))*INVLOG10TWO + 0.5_PR))
     end if

     ! Ionizing radiation
     ! -----------------------------------------------------------------------
#if defined(IONIZING_UV_RADIATION)
     HP_ionize = .true.
     if (HP_ionize) call HP_propagate_UV_radiation(i,isource,level,p,pprev,&
          &fh,hep,oldfh,func1,ray_integral,rep(1:NDIM),rep_prev(1:NDIM),&
          &rhoep,rsource(1:NDIM),rsqd_ep,rtree(1:NDIM),htree,rhotree,&
          &treeep,nextptcl(1:ptot),prevptcl(1:ptot))
#endif

     func1 = intfunc
     if (HPray(i)%UVdone) exit

  end do
! ============================================================================


! Set final properties of ray before terminating or splitting
  HPray(i)%rayend            = p
  HPray(i)%rhoep             = rhoep
  HPray(i)%hep               = hep
  HPray(i)%rep(1:NDIM)       = rep(1:NDIM)
#if defined(IONIZING_UV_RADIATION)
  HPray(i)%integral          = ray_integral
  call distance3(rsource(1:NDIM),rep(1:NDIM),dr(1:NDIM),drsqd)
  rmag_ep = sqrt(drsqd)
#endif
#if defined(DEBUG_HP_WALK_RAY)
  write(6,*) "==========================================================="
  write(6,*) "Walked ray : ",i,level,"    isource : ",isource
  write(6,*) "integral : ",HPray(i)%integral,"   maxres : ",res
  write(6,*) "rep : ",HPray(i)%rep(1:NDIM)
#endif
#if !defined(STELLAR_WIND)
  HPray(i)%done = HPray(i)%UVdone
#endif
  if (HPray(i)%UVdone) return
  if (clevel > level) nnotdone = nnotdone + 1
  

  return
END SUBROUTINE HP_walk_ray
