! HP_IONIZING_RADIATION
! T. Bisbas & D. A. Hubber - 15/9/2008
! Compute ionization integral along ray.  If ionization-front is found, 
! perform a binary-chop to find IF to more accuracy before computing new 
! ionized temperature of particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE HP_propagate_UV_radiation(i,isource,level,p,pprev,fh,hep,oldfh,&
     &oldstep,ray_integral,rep,rep_prev,rhoep,rsource,rsqd_ep,rtree,htree,&
     &rhotree,treeep,nextptcl,prevptcl)
  use interface_module, only : distance2,distance3,IFbinarychop
  use particle_module
  use hydro_module
  use HP_module
  implicit none

  logical, intent(inout) :: treeep              ! ??
  integer, intent(in) :: i                      ! HP ray id
  integer, intent(in) :: isource                ! HEALPix source id
  integer, intent(in) :: level                  ! Current HEALPix level
  integer, intent(inout) :: p                   ! Current particle id
  integer, intent(inout) :: pprev               ! Previous particle id
  real(kind=PR), intent(in) :: rep_prev(1:NDIM) ! Position of previous part.
  real(kind=PR), intent(in) :: rsource(1:NDIM)  ! Source position
  real(kind=PR), intent(inout) :: fh            ! Integration step
  real(kind=PR), intent(inout) :: hep           ! h of evaluation point
  real(kind=PR), intent(inout) :: oldfh         ! Previous integration step
  real(kind=PR), intent(inout) :: oldstep       ! ??
  real(kind=PR), intent(inout) :: ray_integral  ! Ionisation integral
  real(kind=PR), intent(inout) :: rep(1:NDIM)   ! Position of evaluation point
  real(kind=PR), intent(inout) :: rhoep         ! Density at evaluation point
  real(kind=PR), intent(inout) :: rsqd_ep       ! Distance squared from e.p.
  real(kind=PR), intent(inout) :: rtree(1:NDIM) ! Position of tree walk
  real(kind=PR), intent(inout) :: htree         ! h from tree walk
  real(kind=PR), intent(inout) :: rhotree       ! rho from tree walk
  integer, intent(in) :: nextptcl(1:ptot)       ! forward linked list
  integer, intent(in) :: prevptcl(1:ptot)       ! backwards linked list

  real(kind=PR) :: dis_s                        ! Norm. distance from IF
  real(kind=PR) :: dr(1:NDIM)                   ! relative displacement vector
  real(kind=PR) :: drsqd                        ! distance squared
  real(kind=PR) :: integrationstep              ! Int. step var.
  real(kind=PR) :: intstep                      ! Int. step var.
  real(kind=PR) :: raux                         ! Aux. real variable
  real(kind=PR) :: rmag_ep                      ! Distance of e.p.
     
! Sum terms for ionization integral
#if defined(TRAPEZOIDAL_RULE)
  intstep         = rsqd_ep*(rhoep**2)
  integrationstep = 0.5_PR*oldfh*(oldstep + intstep)
#elif defined(SIMPSONS_RULE)
  intstep         = rsqd_ep*(rhoep**2)
  integrationstep = 0.5_PR*oldfh*(oldstep + intstep)
#endif
  oldfh           = fh
  ray_integral    = ray_integral + integrationstep

#if defined(DEBUG_HP_WALK_RAY)
  write(6,*) "ray_integral : ",ray_integral,&
       &"    intmax : ",HPsource(isource)%intmax
  write(6,*) "intstep : ",intstep,"  integrationstep : ",integrationstep
#endif


! Check if we have reached ionization front or not
! =============================================================================
  if (ray_integral < HPsource(isource)%intmax) then

     ! Follow linked-list to find all particles behind eval. point
     do 
        call distance2(rsource(1:NDIM),p,dr(1:NDIM),drsqd)
        if (drsqd > rsqd_ep) exit
        if (sph(p)%newtemp == 0) sph(p)%tempmin = 0.0_PR
        sph(p)%newtemp = sph(p)%newtemp + 1
        sph(p)%tempaux = Tion
#if defined(DEBUG_HP_WALK_ALL_RAYS)
        whichHPlevel(p) = level
#endif
        if (p == HPray(i)%last) then
           HPray(i)%UVdone = .true.
           exit
        end if
        pprev = p
        p = nextptcl(p)
     end do


! If we have reached the IF, locate the exact position
! using a binary chop, and then recalculate rho and h.
! =============================================================================
  else

     ray_integral = ray_integral - integrationstep
     call IFbinarychop(HPsource(isource)%intmax,rep(1:NDIM),rep_prev(1:NDIM),&
          &hep,sph(p)%h,rhoep,rsource(1:NDIM),ray_integral,&
          &oldstep,pprev,p,nextptcl,HPray(i)%last,rtree(1:NDIM),htree,&
          &rhotree,treeep)
     call distance3(rsource(1:NDIM),rep(1:NDIM),dr(1:NDIM),rsqd_ep)
     rmag_ep = sqrt(rsqd_ep)   

#if defined(DEBUG_HP_IF)
     write(6,*) "Found IF : ",isource,i,ray_integral,HPsource(isource)%intmax
     write(6,*) "rmag_ep  : ",rmag_ep,rep(1:NDIM)
     write(6,*) "pprev    : ",pprev,HPray(i)%first,HPray(i)%rayend
#endif


     ! If computing PDR chemistry, call subroutine to compute all quantities
     ! from table and set temperatures of particles.
     ! ------------------------------------------------------------------------
#if defined(PDR_CHEMISTRY)

     call compute_pdr_properties(isource,i,p,nextptcl,HPray(i)%last,&
        &rmag_ep,hep,rhoep)
     HPray(i)%UVdone = .true.


     ! Else, use regular temperature smoothin around ionisation front to 
     ! prevent gap formation problems.
     ! ------------------------------------------------------------------------
#else


     ! First check all particles behind e.p.
     ! -----------------------------------------------------------------------
     do 
        if (pprev == -1) exit
        call distance2(rsource(1:NDIM),pprev,dr(1:NDIM),drsqd)
        dis_s = sqrt(drsqd) - rmag_ep
        raux = -0.5_PR*(1.0_PR/(f3*hep))*dis_s + 0.5_PR
        if (dis_s < -f3*hep) raux = 1.0_PR
        if (dis_s > f3*hep)  raux = 0.0_PR
        if (sph(pprev)%newtemp == 0) sph(pprev)%tempmin = 0.0_PR
        sph(pprev)%newtemp = sph(pprev)%newtemp + 1
        sph(pprev)%tempaux = raux*Tion + (1.0_PR - raux)*Tneut
#if defined(DEBUG_HP_IF)
        write(6,*) "dis_s : ",dis_s,f3*hep,raux,newtemp(pprev)
        write(6,*) "temp : ",temp(pprev),temp_aux(pprev)
        write(6,*) "pprev : ",prevptcl(pprev)
#endif
        !if (raux < 0 .OR. raux > 1) stop 'insane raux'
        if (dis_s < -f3*hep) exit
        if (pprev == HPray(i)%first) exit
        pprev = prevptcl(pprev)
     end do

#if defined(DEBUG_HP_IF)
        write(6,*) "p : ",p
#endif
     
     ! Follow linked-list to find all particles in front of e.p.
     ! -----------------------------------------------------------------------
     do 
        call distance2(rsource(1:NDIM),p,dr(1:NDIM),drsqd)
        dis_s = sqrt(drsqd) - rmag_ep
        raux  = -0.5_PR*(1.0_PR/(f3*hep))*dis_s + 0.5_PR
        if (dis_s < -f3*hep) raux = 1.0_PR
        if (dis_s > f3*hep)  raux = 0.0_PR
        if (sph(p)%newtemp == 0) sph(p)%tempmin = 0.0_PR
        sph(p)%newtemp  = sph(p)%newtemp + 1
        sph(p)%tempaux = raux*Tion + (1.0_PR - raux)*Tneut
#if defined(DEBUG_HP_WALK_ALL_RAYS)
        whichHPlevel(p) = level
#endif
#if defined(DEBUG_HP_IF)
        write(6,*) "dis_s : ",dis_s,f3*hep,raux,sph(p)%newtemp
        write(6,*) "temp : ",temp(p),sph(p)%tempaux
        write(6,*) "p : ",nextptcl(p)
#endif
        
        !if (raux < 0 .OR. raux > 1) stop 'insane raux'
        if (dis_s > f3*hep) exit
        if (p == HPray(i)%last) exit
        p = nextptcl(p)
     end do
     
     HPray(i)%UVdone = .true.
     ! -----------------------------------------------------------------------
#endif

  end if
  ! --------------------------------------------------------------------------
 
  return
END SUBROUTINE HP_propagate_UV_radiation




! ============================================================================
! IFLOCATIONBINARYCHOP
! Locates the position of the ionization front using a binary chop method.
! ============================================================================
SUBROUTINE IFbinarychop(rayint_max,rep,rep_prev,hep,hprev,rhoep,&
     &rc,ray_integral,oldstep,pprev,p,nextptcl,plast,&
     &rtree,htree,rhotree,treeep)
  use interface_module, only : distance3,HP_evaluation_point
  use definitions
  use particle_module
  use hydro_module
  use HP_module
  implicit none

  integer, intent(in) :: p                         ! Particle in linked list
  integer, intent(in) :: plast                     ! Last particle in list
  integer, intent(in) :: pprev                     ! Particle before e.p.
  integer, intent(in) :: nextptcl(1:ptot)          ! Linked list
  real(kind=PR), intent(in)    :: rayint_max       ! Max. integral for source
  real(kind=PR), intent(inout) :: rep(1:NDIM)      ! Position of I.F.
  real(kind=PR), intent(in)    :: rep_prev(1:NDIM) ! Previous e.p.
  real(kind=PR), intent(inout) :: hep              ! h at e.p.
  real(kind=PR), intent(in)    :: hprev            ! h of prev. particle
  real(kind=PR), intent(inout) :: rhoep            ! Density at e.p.
  real(kind=PR), intent(in)    :: rc(1:NDIM)       ! Ionizing source position
  real(kind=PR), intent(inout) :: ray_integral     ! I.F. integral
  real(kind=PR), intent(in)    :: oldstep          ! Previous integration step
  real(kind=PR), intent(inout) :: rtree(1:NDIM)    ! ..
  real(kind=PR), intent(inout) :: htree            ! ..
  real(kind=PR), intent(inout) :: rhotree          ! ..
  logical, intent(inout) :: treeep                 ! ..

  integer :: itcounter               ! Iteration counter
  real(kind=PR) :: dr(1:NDIM)        ! Relative displacement vector
  real(kind=PR) :: drmag             ! Distance
  real(kind=PR) :: intbefore         ! ray_integral before this step
  real(kind=PR) :: intstep           ! Current stepsize
  real(kind=PR) :: integrationstep   ! Average integration step (Simpson's??)
  real(kind=PR) :: rmag_prev         ! Distance of previous step from source
  real(kind=PR) :: rsqd_ep           ! Distance squared from source
  real(kind=PR) :: rhigh(1:NDIM)     ! Upper value of bisection search
  real(kind=PR) :: rlow(1:NDIM)      ! Lower value of bisection search

  itcounter     = 0
  intbefore     = ray_integral
  rlow(1:NDIM)  = rep_prev(1:NDIM)
  rhigh(1:NDIM) = rep(1:NDIM)
  rep(1:NDIM)   = 0.5*(rlow(1:NDIM) + rhigh(1:NDIM))
  call distance3(rc(1:NDIM),rep_prev(1:NDIM),dr(1:NDIM),rsqd_ep)
  rmag_prev = sqrt(rsqd_ep)

#if defined(DEBUG_HII)
  write(6,*) "In IFbinarychop",itcounter,intbefore,rlow,rhigh,rep
#endif


! Main binary chop loop
! ----------------------------------------------------------------------------
  do
     itcounter = itcounter + 1
     call distance3(rc(1:NDIM),rep(1:NDIM),dr(1:NDIM),rsqd_ep)
     drmag = sqrt(rsqd_ep) + SMALL_NUMBER

     call HP_evaluation_point(rep(1:NDIM),drmag,rsqd_ep,rc(1:NDIM),pprev,&
          &p,plast,nextptcl,rhoep,hep,rtree(1:NDIM),htree,rhotree,treeep)

     intstep = (rhoep**2)*rsqd_ep
     integrationstep = 0.5_PR*(oldstep + intstep)*(drmag - rmag_prev)
     ray_integral = intbefore + integrationstep

     ! Adjust binary limits for next iteration
     if (ray_integral >= rayint_max) then
        rhigh(1:NDIM) = rep(1:NDIM)
     else
        rlow(1:NDIM) = rep(1:NDIM)
     end if
     rep(1:NDIM) = 0.5_PR*(rlow(1:NDIM) + rhigh(1:NDIM))
     dr(1:NDIM) = rhigh(1:NDIM) - rlow(1:NDIM)

     ! Terminate iteration if distance between limits is small enough, 
     ! or integral is close to maximum value.
     if (dot_product(dr(1:NDIM),dr(1:NDIM)) <= 0.01_PR*hprev*hprev) return
     if (abs((ray_integral - rayint_max)/rayint_max) < 0.01_PR) return

     ! If not converged after 10 iterations, just use the current value
     if (itcounter >= 10) then
        write(6,*) "Not converged after 10 iterations.  Returning"
        rhoep = 0.5_PR*(sph(pprev)%rho + sph(p)%rho)
        hep   = 0.5_PR*(sph(pprev)%h + sph(p)%h)
        return
     end if

  end do
! ----------------------------------------------------------------------------

  return
END SUBROUTINE IFbinarychop
