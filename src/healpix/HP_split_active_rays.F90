! HP_SPLIT_ACTIVE_RAYS
! D. A. Hubber - 28/02/2010
! Split rays on current level into 4 child rays on next level
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE HP_split_active_rays(isource,itot,level,nliverays,rsource,&
     &livelist,nextptcl,nextptclaux,prevptcl,prevptclaux,raylist)
  use interface_module, only : distance3,HP_calculate_basis_vector
  use particle_module
  use healpix_types
  use HP_module
  implicit none

  integer, intent(in) :: isource                      ! i.d. of HEALPix source
  integer, intent(inout) :: itot                      ! Total no. of rays
  integer, intent(inout) :: level                     ! current level
  integer, intent(inout) :: nliverays                 ! ..
  integer, intent(inout) :: livelist(1:imax)          ! list of live rays
  integer, intent(inout) :: nextptcl(1:ptot)          ! forward linked list
  integer, intent(inout) :: nextptclaux(1:ptot)       ! forward linked list
  integer, intent(inout) :: prevptcl(1:ptot)          ! backwards linked list
  integer, intent(inout) :: prevptclaux(1:ptot)       ! backwards linked list
  integer(kind=I4B), intent(inout) :: raylist(1:ptot) ! ray id of particle
  real(kind=PR), intent(in) :: rsource(1:NDIM)        ! ..

  logical :: afterrayend            ! After end of previous 'rayend'
  integer(kind=I4B) :: ipix         ! HEALPix id
  integer(kind=I4B) :: ipixaux      ! Aux. pixel id
  integer(kind=I4B) :: ipixnew      ! ..
  integer(kind=I4B) :: nside        ! No. of divisions of level 0 ray
  integer :: i                      ! ..
  integer :: iaux                   ! Aux. integer variable
  integer :: ii                     ! Aux. ray counter
  integer :: iii                    ! Aux. level id
  integer :: inew                   ! ..
  integer :: j                      ! Aux. particle counter
  integer :: nactive                ! No. of active particles in ray
  integer :: nrays                  ! No. of rays on current level
  integer :: p                      ! particle id
  integer :: prayend                ! first particle
  integer :: pprev                  ! id of particle from previous e.p.
  real(kind=PR) :: dr(1:NDIM)       ! relative displacement vector
  real(kind=PR) :: drsqd            ! distance squared
  real(kind=PR) :: rbasis(1:NDIM)   ! Rotated HEALPix basis vector for ray i
#if defined(DEBUG_HP_SPLIT_ACTIVE_RAYS)
  integer :: plist                  ! ..
#endif

  debug2("Split rays into child rays on new level [HP_split_active_rays.F90]")

  nextptclaux(1:ptot) = -1
  prevptclaux(1:ptot) = -1

! Create list of rays that are still alive
! ----------------------------------------------------------------------------
  nliverays = 0
  do i=HPlevel(level)%ifirst,HPlevel(level)%ilast
     if (HPray(i)%done) cycle
     nactive = 0
     p = HPray(i)%rayend
     do 
        !if (sph(p)%ionizedo) nactive = nactive + 1
        nactive = nactive + 1
        if (p == HPray(i)%last) exit
        p = nextptcl(p)
     end do
     if (nactive == 0) then 
        HPray(i)%done = .true.
     else 
        HPray(i)%done = .false.
        nliverays = nliverays + 1
        livelist(nliverays) = i
     end if
#if defined(DEBUG_HP_SPLIT_ACTIVE_RAYS)
     write(6,*) "Live ray? : ",i,nactive,HPray(i)
#endif
  end do
! ----------------------------------------------------------------------------
  
! If there are no active particles in remaining rays, exit main loop
#if defined(DEBUG_HP_SPLIT_ACTIVE_RAYS)
  write(6,*) "nliverays : ",nliverays
#endif
  if (nliverays == 0) return
  
! Otherwise, increase level counter and move onto next HEALPix level
  level = level + 1
  nside = 2**(int(level,I4B))
  nrays = 12*(int(nside)**2)
  
! Advance to next level.  If level does not yet exist, create HEALPix 
! vectors and allocate other required arrays
  if (level > ltot_hp) ltot_hp = level

#if defined(DEBUG_HP_SPLIT_ACTIVE_RAYS)
  write(6,*) "new level : ",ltot_hp
#endif
  
! Stop program for now if allocated memory for rays is too small
  if (itot + 4*nliverays > imax) stop 'imax not big enough!!'
  
  HPlevel(level)%ifirst = HPlevel(level-1)%ilast + 1
  HPlevel(level)%ilast  = HPlevel(level-1)%ilast + 4*nliverays
  itot = itot + 4*nliverays
  
  
! Now create linked lists for new level from old level linked lists
! ============================================================================
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(afterrayend,dr,drsqd,i,iaux,iii,inew,ipix,ipixaux,ipixnew) &
!$OMP PRIVATE(j,p,pprev,prayend,rbasis)
  do ii=1,nliverays
     i = livelist(ii)
     ipix = HPray(i)%ipix

#if defined(DEBUG_HP_SPLIT_ACTIVE_RAYS)
     write(6,*) "Creating child rays for ",i,ipix,nliverays
#endif
     
     ! Loop over all child rays from parent ray ipix
     ! -----------------------------------------------------------------------
     do j=0,3
        inew    = HPlevel(level)%ifirst + 4*(ii - 1) + j
        ipixnew = int(4*ipix + j,I4B)
        iaux    = int(ipixnew)
        
        ! Find the start of the remaining linked list
!        p = HPray(i)%rayend
!        pprev = prevptcl(p)
        p = HPray(i)%first
        prayend = HPray(i)%rayend
        afterrayend = .false.
        
        ! Loop through linked list creating child rays from parent
        ! --------------------------------------------------------------------
        do
           ! Find HEALPix pixel id for current level
           ipixaux = int(raylist(p))
           iii = HP_LEVELS
           do
              if (iii == level) exit
              iii = iii - 1
              ipixaux = ipixaux/4
           end do
           if (p == prayend) afterrayend = .true.
           
           ! If pixel id matches, add particle to linked list
           if (ipixaux == ipixnew) then
              if (HPray(inew)%first == -1) then
                 prevptclaux(p) = -1
                 HPray(inew)%first = p
                 !HPray(inew)%done = .false.
              else
                 prevptclaux(p) = HPray(inew)%last
                 nextptclaux(HPray(inew)%last) = p
              end if
              HPray(inew)%last = p
              if (afterrayend .and. HPray(inew)%rayend == -1) then
                 HPray(inew)%done = .false.
                 HPray(inew)%rayend = p
              end if
              !if (p == HPray(inew)%last .and. HPray(inew)%rayend == -1) &
              !     &HPray(inew)%rayend = HPray(inew)%last
           end if
           
           ! Move along to next particle unless we have reached the end
           if (p == HPray(i)%last) exit
           p = nextptcl(p)
        end do
        ! --------------------------------------------------------------------
        
#if defined(DEBUG_HP_SPLIT_ACTIVE_RAYS)
        write(6,*) "Creating child ray ",i,ipix,j,level,inew,ipixnew,&
             &HPray(inew)%first,HPray(inew)%last,HPray(inew)%rayend,&
             &HPray(inew)%done
#endif
        
        if (HPray(inew)%done) cycle
        call distance3(rsource(1:NDIM),HPray(i)%rep(1:NDIM),dr(1:NDIM),drsqd)
        call HP_calculate_basis_vector(level,ipixnew,&
             &rbasis(1:NDIM),isource)
        
        HPray(inew)%rep(1:NDIM) = rsource(1:NDIM) + sqrt(drsqd)*rbasis(1:NDIM)
        HPray(inew)%ipix        = ipixnew
        HPray(inew)%rhoep       = HPray(i)%rhoep
        HPray(inew)%hep         = HPray(i)%hep
#if defined(IONIZING_UV_RADIATION)
        HPray(inew)%integral    = HPray(i)%integral
        HPray(inew)%UVdone      = HPray(i)%UVdone
#endif
#if defined(STELLAR_WIND)
        HPray(inew)%winddone    = HPray(i)%winddone
#endif
     end do
     ! -----------------------------------------------------------------------
     
  end do
!$OMP END PARALLEL DO
! ============================================================================
  

  nextptcl(1:ptot) = nextptclaux(1:ptot)
  prevptcl(1:ptot) = prevptclaux(1:ptot)

#if defined(DEBUG_HP_SPLIT_ACTIVE_RAYS)
  plist = 0
  do i=HPlevel(level)%ifirst,HPlevel(level)%ilast
     p = HPray(i)%first
     do 
        if (p == -1) exit
        plist = plist + 1
        p = nextptcl(p)
     end do
     if ((.not. HPray(i)%done) .and. HPray(i)%first == -1) then
        write(6,*) "Problem with something in split-active-rays : ",&
             &i,HPray(i)%done,HPray(i)%first,HPray(i)%rayend,HPray(i)%last
        stop
     end if
  end do
  if (plist /= ptot) then
     write(6,*) "Problem with linked lists!! : ",level,plist,ptot
     stop
  else
     write(6,*) "Linked lists all fine!! : ",level,plist,ptot
  end if
#endif
  
  return
END SUBROUTINE HP_split_active_rays
