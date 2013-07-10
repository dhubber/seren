! HP_INITIALIZE_SOURCE
! T. Bisbas & D. A. Hubber - 29/10/2009
! Initialize all variables for HEALPix source isource.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE HP_initialize_source(isource,distids,nextptcl,&
     &nextptclaux,prevptcl,prevptclaux,raylist,distsqd)
  use interface_module, only : ang2pix_nest,distance2,&
       &heapsort_real,HP_inverse_positions,vec2ang
  use definitions
  use healpix_types
  use sink_module
  use particle_module
  use time_module
  use HP_module
  implicit none

  integer, intent(in) :: isource                      ! i.d. of HEALPix source
  integer, intent(inout) :: distids(1:ptot)           ! list of particles ids
  integer, intent(inout) :: nextptcl(1:ptot)          ! forward linked list
  integer, intent(inout) :: nextptclaux(1:ptot)       ! forward linked list
  integer, intent(inout) :: prevptcl(1:ptot)          ! backwards linked list
  integer, intent(inout) :: prevptclaux(1:ptot)       ! backwards linked list
  integer(kind=I4B), intent(inout) :: raylist(1:ptot) ! ray id of particle
  real(kind=PR), intent(inout) :: distsqd(1:ptot)     ! distance squared

  integer :: i                      ! ray counter
  integer :: iii                    ! Aux. level id
  integer :: j                      ! Aux. particle counter
  integer :: level                  ! current level
  integer :: nactive                ! No. of active particles in ray
  integer :: nrays                  ! No. of rays on current level
  integer :: p                      ! Particle counter
  integer :: pprev                  ! id of particle from previous e.p.
  integer(kind=I4B) :: ipix         ! HEALPix id
  integer(kind=I4B) :: ipixaux      ! Aux. pixel id
  integer(kind=I4B) :: nside        ! No. of divisions of level 0 ray
  real(kind=PR) :: dr(1:NDIM)       ! relative displacement vector
  real(kind=PR) :: drsqd            ! distance squared
  real(kind=DP) :: phi              ! angle
  real(kind=PR) :: raux             ! Aux. real variable
  real(kind=DP) :: rp(1:NDIM)       ! Position of particle
  real(kind=DP) :: rsource(1:NDIM)  ! Source position
  real(kind=DP) :: theta            ! rotation angle

  debug2("Initialize all variables for healpix-source [HP_initialize_source.F90]")

! Initialise variables and arrays
  level = 0
  nside = 2**(int(level,I4B))
  nrays = int(12*(nside**2))
  nextptcl(1:ptot)    = -1
  prevptcl(1:ptot)    = -1
  nextptclaux(1:ptot) = -1
  prevptclaux(1:ptot) = -1
  rsource(1:NDIM) = HPsource(isource)%r(1:NDIM)

! Zero HEALPix source variables
#if defined(IONIZING_UV_RADIATION)
  HPsource(isource)%HIIangle = 0.0_PR
  HPsource(isource)%HIIvolume = 0.0_PR
#endif
#if defined(STELLAR_WIND)
  HPsource(isource)%windangle = 0.0_PR
  HPsource(isource)%windvolume = 0.0_PR
#endif

! Determine (random) orientation of HEALPix rays in the current timestep 
! using random numbers.`
  call random_number(raux)
  HPsource(isource)%Aangle = real(raux*TWOPI,DP)
  call random_number(raux)
  HPsource(isource)%Bangle = real(raux*PI,DP)
  call random_number(raux)
  HPsource(isource)%Cangle = real(raux*TWOPI,DP)
  
! Now calculate rotation matrix
  HPsource(isource)%arot(1,1) = real(cos(HPsource(isource)%Aangle) * &
       &cos(HPsource(isource)%Cangle) - sin(HPsource(isource)%Aangle) * &
       &cos(HPsource(isource)%Bangle)*sin(HPsource(isource)%Cangle),DP)
  HPsource(isource)%arot(2,1) = real(-sin(HPsource(isource)%Aangle) * &
       &cos(HPsource(isource)%Cangle) - sin(HPsource(isource)%Cangle) * &
       &cos(HPsource(isource)%Aangle)*cos(HPsource(isource)%Bangle),DP)
  HPsource(isource)%arot(3,1) = real(sin(HPsource(isource)%Bangle) * &
       &sin(HPsource(isource)%Cangle),DP)
  HPsource(isource)%arot(1,2) = real(cos(HPsource(isource)%Aangle) * &
       &sin(HPsource(isource)%Cangle) + sin(HPsource(isource)%Aangle) * &
       &cos(HPsource(isource)%Bangle)*cos(HPsource(isource)%Cangle),DP)
  HPsource(isource)%arot(2,2) = real(-sin(HPsource(isource)%Aangle) * &
       &sin(HPsource(isource)%Cangle) + cos(HPsource(isource)%Aangle) * &
       &cos(HPsource(isource)%Bangle)*cos(HPsource(isource)%Cangle),DP)
  HPsource(isource)%arot(3,2) = real(-sin(HPsource(isource)%Bangle) * &
       &cos(HPsource(isource)%Cangle),DP)
  HPsource(isource)%arot(1,3) = real(sin(HPsource(isource)%Bangle) * &
       &sin(HPsource(isource)%Aangle),DP)
  HPsource(isource)%arot(2,3) = real(cos(HPsource(isource)%Aangle) * &
       &sin(HPsource(isource)%Bangle),DP)
  HPsource(isource)%arot(3,3) = real(cos(HPsource(isource)%Bangle),DP)

! Determine rays for each particle on finest HEALPix level
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dr,drsqd,ipix,p,phi,rp,theta)
  do p=1,ptot
     call distance2(real(rsource(1:NDIM),PR),p,dr(1:NDIM),drsqd)
     rp(1:NDIM) = real(dr(1:NDIM),DP)
     call HP_inverse_positions(rp(1:NDIM),isource)
     call vec2ang(rp(1:NDIM),theta,phi)
     call ang2pix_nest(ns_max,theta,phi,ipix,x2pix,y2pix)
     raylist(p) = ipix
  end do
!$OMP END PARALLEL DO

! Calculate distances in order to update particle-source distance ordered list
  do j=1,ptot
     p = j !HPsource(isource)%distorder(j)
     call distance2(real(rsource(1:NDIM),PR),p,dr(1:NDIM),drsqd)
     distsqd(j) = drsqd
     distids(j) = p
  end do

! Sort lists in order of increasing distance from source
#if defined(DEBUG_HP_WALK_ALL_RAYS)
  write(6,*) "Sorting particle lists..."
#endif
  call heapsort_real(ptot,distsqd(1:ptot),distids(1:ptot))
  do p=1,ptot
     HPsource(isource)%distorder(p) = distids(p)
  end do

  level = 0
  nside = 2**(int(level,I4B))
  nrays = int(12*(nside**2))
  HPlevel(level)%ifirst = 1
  HPlevel(level)%ilast = nrays
  do i=1,nrays
     HPray(i)%ipix = int(i - 1,I4B)
     HPray(i)%rep(1:NDIM) = rsource(1:NDIM)
  end do


! Now loop through particles in increasing distance from ionizing source 
! and construct the linked lists for each ray on the first level
! ----------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(iii,ipix,ipixaux,j,nactive,p,pprev)
  do i=HPlevel(level)%ifirst,HPlevel(level)%ilast
     ipix = HPray(i)%ipix
     nactive = 0
     do j=1,ptot
        p = HPsource(isource)%distorder(j)
        ipixaux = int(raylist(p))
        iii = HP_LEVELS
        do
           if (iii == level) exit
           iii = iii - 1
           ipixaux = ipixaux / 4
        end do
        if (ipix .ne. ipixaux) cycle
        nactive = nactive + 1
        !if (sph(p)%ionizedo) nactive = nactive + 1
        if (HPray(i)%first == -1) then
           HPray(i)%done  = .false.
           HPray(i)%first = p
           HPray(i)%last  = p
        else
           pprev           = HPray(i)%last
           prevptcl(p)     = pprev
           nextptcl(pprev) = p
           HPray(i)%last   = p
        end if
     end do
     if (nactive == 0) HPray(i)%done = .true.
     HPray(i)%rayend = HPray(i)%first
#if defined(DEBUG_HP_WALK_ALL_RAYS)
     write(6,*) "i : ",i,"   ipix : ",ipix,"    nactive : ",nactive
     write(6,*) "pfirst : ",HPray(i)%first,"  pnext : ",&
          &nextptcl(HPray(i)%first),"  plast : ",HPray(i)%last
#endif
  end do
!$OMP END PARALLEL DO

  return
END SUBROUTINE HP_initialize_source
