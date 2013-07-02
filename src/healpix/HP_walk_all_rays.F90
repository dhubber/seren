! HP_WALK_ALL_RAYS.F90
! T. Bisbas & D. A. Hubber - 15/9/2008
! Walk all HEALPix rays along HEALPix data structure.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE HP_walk_all_rays(isource)
  use interface_module, only : HP_initialize_source,&
       &HP_walk_ray,HP_split_active_rays,HP_stellar_feedback
  use healpix_types
  use HP_module
  use particle_module
  use hydro_module
  use time_module
  use tree_module
  use type_module
  implicit none

  integer, intent(in) :: isource        ! id of HP source

  integer :: i                          ! ray counter
  integer :: idtot                      ! Total no. of evalulation points
  integer :: itot                       ! Total no. of rays
  integer :: level                      ! current level
  integer :: nep                        ! No, of evaluation points
  integer :: nliverays                  ! No. of live rays
  integer :: nnotdone                   ! No. of nodes not finished
  real(kind=PR) :: drmag                ! Distance from ray
  real(kind=PR) :: maxres               ! ray level resolution  
  real(kind=PR) :: rsource(1:NDIM)      ! position of ionizing source
  real(kind=PR) :: resolution           ! ray level resolution
  real(kind=PR) :: HIIangle             ! Angle subtended by enclosed HII rays
  real(kind=PR) :: HIIvolume            ! Volume enclosed by HII rays
  real(kind=PR) :: windvolume           ! Angle subtended by enclosed wind rays
  real(kind=PR) :: windangle            ! Volume enclosed by wind rays
  integer, allocatable :: distids(:)           ! List of particles (unsorted)
  integer, allocatable :: livelist(:)          ! List of live rays
  integer, allocatable :: nextptcl(:)          ! forward linked list
  integer, allocatable :: nextptclaux(:)       ! forward linked list
  integer, allocatable :: prevptcl(:)          ! backwards linked list
  integer, allocatable :: prevptclaux(:)       ! backwards linked list
  integer(kind=I4B), allocatable :: raylist(:) ! ray id of particle
  real(kind=PR), allocatable :: distsqd(:)     ! Distance from source squared

  debug2("Propagating ionizing radiation [HP_walk_all_rays.F90]")
  debug_timing("HP_SETUP")

! If source contains no UV flux (e.g. for sinks having too little mass), 
! then we simply return here and do not compute anything else
#if defined(IONIZING_UV_RADIATION)
!  if (HPsource(isource)%intmax < SMALL_NUMBER) return
#endif

! Make local copy of position of ionizing source
  rsource(1:NDIM) = HPsource(isource)%r(1:NDIM)

  allocate(nextptcl(1:ptot))
  allocate(nextptclaux(1:ptot))
  allocate(prevptcl(1:ptot))
  allocate(prevptclaux(1:ptot))
  allocate(raylist(1:ptot))
  allocate(livelist(1:imax))
  allocate(distsqd(1:ptot))
  allocate(distids(1:ptot))

! Initialise all HEALPix arrays
!$OMP PARALLEL DO DEFAULT(SHARED)
  do i=1,imax
     HPray(i)%done        = .true.
     HPray(i)%first       = -1
     HPray(i)%last        = -1
     HPray(i)%rayend      = -1
     HPray(i)%rep(1:NDIM) = 0.0_PR
#if defined(IONIZING_UV_RADIATION)
     HPray(i)%integral    = 0.0_PR
#endif
  end do
!$OMP END PARALLEL DO

#if defined(IONIZING_UV_RADIATION)
  sph(1:ptot)%tempaux = 0.0_PR
  HIIangle = 0.0_PR
  HIIvolume = 0.0_PR
#endif
#if defined(STELLAR_WIND)
  windangle = 0.0_PR
  windvolume = 0.0_PR
  windsurface(1:pmax) = -1
#endif

! Initialize all variables for current HEALPix source
  call HP_initialize_source(isource,distids,nextptcl(1:ptot),&
       &nextptclaux(1:ptot),prevptcl(1:ptot),prevptclaux(1:ptot),&
       &raylist(1:ptot),distsqd(1:ptot))

#if defined(IONIZING_UV_RADIATION)
  HPray(1:imax)%UVdone = HPray(1:imax)%done
#endif
#if defined(STELLAR_WIND)
  HPray(1:imax)%winddone = HPray(1:imax)%done
#endif

! Initialize all remaining variables for first level
  idtot   = 0
  itot    = 12
  ltot_hp = 0
  maxres  = 0.0_PR
  resolution = 0.0_PR

#if defined(DEBUG_HP_WALK_ALL_RAYS)
  write(6,*) "Now starting main loop"
#endif
  debug_timing("HP_WALK_RAYS")


! Advance through all levels until all rays are done
! ============================================================================
  do 
     nnotdone = 0
     level    = ltot_hp
#if defined(DEBUG_HP_WALK_ALL_RAYS)
     write(6,*) "On level ",level
#endif

     ! Walk all active rays on current HEALPix level
     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) DEFAULT(SHARED) &
     !$OMP REDUCTION(+ : HIIangle,HIIvolume,idtot,nnotdone) &
     !$OMP REDUCTION(+ : windangle,windvolume) REDUCTION(max : maxres) &
     !$OMP PRIVATE(drmag,resolution,nep)
     do i=HPlevel(level)%ifirst,HPlevel(level)%ilast
        if (HPray(i)%done) cycle
#if defined(DEBUG_HP_WALK_ALL_RAYS)
        write(6,*) "Doing ray : ",i,HPray(i)
#endif
#if defined(IONIZING_UV_RADIATION)
        call HP_walk_ray(i,isource,level,resolution,nep,drmag,&
             &rsource(1:NDIM),nextptcl(1:ptot),prevptcl(1:ptot))
        if (HPray(i)%UVdone) then
           HIIangle = HIIangle + ONETHIRD*PI*4.0_PR**(-level)
           HIIvolume = HIIvolume + &
                & ONETHIRD*ONETHIRD*PI*4.0_PR**(-level)*drmag**3
        end if
#endif
#if defined(STELLAR_WIND)
        call HP_stellar_feedback(i,isource,level,drmag,&
             &rsource(1:NDIM),nextptcl(1:ptot))
        if (HPray(i)%winddone) then
           windangle = windangle + ONETHIRD*PI*4.0_PR**(-level)
           windvolume = windvolume + &
                &ONETHIRD*ONETHIRD*PI*4.0_PR**(-level)*drmag**3
        end if
#endif
#if defined(IONIZING_UV_RADIATION) && defined(STELLAR_WIND)
        if (HPray(i)%winddone .and. HPray(i)%UVdone) HPray(i)%done = .true.
#endif
        maxres = max(maxres,resolution)
        idtot = idtot + nep
        if (.not. HPray(i)%done) nnotdone = nnotdone + 1
     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------

#if defined(DEBUG_HP_WALK_ALL_RAYS)
     write(6,*) "nnotdone : ",nnotdone,ltot_hp
#endif

     ! Exit loop if all rays are now done, or maximum level has been reached.
     if (nnotdone == 0 .or. level == lmax_hp) exit
#if defined(DEBUG_HP_WALK_ALL_RAYS)
     nnotdone = 0
#endif

     ! Split current rays into child rays if there are active particles
     call HP_split_active_rays(isource,itot,level,nliverays,rsource(1:NDIM),&
          &livelist(1:imax),nextptcl(1:ptot),nextptclaux(1:ptot),&
          &prevptcl(1:ptot),prevptclaux(1:ptot),raylist(1:ptot))

     ! Exit loop of there are not live rays remaining
#if defined(DEBUG_HP_WALK_ALL_RAYS)
     write(6,*) "nliverays : ",nliverays
#endif
     if (nliverays == 0) exit

  end do
! ============================================================================


  HPmaxres = maxres

#if defined(IONIZING_UV_RADIATION)
  do i=1,ptot
     if (sph(i)%newtemp > 0) sph(i)%tempmin = &
          & max(sph(i)%tempmin,sph(i)%tempaux)
  end do
  HPsource(isource)%HIIangle = HIIangle
  HPsource(isource)%HIIvolume = HIIvolume
  write(6,*) "HIIangle  : ",HIIangle
  write(6,*) "HIIvolume : ",HIIvolume
#endif
#if defined(STELLAR_WIND)
  HPsource(isource)%windangle = windangle
  HPsource(isource)%windvolume = windvolume
  write(6,*) "windangle  : ",windangle
  write(6,*) "windvolume : ",windvolume
#if defined(PARTICLE_INJECTION_WINDS) && defined(MASS_LOADING)
  call mass_loading_winds(isource)
#elif defined(PARTICLE_INJECTION_WINDS)
  call kamikaze(isource)
#endif
#endif
#if defined(DEBUG_HP_WALK_ALL_RAYS)
  write(6,*) "Finished on level :",level,ltot_hp,lmax_hp,HPlevel(level)%ilast
  write(6,*) "No. eval. points :",idtot,"   No. of rays :",itot
#endif

! Deallocate arrays and free up allocated pointers
  if (allocated(distids)) deallocate(distids)
  if (allocated(distsqd)) deallocate(distsqd)
  if (allocated(livelist)) deallocate(livelist)
  if (allocated(raylist)) deallocate(raylist)
  if (allocated(prevptclaux)) deallocate(prevptclaux)
  if (allocated(prevptcl)) deallocate(prevptcl)
  if (allocated(nextptclaux)) deallocate(nextptclaux)
  if (allocated(nextptcl)) deallocate(nextptcl)


  return
END SUBROUTINE HP_walk_all_rays
