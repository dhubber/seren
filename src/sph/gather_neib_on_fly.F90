! GATHER_NEIB_ON_FLY.F90
! C. P. Batty & D. A. Hubber - 3/6/2007
! Determines the list of neighbouring particle ids of particle p by 
! gather-only within a distance hrange (i.e. |dr| < hrange). 
! Returns no. of neibs (pp_tot) plus list (pp_list).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE gather_neib_on_fly(p,pp_max,pp_tot,pp_list,rsearch,hrange,typemask)
  !use interface_module, only : BHhydrowalk_hgather,distance2
  use particle_module
  use neighbour_module
  use type_module
  use tree_module
  implicit none

  integer, intent(in) :: p                            ! chosen particle
  integer, intent(inout) :: pp_max                    ! length of pp_list
  integer, intent(out) :: pp_tot                      ! number of pot. neibs
  integer, allocatable, intent(out) :: pp_list(:)     ! list of pot. neibs
  real(kind=PR), intent(in) :: rsearch(1:NDIM)        ! r of search sphere
  real(kind=PR), intent(in) :: hrange                 ! size of search sphere
  logical, optional, intent(in) :: typemask(1:ntypes) ! part. types to include?

  integer :: pp                          ! neighbouring particles (p')
  real(kind=PR) :: dr(1:NDIM)            ! vector displacements (p'-p)
  real(kind=PR) :: drsqd                 ! p'-p separation squared
  real(kind=PR) :: hrangesqd             ! particle radius squared
#if defined(BH_TREE)
  integer :: i                           ! counter in neighbour search
  integer :: pp_pot                      ! no. of potential neighbours
#endif
#if defined(DEBUG_GATHER_NEIB)
  integer :: n_gather                    ! Number of neighbours by gather
#if defined(HMASS)
  real(kind=PR) :: mfrac                 ! Mass fraction
  real(kind=PR) :: mgather               ! Mass by gather
  real(kind=PR) :: mp                    ! Mass of particle p
  real(kind=PR) :: mtot                  ! Total mass
#endif
#endif

  debug3("Obtaining neighbour list [gather_neib_on_fly.F90] for particle ", p)

! Make local copies of particle position and initialize other variables
  hrangesqd = hrange*hrange
  pp_tot = 0
#if defined(DEBUG_GATHER_NEIB)
  n_gather = 0
#if defined(HMASS)
  mp = sph(p)%m
  mtot = 0.0_PR
#endif
#endif

! Allocate initial array size
#if defined(GHOST_PARTICLES)
  pp_max = min(LISTSIZE,ptot+pghost)
#else
  pp_max = min(LISTSIZE,ptot)
#endif
  allocate(pp_list(1:pp_max))


! Obtain potential neighbour list by tree walk if using tree.  If list 
! obtained from tree is too big for array (signalled by pp_pot = -1), 
! then reallocate array and walk tree again.
! ----------------------------------------------------------------------------
#if defined(BH_TREE)
  do 
     pp_pot = 0
     call BHhydrowalk_hgather(rsearch(1:NDIM),hrange,pp_pot,&
          & pp_max,pp_list,ctot_hydro,BHhydro(0:ctot_hydro))
#if defined(GHOST_PARTICLES)
     if (pghost > 0 .AND. pp_pot /= -1) &
          & call BHhydrowalk_hgather(rsearch(1:NDIM),hrange,&
          & pp_pot,pp_max,pp_list,ctot_ghost,BHghost(0:ctot_ghost))
#endif
     if (pp_pot < 0) then
#if defined(GHOST_PARTICLES)
        pp_max = ptot + pghost
#else
        pp_max = ptot
#endif
        deallocate(pp_list)
        allocate(pp_list(1:pp_max))
     else
        exit
     end if
  end do
#endif
! ----------------------------------------------------------------------------


! Loop over potential list or all particles (depending on flags) and find 
! the real neighbours.  If array is too small (i.e. pp_tot > pp_max), 
! then reallocate array and repeat.
! ----------------------------------------------------------------------------
  do
#if defined(BH_TREE)
     do i=1,pp_pot
        pp = pp_list(i)
#elif defined(GHOST_PARTICLES)
     do pp=1,ptot+pghost
#else
     do pp=1,ptot
#endif
        ! Cycle to next particle if neighbour is particle itself, particle 
        ! type is not allowed, or particle is outside search radius
        if (pp == p) cycle
        if (present(typemask)) then
           if (.not. typemask(sph(pp)%ptype)) cycle
        end if
#if defined(PERIODIC) && !defined(GHOST_PARTICLES)
        call distance2(rsearch(1:NDIM),pp,dr(1:NDIM),drsqd)
#else
        dr(1:NDIM) = sph(pp)%r(1:NDIM) - rsearch(1:NDIM)
        drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif
        if (drsqd > hrangesqd) cycle

        ! If particle is accepted, record id in array
        pp_tot = pp_tot + 1
        if (pp_tot <= pp_max) pp_list(pp_tot) = pp
        
        ! Record numbers of neighbours for debugging 
#if defined(DEBUG_GATHER_NEIB)
        if (drsqd <= hrangesqd) n_gather = n_gather + 1
#if defined(HMASS)
        if (drsqd <= hrangesqd) mtot = mtot + sph(pp)%m
#endif        
#endif
     end do

     ! If the array size was too small, reallocate and try again.
     if (pp_tot > pp_max) then
#if defined(GHOST_PARTICLES)
        pp_max = ptot+pghost
#else
        pp_max = ptot
#endif
        pp_tot = 0
        deallocate(pp_list)
        allocate(pp_list(1:pp_max))
     else
        exit
     end if

  end do
! ----------------------------------------------------------------------------


#if defined(DEBUG_GATHER_NEIB)
  write(6,*) "h(",p,") : ",hp,"   n_gather : ", n_gather, pp_tot
#if defined(HMASS)
  mgather = mp * real(pp_gather,PR)
  mfrac = mtot / mp
  write(6,*) "mp :",mp," mgather :",mgather,"  mtot :",mtot,"  mfrac :",mfrac
#endif
#endif


  return
END SUBROUTINE gather_neib_on_fly
