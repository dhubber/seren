! GET_NEIB_ON_FLY.F90
! C. P. Batty & D. A. Hubber - 3/6/2007
! Calculates neighbour list for particle p (using gather/scatter method) 
! Returns the neighbour list and no. of neighbours as subroutine parameters 
! rather than storing in main arrays.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE get_neib_on_fly(p,pp_tot,pp_totmax,pp_list,rsearch,hrange,typemask)
  !use interface_module, only : BHhydro_walk,binary_neibfind,distance2
  use particle_module
  use neighbour_module
  use tree_module
  use type_module
  implicit none

  integer, intent(in) :: p                            ! chosen particle id
  real(kind=PR), intent(in) :: rsearch(1:NDIM)        ! Position of particle p
  real(kind=PR), intent(in) :: hrange                 ! h of particle p
  integer, intent(out) :: pp_tot                      ! number of pot. neibs
  integer, intent(in) :: pp_totmax                    ! ..
  integer, allocatable, intent(out) :: pp_list(:)     ! list of pot. neibs
  logical, optional, intent(in) :: typemask(1:ntypes) ! part. types to include?

  integer :: pp                         ! neighbouring particles (p')
  integer :: pp_max                     ! length of pp_list array
  real(kind=PR) :: dr(1:NDIM)           ! vector displacements (p'-p)
  real(kind=PR) :: drsqd                ! p'-p separation squared
  real(kind=PR) :: hrangesqd_p          ! particle radius p squared
  real(kind=PR) :: hrangesqd_pp         ! particle radius pp squared
#if defined(BH_TREE) || defined(BINARY_TREE)
  integer :: i                          ! counter in neighbour search
  integer :: pp_pot                     ! no. of potential neighbours
#endif
#if defined(DEBUG_GET_NEIB)
  integer :: n_gather                   ! No. of neighbours by gather
  integer :: n_scatter                  ! No. of neighbours by scatter
#if defined(HMASS)
  real(kind=PR) :: mfrac                ! Mass fraction
  real(kind=PR) :: mgather              ! Mass by gather
  real(kind=PR) :: mp                   ! Mass of particle p
  real(kind=PR) :: mtot                 ! Total mass
#endif
#endif

  debug3("Obtaining neighbour list [get_neib.F90] for particle ", p)

! Make local copy of particle position and initialize other variables
  hrangesqd_p = hrange*hrange
  pp_tot = 0
#if defined(DEBUG_GET_NEIB)
  n_gather = 0; n_scatter = 0
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

! Obtain potential neighbour list by tree walk if using tree
! ----------------------------------------------------------------------------
#if defined(BH_TREE) || defined(BINARY_TREE)
  pp_max = min(LISTSIZE,ptot)
  do 
     pp_pot = 0
#if defined(BH_TREE)
     call BHhydro_walk(rsearch(1:NDIM),hrange,pp_pot,&
          &pp_max,pp_list,ctot_hydro,BHhydro)
#if defined(GHOST_PARTICLES)
     if (pghost > 0 .and. pp_pot /= -1) &
          &call BHhydro_walk(rsearch(1:NDIM),hrange,pp_pot,&
          &pp_max,pp_list,ctot_ghost,BHghost)
#endif
#elif defined(BINARY_TREE)
     call binary_neibfind(rsearch,hp,pp_pot,pp_max,pp_list)
#endif
     if (pp_pot < 0) then
#if defined(GHOST_PARTICLES)
        pp_max = ptot+pghost
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


! Loop over potential list or all particles (depending on flags) and find 
! the real neighbours.  If array is too small (i.e. pp_tot > pp_max), 
! then reallocate array and repeat.
! ----------------------------------------------------------------------------
  do
#if defined(BH_TREE) || defined(BINARY_TREE)
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
        hrangesqd_pp = KERNRANGESQD*sph(pp)%h*sph(pp)%h
#if defined(PERIODIC) && !defined(GHOST_PARTICLES)
        call distance2(rsearch(1:NDIM),pp,dr(1:NDIM),drsqd)
#else
        dr(1:NDIM) = sph(pp)%r(1:NDIM) - rsearch(1:NDIM)
        drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif
        if (drsqd > hrangesqd_p .and. drsqd > hrangesqd_pp) cycle

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

#if defined(DEBUG_GET_NEIB)
  write(6,*) "h(",p,") : ",hp,"   n_gather : ",  & 
       & n_gather,"  n_scatter : ",n_scatter
#if defined(HMASS)
  mgather = mp * real(pp_gather,PR)
  mfrac = mtot / mp
  write(6,*) "mp :",mp," mgather :",mgather,"  mtot :",mtot,"  mfrac :",mfrac
#endif
#endif

  return
END SUBROUTINE get_neib_on_fly
