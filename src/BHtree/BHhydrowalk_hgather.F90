! BHHYDROWALK_HGATHER.F90
! D. A. Hubber - 24/01/2008
! Walks hydro tree to find potential neighbour list.  Used in h_gather to 
! find all potential neighbours by gather. (in contrast, BHhydro_walk.F90 
! finds potential neighbours by gather and scatter.)  Appends particles to 
! pp_list array (i.e. must be zeroed before calling this routine).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHhydrowalk_hgather(rp,hrange,pp_pot,pp_max,pp_list,ctot,BHtree)
  use interface_module, only : distance3
  use definitions
  use tree_module
  implicit none

  integer, intent(in) :: ctot                      ! Total no. of cells in tree
  integer, intent(in) :: pp_max                    ! Max. no. of pot. neibs
  integer, intent(inout) :: pp_pot                 ! No. of pot. neighbours
  integer, intent(inout) :: pp_list(1:pp_max)      ! Potential neighbour list
  real(kind=PR),intent(in) :: hrange               ! Search range of tree walk
  real(kind=PR),intent(in) :: rp(1:NDIM)           ! Position of particle p
  type(BHhydro_node), intent(in) :: BHtree(0:ctot) ! BH tree array

  integer :: i                                ! Auxilary leaf counter
  integer :: c                                ! Cell counter
  real(kind=PR) :: dr(1:NDIM)                 ! Relative displacement vector
  real(kind=PR) :: drsqd                      ! Distance squared
  real(kind=PR) :: maxdistsqd                 ! Maximum gather dist. squared

! Start on root cell (c=0) and work our way down.
  c = 0

! Main cell-walk loop
! ============================================================================
  do

     ! Compute distance from particle to cell centre, plus the maximum 
     ! radial extent (squared) of all particles within the cell
#if defined(PERIODIC) && !defined(GHOST_PARTICLES)
     call distance3(rp(1:NDIM),BHtree(c)%r(1:NDIM),dr(1:NDIM),drsqd)
#else
     dr(1:NDIM) = BHtree(c)%r(1:NDIM) - rp(1:NDIM)
     drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif
     maxdistsqd = (hrange + BHtree(c)%rmax)*(hrange + BHtree(c)%rmax)

     ! Check if search sphere overlaps cell 'neighbour sphere'
     ! -----------------------------------------------------------------------
     if (drsqd <= maxdistsqd) then 

        ! Check if c is a leaf cell or not
        if (BHtree(c)%leaf > 0) then

           ! If templist has become too long, flag with -1 and return
           if (pp_pot + BHtree(c)%leaf > pp_max) then
              pp_pot = -1
              return
           end if
              
           ! Add all particles in leaf cell to potential neighbour list
           do i=1,BHtree(c)%leaf
              pp_pot = pp_pot + 1
              pp_list(pp_pot) = BHtree(c)%plist(i)
           end do
           c = BHtree(c)%nextcell
           
        else if (BHtree(c)%leaf == 0) then
           c = BHtree(c)%ifopen
        else
           c = BHtree(c)%nextcell
        end if
     else
        c = BHtree(c)%nextcell
     end if
     ! -----------------------------------------------------------------------

     if (c > ctot) exit
  end do
! ============================================================================

#if defined(DEBUG_BHTREEWALK)
  write(6,*) hp,pp_pot,pp_list(1:pp_pot)
#endif

  return
END SUBROUTINE BHhydrowalk_hgather
