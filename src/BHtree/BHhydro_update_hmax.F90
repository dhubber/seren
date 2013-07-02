! BHHYDRO_UPDATE_HMAX.F90
! D. A. Hubber - 22/01/2008
! Updates values of hrangemax in all cells once smoothing lengths have been 
! recomputed in sph_update (h's are computed by gather so updated hmax's 
! are not needed when calculating h's).  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHhydro_update_hmax(cmax,ctot,ltot,first_cell,last_cell,BHtree)
  use particle_module, only : sph
  use tree_module
  implicit none

  integer, intent(in) :: cmax                          ! Max. no of cells
  integer, intent(in) :: ctot                          ! Total no of cells
  integer, intent(in) :: ltot                          ! ..
  integer, intent(in) :: first_cell(0:LMAX)            ! 1st cell in level
  integer, intent(in) :: last_cell(0:LMAX)             ! Last cell in level
  type(BHhydro_node), intent(inout) :: BHtree(0:cmax)  ! Tree array

  integer :: c                  ! Cell counter
  integer :: cc                 ! Child cell counter
  integer :: cend               ! Final child marker
  integer :: i                  ! Auxilary counter
  integer :: l                  ! Level counter
  integer :: npart              ! No. of particles in leaf cell c
  integer :: p                  ! Particle id
  real(kind=PR) :: hrangemax    ! Max. particle extent in cell c

  debug2("Updating values of hmax in BHhydro tree [BHhydro_update_hmax.F90]")

! Loop over all tree levels
! ============================================================================
  do l=ltot,0,-1

     ! Loop over all cells on current level
     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO IF (l > 2) DEFAULT(SHARED) &
     !$OMP PRIVATE(cc,cend,hrangemax,i,npart,p) 
     do c=first_cell(l),last_cell(l)

        ! Zero aux. variable for finding maximum particle range
        hrangemax = 0.0_PR

        ! If a dead leaf cell (e.g. due to accretion)
        ! --------------------------------------------------------------------
        if (BHtree(c)%leaf == -1) then
           BHtree(c)%hrangemax = 0.0_PR

        ! If a leaf cell with only 1 particle
        ! --------------------------------------------------------------------
        else if (BHtree(c)%leaf == 1) then
           p = BHtree(c)%plist(1)
           BHtree(c)%hrangemax = KERNRANGE*sph(p)%h

        ! If a leaf cell, add contributions due to particles
        ! --------------------------------------------------------------------
        else if (BHtree(c)%leaf > 0) then
           npart = BHtree(c)%leaf
           
           ! Loop over all particles in leaf cell
           do i=1,npart
              p = BHtree(c)%plist(i)
              hrangemax = max(hrangemax,KERNRANGE*sph(p)%h)
           end do
           
           ! Record hrangemax in main tree array
           BHtree(c)%hrangemax = hrangemax

        ! If it's a branch cell, add contributions due to child cells
        ! --------------------------------------------------------------------
        else if (BHtree(c)%leaf == 0) then

           ! Set first and final child cell (same nextcell as parent)
           cc = BHtree(c)%ifopen
           cend = BHtree(c)%nextcell
           do
              hrangemax = max(hrangemax,BHtree(cc)%hrangemax)
              cc = BHtree(cc)%nextcell
              if (cc == cend) exit
           end do

           ! Record hrangemax in main tree array
           BHtree(c)%hrangemax = hrangemax
           
        end if
        ! --------------------------------------------------------------------

     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------

  end do
! ============================================================================

  return 
END SUBROUTINE BHhydro_update_hmax
