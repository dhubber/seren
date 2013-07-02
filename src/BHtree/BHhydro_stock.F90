! BHHYDRO_STOCK.F90
! D. A. Hubber - 22/01/2008
! Stocks cells in BH tree with various important properties, i.e. 
! Total mass, position of centre of mass, cell radii, etc..
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHhydro_stock(cmax,ctot,ltot,first_cell,last_cell,BHtree)
  use particle_module
  use tree_module
  implicit none

  integer, intent(in) :: cmax                          ! Max. no of cells
  integer, intent(in) :: ctot                          ! Total no of cells
  integer, intent(in) :: ltot                          ! Total no. of levels
  integer, intent(in) :: first_cell(0:LMAX)            ! 1st cell in level
  integer, intent(in) :: last_cell(0:LMAX)             ! Last cell in level
  type(BHhydro_node), intent(inout) :: BHtree(0:cmax)  ! Tree array

  integer :: c                          ! Cell counter
  integer :: cc                         ! Child cell counter
  integer :: cend                       ! Final child marker
  integer :: i                          ! Auxilary counter
  integer :: k                          ! Dimension counter
  integer :: l                          ! Level counter
  integer :: nalive                     ! No. of 'alive' child cells
  integer :: npart                      ! Number of particles in leaf cell c
  integer :: p                          ! Particle id
  real(kind=PR) :: diagonalsqd          ! Diagonal distance squared
  real(kind=PR) :: halfbox(1:NDIM)      ! Halfbox vector
  real(kind=PR) :: rc(1:NDIM)           ! Centre of mass of cell c
  real(kind=PR) :: rp(1:NDIM)           ! Position of particle p
  real(kind=PR) :: raux                 ! Max. particle dist. from cell centre
#if defined(REORDER_TREE)
  integer :: j                          ! Aux. cell counter
#endif

  debug2("Stocking cells with important quantities [BHhydro_stock.F90]")

! Initialize bounding box array
  do c=0,ctot
     BHstock(c)%bbmin(1:NDIM) = BIG_NUMBER
     BHstock(c)%bbmax(1:NDIM) = -BIG_NUMBER
  end do

! Loop over all levels in tree
! ============================================================================
  do l=ltot,0,-1

#if defined(DEBUG_BHTREESTOCK)
     write(6,*) "Stocking cells on level ",l
     write(6,*) "First cell :",first_cell(l),&
          &"  Last cell :",last_cell(l)
#endif

     ! Loop over all cells on current level
     ! -----------------------------------------------------------------------
#if defined(REORDER_TREE)
     !$OMP PARALLEL DO IF (l > 2) DEFAULT(SHARED) PRIVATE(c,cc,cend) &
     !$OMP PRIVATE(diagonalsqd,halfbox,i,k,nalive,npart,p,rc,rp,raux)
     do j=first_cell(l),last_cell(l)
        c = hydroorder(j)
#else
     !$OMP PARALLEL DO IF (l > 2) DEFAULT(SHARED) PRIVATE(cc,cend) &
     !$OMP PRIVATE(diagonalsqd,halfbox,i,k,nalive,npart,p,rc,rp,raux)
     do c=first_cell(l),last_cell(l)
#endif

#if defined(DEBUG_BHTREESTOCK)
        write(6,*) "Cell : ",c,BHtree(c)%leaf,BHtree(c)%r(1:NDIM)
#endif

        ! Zero variables for current cell
        nalive = 0
        raux = 0.0_PR

        ! If a dead cell
        ! --------------------------------------------------------------------
        if (BHtree(c)%leaf == -1) then
           BHtree(c)%r(1:NDIM) = 0.0_PR
           BHtree(c)%rmax = 0.0_PR
           BHtree(c)%hrangemax = 0.0_PR
           BHstock(c)%bbmin(1:NDIM) = 0.0_PR
           BHstock(c)%bbmax(1:NDIM) = 0.0_PR

        ! If a leaf cell with only 1 particle 
        ! --------------------------------------------------------------------
        else if (BHtree(c)%leaf == 1) then
           p = BHtree(c)%plist(1)
           rc(1:NDIM) = sph(p)%r(1:NDIM)
           BHtree(c)%r(1:NDIM) = rc(1:NDIM)
           BHtree(c)%rmax = 0.0_PR
           BHstock(c)%bbmin(1:NDIM) = rc(1:NDIM)
           BHstock(c)%bbmax(1:NDIM) = rc(1:NDIM)

        ! If a leaf cell with several particle, add contributions
        ! --------------------------------------------------------------------
        else if (BHtree(c)%leaf > 0) then
           npart = BHtree(c)%leaf
           
           ! Loop over all particles in leaf cell
           do i=1,npart
              p = BHtree(c)%plist(i)
              rp(1:NDIM) = sph(p)%r(1:NDIM)
              do k=1,NDIM
                 BHstock(c)%bbmax(k) = max(BHstock(c)%bbmax(k),rp(k))
                 BHstock(c)%bbmin(k) = min(BHstock(c)%bbmin(k),rp(k))
              end do
           end do

           ! Calculate centre of position of neighbour bounding box
           diagonalsqd = 0.0_PR
           do k=1,NDIM
              rc(k) = 0.5_PR*(BHstock(c)%bbmax(k) + BHstock(c)%bbmin(k))
              halfbox(k) = BHstock(c)%bbmax(k) - rc(k)
              diagonalsqd = diagonalsqd + halfbox(k)*halfbox(k)
           end do
           raux = sqrt(diagonalsqd)

           ! Record position, rmax and hrangemax in main tree arrays
           BHtree(c)%r(1:NDIM) = rc(1:NDIM)
           BHtree(c)%rmax = raux

        ! If it's a branch cell, add contributions due to child cells
        ! --------------------------------------------------------------------
        else if (BHtree(c)%leaf == 0) then

           ! Set first and final cell pointers for looping over child cells
           cc = BHtree(c)%ifopen
           cend = BHtree(c)%nextcell
           
           do
              if (BHtree(cc)%leaf >= 0) then
                 nalive = nalive + 1
                 do k=1,NDIM
                    BHstock(c)%bbmax(k) = &
                         &max(BHstock(cc)%bbmax(k),BHstock(c)%bbmax(k))
                    BHstock(c)%bbmin(k) = &
                         &min(BHstock(cc)%bbmin(k),BHstock(c)%bbmin(k))
                 end do 
              end if
              cc = BHtree(cc)%nextcell
              if (cc == cend) exit
           end do

           ! If cell contains all dead child cells, mark as dead also
           if (nalive == 0) then
              BHtree(c)%leaf = -1
              BHtree(c)%r(1:NDIM) = BIG_NUMBER
              BHtree(c)%rmax = 0.0_PR
              BHtree(c)%hrangemax = 0.0_PR
              BHstock(c)%bbmin(1:NDIM) = BIG_NUMBER
              BHstock(c)%bbmax(1:NDIM) = BIG_NUMBER
              cycle

           ! Else calculate centre of position of neighbour bounding box
           else
              diagonalsqd = 0.0_PR
              do k=1,NDIM
                 rc(k) = 0.5_PR*(BHstock(c)%bbmax(k) + BHstock(c)%bbmin(k))
                 halfbox(k) = BHstock(c)%bbmax(k) - rc(k)
                 diagonalsqd = diagonalsqd + halfbox(k)*halfbox(k)
              end do
              raux = sqrt(diagonalsqd)
              
              ! Record position, rmax and hrangemax in main tree arrays
              BHtree(c)%r(1:NDIM) = rc(1:NDIM)
              BHtree(c)%rmax = raux
           end if 

        end if
        ! --------------------------------------------------------------------

     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------

  end do
! ============================================================================


#if defined(DEBUG_BHTREESTOCK)
  do c=0,ctot
     write(6,*) c,BHtree(c)
  end do
  write(6,*) "Verifying level status"
  do l=0,ltot_hydro
     write(6,*) l,first_cell(l),last_cell(l)
  end do
#endif

  return 
END SUBROUTINE BHhydro_stock
