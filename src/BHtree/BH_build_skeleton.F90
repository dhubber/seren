! BHGRAV_BUILD.F90
! D. A. Hubber - 22/01/2008
! Builds Barnes-Hut tree for a given distribution of particles. Only builds 
! tree for self-gravitating particles for use in gravity calculation.
! Some parts borrowed/re-derived from DRAGON code written by S. P. Goodwin.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BH_build_skeleton(nptcls,ids,r)
  use interface_module, only : bounding_box
  use particle_module
  use type_module
  use tree_module
#if defined(OPENMP)
  use omp_lib
#endif
  use debug_tests
  implicit none

  integer, intent(in) :: nptcls                    ! No. of particles
  integer, intent(in) :: ids(1:nptcls)             ! Particle ids
  real(kind=PR), intent(in) :: r(1:NDIM,1:nptcls)  ! Particle positions

  logical :: alldone                   ! Flag for when all particles are done
  integer :: c                         ! Cell counter
  integer :: ckid                      ! Auxillary child id
  integer :: clist(1:NCHILD)           ! List of new children cells
  integer :: i                         ! Auxilary loop counter
  integer :: j                         ! Auxilary loop counter
  integer :: k                         ! Dimension counter
  integer :: l                         ! Tree level counter
  integer :: nkids                     ! No. of new children cells
  integer :: nlist                     ! No. of cells to divide on level l
  integer :: nptclin(1:NCHILD)         ! No. of child cell divisions in BH tree
  integer :: p                         ! Particle counter
  integer, allocatable :: templist(:)  ! List of cells to divide
  real(kind=PR) :: cellsize            ! Size of cell
  real(kind=PR) :: rc(1:NDIM)          ! Local copy of position of cell c
  real(kind=PR) :: rp(1:NDIM)          ! Local copy of position of particle p

  debug2("Building Barnes-Hut tree skeleton [BH_build_skeleton.F90]")

  nlist = 0
  allocate(templist(1:cmax_skeleton))
  cellof(1:ptot) = 0

! Initialize variables before spatial decomposition
  do i=1,nptcls-1
     p = ids(i)
     BHnextptcl(p) = ids(i + 1)
  end do

  do c=0,cmax_skeleton
     BHskeleton(c)%childof(1:NCHILD) = -2
     BHskeleton(c)%pfirst = -1
     BHskeleton(c)%plast = -1
     BHskeleton(c)%leaf = 0
  end do

  BHskeleton(0)%pfirst = ids(1)
  BHskeleton(0)%plast = ids(nptcls)
  first_cell_skeleton(0:LMAX) = 0
  last_cell_skeleton(0:LMAX) = 0

! Start with root cell and work downwards
  c = 0
  l = 0
  ctot_skeleton = 0
  ltot_skeleton = 0
  alldone = .false.
  BHskeleton(0)%leaf = 0
  BHskeleton(0)%nextcell = LARGEST_INT
  BHskeleton(0)%ifopen = LARGEST_INT
  BHskeleton(0)%plist(1:LEAFMAX) = 0
  BHskeleton(0)%r(1:NDIM) = 0.0_PR
  nlist = 1
  templist(1) = 0

! Return if there are no particles in the tree
  if (nptcls == 0) return

! Determine square bounding box of all self-gravitating gas particles and 
! set position of root cell at centre
  call bounding_box(nptcls,r(1:NDIM,1:nptcls),rmax(1:NDIM),rmin(1:NDIM))
  cellsize = 0.5_PR*maxval(rmax(1:NDIM) - rmin(1:NDIM))
  BHskeleton(0)%r(1:NDIM) = 0.5_PR*(rmin(1:NDIM) + rmax(1:NDIM))


! Main loop
! Loop through deeper and deeper levels until all particles are in leafs
! ============================================================================
  do 
     cellsize = 0.5_PR*cellsize

     ! Loop over all cells on parent level
     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(SHARED) &
     !$OMP PRIVATE(c,ckid,k,nkids,nptclin,p,rc,rp)
     do j=1,nlist
        c = templist(j)
        
        ! Loop over all particles in current cell, starting with the first
        p = BHskeleton(c)%pfirst
        rc(1:NDIM) = BHskeleton(c)%r(1:NDIM)
        do
           rp(1:NDIM) = sph(p)%r(1:NDIM)
           
           ! Determine which quadrant/octant particle p is in 
           ckid = 1
           if (rp(1) > rc(1)) ckid = ckid + 1
#if NDIM==2 || NDIM==3
           if (rp(2) > rc(2)) ckid = ckid + 2
#endif
#if NDIM==3
           if (rp(3) > rc(3)) ckid = ckid + 4
#endif
           whichchild(p) = ckid
           BHskeleton(c)%childof(ckid) = -1
!if (l >= LMAX - 1) write(6,*) j,c,p,rp(1:NDIM)
           if (p == BHskeleton(c)%plast) exit
           p = BHnextptcl(p)
        end do
               
     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------


     ! Find all new child cells (Must be done in serial)
     ! -----------------------------------------------------------------------
     do j=1,nlist
        c = templist(j)           
        do k=1,NCHILD
           if (BHskeleton(c)%childof(k) == -1) then
              ctot_skeleton = ctot_skeleton + 1
              if (ctot_skeleton >= cmax_skeleton) then
                 write(6,*) "Tree too big : ",ctot_skeleton,cmax_skeleton
                 stop
              end if
              BHskeleton(c)%childof(k) = ctot_skeleton
           end if
        end do
     end do
     ! ----------------------------------------------------------------------- 


     ! Find number of particles in each new child cell
     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO IF (l > 0) SCHEDULE(GUIDED,1) DEFAULT(SHARED) &
     !$OMP PRIVATE(c,ckid,clist,i,k,nkids,nptclin,p,rc)
     do j=1,nlist
        c = templist(j)
        nptclin(1:NCHILD) = 0
        p = BHskeleton(c)%pfirst
        
        do
           ckid = BHskeleton(c)%childof(whichchild(p))
           cellof(p) = ckid
           if (BHskeleton(ckid)%pfirst == -1) then
              BHskeleton(ckid)%pfirst = p
           else
              BHnextptcl(BHskeleton(ckid)%plast) = p
           end if
           BHskeleton(ckid)%plast = p
           nptclin(whichchild(p)) = nptclin(whichchild(p)) + 1
           if (p == BHskeleton(c)%plast) exit
           p = BHnextptcl(p)
        end do
        
        ! Now create child cells if there are particles in that division
        nkids = 0
        rc(1:NDIM) = BHskeleton(c)%r(1:NDIM)
        
        ! Assign central positions of new cells
        ! --------------------------------------------------------------------
        do k=1,NCHILD
           if (nptclin(k) == 0) cycle
           ckid = BHskeleton(c)%childof(k)

           ! If leafcell, flag it and store children ids
           if (nptclin(k) > 0 .and. nptclin(k) <= LEAFMAX) then
              BHskeleton(ckid)%leaf = nptclin(k)
              i = 0
              p = BHskeleton(ckid)%pfirst
              do
                 i = i + 1
                 BHskeleton(ckid)%plist(i) = p
                 if (p == BHskeleton(ckid)%plast) exit
                 p = BHnextptcl(p)
              end do

           ! If not, record particles within for next level search
           else if (nptclin(k) > LEAFMAX) then
              BHskeleton(ckid)%leaf = 0
           end if
           
           ! Assign positions of child cells
           if (k==1 .or. k==3 .or. k==5 .or. k==7) then
              BHskeleton(ckid)%r(1) = rc(1) - cellsize
           else
              BHskeleton(ckid)%r(1) = rc(1) + cellsize
           end if
#if NDIM==2 || NDIM==3
           if (k==1 .or. k==2 .or. k==5 .or. k==6) then
              BHskeleton(ckid)%r(2) = rc(2) - cellsize
           else
              BHskeleton(ckid)%r(2) = rc(2) + cellsize
           end if          
#endif
#if NDIM==3
           if (k < 5) then
              BHskeleton(ckid)%r(3) = rc(3) - cellsize
           else
              BHskeleton(ckid)%r(3) = rc(3) + cellsize
           end if
#endif

           ! Record id of cells for following linked list construction
           nkids = nkids + 1
           clist(nkids) = ckid

        end do
        ! --------------------------------------------------------------------

        ! Set up linked lists from parent to children cells
        BHskeleton(c)%ifopen = clist(1)
        do k=1,nkids-1
           ckid = clist(k)
           BHskeleton(ckid)%nextcell = clist(k + 1)
        end do
        BHskeleton(clist(nkids))%nextcell = BHskeleton(c)%nextcell

     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------

     ! Record first and last cells in newly created level
     l = l + 1
     first_cell_skeleton(l) = last_cell_skeleton(l-1) + 1
     last_cell_skeleton(l) = ctot_skeleton
     
     ! Check that all new cells are children cells
     nlist = 0
     alldone = .true.
     do c=first_cell_skeleton(l),last_cell_skeleton(l)
        if (BHskeleton(c)%leaf == 0) then
           alldone = .false.
           nlist = nlist + 1
           templist(nlist) = c
        end if
     end do

     ! If all children are leafs, finish building tree
     if (alldone) exit
     alldone = .false.

     ! Check we have not reached maximum level
     if (l >= LMAX) then
        write (6,*) 'Reached maximum level in BH tree build'
        write (6,*) 'There are ', nlist, ' cells left to split'
        do c=first_cell_skeleton(l), last_cell_skeleton(l)
           if (BHskeleton(c)%leaf == 0) then
              write (6,*) 'Cell ', c, ' is still a node cell.'
              write (6,*) "Particles in 'cell': "
              p = BHskeleton(c)%pfirst
              do
                 write (6,*) '   Particle: ', p
                 write (6,*) '   porig = ', sph(p)%porig
                 write (6,*) '   r = ', sph(p)%r
                 write (6,*) '   --------------------------------------'
                 if (p == BHskeleton(c)%plast) exit
                 p = BHnextptcl(p)
              end do
           end if
        end do
        call doublecheck_particle_uniqueness(.TRUE.)
        stop 1
     end if

  end do
! ============================================================================


! Record total number of levels in the tree
  ltot_skeleton  =  l

  deallocate(templist)

  return 
END SUBROUTINE BH_build_skeleton
