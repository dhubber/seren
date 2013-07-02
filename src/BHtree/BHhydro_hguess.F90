! BHHYDRO_HGUESS.F90
! D. A. Hubber - 9/1/2009
! Calculates mass and total contained number of particles in each cell. 
! Then walks the tree to determine reasonable guesses of the smoothing 
! length for all particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHhydro_hguess
  use particle_module
  use tree_module
  use type_module
  use neighbour_module, only : pp_gather
  implicit none

  integer :: c                      ! Cell counter
  integer :: cc                     ! Child cell counter
  integer :: cend                   ! Final child marker
  integer :: i                      ! Auxilary counter
  integer :: l                      ! Level counter
  integer :: laux                   ! Aux. level counter
  integer :: npart                  ! Number of particles in leaf cell c
  integer :: p                      ! Particle id
  integer :: pstack(0:LMAX)         ! Stack of particles numbers for walk
  real(kind=PR) :: hguess           ! Guess of smoothing length
  real(kind=PR) :: rstack(0:LMAX)   ! Stack of cell size for walk
  integer, allocatable :: levelof(:)           ! Level of cells
  integer, allocatable :: pcontained(:)        ! No. of particles in cells
  real(kind=PR), allocatable :: mcontained(:)  ! Mass of particles in cells
#if defined(REORDER_TREE)
  integer :: j                      ! Aux. cell counter
#endif
#if defined(HMASS)
  real(kind=PR) :: mgather          ! Target mass for hguess search
  real(kind=PR) :: mstack(0:LMAX)   ! Stack of particle masses for walk
#endif

  debug2("Walking tree to obtain initial guesses for h [BHhydro_hguess.F90]")

  allocate(levelof(0:ctot_hydro))
  allocate(pcontained(0:ctot_hydro))
  allocate(mcontained(0:ctot_hydro))

  levelof(0:ctot_hydro)    = 0
  pcontained(0:ctot_hydro) = 0
  mcontained(0:ctot_hydro) = 0.0_PR

#if defined(HMASS)
  mgather = real(pp_gather,PR)*sph(p)%m
#endif

! Loop over all levels and calculate pcontained and mcontained for all cells
! ============================================================================
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(c,cc,cend,i,npart,p)
  do l=ltot_hydro,0,-1

     ! Loop over all cells on current level
     ! -----------------------------------------------------------------------
#if defined(REORDER_TREE)
     !$OMP DO SCHEDULE(GUIDED,2) PRIVATE(c)
     do j=first_cell_hydro(l),last_cell_hydro(l)
        c = hydroorder(j)
#else
     !$OMP DO SCHEDULE(GUIDED,2)
     do c=first_cell_hydro(l),last_cell_hydro(l)
#endif

        levelof(c) = l

        ! If a leaf cell, add contributions due to particles
        ! --------------------------------------------------------------------
        if (BHhydro(c)%leaf > 0) then
           npart = BHhydro(c)%leaf
           
           ! Loop over all particles in leaf cell
           pcontained(c) = npart
           do i=1,npart
              p = BHhydro(c)%plist(i)
              mcontained(c) = mcontained(c) + sph(p)%m
           end do
           
        ! If it's a branch cell, add contributions due to child cells
        ! --------------------------------------------------------------------
        else if (BHhydro(c)%leaf == 0) then

           ! Loop over list of child cells
           cc = BHhydro(c)%ifopen
           cend = BHhydro(c)%nextcell
           do
              pcontained(c) = pcontained(c) + pcontained(cc)
              mcontained(c) = mcontained(c) + mcontained(cc)
              cc = BHhydro(cc)%nextcell
              if (cc == cend) exit
           end do
           
        end if
        ! --------------------------------------------------------------------

     end do
     !$OMP END DO
     ! -----------------------------------------------------------------------

  end do
!$OMP END PARALLEL
! ============================================================================


#if defined(DEBUG_H_GUESS)
  write(6,*) "Finished creating temporary tree variables"
#endif


! Now walk through the tree to obtain guesses for smoothing lengths
! ============================================================================
  c = 0

! ----------------------------------------------------------------------------
  do 
     l = levelof(c)
     pstack(l) = pcontained(c)
     rstack(l) = BHhydro(c)%rmax
#if defined(HMASS)
     mstack(l) = mcontained(c)
#endif

     ! If a leaf cell, walk up number/mass stack to find guess of h
     ! -----------------------------------------------------------------------
     if (BHhydro(c)%leaf > 0) then
        laux = l
        
        ! Walk up the stack to find best guess of h
        do 
#if defined(HGATHER) || defined(H_RHO)
           if (pstack(laux) > int(0.5_PR*real(pp_gather,PR)) .or. laux == 0) &
                & exit
#elif defined(HMASS)
           if (mstack(laux) > 0.5_PR*mgather .or. laux == 0) exit
#endif
           laux = laux - 1
        end do
        hguess = INVKERNRANGE*rstack(laux)

#if defined(DEBUG_H_GUESS)
        write(6,*) "Estimated smoothing length =", hguess,c,l
        if (hguess < 0.0_PR) stop "Fatal error: hguess is negative"
#endif

        do i=1,BHhydro(c)%leaf
           p = BHhydro(c)%plist(i)
           sph(p)%h = hguess
        end do
        c = BHhydro(c)%nextcell

     else if (BHhydro(c)%leaf == 0) then
        c = BHhydro(c)%ifopen
     else
        c = BHhydro(c)%nextcell
     end if
     ! -----------------------------------------------------------------------

     if (c > ctot_hydro) exit
  end do
! ----------------------------------------------------------------------------

  deallocate(mcontained)
  deallocate(pcontained)
  deallocate(levelof)

  return
END SUBROUTINE BHhydro_hguess
