! BH_REMOVE_PARTICLES.F90
! D. A. Hubber - 29/8/2008
! Removes all particles from the BH tree that have been accreted by sinks, 
! or removed by some other algorithm.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BH_remove_particles(newid)
  use tree_module
  use particle_module, only : ptot
  implicit none

  integer, intent(in) :: newid(1:ptot)  ! Array of new particle ids

  integer :: c                          ! Cell counter
  integer :: i                          ! Aux. particle counter
  integer :: p                          ! Particle id
  integer :: pnew                       ! New particle id
  integer :: new_leaf                   ! New no. of particles in leaf
  integer :: new_plist(1:LEAFMAX)       ! List of new particle ids

  debug2("Removing dead particles from BHtree [BH_remove_particles.F90]")


! Walk through hydro tree and find all live and dead leaf cells.
! ----------------------------------------------------------------------------
!$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(SHARED) &
!$OMP PRIVATE(i,new_leaf,new_plist,p,pnew)
  do c=0,ctot_hydro
     new_leaf = 0
     if (BHhydro(c)%leaf > 0) then
        do i=1,BHhydro(c)%leaf
           p = BHhydro(c)%plist(i)
           pnew = newid(p)
           if (pnew > 0) then
              new_leaf = new_leaf + 1
              new_plist(new_leaf) = pnew
           end if
        end do

        ! Record new plist and leaf into main tree array
        if (new_leaf > 0) then
           BHhydro(c)%leaf = new_leaf
           BHhydro(c)%plist(1:new_leaf) = new_plist(1:new_leaf)
        else
           BHhydro(c)%r(1:NDIM) = 0.0_PR
           BHhydro(c)%rmax      = 0.0_PR
           BHhydro(c)%hrangemax = 0.0_PR
           BHhydro(c)%leaf      = -1
        end if
     end if
  end do
!$OMP END PARALLEL DO


! Do the same for the gravity tree
! ----------------------------------------------------------------------------
#if defined(SELF_GRAVITY)
!$OMP PARALLEL DO SCHEDULE(GUIDED,1) DEFAULT(SHARED) &
!$OMP PRIVATE(i,new_leaf,new_plist,p,pnew)
  do c=0,ctot_grav
     new_leaf = 0
     if (BHgrav(c)%leaf > 0) then
        do i=1,BHgrav(c)%leaf
           p = BHgrav(c)%plist(i)
           pnew = newid(p)
           if (pnew > 0) then
              new_leaf = new_leaf + 1
              new_plist(new_leaf) = pnew
           end if
        end do
        
        ! Record new plist and leaf into main tree array
        if (new_leaf > 0) then
           BHgrav(c)%leaf = new_leaf
           BHgrav(c)%plist(1:new_leaf) = new_plist(1:new_leaf)
        else
           BHgrav(c)%leaf      = -1
           BHgrav(c)%dminsqd   = BIG_NUMBER
           BHgrav(c)%r(1:NDIM) = 0.0_PR
#ifndef GEOMETRIC_MAC
           BHgrav(c)%mac       = BIG_NUMBER
#endif
        end if
     end if
  end do
!$OMP END PARALLEL DO
#endif


  return
END SUBROUTINE BH_remove_particles
