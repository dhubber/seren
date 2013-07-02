! BHGRAV_BUILD.F90
! D. A. Hubber - 22/01/2008
! Builds Barnes-Hut tree for a given distribution of particles. Only builds 
! tree for self-gravitating particles for use in gravity calculation 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHgrav_build
  use definitions
  use tree_module
  use particle_module
  use type_module
  implicit none

  integer :: c                             ! Cell counter
  integer :: i                             ! Aux. particle counter
  integer :: nptcls                        ! No. of particles in tree
  integer :: p                             ! Particle id
  integer, allocatable :: plist(:)         ! List if particle ids
  real(kind=PR), allocatable :: r(:,:)     ! Particle positions

  allocate(plist(1:ptot))
  allocate(r(1:NDIM,1:ptot))

  ! Create list of self-gravitating particles
  call create_particle_list(nptcls,plist,gravmask)

  ! Create list of particle positions
  do i=1,nptcls
     p = plist(i)
     r(1:NDIM,i) = sph(p)%r(1:NDIM)
  end do

  ! Build skeleton of tree using reduced particle list
  call BH_build_skeleton(nptcls,plist(1:nptcls),r(1:NDIM,1:nptcls))

  ! Now copy skeleton of tree to main gravity tree
  ctot_grav = ctot_skeleton
  ltot_grav = ltot_skeleton

  first_cell_grav(0:ltot_grav) = first_cell_skeleton(0:ltot_grav)
  last_cell_grav(0:ltot_grav) = last_cell_skeleton(0:ltot_grav)

  do c=0,ctot_grav
     BHgrav(c)%leaf = BHskeleton(c)%leaf
     BHgrav(c)%nextcell = BHskeleton(c)%nextcell
     BHgrav(c)%ifopen = BHskeleton(c)%ifopen
     BHgrav(c)%plist(1:LEAFMAX) = BHskeleton(c)%plist(1:LEAFMAX)
  end do

  deallocate(r)
  deallocate(plist)

  return 
END SUBROUTINE BHgrav_build
