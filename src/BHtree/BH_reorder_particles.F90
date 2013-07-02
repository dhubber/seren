! BH_REORDER_PARTICLES.F90
! D. A. Hubber - 19/7/2008
! Reorder the particle arrays such that the particles are contiguous in 
! memory when performing a Barnes-Hut tree walk.  If self-gravity is 
! employed, then the gas particles (pgravitystart -> ptot) are ordered 
! according to the BHgrav tree.  Else, all particles (boundary, IcM and gas) 
! are ordered as according to the BHhydro tree.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BH_reorder_particles
  use interface_module, only : HP_reorder_lists,reorder_particle_arrays
  use particle_module
  use hydro_module
  use time_module
  use type_module
  use tree_module
  use filename_module
#if defined(HEALPIX)
  use HP_module, only : HPtot
#endif
  implicit none
  
  integer :: boundcount                  ! Aux. counter for boundary particles
  integer :: cdmcount                    ! Aux. counter for cdm particles
  integer :: gascount                    ! Aux. counter for gas particles
  integer :: icmcount                    ! Aux. counter for icm particles
  integer :: c                           ! Cell counter
  integer :: i                           ! Aux. counter
  integer :: p                           ! Particle id
  integer :: pstart                      ! First particle id
  integer :: pend                        ! Last particle id
  integer, allocatable :: dummylist(:)   ! List of particles in tree order
#if defined(HEALPIX)
  integer, allocatable :: newid(:)       ! New ids of particles
#endif

  debug2("Order particles according to BH tree [BH_reorder_particles.F90]")

  allocate(dummylist(1:ptot))
  do p=1,ptot
     dummylist(p) = p
  end do

! Make dummy list of particles in tree order
! ============================================================================
  pstart     = 1
  pend       = ptot
  boundcount = pboundarystart - 1
  icmcount   = picmstart - 1
  gascount   = pgasstart - 1
  cdmcount   = pcdmstart - 1
  c = 0
  do
     if (BHhydro(c)%leaf > 0) then
        do i=1,BHhydro(c)%leaf
           p = BHhydro(c)%plist(i)
           if (p >= pboundarystart .and. p <= pboundaryend) then
              boundcount = boundcount + 1
              dummylist(boundcount) = p
              BHhydro(c)%plist(i) = boundcount
           else if (p >= picmstart .and. p <= picmend) then
              icmcount = icmcount + 1
              dummylist(icmcount) = p
              BHhydro(c)%plist(i) = icmcount
           else if (p >= pgasstart .and. p <= pgasend) then
              gascount = gascount + 1
              dummylist(gascount) = p
              BHhydro(c)%plist(i) = gascount
           else if (p >= pcdmstart .and. p <= pcdmend) then
              cdmcount = cdmcount + 1
              dummylist(cdmcount) = p
              BHhydro(c)%plist(i) = cdmcount
           end if
        end do
        c = BHhydro(c)%nextcell
     else if (BHhydro(c)%leaf == 0) then
        c = BHhydro(c)%ifopen
     else
        c = BHhydro(c)%nextcell
     end if
     if (c > ctot_hydro) exit
  end do
! ============================================================================


! Now reorder all arrays based on tree traversal order
  call reorder_particle_arrays(pstart,ptot,dummylist(1:ptot))


! Create inverse id array HEALPix structures.
! ----------------------------------------------------------------------------
#if defined(HEALPIX)
  if (HPtot > 0) then
     allocate(newid(1:ptot))
     newid(1:ptot) = -1
     
     ! Find new ids of all live particles (so dead particles all have id = -1)
     do p=1,ptot
        newid(dummylist(p)) = p
     end do
     !call HP_reorder_lists(newid)
     
     deallocate(newid)
  end if
#endif

! Free allocated memory
  deallocate(dummylist)

  return
END SUBROUTINE BH_reorder_particles
