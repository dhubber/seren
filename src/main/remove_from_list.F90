! REMOVE_FROM_LIST.F90
! D. A. Hubber - 26/7/2007
! Removes all dead particles from data arrays and structures and moves all 
! live particles together in arrays so that they are contiguous in memory.  
! Also removes any dead particles from the tree and the HEALPix ordered 
! lists if they are activated.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE remove_from_list(ndead,deadlist)
  use interface_module, only : BH_remove_particles,reorder_particle_arrays
  use particle_module, only : ptot
  use type_module, only : pgas
#if defined(HEALPIX)
  use HP_module, only : HPtot
#endif
  implicit none

  integer, intent(in) :: ndead             ! No. of dead particles
  integer, intent(in) :: deadlist(1:ptot)  ! List of dead particles

  integer :: idead                         ! Dead particle counter
  integer :: ilive                         ! Live particle counter
  integer :: p                             ! Old particle id counter
  integer, allocatable :: dummylist(:)     ! New order of particles in arrays
#if defined(BH_TREE) || defined(IONIZING_UV_RADIATION)
  integer, allocatable :: newid(:)         ! New ids of particles
#endif

  debug2("Removing dead particles from main arrays [remove_from_list.F90]")

  allocate(dummylist(1:ptot))
#if defined(BH_TREE) || defined(IONIZING_UV_RADIATION)
  allocate(newid(1:ptot))
#endif

  idead = 1
  ilive = 1

! Determine new order of particles in arrays
! ----------------------------------------------------------------------------
  do p=1,ptot
     if (idead <= ndead) then
        if (p == deadlist(idead)) then
           dummylist(ptot - ndead + idead) = p
           idead = idead + 1
        else 
           dummylist(ilive) = p
           ilive = ilive + 1
        end if
     else 
        dummylist(ilive) = p
        ilive = ilive + 1
     end if
  end do


! Now reorder all arrays
  call reorder_particle_arrays(1,ptot,dummylist)


! Create inverse id array for tree or UV HEALPix structures.
! ----------------------------------------------------------------------------
#if defined(BH_TREE) || defined(IONIZING_UV_RADIATION)
  newid(1:ptot) = -1

! Find new ids of all live particles (so dead particles all have id = -1)
  do p=1,ptot-ndead
     newid(dummylist(p)) = p
  end do

! Now change ids of all particles in tree to new re-ordered ids
#if defined(BH_TREE)
  call BH_remove_particles(newid)
#endif

! Re-order distance lists for HEALPix sources
#if defined(IONIZING_UV_RADIATION)
  if (HPtot > 0) call HP_reorder_lists(newid)
#endif
#endif


! Reduce ptot once dead particles have been removed from all arrays
! ----------------------------------------------------------------------------
  ptot = ptot - ndead
  pgas = pgas - ndead
  call types

#if defined(BH_TREE) || defined(IONIZING_UV_RADIATION)
  deallocate(newid)
#endif
  deallocate(dummylist)

  return
END SUBROUTINE remove_from_list
