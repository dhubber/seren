! CREATE_PARTICLE_LIST.F90
! D. A. Hubber & A. McLeod - 17/03/2011
! Create a list of the ids of all particles that are selected in the 
! (optional) mask variable.  If no mask is passed, pass all particle ids.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE create_particle_list(nptcls,plist,typemask)
  use definitions
  use type_module
  use type_module
  use time_module
  use particle_module, only : ptot,sph
  implicit none

  integer, intent(out) :: nptcls            ! No. of active particles
  integer, intent(out) :: plist(1:ptot)     ! List of active particle ids
  logical, optional :: typemask(1:ntypes)   ! Part. types to include?

  integer :: p                              ! Particle counter

  nptcls = 0

  ! --------------------------------------------------------------------------
  do p=1,ptot

     ! If mask is present, include only particles of that type
     if (present(typemask)) then
        if (typemask(sph(p)%ptype)) then
           nptcls = nptcls + 1
           plist(nptcls) = p
        end if

     ! If no mask, include all particles
     else 
        nptcls = nptcls + 1
        plist(nptcls) = p
     end if

  end do
  ! --------------------------------------------------------------------------

  return
END SUBROUTINE create_particle_list
