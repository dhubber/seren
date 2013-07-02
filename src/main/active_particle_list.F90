! ACTIVE_PARTICLE_LIST.F90
! D. A. Hubber & A. McLeod - 15/03/2011
! Create a list of the ids of all active SPH particles.  Can (optionally) 
! create a list only from particle types set in the inputted mask variable.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE active_particle_list(acctot,acclist,typemask)
  use type_module
  use type_module
  use time_module
  use particle_module, only : ptot,sph
  implicit none

  integer, intent(out) :: acctot              ! No. of active particles
  integer, intent(out) :: acclist(1:ptot)     ! List of active particle ids
  logical, optional :: typemask(1:ntypes)     ! part. types to include?

  integer :: p                                ! particle counter

  debug2("Create list of active SPH particles [active_particle_list.F90]")

  acctot = 0

  ! --------------------------------------------------------------------------
  do p=1,ptot

     ! If mask is present, check particle type is selected as well as accdo.
     if (present(typemask)) then
        if (typemask(sph(p)%ptype) .and. sph(p)%accdo) then
           acctot = acctot + 1
           acclist(acctot) = p
        end if

     ! Else if no mask, consider all active particles of all types.
     else if (sph(p)%accdo) then
        acctot = acctot + 1
        acclist(acctot) = p
     end if

  end do
  ! --------------------------------------------------------------------------

  return
END SUBROUTINE active_particle_list
