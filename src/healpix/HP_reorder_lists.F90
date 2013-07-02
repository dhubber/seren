! UV_REORDER_LISTS.F90
! D. A. Hubber - 19/10/2009
! Re-order distance lists for all HP sources.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE HP_reorder_lists(newid)
  use particle_module
  use HP_module
  implicit none

  integer, intent(in) :: newid(1:ptot)   ! List of new particle ids

  integer :: i                           ! HP source counter
  integer :: j                           ! Aux. particle counter
  integer :: nlist                       ! Length of list
  integer :: p                           ! Particle id
  integer :: paux                        ! Aux. particle id

  debug2("Reorder HEALPix distance lists [HP_reorder_lists.F90]")


! Loop over all active HP sources
! ----------------------------------------------------------------------------
  do i=1,HPtot

     nlist = 0

     ! Now loop over list in order and rename particles with new ids.  
     ! Any dead particles (p = -1) are ignored and the list is shortened.
     ! -----------------------------------------------------------------------
     do j=1,ptot
        paux = HPsource(i)%distorder(j)
        p = newid(paux)
        if (p /= -1) then
           nlist = nlist + 1
           HPsource(i)%distorder(nlist) = p
        end if
     end do
     ! -----------------------------------------------------------------------

     HPsource(i)%Nlist = nlist

  end do
! ----------------------------------------------------------------------------


  return
END SUBROUTINE HP_reorder_lists

