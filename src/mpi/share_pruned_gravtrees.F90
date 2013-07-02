! SHARE_PRUNED_GRAVTREES.F90
! A. McLeod - 21/11/2011
! Shares pruned gravity trees with all tasks
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE share_pruned_gravtrees()
  use mpi_communication_module
  use mpi
  implicit none

  integer :: i                     ! loop counter
  integer :: d                     ! domain counter
  integer :: ierr                  ! Return code

  debug_timing("SHARE_LOCAL_GRAVTREES")

#ifdef DEBUG_EXPORT
    write (6,*) "Rank ",rank," posting receives for pruned trees"
#endif

  i = 0

  do d=0,lastrank

    if (d==rank) cycle
#ifdef DEBUG_EXPORT
    write (6,*) "Rank ",rank," posting a receive for task ", d, "; i = ", i
#endif
    ! Post the receive - until this is posted the receive does not actually occur
    ! and may or may not be buffered
    ! This buffer is normally larger than required, but this is OK
    call MPI_IRECV(BHremote_grav(0,d),cmax_remotegrav,MPI_REMOTE_BH_GRAV,&
         &d,GRAVTREES_TAG,MPI_COMM_WORLD,treerecv(i),ierr)
    i = i + 1

  end do

#ifdef DEBUG_EXPORT
    write (6,*) "Rank ",rank," posting sends for pruned trees"
#endif

  i = 0

  do d=0,lastrank
    if (d==rank) cycle

#ifdef DEBUG_EXPORT
    write (6,*) "Rank ",rank," sending pruned gravity tree of ",ctot_localgrav+1,&
    &" cells to task ", d
#endif

    call MPI_ISEND(BHlocal_grav(0), ctot_localgrav+1, MPI_REMOTE_BH_GRAV, &
         & d, GRAVTREES_TAG, MPI_COMM_WORLD,treesend(i),ierr)
    i = i + 1
  end do

  return
END SUBROUTINE share_pruned_gravtrees
