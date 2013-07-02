! MPI_SETUP_TYPES.F90
! A. McLeod - 16/09/08
! Sets up derived types for MPI communication, then commits them.
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine mpi_setup_types
  use mpi
  use mpi_communication_module
  use particle_module
#if defined(SINKS)
  use sink_module
#endif
  use definitions
  implicit none

  character (len=10), parameter :: frm1="(A28,I4,A)"
  integer               :: ierr                ! Return code
                                               ! Number of data, type of data
  integer (kind=MPI_ADDRESS_KIND) :: lb, ub    ! Lower bound, upper bound
  integer                         :: extent    ! Extent of data
  integer (kind=ILP)      :: testintegerILP(1:2)
  type(synchronisetype)   :: testsynchronise(1:2)
  type(sph_particle)      :: testparticle(1:2)
  type(ghosttype)         :: testghosttype(1:2)
! #ifdef HYDRO
!   type(exporthydrotype)   :: testexporthydro(1:2)
!   type(returnhydrotype)   :: testreturnhydro(1:2)
! #endif
#if defined(SELF_GRAVITY)
  type(exportgravitytype) :: testexportgravity(1:2)
  type(returngravitytype) :: testreturngravity(1:2)
#ifdef BH_TREE
  type(BHgrav_node)       :: testBHgrav_node(1:2)
#endif
#endif
#ifdef SINKS
  type(sink_node)         :: testsink_node(1:2)
!   type(sendsink)          :: testsendsink(1:2)
#endif
#ifdef CHECK_NEIGHBOUR_TIMESTEPS
  type(particletimesteptype) :: testparticletimestep(1:2)
#endif

  ! Commit the MPI type MPI_INTEGER_ILP for ILP integer
  ! --------------------------------------------------------------------------

  ! Just in case MPI_INTEGER8 is not supported, do it manually
  call MPI_GET_ADDRESS(testintegerILP(1), lb, ierr)
  call MPI_GET_ADDRESS(testintegerILP(2), ub, ierr)
  extent = int(ub - lb)

  call MPI_TYPE_CONTIGUOUS(extent,MPI_BYTE,MPI_INTEGER_ILP,ierr)
  call MPI_TYPE_COMMIT(MPI_INTEGER_ILP, ierr)
#ifdef DEBUG1
  if (rank==0) write (6,frm1) "MPI_INTEGER_ILP size: ", extent, " bytes"
#endif

  ! Commit the MPI type MPI_SYNCHRONISE for type synchronisetype
  ! --------------------------------------------------------------------------

  ! Check addresses to find correct offsets
  call MPI_GET_ADDRESS(testsynchronise(1), lb, ierr)
  call MPI_GET_ADDRESS(testsynchronise(2), ub, ierr)
  extent = int(ub - lb)

  ! Commit type
  call MPI_TYPE_CONTIGUOUS(extent, MPI_BYTE, MPI_SYNCHRONISE, ierr)
  call MPI_TYPE_COMMIT(MPI_SYNCHRONISE, ierr)
#ifdef DEBUG1
  if (rank==0) write (6,frm1) "MPI_SYNCHRONISE size: ", extent, " bytes"
#endif


  ! Commit the MPI particle type MPI_PARTICLE for type sph_particle
  ! --------------------------------------------------------------------------

  ! Check addresses to find correct offsets
  call MPI_GET_ADDRESS(testparticle(1), lb, ierr)
  call MPI_GET_ADDRESS(testparticle(2), ub, ierr)
  extent = int(ub - lb)

  ! Commit type
  call MPI_TYPE_CONTIGUOUS(extent, MPI_BYTE, MPI_PARTICLE, ierr)
  call MPI_TYPE_COMMIT(MPI_PARTICLE, ierr)
#ifdef DEBUG1
  if (rank==0) write (6,frm1) "MPI_PARTICLE size: ", extent, " bytes"
#endif
  mpi_particle_extent = int(extent,ILP)
  max_message_size = huge(MAX_CHUNKS) / (mpi_particle_extent * safety_factor)


  ! Commit the MPI particle type MPI_GHOST_TYPE for type sph_particle
  ! --------------------------------------------------------------------------

  ! Check addresses to find correct offsets
  call MPI_GET_ADDRESS(testghosttype(1), lb, ierr)
  call MPI_GET_ADDRESS(testghosttype(2), ub, ierr)
  extent = int(ub - lb)

  ! Commit type
  call MPI_TYPE_CONTIGUOUS(extent, MPI_BYTE, MPI_GHOST_TYPE, ierr)
  call MPI_TYPE_COMMIT(MPI_GHOST_TYPE, ierr)
#ifdef DEBUG1
  if (rank==0) write (6,frm1) "MPI_GHOST_TYPE size: ", extent, " bytes"
#endif


! #ifdef HYDRO
!   ! Commit the MPI particle type MPI_HYDROTYPE for type exporthydrotype
!   ! --------------------------------------------------------------------------
! 
!   ! Check addresses to find correct offsets
!   call MPI_GET_ADDRESS(testexporthydro(1), lb, ierr)
!   call MPI_GET_ADDRESS(testexporthydro(2), ub, ierr)
!   extent = int(ub - lb)
!   
!   ! Commit type
!   call MPI_TYPE_CONTIGUOUS(extent, MPI_BYTE, MPI_HYDROTYPE, ierr)
!   call MPI_TYPE_COMMIT(MPI_HYDROTYPE, ierr)
! #ifdef DEBUG1
!   if (rank==0) write (6,frm1) "MPI_HYDROTYPE size: ", extent, " bytes"
! #endif
! 
! 
!   ! Commit the MPI particle type MPI_RETURNHYDRO for type returnhydrotype
!   ! --------------------------------------------------------------------------
! 
!   ! Check addresses to find correct offsets
!   call MPI_GET_ADDRESS(testreturnhydro(1), lb, ierr)
!   call MPI_GET_ADDRESS(testreturnhydro(2), ub, ierr)
!   extent = int(ub - lb)
! 
!   ! Commit type
!   call MPI_TYPE_CONTIGUOUS(extent, MPI_BYTE, MPI_RETURNHYDRO, ierr)
!   call MPI_TYPE_COMMIT(MPI_RETURNHYDRO, ierr)
! #ifdef DEBUG1
!   if (rank==0) write (6,frm1) "MPI_RETURNHYDRO size: ", extent, " bytes"
! #endif
! #endif


#if defined(SELF_GRAVITY)
  ! Commit the MPI particle type MPI_GRAVITYTYPE for type exportgravitytype
  ! ----------------------------------------------------------------------------

  ! Check addresses to find correct offsets
  call MPI_GET_ADDRESS(testexportgravity(1), lb, ierr)
  call MPI_GET_ADDRESS(testexportgravity(2), ub, ierr)
  extent = int(ub - lb)

  ! Commit type
  call MPI_TYPE_CONTIGUOUS(extent, MPI_BYTE, MPI_GRAVITYTYPE, ierr)
  call MPI_TYPE_COMMIT(MPI_GRAVITYTYPE, ierr)
#ifdef DEBUG1
  if (rank==0) write (6,frm1) "MPI_GRAVITYTYPE size: ", extent, " bytes"
#endif


  ! Commit the MPI particle type MPI_RETURNGRAVITY for type returngravitytype
  ! --------------------------------------------------------------------------

  ! Check addresses to find correct offsets
  call MPI_GET_ADDRESS(testreturngravity(1), lb, ierr)
  call MPI_GET_ADDRESS(testreturngravity(2), ub, ierr)
  extent = int(ub - lb)

  ! Commit type
  call MPI_TYPE_CONTIGUOUS(extent, MPI_BYTE, MPI_RETURNGRAVITY, ierr)
  call MPI_TYPE_COMMIT(MPI_RETURNGRAVITY, ierr)
#ifdef DEBUG1
  if (rank==0) write (6,frm1) "MPI_RETURNGRAVITY size: ", extent, " bytes"
#endif


#ifdef BH_TREE
  ! Commit the MPI particle type MPI_REMOTE_BH_GRAV for type BHgrav_node
  ! ----------------------------------------------------------------------------
  if (remotetreedepth > LMAX) then
    write (6,*) "remotetreedepth > LMAX!!"
    stop
  end if

  call MPI_GET_ADDRESS(testBHgrav_node(1), lb, ierr)
  call MPI_GET_ADDRESS(testBHgrav_node(2), ub, ierr)
  extent = int(ub - lb)

  ! Commit type
  call MPI_TYPE_CONTIGUOUS(extent, MPI_BYTE, MPI_REMOTE_BH_GRAV, ierr)
  call MPI_TYPE_COMMIT(MPI_REMOTE_BH_GRAV, ierr)
#ifdef DEBUG1
  if (rank==0) write (6,frm1) "MPI_REMOTE_BH_GRAV size: ", extent, " bytes"
#endif

#endif
  ! End #ifdef BH_TREE

#endif
  ! End #ifdef SELF_GRAVITY

#ifdef SINKS
  ! Commit the MPI particle type MPI_SINK_NODE for type sink_node
  ! --------------------------------------------------------------------------

  call MPI_GET_ADDRESS(testsink_node(1), lb, ierr)
  call MPI_GET_ADDRESS(testsink_node(2), ub, ierr)
  extent = int(ub - lb)

  ! Commit type
  call MPI_TYPE_CONTIGUOUS(extent, MPI_BYTE, MPI_SINK_NODE, ierr)
  call MPI_TYPE_COMMIT(MPI_SINK_NODE, ierr)
#ifdef DEBUG1
  if (rank==0) write (6,frm1) "MPI_SINK_NODE size: ", extent, " bytes"
#endif

!   ! Commit the MPI particle type MPI_SEND_SINK for type sendsink
!   ! --------------------------------------------------------------------------
!   call MPI_GET_ADDRESS(testsendsink(1), lb, ierr)
!   call MPI_GET_ADDRESS(testsendsink(2), ub, ierr)
!   extent = int(ub - lb)
! 
!   ! Commit type
!   call MPI_TYPE_CONTIGUOUS(extent, MPI_BYTE, MPI_SEND_SINK, ierr)
!   call MPI_TYPE_COMMIT(MPI_SEND_SINK, ierr)
! #ifdef DEBUG1
!   if (rank==0) write (6,frm1) "MPI_SEND_SINK size: ", extent, " bytes"
! #endif

#endif


#ifdef CHECK_NEIGHBOUR_TIMESTEPS
  ! Commit the MPI particle type MPI_PARTICLETIMESTEP for type particletimesteptype
  ! --------------------------------------------------------------------------
  ! Check addresses to find correct offsets
  call MPI_GET_ADDRESS(testparticletimestep(1), lb, ierr)
  call MPI_GET_ADDRESS(testparticletimestep(2), ub, ierr)
  extent = int(ub - lb)

  ! Commit type
  call MPI_TYPE_CONTIGUOUS(extent, MPI_BYTE, MPI_PARTICLETIMESTEP, ierr)
  call MPI_TYPE_COMMIT(MPI_PARTICLETIMESTEP, ierr)
#ifdef DEBUG1
  if (rank==0) write (6,frm1) "MPI_PARTICLETIMESTEP size: ", extent, " bytes"
#endif

#endif

  return
end subroutine mpi_setup_types
