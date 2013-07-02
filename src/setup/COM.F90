! COM.F90
! D. A. Hubber & A. McLeod - 24/2/2008
! Calculate the position and velocity of the centre of mass at the 
! beginning of the simulation.  If required, can transform the particle 
! data to the centre of mass frame (logical flag in params.dat).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE COM
  use particle_module
  use diagnostics_module
  use filename_module, only : restart
#if defined(SINKS)
  use sink_module
#endif
#if defined(USE_MPI)
  use mpi_communication_module
  use mpi
#endif
  implicit none

  integer :: p             ! Particle counter
#if defined(SINKS) 
  integer :: s             ! Sink counter
#endif
#if defined(USE_MPI)
  integer :: ierr          ! MPI error value
#if !defined(PERIODIC)
  integer :: d             ! Domain counter
#endif
#endif

  debug2("Calculating the centre of mass of the system [COM.F90]")
  
  rcom0(1:NDIM) = 0.0_DP
  vcom0(1:VDIM) = 0.0_DP
  mtot0 = 0.0_DP

! Calculate position and velocity of centre of mass
  do p=1,ptot
     rcom0(1:NDIM) = rcom0(1:NDIM) + real(sph(p)%m*sph(p)%r(1:NDIM),DP)
     vcom0(1:VDIM) = vcom0(1:VDIM) + real(sph(p)%m*sph(p)%v(1:VDIM),DP)
     mtot0 = mtot0 + real(sph(p)%m,DP)
  end do

! Include sinks
#if defined(SINKS)
  do s=1,stot
     rcom0(1:NDIM) = rcom0(1:NDIM) + real(sink(s)%m*sink(s)%r(1:NDIM),DP)
     vcom0(1:VDIM) = vcom0(1:VDIM) + real(sink(s)%m*sink(s)%v(1:VDIM),DP)
     mtot0 = mtot0 + real(sink(s)%m,DP)
  end do
#endif

#if defined(USE_MPI)
  call MPI_ALLREDUCE(MPI_IN_PLACE, mtot0, 1, MPI_DOUBLE_PRECISION,&
     &MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, rcom0, NDIM, MPI_DOUBLE_PRECISION,&
     &MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, vcom0, VDIM, MPI_DOUBLE_PRECISION,&
     &MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

! Normalise rcom0 and vcom0
  rcom0(1:NDIM) = rcom0(1:NDIM) / mtot0
  vcom0(1:VDIM) = vcom0(1:VDIM) / mtot0
  rcom(1:NDIM) = rcom0(1:NDIM)
  vcom(1:VDIM) = vcom0(1:VDIM)

#if defined(DEBUG2)
  write(6,*) "rcom0 :",rcom0(1:NDIM)
  write(6,*) "vcom0 :",vcom0(1:VDIM)
  write(6,*) "mtot0, com_frame, restart : ",mtot0,com_frame,restart
#endif

! Only convert to centre of mass frame if com_frame flag is on, if this is 
! not a restarted run, and periodic boundary conditions are not employed.  
! ----------------------------------------------------------------------------
#if !defined(PERIODIC)
  if (com_frame  .and. (.not. restart)) then

#if defined(DEBUG2)
     write(6,*) "Changing to COM frame"
     write(6,*) "rcom0 :",rcom0(1:NDIM)
     write(6,*) "vcom0 :",vcom0(1:VDIM)
#endif

     do p=1,ptot
        sph(p)%r(1:NDIM) = sph(p)%r(1:NDIM) - real(rcom0(1:NDIM),PR)
        sph(p)%v(1:VDIM) = sph(p)%v(1:VDIM) - real(vcom0(1:VDIM),PR)
     end do
#if defined(SINKS)
     do s=1,stot
        sink(s)%r(1:NDIM) = &
             &sink(s)%r(1:NDIM) - real(rcom0(1:NDIM),PR)
        sink(s)%v(1:VDIM) = &
             &sink(s)%v(1:VDIM) - real(vcom0(1:VDIM),PR)
     end do
#endif

     rcom0(1:NDIM) = 0.0_DP
     vcom0(1:VDIM) = 0.0_DP
  
#if defined(USE_MPI)
     do d=0,lastrank
        domain_bbmin(1:NDIM,d) = &
             &domain_bbmin(1:NDIM,d) - real(rcom0(1:NDIM),PR)
        domain_bbmax(1:NDIM,d) = &
             &domain_bbmax(1:NDIM,d) - real(rcom0(1:NDIM),PR)
     end do
     local_min = domain_bbmin(1:NDIM,rank)
     local_max = domain_bbmax(1:NDIM,rank)
#endif

  end if
#endif
! ----------------------------------------------------------------------------

  return
END SUBROUTINE COM
