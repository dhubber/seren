! WRITE_DATA_SEREN_UNFORM.F90
! D. A. Hubber & A. McLeod - 25/7/2008; 25/6/2010
! Writes snapshot to file using native Seren format.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_data_seren_unform(out_file)
  use particle_module
  use hydro_module
  use time_module
  use scaling_module
  use type_module
  use sink_module
  use neighbour_module
#if defined(USE_MPI)
  use mpi_communication_module
#endif
  implicit none

  character(len=*), intent(in) :: out_file     ! Name of output file 

  character(len=20) :: format_id               ! File format (for verification)
  character(len=20) :: data_id(1:50)           ! Char ids of arrays written
  character(len=20) :: unit_data(1:21)         ! Unit data
  integer(kind=ILP) :: ilpdata(1:50)           ! Long integer data
  integer :: ndata                             ! Number of arrays written
  integer :: nunit                             ! Number of units
  integer :: p                                 ! Particle counter
  integer :: idata(1:50)                       ! Integer data
  integer :: typedata(1:5,1:50)                ! type data header array
  integer, allocatable :: idummy(:)            ! ..
  real(kind=DP) :: dpdata(1:50)                ! ..
  real(kind=PR) :: rdata(1:50)                 ! Real data array
  real(kind=PR), allocatable :: rdummy1(:)     ! real dummy array
  real(kind=PR), allocatable :: rdummy3(:,:)   ! real vector dummy array
#if defined(SINKS)
  integer, parameter :: sink_data_length=12+NDIM+VDIM+2*DMDT_RANGE
  logical :: ldummy(1:2)                       ! Logical dummy array
  integer :: idummy2(1:2)                      ! Integer dummy array 2
  integer :: s                                 ! Sink counter
  real(kind=PR) :: raux(1:sink_data_length)    ! Aux. variable
#endif

  debug2("Writing snapshot to file [write_data_seren.F90]")

! Append .sf extension indicating seren format and open file
  open(1, file=out_file, status="unknown", form="unformatted")
#if defined(DEBUG1)
#if defined(USE_MPI)
  if (rank==0) write(6,*) "Snapshot file :",trim(out_file)," (SEREN snapshot)"
#else
  write(6,*) "Snapshot file :",trim(out_file)," (SEREN snapshot)"
#endif
#endif

! Allocate main dummy arrays and zero them
  allocate(idummy(1:ptot))
  allocate(rdummy1(1:ptot))
  idata(1:50)     = 0
  ilpdata(1:50)   = 0_ILP
  dpdata(1:50)    = 0.0_DP
  rdata(1:50)     = 0.0_PR
  unit_data(1:21) = ''
  data_id(1:50)   = ''
  ndata           = 0

! Set unit character ids
  if (dimensionless) then
     nunit = 0
  else
     unit_data(1)  = runit
     unit_data(2)  = munit
     unit_data(3)  = tunit
     unit_data(4)  = vunit
     unit_data(5)  = aunit
     unit_data(6)  = rhounit
     unit_data(7)  = sigmaunit
     unit_data(8)  = Punit
     unit_data(9)  = funit
     unit_data(10) = Eunit
     unit_data(11) = momunit
     unit_data(12) = angmomunit
     unit_data(13) = angvelunit
     unit_data(14) = dmdtunit
     unit_data(15) = Lunit
     unit_data(16) = kappaunit
     unit_data(17) = Bunit
     unit_data(18) = Qunit
     unit_data(19) = Junit
     unit_data(20) = uunit
     unit_data(21) = tempunit
     nunit = 21
  end if


! Set array ids and array information data if there are any SPH particles
! ----------------------------------------------------------------------------
  if (ptot > 0) then
     ndata = ndata + 1;    data_id(ndata) = 'porig'
     typedata(1:5,ndata) = (/1,1,ptot,2,0/)
     
     ndata = ndata + 1;    data_id(ndata) = 'r'
     typedata(1:5,ndata) = (/NDIM,1,ptot,4,1/)
     
     ndata = ndata + 1;    data_id(ndata) = 'm'
     typedata(1:5,ndata) = (/1,1,ptot,4,2/)
     
     ndata = ndata + 1;    data_id(ndata) = 'h'
     typedata(1:5,ndata) = (/1,1,ptot,4,1/)
     
     ndata = ndata + 1;    data_id(ndata) = 'v' 
     typedata(1:5,ndata) = (/VDIM,1,ptot,4,4/)

#if defined(MHD)
     ndata = ndata + 1;    data_id(ndata) = 'B' 
     typedata(1:5,ndata) = (/BDIM,1,ptot,4,4/)
#endif
     
     ndata = ndata + 1;    data_id(ndata) = 'rho'
     typedata(1:5,ndata) = (/1,1,ptot,4,6/)
#if defined(HYDRO)
     ndata = ndata + 1;    data_id(ndata) = 'temp'
     typedata(1:5,ndata) = (/1,1,ptot,4,21/)

     ndata = ndata + 1;    data_id(ndata) = 'u'
     typedata(1:5,ndata) = (/1,1,ptot,4,20/)
#endif
#if defined(DEBUG_WRITE_MPI_TASK)
     ndata = ndata + 1;    data_id(ndata) = 'mpitask'
     typedata(1:5,ndata) = (/1,1,ptot,2,0/)
#endif

     if (dimensionless) typedata(5,:) = 0
  end if

#if defined(SINKS)
#if defined(USE_MPI)
  if (stot > 0 .AND. rank == 0) then
#else
  if (stot > 0) then
#endif
     ndata = ndata + 1;    data_id(ndata) = 'sink_v1'
        typedata(1:5,ndata) = (/0,1,stot,7,0/)
  end if
#endif

! Set important header information
  idata(1)    = ptot
  idata(2)    = stot
  idata(3)    = pboundary
  idata(4)    = picm
  idata(5)    = pgas
  idata(6)    = pcdm
  idata(7)    = pdust
  idata(8)    = pion
  idata(20)   = nunit
  idata(21)   = ndata
  idata(30)   = DMDT_RANGE
  idata(31)   = pgas_orig
  idata(32)   = pp_gather
#if defined(USE_MPI)
  idata(40)   = rank
  idata(41)   = numtasks
#else
  idata(40)   = 0            ! MPI rank, zero for non-MPI
  idata(41)   = 0            ! MPI numtasks, zero for non-MPI
#endif
  ilpdata(1)  = snapshot
  ilpdata(2)  = nsteps
  ilpdata(3)  = ntempnext
  ilpdata(4)  = ndiagnext
  ilpdata(5)  = nsnapnext
  ilpdata(6)  = nsinknext
  rdata(1)    = h_fac
  rdata(2)    = gamma
  rdata(3)    = mu_bar
  rdata(4)    = hmin*rscale
  dpdata(1)   = time*tscale
  dpdata(2)   = lastsnap*tscale
  dpdata(3)   = mgas_orig*mscale

! Write information identifying format and precision of file
  format_id = 'SERENBINARYDUMPV2'
  write(1) format_id
  write(1) PR
  write(1) NDIM
  write(1) VDIM
  write(1) BDIM

! Write header information to file
  write(1) idata
  write(1) ilpdata
  write(1) rdata
  write(1) dpdata
  if (nunit > 0) write(1) unit_data(1:nunit)
  if (ndata > 0) write(1) data_id(1:ndata)
  if (ndata > 0) write(1) typedata(1:5,1:ndata)


! Write arrays for SPH particles
! ============================================================================
  if (ptot > 0) then

     ! Original ids
     ! -----------------------------------------------------------------------
     do p=1,ptot
        idummy(p) = sph(p)%porig
     end do
     write(1) idummy

     ! Positions
     ! -----------------------------------------------------------------------
     allocate(rdummy3(1:NDIM,1:ptot))
     do p=1,ptot
        rdummy3(1:NDIM,p) = sph(p)%r(1:NDIM)*real(rscale,PR)
#if defined(USE_MPI) && defined(PERIODIC)
        call unwrap_particle_position(rdummy3(1:NDIM,p))
#endif
     end do
     write(1) rdummy3
     deallocate(rdummy3)
     
     ! Mass
     ! -----------------------------------------------------------------------
     do p=1,ptot
        rdummy1(p) = sph(p)%m*real(mscale,PR)
     end do
     write(1) rdummy1
     
     ! Smoothing length
     ! -----------------------------------------------------------------------
     do p=1,ptot
        rdummy1(p) = sph(p)%h*real(rscale,PR)
     end do
     write(1) rdummy1
     
     ! Velocities
     ! -----------------------------------------------------------------------
     allocate(rdummy3(1:VDIM,1:ptot))
     do p=1,ptot
        rdummy3(1:VDIM,p) = sph(p)%v(1:VDIM)*real(vscale,PR)
     end do
     write(1) rdummy3
     deallocate(rdummy3)
     
#if defined(MHD)
     ! Magnetic field
     ! -----------------------------------------------------------------------
     allocate(rdummy3(1:BDIM,1:ptot))
     do p=1,ptot
        rdummy3(1:BDIM,p) = sph(p)%B(1:BDIM)*real(Bscale,PR)
     end do
     write(1) rdummy3
     deallocate(rdummy3)
#endif
     
     ! Density
     ! -----------------------------------------------------------------------
     do p=1,ptot
        rdummy1(p) = sph(p)%rho*real(rhoscale,PR)
     end do
     write(1) rdummy1
     
     ! Temperature
     ! -----------------------------------------------------------------------
#if defined(HYDRO)
     do p=1,ptot
        rdummy1(p) = sph(p)%temp
     end do
     write(1) rdummy1
     
     ! Internal energy
     ! -----------------------------------------------------------------------
     do p=1,ptot
        rdummy1(p) = sph(p)%u*real(uscale,PR)
     end do
     write(1) rdummy1
#endif

#if defined(DEBUG_WRITE_MPI_TASK)
     ! MPI tasks
     ! -----------------------------------------------------------------------
#if defined(USE_MPI)
     idummy(1:ptot) = rank
#else
     idummy(1:ptot) = 0
#endif
     write(1) idummy
#endif
     
  end if
! ============================================================================


! Sinks
! ----------------------------------------------------------------------------
#if defined(SINKS)
#if defined(USE_MPI)
  if (stot > 0 .AND. rank == 0) then
#else
  if (stot > 0) then
#endif
     write(1) 2,2,0,sink_data_length,0,0
     do s=1,stot
        ldummy(1)                     = sink(s)%accrete
        ldummy(2)                     = sink(s)%static
        idummy2(1)                    = sink(s)%id
        idummy2(2)                    = sink(s)%ncreate
        raux(1:sink_data_length)      = 0.0_PR
        raux(1)                       = real(sink(s)%tcreate*tscale,PR)
        raux(2:NDIM+1)                = sink(s)%r(1:NDIM)*real(rscale,PR)
        raux(NDIM+2:NDIM+VDIM+1)      = sink(s)%v(1:NDIM)*real(vscale,PR)
        raux(NDIM+VDIM+2)             = sink(s)%m*real(mscale,PR)
        raux(NDIM+VDIM+3)             = sink(s)%h*real(rscale,PR)
        raux(NDIM+VDIM+4)             = sink(s)%radius*real(rscale,PR)
        raux(NDIM+VDIM+5:NDIM+VDIM+7) = real(sink(s)%angmom(1:3)*angmomscale,PR)
        raux(NDIM+VDIM+8)             = real(sink(s)%dmdt*dmdtscale,PR)
        raux(NDIM+VDIM+9)             = real(sink(s)%star_radius*rscale,PR)
        raux(NDIM+VDIM+10)            = real(sink(s)%luminosity*Lscale,PR)
        raux(NDIM+VDIM+11)            = real(sink(s)%temperature,PR)
        if (DMDT_RANGE > 0) then
           raux(NDIM+VDIM+12:NDIM+VDIM+11+DMDT_RANGE) = &
              & real(sink(s)%macc(1:DMDT_RANGE)*mscale,PR)
           raux(NDIM+VDIM+12+DMDT_RANGE:NDIM+VDIM+11+2*DMDT_RANGE) = &
              & real(sink(s)%tacc(1:DMDT_RANGE)*tscale,PR)
        end if
        raux(NDIM+VDIM+12+2*DMDT_RANGE) = real(sink(s)%mmax*mscale,PR)
        write(1) ldummy
        write(1) idummy2
        write(1) raux
     end do
  end if
#endif


! Close file once finished
! ----------------------------------------------------------------------------
  close(1)

  deallocate(rdummy1)
  deallocate(idummy)

  return
END SUBROUTINE write_data_seren_unform
