! WRITE_DATA_SEREN_FORM.F90
! D. A. Hubber & A. McLeod - 25/7/2008; 25/6/2010
! Writes snapshot to file using native Seren format. 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_data_seren_form(out_file)
  use particle_module
  use hydro_module
  use neighbour_module
  use time_module
  use scaling_module
  use type_module
  use sink_module
#if defined(USE_MPI)
  use mpi_communication_module
#endif
  implicit none

  character(len=*), intent(in) :: out_file   ! Name of output file 
 
  character(len=20) :: data_id(1:50)         ! Char ids of arrays written
  character(len=20) :: unit_data(1:50)       ! Unit data
  integer :: i                               ! Aux. loop counter
  integer :: ndata                           ! Number of arrays written
  integer :: nunit                           ! Number of units
  integer :: p                               ! Particle counter
  integer :: idata(1:50)                     ! Integer data
  integer(kind=ILP) :: ilpdata(1:50)         ! Long integer data
  real(kind=DP) :: dpdata(1:50)              ! Double precision data array
  real(kind=PR) :: rdata(1:50)               ! Real data array
  integer :: typedata(1:5,1:50)              ! type data header array
  real(kind=PR) :: rp(1:NDIM)                ! Store for particle positions
#if defined(SINKS)
  integer, parameter :: sink_data_length=12+NDIM+VDIM+2*DMDT_RANGE
  character(len=30)  :: sink_format_string
  integer :: s                               ! Sink counter
#endif

  debug2("Writing snapshot to file [write_data_seren.F90]")

! Append .sf extension indicating seren format and open file
  open(1, file=out_file, status="unknown", form="formatted")
#if defined(DEBUG1)
#if defined(USE_MPI)
  if (rank==0) write(6,*) "Snapshot file : ",trim(out_file),&
                         &"   (SEREN snapshot)"
#else
  write(6,*) "Snapshot file : ",trim(out_file),"   (SEREN snapshot)"
#endif
#endif

#if defined(SINKS)
  write (sink_format_string,'(A,I0,A)') "(", sink_data_length, "E18.10)"
#endif

! Zero main arrays
  idata(1:50)     = 0
  ilpdata(1:50)   = 0_ILP
  dpdata(1:50)    = 0.0_DP
  rdata(1:50)     = 0.0_PR
  unit_data(1:50) = ''
  data_id(1:50)   = ''
  ndata           = 0

! Set unit character ids
#if defined(DIMENSIONLESS)
  nunit = 0
#else
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
#endif

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
     
     ndata = ndata + 1;    data_id(ndata) = 'rho'
     typedata(1:5,ndata) = (/1,1,ptot,4,6/)
#if defined(HYDRO)
     ndata = ndata + 1;    data_id(ndata) = 'temp'
     typedata(1:5,ndata) = (/1,1,ptot,4,21/)

     ndata = ndata + 1;    data_id(ndata) = 'u'
     typedata(1:5,ndata) = (/1,1,ptot,4,20/)
#endif
#if defined(ENTROPIC_FUNCTION)
     ndata = ndata + 1;    data_id(ndata) = 'Aent'
     typedata(1:5,ndata) = (/1,1,ptot,4,21/)
#endif
#if defined(DEBUG_WRITE_MPI_TASK)
     ndata = ndata + 1;    data_id(ndata) = 'mpitask'
     typedata(1:5,ndata) = (/1,1,ptot,2,0/)
#endif

#if defined(DIMENSIONLESS)
     typedata(5,:) = 0
#endif
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
  write(1,'(A16)') "SERENASCIIDUMPV2"
  write(1,'(I2)') PR
  write(1,'(I2)') NDIM
  write(1,'(I2)') VDIM
  write(1,'(I2)') BDIM

! Write the rest of the header information to file
  do i=1,50
     write(1,'(I10)') idata(i)
  end do
  do i=1,50
     write(1,'(I10)') ilpdata(i)
  end do
  do i=1,50
     write(1,'(E18.10)') rdata(i)
  end do
  do i=1,50
     write(1,'(E18.10)') dpdata(i)
  end do
  if (nunit > 0) then
     do i=1,nunit
        write(1,'(a)') unit_data(i)
     end do
  end if
  if (ndata > 0) then
     do i=1,ndata
        write(1,'(a)') data_id(i)
     end do
     do i=1,ndata
        write(1,'(5I10)') typedata(1:5,i)
     end do
  end if
  i = 0

! Write arrays for SPH particles
! ============================================================================
  if (ptot > 0) then

     ! Original particle ids
     ! -----------------------------------------------------------------------
     i = i + 1
     do p=1,ptot
        write(1,'(I10)') sph(p)%porig
     end do
     
     ! Positions
     ! -----------------------------------------------------------------------
     i = i + 1

     do p=1,ptot
        rp = sph(p)%r
#if defined(USE_MPI) && defined(PERIODIC)
        call unwrap_particle_position(rp)
#endif
#if NDIM==1
        write(1,'(E18.10)') rp(1)*rscale
#elif NDIM==2
        write(1,'(2E18.10)') rp(1)*rscale, rp(2)*rscale
#else
        write(1,'(3E18.10)') rp(1)*rscale, rp(2)*rscale, rp(3)*rscale
#endif
     end do

     ! Mass
     ! -----------------------------------------------------------------------
     i = i + 1
     do p=1,ptot
        write(1,'(E18.10)') sph(p)%m*mscale
     end do
     
     ! Smoothing length
     ! -----------------------------------------------------------------------
     i = i + 1
     do p=1,ptot
        write(1,'(E18.10)') sph(p)%h*rscale
     end do
     
     ! Velocities
     ! -----------------------------------------------------------------------
     i = i + 1
     do p=1,ptot
#if VDIM==1
        write(1,'(E18.10)') sph(p)%v(1)*vscale
#elif VDIM==2
        write(1,'(2E18.10)') sph(p)%v(1)*vscale, sph(p)%v(2)*vscale
#else
        write(1,'(3E18.10)') sph(p)%v(1)*vscale, sph(p)%v(2)*vscale, &
             &sph(p)%v(3)*vscale
#endif
     end do

     ! Density
     ! -----------------------------------------------------------------------
     i = i + 1
     do p=1,ptot
        write(1,'(E18.10)') sph(p)%rho*rhoscale
     end do

     ! Temperature
     ! -----------------------------------------------------------------------
#if defined(HYDRO)
     i = i + 1
     do p=1,ptot
        write(1,'(E18.10)') sph(p)%temp
     end do

     ! Internal energy
     ! -----------------------------------------------------------------------
     i = i + 1
     do p=1,ptot
        write(1,'(E18.10)') sph(p)%u*uscale
     end do
#endif

     ! Entropic function
     ! -----------------------------------------------------------------------
#if defined(ENTROPIC_FUNCTION)
     i = i + 1
     do p=1,ptot
        write(1,'(E18.10)') sph(p)%Aent
     end do
#endif

#if defined(DEBUG_WRITE_MPI_TASK)
     ! MPI tasks
     ! -----------------------------------------------------------------------
     i = i + 1
     do p=1,ptot
#if defined(USE_MPI)
        write(1,'(I10)') rank
#else
        write(1,'(I10)') 0
#endif
     end do
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
     i = i + 1
     write(1,'(6I10)') 2,2,0,sink_data_length,0,0
     do s=1,stot
        write(1,'(2L1)') sink(s)%accrete,sink(s)%static
        write(1,'(2I8)') sink(s)%id,sink(s)%ncreate
        write(1,sink_format_string) &
             &real(sink(s)%tcreate*tscale,PR),&
             &sink(s)%r(1:NDIM)*rscale,&
             &sink(s)%v(1:VDIM)*vscale,&
             &sink(s)%m*mscale,&
             &sink(s)%h*rscale,&
             &sink(s)%radius*rscale,&
             &sink(s)%angmom(1:3)*angmomscale,&
             &sink(s)%dmdt*dmdtscale,&
             &sink(s)%star_radius*rscale,&
             &sink(s)%luminosity*Lscale,&
             &sink(s)%temperature,&
             &sink(s)%macc(1:DMDT_RANGE)*mscale,&
             &sink(s)%tacc(1:DMDT_RANGE)*tscale,&
             &sink(s)%mmax*mscale
     end do
  end if
#endif


! Close file once finished
  close(1)

  return
END SUBROUTINE write_data_seren_form
