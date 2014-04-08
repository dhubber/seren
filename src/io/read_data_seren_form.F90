! READ_DATA_SEREN_FORM.F90
! D. A. Hubber & A. McLeod - 25/7/2008; 1/9/2010
! Read snapshot in Seren format.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_data_seren_form(out_file,decomp_read)
  use particle_module
  use hydro_module
  use time_module
  use type_module
  use sink_module
#if defined(USE_MPI)
  use mpi_communication_module
#endif
  implicit none

  character(len=*), intent(in) :: out_file   ! Name of output file
  logical, intent(in)          :: decomp_read ! MPI decomposition read

  character(len=20) :: unit_data(1:50)       ! Char ids of arrays written
  character(len=20) :: data_id(1:50)         ! Char ids of arrays written
  character(len=20) :: format_id             ! File format (for verification)
  logical :: ldata(1:50)                     ! ..
  integer :: dim_check                       ! Dimension check
  integer :: dmdt_range_aux                  ! Accretion history array size
  integer :: i                               ! Aux. counter
  integer :: itemp                           ! Aux. integer variable
  integer :: ndata                           ! No. of arrays written
  integer :: nvartype(1:6)                   ! ..
  integer :: nunit                           ! No. of unit strings
  integer :: p                               ! Particle counter
  integer :: pfirst                          ! ..
  integer :: plast                           ! ..
  integer :: idata(1:50)                     ! Integer data
  integer (kind=ILP) :: ilpdata(1:50)        ! Long integer data
  real(kind=DP) :: dpdata(1:50)              ! Double precision data array
  real(kind=PR) :: rdata(1:50)               ! Real data array
  real(kind=PR) :: rtemp(1:3)                ! parray temp variable
  integer :: typedata(1:5,1:50)              ! type data header array
#if defined(SINKS)
  integer, parameter :: sink_data_length=11+NDIM+VDIM+2*DMDT_RANGE
  character(len=30)  :: sink_format_string
  integer :: sink_data_length_aux            ! ..
  integer :: j                               ! Aux. counter
  integer :: idata3(1:6)                     ! ..
  !integer :: nl,ni,nli,npr,ndp,nchar         ! ..
  integer :: s                               ! Sink counter
  real(kind=PR) :: raux(1:sink_data_length)  ! Aux. time variable
#endif

  debug1("Reading in formatted data file : "//trim(out_file)//" [read_data_seren_form.F90]")

  open(1, file=out_file, status="old", form="formatted")
  write(6,*) "Snapshot file : ",trim(out_file),"   (SEREN ASCII snapshot)"

! Read information identifying format and precision of file.
! Then check if each value corresponds to current Seren values.
  read(1,*) format_id
  format_id = trim(adjustl(format_id))
  if (trim(adjustl(format_id)) /= "SERENASCIIDUMPV2" .AND.&
     &trim(adjustl(format_id)) /= "SERENASCIIDUMPV3") then
     stop 'Incorrect format of IC file'
  end if 
  read(1,*) dim_check
  read(1,*) dim_check
  if (dim_check /= NDIM) then
     stop 'Incorrect NDIM of IC file'
  end if
  read(1,*) dim_check
  if (dim_check /= VDIM) then
     stop 'Incorrect VDIM of IC file'
  end if
  read(1,*) dim_check
  if (dim_check /= BDIM) then
     stop 'Incorrect BDIM of IC file'
  end if

#if defined(SINKS)
  write (sink_format_string,'(A,I0,A)') "(", sink_data_length, "E18.10)"
#endif

! First, read in header information 
  do p=1,50
     read(1,'(I10)') idata(p)
  end do
  do p=1,50
     read(1,'(I10)') ilpdata(p)
  end do
  ptot           = idata(1)
  stot           = idata(2)
  pboundary      = idata(3)
  picm           = idata(4)
  pgas           = idata(5)
  pcdm           = idata(6)
  pdust          = idata(7)
  pion           = idata(8)
  nunit          = idata(20)
  ndata          = idata(21)
  dmdt_range_aux = idata(30)
  pgas_orig      = idata(31)
#if defined(USE_MPI)
  if (rank/=0) then
     ! Check rank and domain (but ignore if rank==0 as we could be doing
     ! decomposition or a single MPI task run
     if (idata(40) /= rank .OR. idata(41) /= numtasks) then
        write (6,*) "File for rank ", rank, " of ", numtasks
        write (6,*) "But this is rank ", rank, " of ", numtasks, "!"
        stop
     end if
  end if
#else
  if (idata(40) /= 0) write (6,*) "WARNING - MPI rank of ", idata(40), &
     & " detected, but this is standard Seren!"
  if (idata(41) /= 0) write (6,*) "WARNING - MPI numtasks of ", idata(41), &
     & " detected, but this is standard Seren!"
#endif
  snapshot       = ilpdata(1)
  nsteps         = ilpdata(2)
  ntempnext      = ilpdata(3)
  ndiagnext      = ilpdata(4)
  nsnapnext      = ilpdata(5)
  nsinknext      = ilpdata(6)

#if !defined(USE_MPI)
  write(6,*) "SPH Particles  = ", ptot,"    Sink Particles = ", stot
  write(6,*) "Gas            = ", pgas
  write(6,*) "Boundary       = ", pboundary
  write(6,*) "Intercloud     = ", picm
  write(6,*) "Dark matter    = ", pcdm
  write(6,*) "Dust           = ", pdust
  write(6,*) "Ions           = ", pion
#endif
  if (ptot /= (pgas + picm + pboundary + pcdm + pdust + pion)) &
       &stop "Fatal error: particles do not add up"

! Read in real information
  do p=1,50
     read(1,*) rdata(p)
  end do
  do p=1,50
     read(1,*) dpdata(p)
  end do
  time      = dpdata(1)
  lastsnap  = dpdata(2)
  mgas_orig = dpdata(3)

! Read in unit data
  do p=1,nunit
     read(1,'(20A)') unit_data(p)
  end do

! Read in ids of arrays contained in file
  do p=1,ndata
     read(1,*) data_id(p)
  end do

! Read in array of typedata information
  do p=1,ndata
     read(1,*) typedata(1:5,p)
  end do

! Allocate memory now we know ptot
! If MPI decomposition read, allocate 'minimal = TRUE'
  call allocate_memory(decomp_read)


! Now loop through array ids and read each array in turn
! ============================================================================
  do i=1,ndata

     ! Find pfirst, plast from typedata for this data set
     pfirst = typedata(2,i); plast = typedata(3,i)
     write(6,*) i,pfirst,plast,data_id(i)

     ! porig
     ! -----------------------------------------------------------------------
     if (data_id(i)=='porig') then
        do p=pfirst,plast
           read(1,*) itemp
           if (decomp_read) then
              minimal_sph(p)%porig = itemp
           else
              sph(p)%porig = itemp
           end if
        end do

     ! Positions
     ! -----------------------------------------------------------------------
     else if (data_id(i)=='r') then
        do p=pfirst,plast
#if NDIM==1
           read(1,*) rtemp(1)
#elif NDIM==2
           read(1,*) rtemp(1:2)
#elif NDIM==3
           read(1,*) rtemp(1:3)
#endif
           if (decomp_read) then
              minimal_sph(p)%r(1:NDIM) = real(rtemp(1:NDIM),PR)
           else
              sph(p)%r(1:NDIM) = real(rtemp(1:NDIM),PR)
           end if
        end do

     ! Masses
     ! ----------------------------------------------------------------------- 
     else if (data_id(i)=='m') then
        do p=pfirst,plast
           read(1,*) rtemp(1)
           if (decomp_read) then
              minimal_sph(p)%m = real(rtemp(1),PR)
           else
              sph(p)%m = real(rtemp(1),PR)
           end if
        end do

     ! Smoothing lengths
     ! ----------------------------------------------------------------------- 
     else if (data_id(i)=='h') then
        do p=pfirst,plast
           read(1,*) rtemp(1)
           if (decomp_read) then
              minimal_sph(p)%h = real(rtemp(1),PR)
           else
              sph(p)%h = real(rtemp(1),PR)
           end if
        end do

     ! Velocities
     ! -----------------------------------------------------------------------
     else if (data_id(i)=='v') then
        do p=pfirst,plast
#if VDIM==1
           read(1,*) rtemp(1)
#elif VDIM==2
           read(1,*) rtemp(1:2)
#elif VDIM==3
           read(1,*) rtemp(1:3)
#endif
           if (decomp_read) then
              minimal_sph(p)%v(1:VDIM) = real(rtemp(1:VDIM),PR)
           else
              sph(p)%v(1:VDIM) = real(rtemp(1:VDIM),PR)
           end if
        end do

     ! Density
     ! ----------------------------------------------------------------------- 
     else if (data_id(i)=='rho') then
        do p=pfirst,plast
           read(1,*) rtemp(1)
           if (decomp_read) then
              minimal_sph(p)%rho = real(rtemp(1),PR)
           else
              sph(p)%rho = real(rtemp(1),PR)
           end if
        end do

     ! Temperature
     ! ----------------------------------------------------------------------- 
     else if (data_id(i)=='temp') then
        do p=pfirst,plast
           read(1,*) rtemp(1)
#if defined(HYDRO)
           if (decomp_read) then
              minimal_sph(p)%temp = real(rtemp(1),PR)
           else
              sph(p)%temp = real(rtemp(1),PR)
           end if
#endif
        end do

     ! Internal energy
     ! ----------------------------------------------------------------------- 
     else if (data_id(i)=='u') then
        do p=pfirst,plast
           read(1,*) rtemp(1)
#if defined(HYDRO)
           if (decomp_read) then
              minimal_sph(p)%u = real(rtemp(1),PR)
           else
              sph(p)%u = real(rtemp(1),PR)
           end if
#endif
        end do

     ! Entropic function
     ! ----------------------------------------------------------------------- 
     else if (data_id(i)=='Aent') then
        do p=pfirst,plast
           read(1,*) rtemp(1)
#if defined(ENTROPIC_FUNCTION)
           if (decomp_read) then
              minimal_sph(p)%Aent = real(rtemp(1),PR)
           else
              sph(p)%Aent = real(rtemp(1),PR)
           end if
#endif
        end do

     ! B-field
     ! ----------------------------------------------------------------------- 
     else if (data_id(i)=='B') then
        do p=pfirst,plast
#if BDIM==1
           read(1,*) rtemp(1)
#elif BDIM==2
           read(1,*) rtemp(1:2)
#elif BDIM==3
           read(1,*) rtemp(1:3)
#endif
        end do

     ! Sinks
     ! ----------------------------------------------------------------------- 
     else if (data_id(i)=='sink_v1') then
#if defined(SINKS)
        read(1,'(6I19)') idata3(1:6)
        sink_data_length_aux = idata3(4)
        if (stot > 0) then
           do j=pfirst,plast
              s = j
              read(1,'(2L1)') ldata(1:2)
              read(1,'(2I8)') idata(1:2)
              read(1,sink_format_string) raux(1:sink_data_length)
              sink(s)%id          = idata(1)
              sink(s)%ncreate     = idata(2)
              sink(s)%accrete     = ldata(1)
              sink(s)%static      = ldata(2)              
              sink(s)%tcreate     = real(raux(1),DP)
              sink(s)%r(1:NDIM)   = raux(2:NDIM+1)
              sink(s)%v(1:NDIM)   = raux(NDIM+2:NDIM+VDIM+1)
              sink(s)%m           = raux(NDIM+VDIM+2)
              sink(s)%h           = raux(NDIM+VDIM+3)
              sink(s)%radius      = raux(NDIM+VDIM+4)
              sink(s)%angmom(1:3) = raux(NDIM+VDIM+5:NDIM+VDIM+7)
              sink(s)%dmdt        = raux(NDIM+VDIM+8)
              sink(s)%star_radius = raux(NDIM+VDIM+9)
              sink(s)%luminosity  = raux(NDIM+VDIM+10)
              sink(s)%temperature = raux(NDIM+VDIM+11)
              if (dmdt_range_aux > 0) then
                 dmdt_range_aux = min(dmdt_range_aux,DMDT_RANGE)
                 sink(s)%macc(1:dmdt_range_aux) = &
                      &real(raux(NDIM+VDIM+12:&
                      &NDIM+VDIM+11+dmdt_range_aux),DP)
                 sink(s)%tacc(1:dmdt_range_aux) = &
                      &real(raux(NDIM+VDIM+12+dmdt_range_aux:&
                      &NDIM+VDIM+11+2*dmdt_range_aux),DP)
              end if
              if (sink_data_length_aux > sink_data_length) &
                   &sink(s)%mmax = raux(NDIM+VDIM+12+2*dmdt_range_aux)
           end do
        end if
#else
        stop 'SEREN file contains sinks; not activated in code'
#endif

     ! Skip through arbitrary 1-D or 2-D array
     ! -----------------------------------------------------------------------
     else if (typedata(1,i) >= 1) then
        do p=pfirst,plast
           if (typedata(4,i)==1) read(1,'(999L1)') ldata(1:typedata(1,i))
           if (typedata(4,i)==2) read(1,*) idata(1:typedata(1,i))
           if (typedata(4,i)==3) read(1,*) ilpdata(1:typedata(1,i))
           if (typedata(4,i)==4) read(1,*) rdata(1:typedata(1,i))
           if (typedata(4,i)==5) read(1,*) dpdata(1:typedata(1,i))
        end do

     ! Skip through arbitrary data structure
     ! -----------------------------------------------------------------------
     else if (typedata(1,i) == 0 .and. typedata(4,i) == 7) then
        read(1,*) nvartype(1:6)
        do p=pfirst,plast
           if (nvartype(1) > 0) read(1,'(999L1)') ldata(1:nvartype(1))
           if (nvartype(2) > 0) read(1,*) idata(1:nvartype(2))
           if (nvartype(3) > 0) read(1,*) ilpdata(1:nvartype(3))
           if (nvartype(4) > 0) read(1,*) rdata(1:nvartype(4))
           if (nvartype(5) > 0) read(1,*) dpdata(1:nvartype(5))
        end do

     end if
     ! -----------------------------------------------------------------------

  end do
! ============================================================================

#if defined(RTSPH) && defined(HYDRO) && !defined(INTERNAL_ENERGY)
  temp(1:ptot) = 1.0_PR
#endif

! Close file once finished
  close(1)

#if defined(USE_MPI) && defined(SINKS)
! If this is MPI and we have sinks, load them now from the root task
! unless this is a root-only read (for decomposition)
  if (stot > 0) then
     sink(1:stot)%domain = 0 ! Assign all sinks to the root temporarily
     if (.NOT. decomp_read) call sink_share
  end if
#endif

  return
END SUBROUTINE read_data_seren_form
