! READ_DATA_SEREN_UNFORM.F90
! D. A. Hubber & A. McLeod - 25/7/2008; 25/6/2010
! Read binary snapshot in Seren format
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_data_seren_unform(out_file,decomp_read)
  use particle_module
  use hydro_module
  use time_module
  use type_module
  use sink_module
#if defined(USE_MPI)
  use mpi_communication_module
#endif
  implicit none

  character(len=*), intent(in) :: out_file    ! Name of output file
  logical, intent(in)          :: decomp_read ! MPI decomposition read

  character(len=20) :: unit_data(1:50)        ! Char ids of arrays written
  character(len=20) :: data_id(1:50)          ! Char ids of arrays written
  character(len=20) :: format_id              ! File format (for verification)
  logical, allocatable :: ldummy(:)           ! Logical dummy array
  integer :: dim_check                        ! Dimension check
  integer :: dmdt_range_aux                   ! Accretion history array size
  integer :: i                                ! Aux. loop counter
  integer :: idata(1:50)                      ! Integer data
  integer :: j                                ! Aux. particle id
  integer :: jtot                             ! No. of particles in array
  integer :: ndata                            ! Number of arrays written
  integer :: nunit                            ! Number of units
  integer :: nvartype(1:6)                    ! ..
  integer :: p                                ! Particle counter
  integer :: pr_check                         ! Precision check
  integer :: pfirst                           ! ..
  integer :: plast                            ! ..
  integer :: typedata(1:5,1:50)               ! type data header array
  integer, allocatable :: idummy(:)           ! ..
  integer(kind=ILP) :: ilpdata(1:50)          ! Long integer data
  integer, allocatable :: ilpdummy(:)         ! ..
  real(kind=DP) :: dpdata(1:50)               ! Double precision data array
  real(kind=PR) :: rdata(1:50)                ! Real data array
  real(kind=PR), allocatable :: rdummy1(:)    ! real dummy array
  real(kind=PR), allocatable :: rdummy2(:,:)  ! real vector dummy array
  real(kind=DP), allocatable :: dpdummy1(:)   ! real dummy array
  real(kind=DP), allocatable :: dpdummy2(:,:) ! real vector dummy array
#if defined(SINKS)
  integer :: dmdt_range_min                   ! ..
  integer :: sink_data_length                 ! ..
  integer :: sink_data_length_aux             ! ..
  integer :: idummy2(1:2)                     ! Integer dummy array 2
  integer :: idummy3(1:6)                     ! ..
  !integer :: nl,ni,nli,npr,ndp,nchar          ! ..
  integer :: s                                ! Sink counter
  real(kind=PR), allocatable :: raux(:)       ! Aux. variable
#endif
#if defined(MHD)
  logical :: has_B                            ! File has magnetic fields
#endif

  debug2("Reading snapshot [read_data_seren_unform.F90]")

  open(1, file=out_file, status="old", form="unformatted")
  write(6,*) "Snapshot file : ",trim(out_file),"   (SEREN binary snapshot)"

! Read information identifying format and precision of file
! Then check if each value corresponds to current Seren values
  read(1) format_id
  format_id = trim(adjustl(format_id))
  write(6,*) "Snapshot file :",trim(out_file),"  ",trim(format_id)
  if (trim(adjustl(format_id)) /= "SERENBINARYDUMPV2") then
     stop 'Incorrect format of IC file'
  end if 
  read(1) pr_check
  if (PR == SP .and. pr_check /= 1 .and. pr_check /=4) then 
     stop 'Incorrect PR of IC file'
  else if (PR == DP .and. pr_check /= 2 .and. pr_check /=8) then 
     stop 'Incorrect PR of IC file'
  end if
  read(1) dim_check
  if (dim_check /= NDIM) then
     stop 'Incorrect NDIM of IC file'
  end if
  read(1) dim_check
  if (dim_check /= VDIM) then
     stop 'Incorrect VDIM of IC file'
  end if
  read(1) dim_check
  if (dim_check /= BDIM) then
     stop 'Incorrect BDIM of IC file'
  end if

! Read in all important header information
  read(1) idata
  read(1) ilpdata
  read(1) rdata
  read(1) dpdata
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
  ! Check rank and domain (but ignore if rank==0 as we could be doing
  ! decomposition or a single MPI task run
  if (rank /= 0) then
     if (idata(40) /= rank .or. idata(41) /= numtasks) then
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
  time           = dpdata(1)
  lastsnap       = dpdata(2)
  mgas_orig      = dpdata(3)

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

  if (nunit > 0) read(1) unit_data(1:nunit)
  if (ndata > 0) read(1) data_id(1:ndata)
  if (ndata > 0) read(1) typedata(1:5,1:ndata)

  write(6,*) "unit_data : ",ndata,unit_data(1:nunit)
  write(6,*) "data_id : ",data_id(1:ndata)

! Allocate memory now we know ptot
! If MPI decomposition read, allocate 'minimal = TRUE'
  call allocate_memory(decomp_read)

#if defined(MHD)
  has_B = .FALSE.
#endif


! Now loop through array ids and read each array in turn
! ============================================================================
  do i=1,ndata

     ! Find pfirst, plast from typedata for this data set
     pfirst = typedata(2,i); plast = typedata(3,i)
     jtot = plast - pfirst + 1

     ! porig
     ! -----------------------------------------------------------------------
     if (data_id(i)=="porig") then
        allocate(idummy(1:jtot))
        read(1) idummy
        if (decomp_read) then
           minimal_sph(pfirst:plast)%porig = idummy(1:jtot)
        else
           sph(pfirst:plast)%porig = idummy(1:jtot)
        end if
        deallocate(idummy)

     ! Positions
     ! -----------------------------------------------------------------------
     else if (data_id(i)=="r") then
        allocate(rdummy2(1:NDIM,1:jtot))
        read(1) rdummy2
        if (decomp_read) then
           do j=1,jtot
              p = pfirst + j - 1
              minimal_sph(p)%r(1:NDIM) = rdummy2(1:NDIM,j)
           end do
        else
           do j=1,jtot
              p = pfirst + j - 1
              sph(p)%r(1:NDIM) = rdummy2(1:NDIM,j)
           end do
        end if
        deallocate(rdummy2)

     ! Masses
     ! -----------------------------------------------------------------------
     else if (data_id(i)=="m") then
        allocate(rdummy1(1:jtot))
        read(1) rdummy1
        if (decomp_read) then
           minimal_sph(pfirst:plast)%m = rdummy1(1:jtot)
        else
           sph(pfirst:plast)%m = rdummy1(1:jtot)
        end if
        deallocate(rdummy1)

     ! Smoothing length
     ! -----------------------------------------------------------------------
     else if (data_id(i)=="h") then
        allocate(rdummy1(1:jtot))
        read(1) rdummy1
        if (decomp_read) then
           minimal_sph(pfirst:plast)%h = rdummy1(1:jtot)
        else
           sph(pfirst:plast)%h = rdummy1(1:jtot)
        end if
        deallocate(rdummy1)

     ! Velocities
     ! -----------------------------------------------------------------------
     else if (data_id(i)=="v") then
        allocate(rdummy2(1:VDIM,1:jtot))
        read(1) rdummy2
        if (decomp_read) then
           do j=1,jtot
              p = pfirst + j - 1
              minimal_sph(p)%v(1:VDIM) = rdummy2(1:VDIM,j)
           end do
        else
           do j=1,jtot
              p = pfirst + j - 1
              sph(p)%v(1:VDIM) = rdummy2(1:VDIM,j)
           end do
        end if
        deallocate(rdummy2)

#if defined(MHD)
     ! Magnetic fields
     ! -----------------------------------------------------------------------
     else if (data_id(i)=="B") then
        has_B = .TRUE.
        allocate(rdummy2(1:BDIM,1:jtot))
        read(1) rdummy2
        if (decomp_read) then
           do j=1,jtot
              p = pfirst + j - 1
              minimal_sph(p)%B(1:BDIM) = rdummy2(1:BDIM,j)
           end do
        else
           do j=1,jtot
              p = pfirst + j - 1
              sph(p)%B(1:BDIM) = rdummy2(1:BDIM,j)
           end do
        end if
        deallocate(rdummy2)
#endif

     ! Density
     ! -----------------------------------------------------------------------
     else if (data_id(i)=="rho") then
        allocate(rdummy1(1:jtot))
        read(1) rdummy1
        if (decomp_read) then
           minimal_sph(pfirst:plast)%rho = rdummy1(1:jtot)
        else
           sph(pfirst:plast)%rho = rdummy1(1:jtot)
        end if
        deallocate(rdummy1)

     ! Temperature
     ! -----------------------------------------------------------------------
     else if (data_id(i)=='temp') then
        allocate(rdummy1(1:jtot))
        read(1) rdummy1
#if defined(HYDRO)
        if (decomp_read) then
           minimal_sph(pfirst:plast)%temp = rdummy1(1:jtot)
        else
           sph(pfirst:plast)%temp = rdummy1(1:jtot)
        end if
#endif
        deallocate(rdummy1)

     ! Internal energy
     ! -----------------------------------------------------------------------
     else if (data_id(i)=='u') then
        allocate(rdummy1(1:jtot))
        read(1) rdummy1
#if defined(HYDRO)
        if (decomp_read) then
           minimal_sph(pfirst:plast)%u = rdummy1(1:jtot)
        else
           sph(pfirst:plast)%u = rdummy1(1:jtot)
        end if
#endif
        deallocate(rdummy1)

     ! Sinks
     ! -----------------------------------------------------------------------
     else if (data_id(i)=='sink_v1') then
#if defined(SINKS)
        read(1) idummy3(1:6)
        sink_data_length = 11+NDIM+VDIM+2*dmdt_range_aux
        sink_data_length_aux = idummy3(4)
        allocate(ldummy(1:2))
        allocate(raux(1:sink_data_length_aux))
        if (stot > 0) then
           do j=pfirst,plast
              s = j
              read(1) ldummy
              read(1) idummy2
              read(1) raux
              sink(s)%id          = idummy2(1)
              sink(s)%ncreate     = idummy2(2)
              sink(s)%accrete     = ldummy(1)
              sink(s)%static      = ldummy(2)
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
                 dmdt_range_min = min(dmdt_range_aux,DMDT_RANGE)
                 sink(s)%macc = 0.0_PR
                 sink(s)%macc(1:dmdt_range_min) = &
                      &real(raux(NDIM+VDIM+12:&
                      &NDIM+VDIM+11+dmdt_range_min),DP)
                 sink(s)%tacc = 0.0_PR
                 sink(s)%tacc(1:dmdt_range_min) = &
                      &real(raux(NDIM+VDIM+12+dmdt_range_aux:&
                      &NDIM+VDIM+11+dmdt_range_min+dmdt_range_aux),DP)
              end if
              if (sink_data_length_aux > sink_data_length) &
                   &sink(s)%mmax = raux(NDIM+VDIM+12+2*dmdt_range_aux)
           end do
        end if
        deallocate(ldummy)
#else
        stop 'SEREN file contains sinks; not activated in code'
#endif

     ! Skip through arbitrary 1-D or 2-D array
     ! -----------------------------------------------------------------------
     else if (typedata(1,i) >= 1) then
        if (typedata(4,i)==1) then
           allocate(ldummy(1:plast-pfirst+1))
           read(1) ldummy
           deallocate(ldummy)
        else if (typedata(4,i)==2) then
           allocate(idummy(1:plast-pfirst+1))
           read(1) idummy
           deallocate(idummy)
        else if (typedata(4,i)==3) then
           allocate(ilpdummy(1:plast-pfirst+1))
           read(1) ilpdummy
           deallocate(ilpdummy)
        else if (typedata(4,i)==4 .and. typedata(1,i)==1) then
           allocate(rdummy1(1:plast-pfirst+1))
           read(1) rdummy1
           deallocate(rdummy1)
        else if (typedata(4,i)==4 .and. typedata(1,i)>1) then
           allocate(rdummy2(1:typedata(1,i),1:plast-pfirst+1))
           read(1) rdummy2
           deallocate(rdummy2)
        else if (typedata(4,i)==5 .and. typedata(1,i)==1) then
           allocate(dpdummy1(1:plast-pfirst+1))
           read(1) dpdummy1
           deallocate(dpdummy1)
        else if (typedata(4,i)==5 .and. typedata(1,i)>1) then
           allocate(dpdummy2(1:typedata(1,i),1:plast-pfirst+1))
           read(1) dpdummy2
           deallocate(dpdummy2)
        end if

     ! Skip through arbitrary data structure
     ! -----------------------------------------------------------------------
     else if (typedata(1,i) == 0 .and. typedata(4,i) == 7) then
        read(1) nvartype
        if (nvartype(1) > 0) allocate(ldummy(1:nvartype(1)))
        if (nvartype(2) > 0) allocate(idummy(1:nvartype(2)))
        if (nvartype(3) > 0) allocate(ilpdummy(1:nvartype(3)))
        if (nvartype(4) > 0) allocate(rdummy1(1:nvartype(4)))
        if (nvartype(5) > 0) allocate(dpdummy1(1:nvartype(5)))
        do p=pfirst,plast
           if (nvartype(1) > 0) read(1) ldummy
           if (nvartype(2) > 0) read(1) idummy
           if (nvartype(3) > 0) read(1) ilpdummy
           if (nvartype(4) > 0) read(1) rdummy1
           if (nvartype(5) > 0) read(1) dpdummy1
        end do
        if (nvartype(5) > 0) deallocate(dpdummy1)
        if (nvartype(4) > 0) deallocate(rdummy1)
        if (nvartype(3) > 0) deallocate(ilpdummy)
        if (nvartype(2) > 0) deallocate(idummy)
        if (nvartype(1) > 0) deallocate(ldummy)

     end if
  end do
! ============================================================================

#if defined(MHD)
  if (.NOT. has_B) then
     do p=1,ptot
        if (decomp_read) then
           minimal_sph(p)%B(1:BDIM) = 0.0_PR
        else
           sph(p)%B(1:BDIM) = 0.0_PR
        end if
     end do
  end if
#endif

#if defined(RTSPH) && defined(HYDRO) && !defined(INTERNAL_ENERGY)
  temp(1:ptot) = 1.0_PR
#endif

! Close file once finished
! ----------------------------------------------------------------------------
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
END SUBROUTINE read_data_seren_unform
