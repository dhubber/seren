! Program to join a simulation from multiple domains
! Andrew McLeod

program join
  use systemargs
  use datamodule
  use progmodule
  use seren_data_store
  implicit none
  character (LEN=3)    :: yn, yn2            ! Yes/No input
  integer              :: p, k               ! Loop variables
  integer              :: ierr               ! error variable
  integer              :: idata(1:50)        ! Dragon header integers
  integer(kind=ILP)    :: ilpdata(1:50)      ! Seren header ILP integers
  real (kind=PR)       :: rdata(1:50)        ! Dragon header reals
  real (kind=DP)       :: dpdata(1:50)       ! Seren header DP reals
  logical              :: commandline        ! Use the command line arguments
  logical              :: commandyn          ! We have the overwrite flag from the command line
  logical              :: remove_old         ! Remove old files? (implies overwrite on)
  integer              :: nargs              ! Number of command line arguments
  character (LEN=200)  :: runid              ! Command line argument 1
  character (LEN=200)  :: outputfile         ! Command line argument 2
  character(len=8)     :: file_ext           ! filename extension for data output
  character(len=200)    :: tempfilename       ! filename for data output
  integer, allocatable :: idummy(:)            ! ..
  real(kind=PR), allocatable :: pdummy(:,:)    ! parray dummy array
  real(kind=PR), allocatable :: rdummy1(:)     ! real dummy array
  real(kind=PR), allocatable :: rdummy3(:,:)   ! real vector dummy array
  integer :: idummy2(1:2)                      ! Integer dummy array
  logical :: ldummy(1:2)                       ! Logical dummy array
  integer :: s                                 ! Sink counter
  real(kind=PR), allocatable :: macc(:)        ! Accreted mass aux. array
  real(kind=DP), allocatable :: tacc(:)        ! Accretion time aux. array
  real(kind=PR) :: sdummy(1:20)                ! Dummy sink data array
  real(kind=PR), allocatable :: raux(:)        ! Dummy sink data array (new_format)
  integer              :: unknown_id           ! ID of unknown data type
  integer              :: type_id, width       ! type_id and width of unknown type
  integer              :: sink_data_length     ! length of sink data
  character(len=30)    :: sink_format_string   ! format string for reading sinks (new_format)

  ! Read data from screen

  nargs = num_of_args()
  commandline = .FALSE.
  commandyn = .FALSE.
  remove_old = .FALSE.

  if (nargs == 3 .OR. nargs == 4) then
    commandline = .TRUE.
    joinmpi = .TRUE.
    call get_cmd_arg(1,runid)
    call get_cmd_arg(2,MPInumthreads)
    call get_cmd_arg(3,outputfile)
    runid = trim(adjustl(runid))
    outputfile = trim(adjustl(outputfile))
    filename = trim(runid)
    if (nargs == 4) then
      commandyn = .TRUE.
      call get_cmd_arg(4,yn)
      if (.NOT.(yn=="yes".OR.yn=="y".OR.yn=="no".OR.yn=="n".OR.yn=="r".OR.yn=="rm")) then
        write (6,*) "Invalid command line flag! (third argument should be 'y', 'n' or 'r')"
        stop
      end if
      if (yn(1:1)=="r") then
         remove_old = .TRUE.
         yn="y"
      end if
    end if
  else if (nargs /= 0) then
    write (6,*) "Invalid number (",nargs,") of command line arguments!"
    stop
  else

    write (6,'(A)')   "   **********************************************************"
    write (6,'(A)')   "   *                     Join Simulations                   *"
    write (6,'(A//)') "   **********************************************************"

    write (6,'(A)',ADVANCE="NO") "Are we joining MPI simulations? (y/n): "
    read (5,*) yn2
    if (.NOT.(yn2=="yes".OR.yn2=="y".OR.yn2=="no".OR.yn2=="n")) then
      write (6,*) "Invalid response! (should be 'y' or 'n')"
      stop
    else if (yn2(1:1)=="y") then
      joinmpi = .TRUE.
      write (6,'(A)',ADVANCE="NO") "Enter input basename (e.g. runid.0001): "
      read (5,*) filename

      write (6,'(A)',ADVANCE="NO") "Enter number of MPI threads used: "
      read (5,*) MPInumthreads

    else
      joinmpi = .FALSE.
      write (6,'(A)',ADVANCE="NO") "How many files are we joining?: "
      read (5,*) MPInumthreads
    end if
  end if

!   allocate(MPIgeometry(0:MPInumthreads-1))
!
!   open(unit=10,file=trim(adjustl(filename))//".MPI",status="OLD",form="FORMATTED",iostat=ierr)
!   if (ierr /= 0) then
!     write (6,*) "File error opening MPI data file!"
!     close(10)
!     stop
!   end if
!
!   read(10,*) MPInumthreads, MPItreedepth, ptot, nsteps, snapshot, time
!   read(10,*) numtypes(-1:9)
!
!   allocate(MPItreeptot(0:0,0:MPInumthreads-1))
!
!   do i=0,MPInumthreads-1
!     read(10,*) idata(1:3), rdata(1:6)
!     MPItreeptot(0,idata(1)) = idata(3)
!   end do
!
!   close(10)

  call joinsimulation()

!   do i=1,20
!     read (10,*) idata(i)     ! Read dragon integer header
!   end do
!
!   ptot = idata(1)            ! Set particle number
!
!
!   allocate(r(1:3,1:ptot), v(1:3,1:ptot), h(1:ptot), m(1:ptot), temp(1:ptot), rho(1:ptot), ptype(1:ptot))
!   allocate(scratchspace(1:3,1:ptot),scratchscalar(1:ptot),tempids(1:ptot))
!
!   do i=1,50
!     read (10,*) rdata(i)     ! Read dragon real header
!   end do
!
!   rdata(1) = 0               ! Set time=0
!
!   do i=1,ptot
!     read(10,*) r(1:3,i)        ! Read positions
!   end do
!
!   do i=1,ptot
!     read(10,*) v(1:3,i)        ! Read velocities
!   end do
!
!   do i=1,ptot
!     read(10,*) temp(i)         ! Read temperatures
!   end do
!
!   do i=1,ptot
!     read(10,*) h(i)            ! Read smoothing lengths
!   end do
!
!   do i=1,ptot
!     read(10,*) rho(i)          ! Read densities
!   end do
!
!   do i=1,ptot
!     read(10,*) m(i)            ! Read particle mass
!   end do
!
!   do i=1,ptot
!     read(10,*) ptype(i)         ! Read particle type
!   end do
!
!   close(10)
!
!   call splitsimulation()               ! Split the simulation into separate files
!
!   call writesplitfiles()               ! Write these split files to disk

  if (commandline) then
    filename = outputfile
  else
    write (6,*) "Enter output file name:"
    read (5,*) filename
  end if

  ! Open output file
  if (formatted_files) then
    open(unit=10,file=trim(filename),status="NEW",form="FORMATTED",iostat=ierr)
  else
    open(unit=10,file=trim(filename),status="NEW",form="UNFORMATTED",iostat=ierr)
  end if
  if (ierr /= 0) then
    if (.NOT.commandyn) then
      write (6,*) "File already exists. Overwrite? (y/n)"
      read (5,*) yn
    end if
    if (yn(1:1) == "y") then
      if (formatted_files) then
        open(unit=10,file=trim(filename),status="REPLACE",form="FORMATTED",iostat=ierr)
      else
        open(unit=10,file=trim(filename),status="REPLACE",form="UNFORMATTED",iostat=ierr)
      end if
      if (ierr /= 0) then
        write (6,*) "File error, iostat=", ierr
        stop
      end if
    else
      close (10)
      stop
    end if
  end if

  if (seren_format) then

     if (new_format) then
        sink_data_length = 11+NDIM+VDIM+2*dmdt_range
        write (sink_format_string,'(A,I0,A)') "(", sink_data_length, "E18.10)"
        allocate(raux(1:sink_data_length))
     end if

     allocate(macc(1:dmdt_range))
     allocate(tacc(1:dmdt_range))
     idata(1:50)     = 0
     ilpdata(1:50)   = 0
     dpdata(1:50)    = 0.0_DP
     rdata(1:50)     = 0.0_PR
     !unit_data(1:50) = ''
     data_id(1:50)   = ''
     ndata           = 0

     if (is_porig) then
        ndata = ndata + 1
        data_id(ndata) = 'porig'
        if (new_format) typedata(1:5,ndata) = (/1,1,ptot,2,0/)
     end if
     if (is_parray) then
        ndata = ndata + 1
        data_id(ndata) = 'parray'
        if (new_format) typedata(1:5,ndata) = (/NDIM+2,1,ptot,4,0/)
     end if
     if (is_r) then
        ndata = ndata + 1
        data_id(ndata) = 'r'
        if (new_format) typedata(1:5,ndata) = (/NDIM,1,ptot,4,unit_r/)
     end if
     if (is_m) then
        ndata = ndata + 1
        data_id(ndata) = 'm'
        if (new_format) typedata(1:5,ndata) = (/1,1,ptot,4,unit_m/)
     end if
     if (is_h) then
        ndata = ndata + 1
        data_id(ndata) = 'h'
        if (new_format) typedata(1:5,ndata) = (/1,1,ptot,4,unit_h/)
     end if
     if (is_v) then
        ndata = ndata + 1
        data_id(ndata) = 'v'
        if (new_format) typedata(1:5,ndata) = (/VDIM,1,ptot,4,unit_v/)
     end if
     if (is_rho) then
        ndata = ndata + 1
        data_id(ndata) = 'rho'
        if (new_format) typedata(1:5,ndata) = (/1,1,ptot,4,unit_rho/)
     end if
     if (is_temp) then
        ndata = ndata + 1
        data_id(ndata) = 'temp'
        if (new_format) typedata(1:5,ndata) = (/1,1,ptot,4,unit_temp/)
     end if
     if (is_u) then
        ndata = ndata + 1
        data_id(ndata) = 'u'
        if (new_format) typedata(1:5,ndata) = (/1,1,ptot,4,unit_u/)
     end if
     if (is_B) then
        ndata = ndata + 1
        data_id(ndata) = 'B'
        if (new_format) typedata(1:5,ndata) = (/BDIM,1,ptot,4,unit_B/)
     end if
     if (is_sink_v1) then
        ndata = ndata + 1
        data_id(ndata) = 'sink_v1'
        if (new_format) typedata(1:5,ndata) = (/0,1,stot,7,0/)
     end if
     do i=1,nunknown
        ndata = ndata + 1
        data_id(ndata) = trim(unknown_data(i)%data_id)
        if (new_format) then
           width = unknown_data(i)%width
           type_id = unknown_data(i)%type_id
           typedata(1:5,ndata) = (/width,1,ptot,type_id,0/)
        end if
     end do

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
     idata(30)   = dmdt_range
     idata(31)   = pgas_orig
     idata(32)   = pp_gather
     ilpdata(1)  = snapshot
     ilpdata(2)  = nsteps
     ilpdata(3)  = ntempnext
     ilpdata(4)  = ndiagnext
     ilpdata(5)  = nsnapnext
     ilpdata(6)  = nsinknext
     rdata(1)    = h_fac
     dpdata(1)   = time
     dpdata(2)   = lastsnap
     dpdata(3)   = mgas_orig

     if (formatted_files) then

        ! Write information identifying format and precision of file
        if (new_format) then
           write(10,'(A16)') "SERENASCIIDUMPV2"
        else
           write(10,'(A16)') "SERENASCIIDUMPV1"
        end if
        write(10,'(I2)') PR
        write(10,'(I2)') NDIM
        write(10,'(I2)') VDIM
        write(10,'(I2)') BDIM

        ! Write the rest of the header information to file
        do i=1,50
           write(10,'(I10)') idata(i)
        end do
        do i=1,50
           write(10,'(I10)') ilpdata(i)
        end do
        do i=1,50
           write(10,'(E18.10)') rdata(i)
        end do
        do i=1,50
           write(10,'(E18.10)') dpdata(i)
        end do
        if (nunit > 0) then
           do i=1,nunit
              write(10,'(a)') unit_data(i)
           end do
        end if
        if (ndata > 0) then
           if (new_format) then
              do i=1,ndata
                 write(10,'(a)') data_id(i)
              end do
              do i=1,ndata
                 write(1,'(5I10)') typedata(1:5,i)
              end do
           else
              do i=1,ndata
                 write(10,'(a)') data_id(i)
              end do
           end if
        end if

     else
        ! Write information identifying format and precision of file
        if (new_format) then
           format_id = 'SERENBINARYDUMPV2'
        else
           format_id = 'SERENBINARYDUMPV1'
        end if
        write(10) format_id
        write(10) PR
        write(10) NDIM
        write(10) VDIM
        write(10) BDIM

        ! Write header information to file
        write(10) idata
        write(10) ilpdata
        write(10) rdata
        write(10) dpdata
        if (new_format) then
           if (nunit > 0) write(10) unit_data(1:nunit)
           if (ndata > 0) then
              write(10) data_id(1:ndata)
              write(10) typedata(1:5,1:ndata)
           end if
        else
           if (nunit > 0) write(10) unit_data
           if (ndata > 0) write(10) data_id
        end if
     end if

     do i=1,ndata
        select case (trim(data_id(i)))
           case ("porig")
              ! Original particle number
              if (allocated(rdummy1)) deallocate(rdummy1)
              allocate(idummy(1:ptot))
              if (formatted_files) then
                 do p=1,ptot
                    write(10,'(I10)') porig(p)
                 end do
              else
                 do p=1,ptot
                    idummy(p) = porig(p)
                 end do
                 write(10) idummy
              end if
              deallocate(idummy)
           case ("parray")
              ! Positions, masses and smoothing lengths
              if (allocated(rdummy1)) deallocate(rdummy1)
              allocate(pdummy(1:NDIM+2,1:ptot))
              if (formatted_files) then
                 do p=1,ptot
#if NDIM==1
                    write(10,'(3E18.10)') parray(1:NDIM+2,p)
#elif NDIM==2
                    write(10,'(4E18.10)') parray(1:NDIM+2,p)
#else
                    write(10,'(5E18.10)') parray(1:NDIM+2,p)
#endif
                 end do
              else
                 do p=1,ptot
                    pdummy(1:NDIM+2,p) = parray(1:NDIM+2,p)
                 end do
                 write(10) pdummy
              end if
              deallocate(pdummy)
           case ("r")
              ! Position (output by SEREN in parray instead)
              if (allocated(rdummy1)) deallocate(rdummy1)
              if (formatted_files) then
                 do p=1,ptot
#if NDIM==1
                    write(10,'(E18.10)') parray(1:NDIM,p)
#elif NDIM==2
                    write(10,'(2E18.10)') parray(1:NDIM,p)
#else
                    write(10,'(3E18.10)') parray(1:NDIM,p)
#endif
                 end do
              else
                 allocate(rdummy3(1:NDIM,1:ptot))
                 do p=1,ptot
                    rdummy3(1:NDIM,p) = parray(1:NDIM,p)
                 end do
                 write(10) rdummy3
                 deallocate(rdummy3)
              end if
           case ("h")
              if (.NOT. allocated(rdummy1)) allocate(rdummy1(1:ptot))
              ! Smoothing lengths (output by SEREN in parray instead)
              if (formatted_files) then
                 do p=1,ptot
                    write(10,'(E18.10)') parray(NDIM+2,p)
                 end do
              else
                 do p=1,ptot
                    rdummy1(p) = parray(NDIM+2,p)
                 end do
                 write(10) rdummy1
              end if
           case ("m")
              if (.NOT. allocated(rdummy1)) allocate(rdummy1(1:ptot))
              ! Mass (output by SEREN in parray instead)
              if (formatted_files) then
                 do p=1,ptot
                    write(10,'(E18.10)') parray(NDIM+1,p)
                 end do
              else
                 do p=1,ptot
                    rdummy1(p) = parray(NDIM+1,p)
                 end do
                 write(10) rdummy1
              end if
           case ("v")
              ! Velocities
              if (allocated(rdummy1)) deallocate(rdummy1)
              if (formatted_files) then
                 do p=1,ptot
#if VDIM==1
                    write(10,'(E18.10)') v(1:VDIM,p)
#elif VDIM==2
                    write(10,'(2E18.10)') v(1:VDIM,p)
#else
                    write(10,'(3E18.10)') v(1:VDIM,p)
#endif
                 end do
              else
                 allocate(rdummy3(1:VDIM,1:ptot))
                 do p=1,ptot
                    rdummy3(1:VDIM,p) = v(1:VDIM,p)
                 end do
                 write(10) rdummy3
                 deallocate(rdummy3)
              end if
           case ("rho")
              if (.NOT. allocated(rdummy1)) allocate(rdummy1(1:ptot))
              ! Densities
              if (formatted_files) then
                 do p=1,ptot
                    write(10,'(E18.10)') rho(p)
                 end do
              else
                 do p=1,ptot
                    rdummy1(p) = rho(p)
                 end do
                 write(10) rdummy1
              end if
           case ("temp")
              if (.NOT. allocated(rdummy1)) allocate(rdummy1(1:ptot))
              ! Temperatures
              if (formatted_files) then
                 do p=1,ptot
                    write(10,'(E18.10)') temp(p)
                 end do
              else
                 do p=1,ptot
                    rdummy1(p) = temp(p)
                 end do
                 write(10) rdummy1
              end if
           case ("u")
              if (.NOT. allocated(rdummy1)) allocate(rdummy1(1:ptot))
              ! Internal energy
              if (formatted_files) then
                 do p=1,ptot
                    write(10,'(E18.10)') u(p)
                 end do
              else
                 do p=1,ptot
                    rdummy1(p) = u(p)
                 end do
                 write(10) rdummy1
              end if
           case ("B")
              ! Magnetic fields
              if (allocated(rdummy1)) deallocate(rdummy1)
              if (formatted_files) then
                 do p=1,ptot
#if VDIM==1
                    write(10,'(E18.10)') B(1:BDIM,p)
#elif VDIM==2
                    write(10,'(2E18.10)') B(1:BDIM,p)
#else
                    write(10,'(3E18.10)') B(1:BDIM,p)
#endif
                 end do
              else
                 allocate(rdummy3(1:BDIM,1:ptot))
                 do p=1,ptot
                    rdummy3(1:BDIM,p) = B(1:BDIM,p)
                 end do
                 write(10) rdummy3
                 deallocate(rdummy3)
              end if
           case ("sink_v1")
              if (allocated(rdummy1)) deallocate(rdummy1)
              ! Sinks
              if (formatted_files) then
                 if (new_format) then
                    if (stot > 0) then
                       write(10,'(6I10)') 2,2,0,sink_data_length,0,0
                       do s=1,stot
                          write(10,'(2L1)') sink_array(s)%accrete,sink_array(s)%static
                          write(10,'(2I8)') sink_array(s)%id, sink_array(s)%ncreate
                          write(10,sink_format_string) &
                               &sink_array(s)%tcreate,&
                               &sink_array(s)%r(1:NDIM),&
                               &sink_array(s)%v(1:VDIM),&
                               &sink_array(s)%m,&
                               &sink_array(s)%h,&
                               &sink_array(s)%radius,&
                               &sink_array(s)%angmom(1:3),&
                               &sink_array(s)%dmdt,&
                               &sink_array(s)%star_radius,&
                               &sink_array(s)%luminosity,&
                               &sink_array(s)%temperature,&
                               &sink_array(s)%macc(1:dmdt_range),&
                               &sink_array(s)%tacc(1:dmdt_range)
                       end do
                    end if
                 else
                    if (stot > 0) then
                       do s=1,stot
                          write(10,'(2I8)')     sink_array(s)%id, sink_array(s)%ncreate
                          write(10,'(2L1)')     sink_array(s)%accrete,sink_array(s)%static
                          write(10,'(E18.10)')  real(sink_array(s)%tcreate,PR)
                          write(10,'(3E18.10)') sink_array(s)%r(1:NDIM)
                          write(10,'(3E18.10)') sink_array(s)%v(1:VDIM)
                          write(10,'(E18.10)')  sink_array(s)%m
                          write(10,'(E18.10)')  sink_array(s)%h
                          write(10,'(E18.10)')  sink_array(s)%radius
                          write(10,'(3E18.10)') sink_array(s)%angmom(1:3)
                          write(10,'(E18.10)')  sink_array(s)%dmdt
                          write(10,'(E18.10)')  sink_array(s)%star_radius
                          write(10,'(E18.10)')  sink_array(s)%luminosity
                          write(10,'(E18.10)')  sink_array(s)%temperature
                          if (dmdt_range > 0) then
                             do k=1,dmdt_range
                                write(10,'(2E18.10)') sink_array(s)%macc(k),&
                                     &sink_array(s)%tacc(k)
                             end do
                          end if
                       end do
                    end if
                 end if
              else
                 if (new_format) then
                    if (stot > 0) then
                       write(10) 2,2,0,sink_data_length,0,0
                       do s=1,stot
                          ldummy(1)                     = sink_array(s)%accrete
                          ldummy(2)                     = sink_array(s)%static
                          idummy2(1)                    = sink_array(s)%id
                          idummy2(2)                    = sink_array(s)%ncreate
                          raux                          = 0.0_PR
                          raux(1)                       = real(sink_array(s)%tcreate,PR)
                          raux(2:NDIM+1)                = sink_array(s)%r(1:NDIM)
                          raux(NDIM+2:NDIM+VDIM+1)      = sink_array(s)%v(1:NDIM)
                          raux(NDIM+VDIM+2)             = sink_array(s)%m
                          raux(NDIM+VDIM+3)             = sink_array(s)%h
                          raux(NDIM+VDIM+4)             = sink_array(s)%radius
                          raux(NDIM+VDIM+5:NDIM+VDIM+7) = sink_array(s)%angmom(1:3)
                          raux(NDIM+VDIM+8)             = sink_array(s)%dmdt
                          raux(NDIM+VDIM+8)             = sink_array(s)%star_radius
                          raux(NDIM+VDIM+8)             = sink_array(s)%luminosity
                          raux(NDIM+VDIM+8)             = sink_array(s)%temperature
                          if (DMDT_RANGE > 0) then
                             raux(NDIM+VDIM+12:NDIM+VDIM+11+dmdt_range) = &
                                & sink_array(s)%macc(1:dmdt_range)
                             raux(NDIM+VDIM+12+dmdt_range:NDIM+VDIM+11+2*dmdt_range) = &
                                & sink_array(s)%tacc(1:dmdt_range)
                          end if
                          write(10) ldummy
                          write(10) idummy2
                          write(10) raux
                       end do
                    end if
                 else
                    if (stot > 0) then
                       do s=1,stot
                          idummy2(1)       = sink_array(s)%id
                          idummy2(2)       = sink_array(s)%ncreate
                          ldummy(1)        = sink_array(s)%accrete
                          ldummy(2)        = sink_array(s)%static
                          sdummy           = 0.0_PR
                          sdummy(1)        = real(sink_array(s)%tcreate,PR)
                          sdummy(2:1+NDIM) = sink_array(s)%r(1:NDIM)
                          sdummy(5:4+NDIM) = sink_array(s)%v(1:NDIM)
                          sdummy(8)        = sink_array(s)%m
                          sdummy(9)        = sink_array(s)%h
                          sdummy(10)       = sink_array(s)%radius
                          sdummy(11:13)    = sink_array(s)%angmom(1:3)
                          sdummy(14)       = sink_array(s)%dmdt
                          sdummy(15)       = sink_array(s)%star_radius
                          sdummy(16)       = sink_array(s)%luminosity
                          sdummy(17)       = sink_array(s)%temperature
                          write(10) idummy2
                          write(10) ldummy
                          write(10) sdummy
                          if (dmdt_range > 0) then
                             macc(1:dmdt_range) = sink_array(s)%macc(1:dmdt_range)
                             tacc(1:dmdt_range) = sink_array(s)%tacc(1:dmdt_range)
                             write(10) macc
                             write(10) tacc
                          end if
                       end do
                    end if
                 end if
              end if
           case default
              ! Unknown data types
              unknown_id = 0
              do k=1,nunknown
                 if (trim(unknown_data(k)%data_id) == trim(data_id(i))) then
                    unknown_id = k
                    exit
                 end if
              end do
              if (unknown_id == 0) stop "Could not find unknown data type?"
              width = unknown_data(k)%width
              type_id = unknown_data(k)%type_id
              if (type_id == 4) then
                 ! PR data
                 if (width==1) then
                    ! PR scalar
                    if (.NOT. allocated(rdummy1)) allocate(rdummy1(1:ptot))
                    if (formatted_files) then
                       do p=1,ptot
                          write(10,'(E18.10)') unknown_data(unknown_id)%data(p)
                       end do
                    else
                       do p=1,ptot
                          rdummy1(p) = unknown_data(unknown_id)%data(p)
                       end do
                       write(10) rdummy1
                    end if
                 else
                    ! PR vector
                    if (formatted_files) then
                       do p=1,ptot
                          write(10,'(9999E18.10)') unknown_data(unknown_id)%vector(1:width,p)
                       end do
                    else
                       allocate(pdummy(1:width,1:ptot))
                       do p=1,ptot
                          pdummy(1:width,p) = unknown_data(unknown_id)%vector(1:width,p)
                       end do
                       write(10) pdummy
                       deallocate(pdummy)
                    end if
                 end if
              else if (type_id == 2) then
                 if (allocated(rdummy1)) deallocate(rdummy1)
                 if (.NOT. allocated(idummy)) allocate(idummy(1:ptot))
                 if (formatted_files) then
                    do p=1,ptot
                       write(10,'(I10)') unknown_data(unknown_id)%idata(p)
                    end do
                 else
                    do p=1,ptot
                       idummy(p) = unknown_data(unknown_id)%idata(p)
                    end do
                    write(10) idummy
                 end if
              else
                 stop "Haven't written non-PR or integer data yet..."
              end if
        end select
     end do

  else
     idata = 0
     rdata = real(0,PR)

     idata(2) = nsteps
     idata(4) = snapshot
     rdata(1) = time

     idata(1) = ptot            ! Set particle number
     idata(3) = ptot            ! and again

     if (formatted_files) then
       ! Write DRAGON header
       write(10,'(19(I0/),I0)') idata(1:20)
       write(10,'(49(E16.7/),E16.7)') rdata

       do i=1,ptot
         write(10,*) r(1:3,i)        ! Write positions
       end do

       do i=1,ptot
         write(10,*) v(1:3,i)        ! Write velocities
       end do

       do i=1,ptot
         write(10,*) temp(i)         ! Write temperatures
       end do

       do i=1,ptot
         write(10,*) h(i)            ! Write smoothing lengths
       end do

       do i=1,ptot
         write(10,*) rho(i)          ! Write densities
       end do

       do i=1,ptot
         write(10,*) m(i)            ! Write particle mass
       end do

       do i=1,ptot
         write(10,*) ptype(i)         ! Write particle type
       end do

     else
       write (10) idata(1:20)
       write (10) rdata
       write (10) r(1:3,1:ptot)
       write (10) v(1:3,1:ptot)
       write (10) temp(1:ptot)
       write (10) h(1:ptot)
       write (10) rho(1:ptot)
       write (10) m(1:ptot)
       write (10) ptype(1:ptot)
     end if

  end if

  close(10)

  if (remove_old) then
    do i=0,MPInumthreads-1

      write(file_ext,"(I0)") i
      tempfilename = trim(adjustl(filename))//".MPI."//trim(adjustl(file_ext))

      open(10, file=tempfilename, status="unknown")
      close(10,status="delete")

    end do
  end if

  stop

end program join
