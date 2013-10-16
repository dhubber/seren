! Program to modify a simulation
! Andrew McLeod
! ----------------------------------------------------------------------------

module progmodule
  use datamodule
  use seren_data_store
  implicit none

  contains

   subroutine try_read_seren_header(unitno,filename,iheader,ilpheader,rheader,dpheader,formatted,is_seren_form)
      integer, intent(in)               :: unitno        ! File IO unit number
      character (len=*), intent(in)     :: filename      ! File name
      integer, intent(out)              :: iheader(1:50) ! Dragon INTEGER header
      integer (kind=ILP), intent(out)   :: ilpheader(1:50) ! Seren ILP INTEGER header
      real (kind=PR), intent(out)       :: rheader(1:50) ! Dragon REAL header
      real (kind=DP), intent(out)       :: dpheader(1:50)! Seren DP REAL header
      real                              :: rheader_sp(1:50) ! Seren SP REAL header
      real                              :: rheader_dp(1:50) ! Seren SP REAL header
      logical, intent(out)              :: formatted     ! Is file formatted?
      logical, intent(out)              :: is_seren_form ! Is file Seren format?
      logical                           :: exists        ! Does file exist?
      integer                           :: ierr          ! Error value
      integer                           :: i             ! Loop counter
   
      new_format = .FALSE.

      inquire (file=filename,exist=exists)
      if (.NOT.exists) then
         write (6,*) "File ", filename, " does not exist!"
         stop
      endif

      ! Open file as unformatted file first
      open(unit=unitno,file=filename,status="OLD",form="UNFORMATTED",iostat=ierr)
      if (ierr /= 0) then
         write (6,*) "File error opening file ", filename, "!"
         close(unitno)
         stop
      end if
      rewind(unit=unitno)
   
      read (unit=unitno,iostat=ierr) format_id
   
      if (ierr /= 0 .OR. (trim(adjustl(format_id)) /= "SERENBINARYDUMPV1" .AND. &
            & trim(adjustl(format_id)) /= "SERENBINARYDUMPV2")) then
         ! Try formatted read
         formatted = .TRUE.
         close (unit=unitno)
         open(unit=unitno,file=filename,status='old',form='formatted',iostat=ierr)
         if (ierr /= 0) then
            write (6,*) "File error opening file ", filename, " as formatted!"
            stop
         endif
         rewind(unit=unitno)
         read (unit=unitno,fmt=*,iostat=ierr) format_id
         if (ierr /= 0 .OR. (trim(adjustl(format_id)) /= "SERENASCIIDUMPV1" .AND. &
            & trim(adjustl(format_id)) /= "SERENASCIIDUMPV2")) then
            is_seren_form = .FALSE.
            return
         endif
         if (trim(adjustl(format_id))=="SERENASCIIDUMPV2") new_format = .TRUE.
      else
         formatted = .FALSE.
         if (trim(adjustl(format_id))=="SERENBINARYDUMPV2") new_format = .TRUE.
      end if

      is_seren_form = .TRUE.

      if (formatted) then
         read (unitno,*) PRtemp
         read (unitno,*) NDIMtemp
         read (unitno,*) VDIMtemp
         read (unitno,*) BDIMtemp
      else
         read (unitno) PRtemp
         read (unitno) NDIMtemp
         read (unitno) VDIMtemp
         read (unitno) BDIMtemp
      end if
   
      if (PRtemp == 8 .OR. PRtemp == 2) then
         ! Double precision file
         doubleprec = .TRUE.
      else if (PR == 4 .OR. PR == 1) then
         ! Single precision file
         doubleprec = .FALSE.
      else
         write (6,*) "Unknown kind value of ", PR, "!"
         stop
      end if
   
      if (formatted) then
         do i=1,size(iheader)
            read (unitno,*,end=55) iheader(i)
         end do
      
         do i=1,size(ilpheader)
            read (unitno,*,end=55) ilpheader(i)
         end do
      
         do i=1,size(rheader)
            read (unitno,*,end=55) rheader(i)
         end do

         do i=1,size(dpheader)
            read (unitno,*,end=55) dpheader(i)
         end do

         ndata = iheader(21)
         nunit = iheader(20)

         do i=1,nunit
            read (unitno,'(A)') unit_data(i)
         end do
         do i=1,ndata
            read (unitno,'(A)') data_id(i)
         end do
         if (new_format) then
            do i=1,ndata
               read (unitno,'(5I10)') typedata(1:5,i)
            end do
         end if
      else
         read (unitno,end=55) iheader
         read (unitno,end=55) ilpheader
         if (doubleprec) then
            read (unitno,end=55) rheader_dp
            rheader = real(rheader_dp,PR)
         else
            read (unitno,end=55) rheader_sp
            rheader = real(rheader_sp,PR)
         endif
         read (unitno,end=55) dpheader

         ndata = iheader(21)
         nunit = iheader(20)

         if (nunit > 0) then
            if (.NOT. new_format) read (unitno) unit_data
            if (new_format) read (unitno) unit_data(1:nunit)
         end if
         if (ndata > 0) then
            if (.NOT. new_format) read (unitno) data_id
            if (new_format) then
               read (unitno) data_id(1:ndata)
               read (unitno) typedata(1:5,1:ndata)
            end if
         end if
      end if

      return

55    continue
      write (6,*) "End of file reading header!!"
      stop
   
      end subroutine try_read_seren_header

  subroutine openandreadheader(unitno,filename,iheader,rheader,formatted)
    integer, intent(in)               :: unitno        ! File IO unit number
    character (len=*), intent(in)     :: filename      ! File name
    integer, intent(out)              :: iheader(1:20) ! Dragon INTEGER header
    real (kind=PR), intent(out)       :: rheader(1:50) ! Dragon REAL header
    logical, intent(out)              :: formatted     ! Is file formatted?
    logical                           :: exists        ! Does file exist?
    integer                           :: ierr          ! Error value
    integer                           :: i             ! Loop counter

    inquire (file=filename,exist=exists)
    if (.NOT.exists) then
      write (6,*) "File ", filename, " does not exist!"
      stop
    endif

    ! Open file as formatted file first
    open(unit=unitno,file=filename,status="OLD",form="FORMATTED",iostat=ierr)
    if (ierr /= 0) then
      write (6,*) "File error opening file ", filename, "!"
      close(unitno)
      stop
    end if
    rewind(unit=unitno)

    do i=1,20
      read(unit=unitno,fmt=*,err=10,iostat=ierr) iheader(i)
    end do
    do i=1,50
      read(unit=unitno,fmt=*,err=10,iostat=ierr) rheader(i)
    end do

    formatted = .TRUE.

    return
    10 continue

    ! Formatted read has failed - try unformatted read
    close(unitno)
    formatted = .FALSE.
    open(unit=unitno,file=filename,status="OLD",form="UNFORMATTED",iostat=ierr)
    rewind(unit=unitno)
    read(unitno) iheader(1:20)
    read(unitno) rheader

    return
  end subroutine openandreadheader

  subroutine joinsimulation
    ! Join multiple files together
    integer              :: i, p               ! Loop variables
    character(len=8)     :: file_ext           ! filename extension for data output
    character(len=200)    :: tempfilename       ! filename for data output
    character(len=200), allocatable :: nonmpifiles(:) ! Array of filenames for non-mpi adding
    integer              :: threadptot         ! Number of particles in this file
    integer              :: iheader(1:50)      ! Dragon integer header
    integer (kind=ILP)   :: ilpheader(1:50)    ! Seren ILP integer header
    real (kind=PR)       :: rheader(1:50)      ! Dragon real header
    real (kind=DP)       :: dpheader(1:50)     ! Seren DP real header
    real (kind=PR)       :: rdummy(1:NDIM)     ! Dragon dummy real
    logical              :: formatted          ! Is file formatted?
    integer, allocatable :: idummy1(:)         ! dummy array for integers
    real(kind=PR), allocatable :: rdummy1(:)   ! dummy array for reals
    real(kind=PR), allocatable :: rdummy3(:,:) ! dummy array for 3D reals

    ptot = 0
    numtypes = 0
    pboundary = 0
    picm = 0
    pgas = 0
    pcdm = 0
    pdust = 0
    pion = 0
    mgas_orig = 0._DP
    pgas_orig = 0

    allocate(nonmpifiles(0:MPInumthreads-1))

    ! Loop over all the files, reading in the ptot and number of particle types
    do i=0,MPInumthreads-1
!       file_ext = ""
!       ! Data output needs extension
!       if (i < 10) then
!         write(file_ext,"(I1)") i
!         file_ext="000"//file_ext
!       else if (i < 100) then
!         write(file_ext,"(I2)") i
!         file_ext="00"//file_ext
!       else if (i < 1000) then
!         write(file_ext,"(I3)") i
!         file_ext="0"//file_ext
!       else if (i < 10000) then
!         write(file_ext,"(I4)") i
!       else if (i < 100000) then
!         write(file_ext,"(A)") "too_big"  ! *** too many threads ***
!       end if

      if (joinmpi) then
        write(file_ext,"(I0)") i

        !tempfilename = trim(adjustl(filename))//"."//trim(adjustl(file_ext))
        tempfilename = trim(adjustl(filename))//".MPI."//trim(adjustl(file_ext))
      else
        write (6,'(A)',ADVANCE="NO") "Enter filename of file to be added: "
        read (5,*) nonmpifiles(i)
        tempfilename = nonmpifiles(i)
      end if

      !write(6,*) "tempfilename is ", tempfilename, i, file_ext

      ! First pass to scan ptot and particle type numbers
      ! First, read in header information

      call try_read_seren_header(10,tempfilename,iheader,ilpheader,rheader,dpheader,formatted,seren_format)

      if (.NOT. seren_format) call openandreadheader(10,tempfilename,iheader(1:20),rheader,formatted)

      if (i==0) formatted_files = formatted

      ! Assign variables for important information
      if (seren_format) then
         close(10)
         threadptot = iheader(1)
         t_stot = iheader(2)
         t_pboundary = iheader(3)
         t_picm = iheader(4)
         t_pgas = iheader(5)
         t_pcdm = iheader(6)
         t_pdust = iheader(7)
         t_pion = iheader(8)
         dmdt_range = iheader(30)
         if (joinmpi) then
            pgas_orig = iheader(31)
         else
            pgas_orig = pgas_orig + iheader(31)
         end if
         pp_gather = iheader(32)
         snapshot = ilpheader(1)
         nsteps = ilpheader(2)
         ntempnext = ilpheader(3)
         ndiagnext = ilpheader(4)
         nsnapnext = ilpheader(5)
         nsinknext = ilpheader(6)
         h_fac = rheader(1)
         time = dpheader(1)
         lastsnap = dpheader(2)
         if (joinmpi) then
            mgas_orig = dpheader(3)
         else
            mgas_orig = mgas_orig + dpheader(3)
         end if
         ptot = ptot + threadptot
         pboundary = pboundary + t_pboundary
         picm = picm + t_picm
         pgas = pgas + t_pgas
         pcdm = pcdm + t_pcdm
         pdust = pdust + t_pdust
         pion = pion + t_pion
         if (.NOT. joinmpi) stot = stot + t_stot
         if (joinmpi) stot = t_stot
         is_porig = .TRUE.
         if (new_format) then
            is_parray = .FALSE.
            is_r = .TRUE.
            is_m = .TRUE.
            is_h = .TRUE.
         else
            is_parray = .TRUE.
            is_r = .FALSE.
            is_m = .FALSE.
            is_h = .FALSE.
         end if
         is_v = .TRUE.
         is_temp = .TRUE.
         is_rho = .TRUE.
         is_u = .TRUE.
         is_sink_v1 = .FALSE.
      else
         threadptot = iheader(1)
         nsteps     = iheader(2)
         snapshot   = iheader(4)
         time       = rheader(1)
         ptot = ptot + threadptot
      end if

      if (.NOT. seren_format) then
         allocate(ptype(1:threadptot))
         allocate(rdummy1(1:threadptot),rdummy3(1:NDIM,1:threadptot),idummy1(1:threadptot))
   
         if (formatted) then
           do p=1,threadptot*2
             read(10,*) rdummy(1:NDIM)
           end do
           do p=1,threadptot*4
             read(10,*) rdummy(1)
           end do
           do p=1,threadptot
             read(10,*) ptype(p)
             numtypes(ptype(p)) = numtypes(ptype(p)) + 1
           end do
         else
           do p=1,2
             read(10) rdummy3
           end do
           do p=1,4
             read(10) rdummy1
           end do
           read(10) ptype
           do p=1,threadptot
             numtypes(ptype(p)) = numtypes(ptype(p)) + 1
           end do
         end if
   
         close(10)
   
         deallocate(ptype)
         deallocate(rdummy1, rdummy3, idummy1)
      end if

    end do

    if (seren_format) then
      numtypes(9) = numtypes(9) + picm
      numtypes(6) = numtypes(6) + pboundary
      numtypes(1) = numtypes(1) + pgas
      numtypes(10) = numtypes(10) + pcdm
      numtypes(11) = numtypes(11) + pdust
      numtypes(12) = numtypes(12) + pion
    end if

    ! If there are any other types of particles, choke at this point

    if ((numtypes(1) + numtypes(4) + numtypes(6) + numtypes(9) + numtypes(10) +&
        &numtypes(11) + numtypes(12) + numtypes(-1)) /= ptot) then
      write (6,*) "Particles of non-permitted types!"
      stop
    end if

    ! For now, forbid split particles
    if (numtypes(4) > 0) then
      write (6,*) "Split particles not yet implemented!"
      stop
    end if

    ! Order of particles is: boundary (6), icm (9), gas (1), cdm(10),
    !                        dust(11), ion(12), sinks(-1)
    slot(6)  = 1
    slot(9)  = 1 + numtypes(6)
    slot(1)  = slot(9) + numtypes(9)
    slot(10) = slot(1) + numtypes(1)
    slot(11) = slot(10) + numtypes(10)
    slot(12) = slot(11) + numtypes(11)
    slot(-1) = slot(12) + numtypes(12)

    !write (6,*) "slots(6,9,1,4,-1) = ", slot(6), slot(9), slot(1), slot(4), slot(-1)
    !write (6,*) "numtypes(6,9,1,4,-1) = ", numtypes(6), numtypes(9), numtypes(1), numtypes(4), numtypes(-1)

    if (seren_format) then
       allocate(porig(1:ptot), parray(1:NDIM+2,1:ptot), v(1:NDIM,1:ptot), &
          & temp(1:ptot), rho(1:ptot), u(1:ptot))
       allocate(sink_array(1:stot))
       do i=1,stot
          allocate(sink_array(i)%tacc(1:dmdt_range))
          allocate(sink_array(i)%macc(1:dmdt_range))
       end do
       sink_slot = 1
    else
       allocate(r(1:NDIM,1:ptot), v(1:NDIM,1:ptot), h(1:ptot), m(1:ptot), temp(1:ptot), rho(1:ptot), ptype(1:ptot))
       ! Allocate ptype array on assumption particles are correctly ordered
       if (numtypes(6) > 0) ptype(1:numtypes(6)) = 6
       if (numtypes(9) > 0) ptype(slot(9):slot(9)+numtypes(9)-1) = 9
       if (numtypes(1) > 0) ptype(slot(1):slot(1)+numtypes(1)-1) = 1
       if (numtypes(4) > 0) ptype(slot(4):slot(4)+numtypes(4)-1) = 4
       if (numtypes(10) > 0) ptype(slot(10):slot(10)+numtypes(10)-1) = 10
       if (numtypes(11) > 0) ptype(slot(11):slot(11)+numtypes(11)-1) = 11
       if (numtypes(12) > 0) ptype(slot(12):slot(12)+numtypes(12)-1) = 12
       if (numtypes(-1) > 0) ptype(slot(-1):ptot) = -1
    end if

    ! Now begin reading files into main arrays

    do i=0,MPInumthreads-1
!       file_ext = ""
!       ! Data output needs extension
!       if (i < 10) then
!         write(file_ext,"(I1)") i
!         file_ext="000"//file_ext
!       else if (i < 100) then
!         write(file_ext,"(I2)") i
!         file_ext="00"//file_ext
!       else if (i < 1000) then
!         write(file_ext,"(I3)") i
!         file_ext="0"//file_ext
!       else if (i < 10000) then
!         write(file_ext,"(I4)") i
!       else if (i < 100000) then
!         write(file_ext,"(A)") "too_big"  ! *** too many threads ***
!       end if

      if (joinmpi) then
        write(file_ext,"(I0)") i

        !tempfilename = trim(adjustl(filename))//"."//trim(adjustl(file_ext))
        tempfilename = trim(adjustl(filename))//".MPI."//trim(adjustl(file_ext))
      else
        tempfilename = nonmpifiles(i)
      end if

      !write(6,*) "tempfilename is ", tempfilename, i, file_ext

      ! Load data from file into main arrays
      if (seren_format) then
         call loadmpifile_serenform(tempfilename)
      else
         call loadmpifile(tempfilename)
      end if
      
    end do

    return
  end subroutine joinsimulation

  subroutine loadmpifile_serenform(file)
    ! Load a MPI thread data file into main arrays
    integer                                   :: j, k               ! Loop counters
    integer                                   :: ierr,ierr_out      ! Error variables
    character (len=200), intent(in)           :: file               ! Filename of MPI thread file
    integer                                   :: threadptot         ! Number of particles in MPI file
    integer                                   :: iheader(1:50)      ! Dragon header integers
    integer (kind=ILP)                        :: ilpheader(1:50)    ! Seren ILP header
    real (kind=PR)                            :: rheader(1:50)      ! Dragon header reals
    real (kind=DP)                            :: dpheader(1:50)     ! Dragon header reals
    integer                                   :: tempslot(-1:12)    ! Temporary slot variable
    logical                                   :: formatted          ! Is file formatted?
    logical                                   :: seren_format       ! Is file formatted?
    integer                                   :: idummy2(1:2)
    logical                                   :: ldummy(1:2)
    real(kind=PR)                             :: sdummy(1:20)
    real(kind=SP)                             :: sdummy_sp(1:20)
    real(kind=DP)                             :: sdummy_dp(1:20)
    integer, allocatable                      :: dummy_int(:)
    real(kind=SP), allocatable                :: dummy_sp(:,:)
    real(kind=DP), allocatable                :: dummy_dp(:,:)
    real(kind=PR), allocatable                :: dummy(:,:)
    real(kind=SP), allocatable                :: dummy_sp_scalar(:)
    real(kind=DP), allocatable                :: dummy_dp_scalar(:)
    real(kind=PR), allocatable                :: dummy_scalar(:)
    real(kind=SP), allocatable                :: macc_sp(:)
    real(kind=DP), allocatable                :: macc_dp(:)
    real(kind=PR), allocatable                :: macc(:)
    real(kind=DP), allocatable                :: tacc(:)
    real(kind=SP), allocatable :: raux(:)
    real(kind=DP), allocatable :: raux_dp(:)
    integer :: sink_data_length
    character(len=30)  :: sink_format_string
    integer :: nl,ni,nli,npr,ndp,nchar
    integer :: width, type_id
    integer :: pfirst
    integer :: plast
    type(sink_node)                           :: sink_dummy

    integer                                   :: unknown_id

    call try_read_seren_header(10,file,iheader,ilpheader,rheader,dpheader,formatted,seren_format)

    if (.NOT. seren_format) then
       write (6,*) "Expecting SEREN format!"
       stop
    end if

    if (new_format) then
       unit_r = 0
       unit_m = 0
       unit_h = 0
    end if

    if (stot>0 .AND. new_format) then
        sink_data_length = 11+NDIMtemp+VDIMtemp+2*dmdt_range
        allocate(raux(1:sink_data_length))
        if (doubleprec) allocate(raux_dp(1:sink_data_length))
        write (sink_format_string,'(A,I0,A)') "(", sink_data_length, "E18.10)"
     end if

    threadptot = iheader(1)

    threadptot = iheader(1)
    t_stot = iheader(2)
    t_pboundary = iheader(3)
    t_picm = iheader(4)
    t_pgas = iheader(5)
    t_pcdm = iheader(6)
    t_pdust = iheader(7)
    t_pion = iheader(8)

    tempslot = slot

    if (stot>0 .AND. new_format) then
       if (dmdt_range>0) then
          allocate(macc_sp(1:dmdt_range))
          allocate(macc_dp(1:dmdt_range))
          allocate(macc(1:dmdt_range))
          allocate(tacc(1:dmdt_range))
       end if
    end if
    allocate(sink_dummy%macc(1:dmdt_range))
    allocate(sink_dummy%tacc(1:dmdt_range))

    ierr_out = 0

    t_is_porig = .FALSE.
    t_is_parray = .FALSE.
    t_is_r = .FALSE.
    t_is_m = .FALSE.
    t_is_h = .FALSE.
    t_is_v = .FALSE.
    t_is_temp = .FALSE.
    t_is_rho = .FALSE.
    t_is_u = .FALSE.
    t_is_sink_v1 = .FALSE.

    do i=1,ndata
        select case (trim(data_id(i)))
           case ("porig")
              if (new_format) then
                 pfirst = typedata(2,i); plast = typedata(3,i)
                 if (pfirst /= 1) stop "Invalid range of particle data for MPI joining"
                 if (plast /= threadptot) stop "Invalid range of particle data for MPI joining"
              end if
              if (allocated(dummy_scalar)) then
                 deallocate(dummy_scalar)
                 if (doubleprec) deallocate(dummy_dp_scalar)
                 if (.NOT. doubleprec) deallocate(dummy_sp_scalar)
              end if
              ! Original particle number
              ! Read through porig numbers
              allocate(dummy_int(1:threadptot))
              t_is_porig = .TRUE.
              if (.NOT.formatted) then
                 read(10,iostat=ierr) dummy_int(1:threadptot)
              else
                 do k=1,threadptot
                    read(10,*,iostat=ierr) dummy_int(k)
                 end do
              end if
              if (ierr /= 0) then
                 print*,' WARNING: errors reading through porig '
                 ierr_out = -1
              end if
              porig(tempslot(6):tempslot(6)+t_pboundary-1) = &
                 & dummy_int(1:t_pboundary)
              porig(tempslot(9):tempslot(9)+t_picm-1) = &
                 & dummy_int(t_pboundary+1:t_pboundary+t_picm)
              porig(tempslot(1):tempslot(1)+t_pgas-1) = &
                 & dummy_int(t_pboundary+t_picm+1:t_pboundary+t_picm+t_pgas)
              porig(tempslot(10):tempslot(10)+t_pcdm-1) = &
                 & dummy_int(t_pboundary+t_picm+t_pgas+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm)
              porig(tempslot(11):tempslot(11)+t_pdust-1) = &
                 & dummy_int(t_pboundary+t_picm+t_pgas+t_pcdm+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust)
              porig(tempslot(12):tempslot(12)+t_pion-1) = &
                 & dummy_int(t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+t_pion)
              deallocate(dummy_int)
           case ("parray")
              if (new_format) then
                 pfirst = typedata(2,i); plast = typedata(3,i)
                 if (pfirst /= 1) stop "Invalid range of particle data for MPI joining"
                 if (plast /= threadptot) stop "Invalid range of particle data for MPI joining"
              end if
              allocate(dummy(1:NDIMtemp+2,1:threadptot))
              if (allocated(dummy_scalar)) then
                 deallocate(dummy_scalar)
                 if (doubleprec) deallocate(dummy_dp_scalar)
                 if (.NOT. doubleprec) deallocate(dummy_sp_scalar)
              end if
              t_is_parray = .TRUE.
              ! Positions, masses and smoothing lengths
              if (doubleprec) allocate(dummy_dp(1:NDIMtemp+2,1:threadptot))
              if (.NOT. doubleprec) allocate(dummy_sp(1:NDIMtemp+2,1:threadptot))
              if (.NOT.formatted) then
                 if (doubleprec) then
                    read(10,iostat=ierr) dummy_dp(1:NDIMtemp+2,1:threadptot)
                    dummy(1:NDIMtemp+2,1:threadptot) = real(dummy_dp(1:NDIMtemp+2,1:threadptot),PR)
                 else
                    read(10,iostat=ierr) dummy_sp(1:NDIMtemp+2,1:threadptot)
                    dummy(1:NDIMtemp+2,1:threadptot) = real(dummy_sp(1:NDIMtemp+2,1:threadptot),PR)
                 end if
              else
                 do k=1,threadptot
                    read(10,*,iostat=ierr) dummy(1:NDIMtemp+2,k)
                 end do
              end if
              parray(1:NDIMtemp+2,tempslot(6):tempslot(6)+t_pboundary-1) = &
                 & dummy(1:NDIMtemp+2,1:t_pboundary)
              parray(1:NDIMtemp+2,tempslot(9):tempslot(9)+t_picm-1) = &
                 & dummy(1:NDIMtemp+2,t_pboundary+1:t_pboundary+t_picm)
              parray(1:NDIMtemp+2,tempslot(1):tempslot(1)+t_pgas-1) = &
                 & dummy(1:NDIMtemp+2,t_pboundary+t_picm+1:t_pboundary+t_picm+t_pgas)
              parray(1:NDIMtemp+2,tempslot(10):tempslot(10)+t_pcdm-1) = &
                 & dummy(1:NDIMtemp+2,t_pboundary+t_picm+t_pgas+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm)
              parray(1:NDIMtemp+2,tempslot(11):tempslot(11)+t_pdust-1) = &
                 & dummy(1:NDIMtemp+2,t_pboundary+t_picm+t_pgas+t_pcdm+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust)
              parray(1:NDIMtemp+2,tempslot(12):tempslot(12)+t_pion-1) = &
                 & dummy(1:NDIMtemp+2,t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+t_pion)
              if (ierr /= 0) then
                 print*,' WARNING: errors reading positions, masses, smoothing lengths '
                 ierr_out = -1
              end if
              if (doubleprec) deallocate(dummy_dp)
              if (.NOT. doubleprec) deallocate(dummy_sp)
              deallocate(dummy)
           case ("r")
              if (new_format) then
                 pfirst = typedata(2,i); plast = typedata(3,i)
                 if (pfirst /= 1) stop "Invalid range of particle data for MPI joining"
                 if (plast /= threadptot) stop "Invalid range of particle data for MPI joining"
                 unit_r = typedata(5,i)
              end if
              allocate(dummy(1:NDIMtemp,1:threadptot))
              if (allocated(dummy_scalar)) then
                 deallocate(dummy_scalar)
                 if (doubleprec) deallocate(dummy_dp_scalar)
                 if (.NOT. doubleprec) deallocate(dummy_sp_scalar)
              end if
              t_is_r = .TRUE.
              if (.NOT. new_format .AND. t_is_m .AND. t_is_h) t_is_parray = .TRUE.
              ! Position (output by SEREN in parray instead)
              if (doubleprec) allocate(dummy_dp(1:NDIMtemp,1:threadptot))
              if (.NOT. doubleprec) allocate(dummy_sp(1:NDIMtemp,1:threadptot))
              if (.NOT.formatted) then
                 if (doubleprec) then
                    read(10,iostat=ierr) dummy_dp(1:NDIMtemp,1:threadptot)
                    dummy(1:NDIMtemp,1:threadptot) = real(dummy_dp(1:NDIMtemp,1:threadptot),PR)
                 else
                    read(10,iostat=ierr) dummy_sp(1:NDIMtemp,1:threadptot)
                    dummy(1:NDIMtemp,1:threadptot) = real(dummy_sp(1:NDIMtemp,1:threadptot),PR)
                 end if
              else
                 do k=1,threadptot
                    read(10,*,iostat=ierr) dummy(1:NDIMtemp,k)
                 end do
              end if
              parray(1:NDIMtemp,tempslot(6):tempslot(6)+t_pboundary-1) = &
                 & dummy(1:NDIMtemp,1:t_pboundary)
              parray(1:NDIMtemp,tempslot(9):tempslot(9)+t_picm-1) = &
                 & dummy(1:NDIMtemp,t_pboundary+1:t_pboundary+t_picm)
              parray(1:NDIMtemp,tempslot(1):tempslot(1)+t_pgas-1) = &
                 & dummy(1:NDIMtemp,t_pboundary+t_picm+1:t_pboundary+t_picm+t_pgas)
              parray(1:NDIMtemp,tempslot(10):tempslot(10)+t_pcdm-1) = &
                 & dummy(1:NDIMtemp,t_pboundary+t_picm+t_pgas+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm)
              parray(1:NDIMtemp,tempslot(11):tempslot(11)+t_pdust-1) = &
                 & dummy(1:NDIMtemp,t_pboundary+t_picm+t_pgas+t_pcdm+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust)
              parray(1:NDIMtemp,tempslot(12):tempslot(12)+t_pion-1) = &
                 & dummy(1:NDIMtemp,t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+t_pion)
              if (ierr /= 0) then
                 print*,' WARNING: errors reading positions '
                 ierr_out = -1
              end if
              if (doubleprec) deallocate(dummy_dp)
              if (.NOT. doubleprec) deallocate(dummy_sp)
              deallocate(dummy)
           case ("h")
              if (new_format) then
                 pfirst = typedata(2,i); plast = typedata(3,i)
                 if (pfirst /= 1) stop "Invalid range of particle data for MPI joining"
                 if (plast /= threadptot) stop "Invalid range of particle data for MPI joining"
                 unit_h = typedata(5,i)
              end if
              if (.NOT. allocated(dummy_scalar)) then
                 allocate(dummy_scalar(1:threadptot))
                 if (doubleprec) allocate(dummy_dp_scalar(1:threadptot))
                 if (.NOT. doubleprec) allocate(dummy_sp_scalar(1:threadptot))
              end if
              t_is_h = .TRUE.
              if (.NOT. new_format .AND. t_is_r .AND. t_is_m) is_parray = .TRUE.
              ! Smoothing lengths (output by SEREN in parray instead)
              if (.NOT.formatted) then
                 if (doubleprec) then
                    read(10,iostat=ierr) dummy_dp_scalar(1:threadptot)
                    dummy_scalar(1:threadptot) = real(dummy_dp_scalar(1:threadptot),PR)
                 else
                    read(10,iostat=ierr) dummy_sp_scalar(1:threadptot)
                    dummy_scalar(1:threadptot) = real(dummy_sp_scalar(1:threadptot),PR)
                 end if
              else
                 do k=1,threadptot
                    read(10,*,iostat=ierr) dummy_scalar(k)
                 end do
              end if
              parray(NDIMtemp+2,tempslot(6):tempslot(6)+t_pboundary-1) = &
                 & dummy_scalar(1:t_pboundary)
              parray(NDIMtemp+2,tempslot(9):tempslot(9)+t_picm-1) = &
                 & dummy_scalar(t_pboundary+1:t_pboundary+t_picm)
              parray(NDIMtemp+2,tempslot(1):tempslot(1)+t_pgas-1) = &
                 & dummy_scalar(t_pboundary+t_picm+1:t_pboundary+t_picm+t_pgas)
              parray(NDIMtemp+2,tempslot(10):tempslot(10)+t_pcdm-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm)
              parray(NDIMtemp+2,tempslot(11):tempslot(11)+t_pdust-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+t_pcdm+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust)
              parray(NDIMtemp+2,tempslot(12):tempslot(12)+t_pion-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+t_pion)
              if (ierr /= 0) then
                 print*,' WARNING: errors reading smoothing lengths'
                 ierr_out = -1
              end if
           case ("m")
              if (new_format) then
                 pfirst = typedata(2,i); plast = typedata(3,i)
                 if (pfirst /= 1) stop "Invalid range of particle data for MPI joining"
                 if (plast /= threadptot) stop "Invalid range of particle data for MPI joining"
                 unit_m = typedata(5,i)
              end if
              if (.NOT. allocated(dummy_scalar)) then
                 allocate(dummy_scalar(1:threadptot))
                 if (doubleprec) allocate(dummy_dp_scalar(1:threadptot))
                 if (.NOT. doubleprec) allocate(dummy_sp_scalar(1:threadptot))
              end if
              t_is_m = .TRUE.
              if (.NOT. new_format .AND. t_is_r .AND. t_is_h) is_parray = .TRUE.
              ! Mass (output by SEREN in parray instead)
              if (.NOT.formatted) then
                 if (doubleprec) then
                    read(10,iostat=ierr) dummy_dp_scalar(1:threadptot)
                    dummy_scalar(1:threadptot) = real(dummy_dp_scalar(1:threadptot),PR)
                 else
                    read(10,iostat=ierr) dummy_sp_scalar(1:threadptot)
                    dummy_scalar(1:threadptot) = real(dummy_sp_scalar(1:threadptot),PR)
                 end if
              else
                 do k=1,threadptot
                    read(10,*,iostat=ierr) dummy_scalar(k)
                 end do
              end if
              parray(NDIMtemp+1,tempslot(6):tempslot(6)+t_pboundary-1) = &
                 & dummy_scalar(1:t_pboundary)
              parray(NDIMtemp+1,tempslot(9):tempslot(9)+t_picm-1) = &
                 & dummy_scalar(t_pboundary+1:t_pboundary+t_picm)
              parray(NDIMtemp+1,tempslot(1):tempslot(1)+t_pgas-1) = &
                 & dummy_scalar(t_pboundary+t_picm+1:t_pboundary+t_picm+t_pgas)
              parray(NDIMtemp+1,tempslot(10):tempslot(10)+t_pcdm-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm)
              parray(NDIMtemp+1,tempslot(11):tempslot(11)+t_pdust-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+t_pcdm+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust)
              parray(NDIMtemp+1,tempslot(12):tempslot(12)+t_pion-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+t_pion)
              if (ierr /= 0) then
                 print*,' WARNING: errors reading masses'
                 ierr_out = -1
              end if
           case ("v")
              if (new_format) then
                 pfirst = typedata(2,i); plast = typedata(3,i)
                 if (pfirst /= 1) stop "Invalid range of particle data for MPI joining"
                 if (plast /= threadptot) stop "Invalid range of particle data for MPI joining"
                 unit_v = typedata(5,i)
              end if
              allocate(dummy(1:VDIMtemp,1:threadptot))
              if (allocated(dummy_scalar)) then
                 deallocate(dummy_scalar)
                 if (doubleprec) deallocate(dummy_dp_scalar)
                 if (.NOT. doubleprec) deallocate(dummy_sp_scalar)
              end if
              t_is_v = .TRUE.
              ! Velocities
              if (doubleprec) allocate(dummy_dp(1:VDIMtemp,1:threadptot))
              if (.NOT. doubleprec) allocate(dummy_sp(1:VDIMtemp,1:threadptot))
              if (.NOT.formatted) then
                 if (doubleprec) then
                    read(10,iostat=ierr) dummy_dp(1:VDIMtemp,1:threadptot)
                    dummy(1:VDIMtemp,1:threadptot) = real(dummy_dp(1:VDIMtemp,1:threadptot),PR)
                 else
                    read(10,iostat=ierr) dummy_sp(1:VDIMtemp,1:threadptot)
                    dummy(1:VDIMtemp,1:threadptot) = real(dummy_sp(1:VDIMtemp,1:threadptot),PR)
                 end if
              else
                 do k=1,threadptot
                    read(10,*,iostat=ierr) dummy(1:VDIMtemp,k)
                 end do
              end if
              v(1:VDIMtemp,tempslot(6):tempslot(6)+t_pboundary-1) = &
                 & dummy(1:VDIMtemp,1:t_pboundary)
              v(1:VDIMtemp,tempslot(9):tempslot(9)+t_picm-1) = &
                 & dummy(1:VDIMtemp,t_pboundary+1:t_pboundary+t_picm)
              v(1:VDIMtemp,tempslot(1):tempslot(1)+t_pgas-1) = &
                 & dummy(1:VDIMtemp,t_pboundary+t_picm+1:t_pboundary+t_picm+t_pgas)
              v(1:NDIMtemp,tempslot(10):tempslot(10)+t_pcdm-1) = &
                 & dummy(1:VDIMtemp,t_pboundary+t_picm+t_pgas+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm)
              v(1:NDIMtemp,tempslot(11):tempslot(11)+t_pdust-1) = &
                 & dummy(1:VDIMtemp,t_pboundary+t_picm+t_pgas+t_pcdm+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust)
              v(1:NDIMtemp,tempslot(12):tempslot(12)+t_pion-1) = &
                 & dummy(1:VDIMtemp,t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+t_pion)
              if (ierr /= 0) then
                 print*,' WARNING: errors reading velocities '
                 ierr_out = -1
              end if
              if (doubleprec) deallocate(dummy_dp)
              if (.NOT. doubleprec) deallocate(dummy_sp)
              deallocate(dummy)
           case ("rho")
              if (new_format) then
                 pfirst = typedata(2,i); plast = typedata(3,i)
                 if (pfirst /= 1) stop "Invalid range of particle data for MPI joining"
                 if (plast /= threadptot) stop "Invalid range of particle data for MPI joining"
                 unit_rho = typedata(5,i)
              end if
              if (.NOT. allocated(dummy_scalar)) then
                 allocate(dummy_scalar(1:threadptot))
                 if (doubleprec) allocate(dummy_dp_scalar(1:threadptot))
                 if (.NOT. doubleprec) allocate(dummy_sp_scalar(1:threadptot))
              end if
              t_is_rho = .TRUE.
              ! Densities
              if (.NOT.formatted) then
                 if (doubleprec) then
                    read(10,iostat=ierr) dummy_dp_scalar(1:threadptot)
                    dummy_scalar(1:threadptot) = real(dummy_dp_scalar(1:threadptot),PR)
                 else
                    read(10,iostat=ierr) dummy_sp_scalar(1:threadptot)
                    dummy_scalar(1:threadptot) = real(dummy_sp_scalar(1:threadptot),PR)
                 end if
              else
                 do k=1,threadptot
                    read(10,*,iostat=ierr) dummy_scalar(k)
                 end do
              end if
              rho(tempslot(6):tempslot(6)+t_pboundary-1) = &
                 & dummy_scalar(1:t_pboundary)
              rho(tempslot(9):tempslot(9)+t_picm-1) = &
                 & dummy_scalar(t_pboundary+1:t_pboundary+t_picm)
              rho(tempslot(1):tempslot(1)+t_pgas-1) = &
                 & dummy_scalar(t_pboundary+t_picm+1:t_pboundary+t_picm+t_pgas)
              rho(tempslot(10):tempslot(10)+t_pcdm-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm)
              rho(tempslot(11):tempslot(11)+t_pdust-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+t_pcdm+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust)
              rho(tempslot(12):tempslot(12)+t_pion-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+t_pion)
              if (ierr /= 0) then
                 print*,' WARNING: errors reading densities'
                 ierr_out = -1
              end if
           case ("temp")
              if (new_format) then
                 pfirst = typedata(2,i); plast = typedata(3,i)
                 if (pfirst /= 1) stop "Invalid range of particle data for MPI joining"
                 if (plast /= threadptot) stop "Invalid range of particle data for MPI joining"
                 unit_temp = typedata(5,i)
              end if
              if (.NOT. allocated(dummy_scalar)) then
                 allocate(dummy_scalar(1:threadptot))
                 if (doubleprec) allocate(dummy_dp_scalar(1:threadptot))
                 if (.NOT. doubleprec) allocate(dummy_sp_scalar(1:threadptot))
              end if
              t_is_temp = .TRUE.
              ! Temperatures
              if (.NOT.formatted) then
                 if (doubleprec) then
                    read(10,iostat=ierr) dummy_dp_scalar(1:threadptot)
                    dummy_scalar(1:threadptot) = real(dummy_dp_scalar(1:threadptot),PR)
                 else
                    read(10,iostat=ierr) dummy_sp_scalar(1:threadptot)
                    dummy_scalar(1:threadptot) = real(dummy_sp_scalar(1:threadptot),PR)
                 end if
              else
                 do k=1,threadptot
                    read(10,*,iostat=ierr) dummy_scalar(k)
                 end do
              end if
              temp(tempslot(6):tempslot(6)+t_pboundary-1) = &
                 & dummy_scalar(1:t_pboundary)
              temp(tempslot(9):tempslot(9)+t_picm-1) = &
                 & dummy_scalar(t_pboundary+1:t_pboundary+t_picm)
              temp(tempslot(1):tempslot(1)+t_pgas-1) = &
                 & dummy_scalar(t_pboundary+t_picm+1:t_pboundary+t_picm+t_pgas)
              temp(tempslot(10):tempslot(10)+t_pcdm-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm)
              temp(tempslot(11):tempslot(11)+t_pdust-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+t_pcdm+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust)
              temp(tempslot(12):tempslot(12)+t_pion-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+t_pion)
              if (ierr /= 0) then
                 print*,' WARNING: errors reading temperatures'
                 ierr_out = -1
              end if
           case ("u")
              if (new_format) then
                 pfirst = typedata(2,i); plast = typedata(3,i)
                 if (pfirst /= 1) stop "Invalid range of particle data for MPI joining"
                 if (plast /= threadptot) stop "Invalid range of particle data for MPI joining"
                 unit_u = typedata(5,i)
              end if
              if (.NOT. allocated(dummy_scalar)) then
                 allocate(dummy_scalar(1:threadptot))
                 if (doubleprec) allocate(dummy_dp_scalar(1:threadptot))
                 if (.NOT. doubleprec) allocate(dummy_sp_scalar(1:threadptot))
              end if
              t_is_u = .TRUE.
              ! Internal energy
              if (.NOT.formatted) then
                 if (doubleprec) then
                    read(10,iostat=ierr) dummy_dp_scalar(1:threadptot)
                    dummy_scalar(1:threadptot) = real(dummy_dp_scalar(1:threadptot),PR)
                 else
                    read(10,iostat=ierr) dummy_sp_scalar(1:threadptot)
                    dummy_scalar(1:threadptot) = real(dummy_sp_scalar(1:threadptot),PR)
                 end if
              else
                 do k=1,threadptot
                    read(10,*,iostat=ierr) dummy_scalar(k)
                 end do
              end if
              u(tempslot(6):tempslot(6)+t_pboundary-1) = &
                 & dummy_scalar(1:t_pboundary)
              u(tempslot(9):tempslot(9)+t_picm-1) = &
                 & dummy_scalar(t_pboundary+1:t_pboundary+t_picm)
              u(tempslot(1):tempslot(1)+t_pgas-1) = &
                 & dummy_scalar(t_pboundary+t_picm+1:t_pboundary+t_picm+t_pgas)
              u(tempslot(10):tempslot(10)+t_pcdm-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm)
              u(tempslot(11):tempslot(11)+t_pdust-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+t_pcdm+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust)
              u(tempslot(12):tempslot(12)+t_pion-1) = &
                 & dummy_scalar(t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+1:&
                 &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+t_pion)
              if (ierr /= 0) then
                 print*,' WARNING: errors reading internal energy'
                 ierr_out = -1
              end if
           case ("B")
              ! Magnetic fields
              stop "No support for magnetic fields yet!"
           case ("sink_v1")
              if (allocated(dummy_scalar)) then
                 deallocate(dummy_scalar)
                 if (doubleprec) deallocate(dummy_dp_scalar)
                 if (.NOT. doubleprec) deallocate(dummy_sp_scalar)
              end if
              if (new_format) then
                 pfirst = typedata(2,i); plast = typedata(3,i)
                 if (pfirst /= 1) stop "didn't write this properly; all sinks should be in one file!"
                 if (plast /= t_stot) stop "didn't write this properly; all sinks should be in one file!"
                 if (.NOT.formatted) then
                    read(10,iostat=ierr) nl,ni,nli,npr,ndp,nchar
                 else
                    read(10,'(6I10)',iostat=ierr) nl,ni,nli,npr,ndp,nchar
                 end if
              end if
              t_is_sink_v1 = .TRUE.
              ! Load sinks in sink data storage, will add them later
              do j=1,t_stot
                 if (new_format) then
                    if (.NOT. formatted) then
                       read(10,iostat=ierr) ldummy
                       read(10,iostat=ierr) idummy2
                       if (doubleprec) then
                          read(10,iostat=ierr) raux_dp
                          raux = real(raux_dp)
                       else
                          read(10,iostat=ierr) raux
                       end if
                    else
                       read(10,'(2L1)',iostat=ierr) ldummy
                       read(10,'(2I8)',iostat=ierr) idummy2
                       read(10,sink_format_string) raux(1:sink_data_length)
                    end if
                    sink_dummy%id          = idummy2(1)
                    sink_dummy%ncreate     = idummy2(2)
                    sink_dummy%accrete     = ldummy(1)
                    sink_dummy%static      = ldummy(2)
                    sink_dummy%tcreate     = real(raux(1),DP)
                    sink_dummy%r(1:NDIMtemp)   = raux(2:NDIMtemp+1)
                    sink_dummy%v(1:NDIMtemp)   = raux(NDIMtemp+2:NDIMtemp+VDIMtemp+1)
                    sink_dummy%m           = raux(NDIMtemp+VDIMtemp+2)
                    sink_dummy%h           = raux(NDIMtemp+VDIMtemp+3)
                    sink_dummy%radius      = raux(NDIMtemp+VDIMtemp+4)
                    sink_dummy%angmom(1:3) = raux(NDIMtemp+VDIMtemp+5:NDIMtemp+VDIMtemp+7)
                    sink_dummy%dmdt        = raux(NDIMtemp+VDIMtemp+8)
                    sink_dummy%star_radius = raux(NDIMtemp+VDIMtemp+9)
                    sink_dummy%luminosity  = raux(NDIMtemp+VDIMtemp+10)
                    sink_dummy%temperature = raux(NDIMtemp+VDIMtemp+11)
                    if (dmdt_range > 0) then
                       sink_dummy%macc(1:dmdt_range) = &
                            &real(raux(NDIMtemp+VDIMtemp+12:&
                            &NDIMtemp+VDIMtemp+11+dmdt_range),DP)
                       sink_dummy%tacc(1:dmdt_range) = &
                            &real(raux(NDIMtemp+VDIMtemp+12+dmdt_range:&
                            &NDIMtemp+VDIMtemp+11+2*dmdt_range),DP)
                    end if
                 else
                    if (.NOT.formatted) then
                       read(10,iostat=ierr) idummy2
                       read(10,iostat=ierr) ldummy
                       if (doubleprec) then
                          read(10,iostat=ierr) sdummy_dp
                          sdummy = real(sdummy_dp,PR)
                          if (dmdt_range > 0) read(10,iostat=ierr) macc_dp
                          macc = real(macc_dp,PR)
                       else
                          read(10,iostat=ierr) sdummy_sp
                          sdummy = real(sdummy_sp,PR)
                          if (dmdt_range > 0) read(10,iostat=ierr) macc_sp
                          macc = real(macc_sp,PR)
                       end if
                       if (dmdt_range > 0) read(10,iostat=ierr) tacc
                       sink_dummy%id          = idummy2(1)
                       sink_dummy%ncreate     = idummy2(2)
                       sink_dummy%accrete     = ldummy(1)
                       sink_dummy%static      = ldummy(2)
                       sink_dummy%tcreate     = real(sdummy(1),DP)
                       sink_dummy%r(1:NDIM)   = sdummy(2:1+NDIM)
                       sink_dummy%v(1:NDIM)   = sdummy(5:4+NDIM)
                       sink_dummy%m           = sdummy(8)
                       sink_dummy%h           = sdummy(9)
                       sink_dummy%radius      = sdummy(10)
                       sink_dummy%angmom(1:3) = real(sdummy(11:13),DP)
                       sink_dummy%dmdt        = real(sdummy(14),DP)
                       sink_dummy%star_radius = real(sdummy(15),DP)
                       sink_dummy%luminosity  = real(sdummy(16),DP)
                       sink_dummy%temperature = real(sdummy(17),DP)
                       do k=1,dmdt_range
                          sink_dummy%macc(k) = macc(k)
                          sink_dummy%tacc(k) = tacc(k)
                       end do
                    else
                       read(10,'(2I8)') sink_dummy%id, sink_dummy%ncreate
                       read(10,'(2L1)') sink_dummy%accrete, sink_dummy%static
                       read(10,'(E18.10)') sink_dummy%tcreate
                       read(10,'(3E18.10)') sink_dummy%r(1:NDIMtemp)
                       read(10,'(3E18.10)') sink_dummy%v(1:VDIMtemp)
                       read(10,'(E18.10)') sink_dummy%m
                       read(10,'(E18.10)') sink_dummy%h
                       read(10,'(E18.10)') sink_dummy%radius
                       read(10,'(3E18.10)') sink_dummy%angmom(1:3)
                       read(10,'(E18.10)') sink_dummy%dmdt
                       read(10,'(E18.10)') sink_dummy%star_radius
                       read(10,'(E18.10)') sink_dummy%luminosity
                       read(10,'(E18.10)') sink_dummy%temperature
                       do k=1,dmdt_range
                          read(10,'(2E18.10)') sink_dummy%macc(k), sink_dummy%tacc(k)
                       end do
                    end if
                 end if
!                  if (joinmpi) then
!                     sink_array(sink_dummy%id) = sink_dummy
!                  else
                    sink_array(sink_slot) = sink_dummy
                    sink_slot = sink_slot + 1
!                  end if
              end do
           case default
              if (new_format) then
                 pfirst = typedata(2,i); plast = typedata(3,i)
                 if (pfirst /= 1) stop "Invalid range of particle data for MPI joining"
                 if (plast /= threadptot) stop "Invalid range of particle data for MPI joining"
                 if (typedata(4,i) == 7) stop "Too lazy to have written structure skipping code yet!"
                 if (typedata(4,i) /= 2 .AND. typedata(4,i) /= 4) &
                    stop "Too lazy to have written code for non-PR or integer data yet!"
                 width = typedata(1,i)
                 type_id = typedata(4,i)
              else
                 type_id = 4
                 width = 1
              end if
              ! Unknown data type; find if we already have it
              unknown_id = 0
              do k=1,nunknown
                 if (trim(unknown_data(k)%data_id) == trim(data_id(i))) then
                    unknown_id = k
                    exit
                 end if
              end do
              if (unknown_id == 0) then
                 nunknown = nunknown + 1
                 unknown_id = nunknown
                 if (type_id==4) then
                    if (width==1) then
                        allocate(unknown_data(unknown_id)%data(1:ptot))
                        unknown_data(unknown_id)%data = 0._PR
                    else
                        allocate(unknown_data(unknown_id)%vector(1:width,1:ptot))
                        unknown_data(unknown_id)%vector = 0._PR
                    end if
                 else if (type_id==2) then
                    allocate(unknown_data(unknown_id)%idata(1:ptot))
                    unknown_data(unknown_id)%idata = 0
                 else
                    stop "Haven't done this yet!"
                 end if
                 unknown_data(unknown_id)%data_id = trim(data_id(i))
                 unknown_data(unknown_id)%width = width
                 unknown_data(unknown_id)%type_id = type_id
              end if
              ! Read unknown data type
              if (type_id==4) then
                 if (width==1) then
                    if (.NOT. allocated(dummy_scalar)) then
                       allocate(dummy_scalar(1:threadptot))
                       if (doubleprec) allocate(dummy_dp_scalar(1:threadptot))
                       if (.NOT. doubleprec) allocate(dummy_sp_scalar(1:threadptot))
                    end if
                    if (.NOT.formatted) then
                       if (doubleprec) then
                          read(10,iostat=ierr) dummy_dp_scalar(1:threadptot)
                          dummy_scalar(1:threadptot) = real(dummy_dp_scalar(1:threadptot),PR)
                       else
                          read(10,iostat=ierr) dummy_sp_scalar(1:threadptot)
                          dummy_scalar(1:threadptot) = real(dummy_sp_scalar(1:threadptot),PR)
                       end if
                    else
                       do k=1,threadptot
                       read(10,*,iostat=ierr) dummy_scalar(k)
                       end do
                    end if
                    unknown_data(unknown_id)%data(tempslot(6):tempslot(6)+t_pboundary-1) = &
                       & dummy_scalar(1:t_pboundary)
                    unknown_data(unknown_id)%data(tempslot(9):tempslot(9)+t_picm-1) = &
                       & dummy_scalar(t_pboundary+1:t_pboundary+t_picm)
                    unknown_data(unknown_id)%data(tempslot(1):tempslot(1)+t_pgas-1) = &
                       & dummy_scalar(t_pboundary+t_picm+1:t_pboundary+t_picm+t_pgas)
                    unknown_data(unknown_id)%data(tempslot(10):tempslot(10)+t_pcdm-1) = &
                       & dummy_scalar(t_pboundary+t_picm+t_pgas+1:&
                       &t_pboundary+t_picm+t_pgas+t_pcdm)
                    unknown_data(unknown_id)%data(tempslot(11):tempslot(11)+t_pdust-1) = &
                       & dummy_scalar(t_pboundary+t_picm+t_pgas+t_pcdm+1:&
                       &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust)
                    unknown_data(unknown_id)%data(tempslot(12):tempslot(12)+t_pion-1) = &
                       & dummy_scalar(t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+1:&
                       &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+t_pion)
                 else
                    if (.NOT. allocated(dummy)) then
                       allocate(dummy(1:width,1:threadptot))
                       if (doubleprec) allocate(dummy_dp(1:width,1:threadptot))
                       if (.NOT. doubleprec) allocate(dummy_sp(1:width,1:threadptot))
                    end if
                    if (.NOT.formatted) then
                       if (doubleprec) then
                          read(10,iostat=ierr) dummy_dp(1:width,1:threadptot)
                          dummy(1:width,1:threadptot) = real(dummy_dp(1:width,1:threadptot),PR)
                       else
                          read(10,iostat=ierr) dummy_sp(1:width,1:threadptot)
                          dummy(1:width,1:threadptot) = real(dummy_sp(1:width,1:threadptot),PR)
                       end if
                    else
                       do k=1,threadptot
                          read(10,*,iostat=ierr) dummy(1:width,k)
                       end do
                    end if
                    unknown_data(unknown_id)%vector(1:width,tempslot(6):tempslot(6)+t_pboundary-1) = &
                       & dummy(1:width,1:t_pboundary)
                    unknown_data(unknown_id)%vector(1:width,tempslot(9):tempslot(9)+t_picm-1) = &
                       & dummy(1:width,t_pboundary+1:t_pboundary+t_picm)
                    unknown_data(unknown_id)%vector(1:width,tempslot(1):tempslot(1)+t_pgas-1) = &
                       & dummy(1:width,t_pboundary+t_picm+1:t_pboundary+t_picm+t_pgas)
                    unknown_data(unknown_id)%vector(1:width,tempslot(10):tempslot(10)+t_pcdm-1) = &
                       & dummy(1:width,t_pboundary+t_picm+t_pgas+1:&
                       &t_pboundary+t_picm+t_pgas+t_pcdm)
                    unknown_data(unknown_id)%vector(1:width,tempslot(11):tempslot(11)+t_pdust-1) = &
                       & dummy(1:width,t_pboundary+t_picm+t_pgas+t_pcdm+1:&
                       &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust)
                    unknown_data(unknown_id)%vector(1:width,tempslot(12):tempslot(12)+t_pion-1) = &
                       & dummy(1:width,t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+1:&
                       &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+t_pion)
                    if (doubleprec) deallocate(dummy_dp)
                    if (.NOT. doubleprec) deallocate(dummy_sp)
                    deallocate(dummy)
                 end if
              else if (type_id==2) then
                 if (.NOT. allocated(dummy_int)) then
                    allocate(dummy_int(1:threadptot))
                 end if
                 if (.NOT.formatted) then
                    read(10,iostat=ierr) dummy_int(1:threadptot)
                 else
                    do k=1,threadptot
                       read(10,*,iostat=ierr) dummy_int(k)
                    end do
                 end if
                 unknown_data(unknown_id)%idata(tempslot(6):tempslot(6)+t_pboundary-1) = &
                    & dummy_int(1:t_pboundary)
                 unknown_data(unknown_id)%idata(tempslot(9):tempslot(9)+t_picm-1) = &
                    & dummy_int(t_pboundary+1:t_pboundary+t_picm)
                 unknown_data(unknown_id)%idata(tempslot(1):tempslot(1)+t_pgas-1) = &
                    & dummy_int(t_pboundary+t_picm+1:t_pboundary+t_picm+t_pgas)
                 unknown_data(unknown_id)%idata(tempslot(10):tempslot(10)+t_pcdm-1) = &
                    & dummy_int(t_pboundary+t_picm+t_pgas+1:&
                    &t_pboundary+t_picm+t_pgas+t_pcdm)
                 unknown_data(unknown_id)%idata(tempslot(11):tempslot(11)+t_pdust-1) = &
                    & dummy_int(t_pboundary+t_picm+t_pgas+t_pcdm+1:&
                    &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust)
                 unknown_data(unknown_id)%idata(tempslot(12):tempslot(12)+t_pion-1) = &
                    & dummy_int(t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+1:&
                    &t_pboundary+t_picm+t_pgas+t_pcdm+t_pdust+t_pion)
                 deallocate(dummy_int)
              else
                 stop "Haven't done this bit yet!"
              end if
              if (ierr /= 0) then
                 print*,' WARNING: errors reading unknown data type ', trim(data_id(i))
                 ierr_out = -1
              end if
        end select
     end do

     if (.NOT. t_is_porig) is_porig = .FALSE.
     if (.NOT.t_is_parray) is_parray = .FALSE.
     if (.NOT.t_is_v) is_v = .FALSE.
     if (.NOT.t_is_temp) is_temp = .FALSE.
     if (.NOT.t_is_rho) is_rho = .FALSE.
     if (.NOT.t_is_u) is_u = .FALSE.
     if (t_is_sink_v1) is_sink_v1 = .TRUE.

     if (doubleprec.AND.allocated(dummy_dp_scalar)) deallocate(dummy_dp_scalar)
     if ((.NOT. doubleprec).AND.allocated(dummy_dp_scalar)) deallocate(dummy_sp_scalar)
     if (allocated(dummy_scalar)) deallocate(dummy_scalar)
     tempslot(1) = tempslot(1) + t_pgas
     tempslot(6) = tempslot(6) + t_pboundary
     tempslot(9) = tempslot(9) + t_picm
     tempslot(10) = tempslot(10) + t_pcdm
     tempslot(11) = tempslot(11) + t_pdust
     tempslot(12) = tempslot(12) + t_pion

     ! Don't need to read particle types

     slot = tempslot

     close (10)

     return
  end subroutine loadmpifile_serenform

  subroutine loadmpifile(file)
    ! Load a MPI thread data file into main arrays
    integer                                   :: p                  ! Loop counters
    character (len=200), intent(in)           :: file               ! Filename of MPI thread file
    integer                                   :: threadptot         ! Number of particles in MPI file
    integer                                   :: iheader(1:50)      ! Dragon header integers
    real (kind=PR)                            :: rheader(1:50)      ! Dragon header reals
    real (kind=PR)                            :: rdummy(1:NDIM)     ! Dragon dummy real
    integer                                   :: tempslot(-1:12)    ! Temporary slot variable
    integer, allocatable                      :: threadptype(:)     ! ptype for this thread only
    logical                                   :: formatted          ! Is file formatted?
    integer, allocatable :: idummy1(:)         ! dummy array for integers
    real(kind=PR), allocatable :: rdummy1(:)   ! dummy array for reals
    real(kind=PR), allocatable :: rdummy3(:,:) ! dummy array for 3D reals

    call openandreadheader(10,file,iheader(1:20),rheader,formatted)

    threadptot = iheader(1)

    allocate(threadptype(1:threadptot))

    if (formatted) then

      do p=1,threadptot*2
         read(10,*) rdummy(1:NDIM)
      end do
      do p=1,threadptot*4
         read(10,*) rdummy(1)
      end do
      do p=1,threadptot
         read(10,*) threadptype(p)
      end do

      rewind(10)

      do i=1,20
         read (10,*) iheader(i)     ! Read dragon integer header
      end do

      do i=1,50
         read (10,*) rheader(i)     ! Read dragon real header
      end do

      ! Read positions
      tempslot = slot
      do p=1,threadptot
         read(10,*) r(1:NDIM,tempslot(threadptype(p)))
         tempslot(threadptype(p)) = tempslot(threadptype(p)) + 1
      end do

      ! Read velocities
      tempslot = slot
      do p=1,threadptot
         read(10,*) v(1:NDIM,tempslot(threadptype(p)))
         tempslot(threadptype(p)) = tempslot(threadptype(p)) + 1
      end do

      ! Read temperatures
      tempslot = slot
      do p=1,threadptot
         read(10,*) temp(tempslot(threadptype(p)))    ! Load all particles of this type
         tempslot(threadptype(p)) = tempslot(threadptype(p)) + 1
      end do

      ! Read smoothing lengths
      tempslot = slot
      do p=1,threadptot
         read(10,*) h(tempslot(threadptype(p)))    ! Load all particles of this type
         tempslot(threadptype(p)) = tempslot(threadptype(p)) + 1
      end do

      ! Read densities
      tempslot = slot
      do p=1,threadptot
         read(10,*) rho(tempslot(threadptype(p)))    ! Load all particles of this type
         tempslot(threadptype(p)) = tempslot(threadptype(p)) + 1
      end do

      ! Read particle masses
      tempslot = slot
      do p=1,threadptot
         read(10,*) m(tempslot(threadptype(p)))    ! Load all particles of this type
         tempslot(threadptype(p)) = tempslot(threadptype(p)) + 1
      end do

    else

      allocate(rdummy1(1:threadptot),rdummy3(1:NDIM,1:threadptot),idummy1(1:threadptot))

! First pass to get numbers of particles of each type
! --------------------------------------------------------------------------
      read(10) rdummy3  ! Positions
      read(10) rdummy3  ! Velocities
      read(10) rdummy1  ! Temperatures
      read(10) rdummy1  ! Smoothing lengths
      read(10) rdummy1  ! Densities
      read(10) rdummy1  ! Masses
      read(10) idummy1  ! Particle types

      do p=1,threadptot
         threadptype(p) = idummy1(p)
      end do

      rewind(10)

      read (10) iheader(1:20) ! Read dragon integer header
      read (10) rheader       ! Read dragon real header

      ! Read positions
      read(10) rdummy3
      tempslot = slot
      do p=1,threadptot
         r(1:NDIM,tempslot(threadptype(p))) = rdummy3(1:NDIM,p)
         tempslot(threadptype(p)) = tempslot(threadptype(p)) + 1
      end do

      ! Read velocities
      read(10) rdummy3
      tempslot = slot
      do p=1,threadptot
         v(1:NDIM,tempslot(threadptype(p))) = rdummy3(1:NDIM,p)
         tempslot(threadptype(p)) = tempslot(threadptype(p)) + 1
      end do

      ! Read temperatures
      read(10) rdummy1
      tempslot = slot
      do p=1,threadptot
         temp(tempslot(threadptype(p))) = rdummy1(p)    ! Load all particles of this type
         tempslot(threadptype(p)) = tempslot(threadptype(p)) + 1
      end do

      ! Read smoothing lengths
      read(10) rdummy1
      tempslot = slot
      do p=1,threadptot
         h(tempslot(threadptype(p))) = rdummy1(p)    ! Load all particles of this type
         tempslot(threadptype(p)) = tempslot(threadptype(p)) + 1
      end do

      ! Read densities
      read(10) rdummy1
      tempslot = slot
      do p=1,threadptot
         rho(tempslot(threadptype(p))) = rdummy1(p)    ! Load all particles of this type
         tempslot(threadptype(p)) = tempslot(threadptype(p)) + 1
      end do

      ! Read particle masses
      read(10) rdummy1
      tempslot = slot
      do p=1,threadptot
         m(tempslot(threadptype(p))) = rdummy1(p)    ! Load all particles of this type
         tempslot(threadptype(p)) = tempslot(threadptype(p)) + 1
      end do


      deallocate(rdummy1,rdummy3,idummy1)

    end if

    ! Don't need to read particle types

    slot = tempslot

    close (10)

    deallocate(threadptype)

    return
  end subroutine loadmpifile

end module progmodule