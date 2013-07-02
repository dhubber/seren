! READ_DATA_ASCII.F90
! D. A. Hubber - 17/06/2010
! Reads in initial conditions file in simple column ASCII format.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_data_ascii(in_file,decomp_read)
  use particle_module
  use hydro_module
  use type_module
  use sink_module
  use time_module
  use scaling_module
  implicit none

  character(len=*), intent(in) :: in_file   ! Name of ascii file
  logical, intent(in)          :: decomp_read ! MPI decomposition read

  character(len=20) :: auxdata         ! aux. data id variable
  character(len=20) :: data_id(1:100)  ! list of data identifiers 
  logical :: tflag                     ! temperature flag
  logical :: uflag                     ! internal energy flag
  integer :: boundaryslot              ! point in list to boundary particle
  integer :: cdmslot                   ! point in list to cdm particle
  integer :: gasslot                   ! point in list to insert particle
  integer :: i                         ! aux. loop variable
  integer :: icmslot                   ! point in list to icm particle
  integer :: j                         ! aux. loop variable
  integer :: ndata                     ! no. of data columns
  integer :: p                         ! counter to loop over particles
  integer :: psplit                    ! counter for particle splitting
  integer :: pdead                     ! counter for accreted particles
  integer :: auxtype                   ! particle type
  integer :: s                         ! sink id
  integer :: sinkslot                  ! point in list sink particle
  real(kind=PR) :: pdata(1:100)        ! aux. particle data array

  debug1("Reading in data file : "//trim(in_file)//" [read_data_ascii.F90]")

! Initialise counters
  ndata     = 0
  pboundary = 0
  psplit    = 0
  pdead     = 0
  picm      = 0
  pgas      = 0
  pcdm      = 0
  pdust     = 0
  pion      = 0
  stot      = 0

! First, read-in ascii column data
  open(1, file="asciicolumns.dat", status="old", form="formatted")
  do
     read(1,*,err=5,end=5) auxdata
     ndata = ndata + 1
     data_id(ndata) = auxdata
  end do
5 close(5)

! Check there are a sufficient no. of data columns, or ptype is not 1st column
  if (ndata <= 1) stop "Invalid no. of data columns"
  if (data_id(1) /= "ptype") stop "First column is not particle type id"
  open(1, file=in_file, status="old", form="formatted")

! First, count the total no. of particles and the no. of each particle type
! ----------------------------------------------------------------------------
  do
     read(1,*,err=10,end=10) auxtype,pdata(2:ndata)
     if (auxtype == -1) stot = stot + 1
     if (auxtype == 0)  pdead = pdead + 1
     if (auxtype == 1)  pgas = pgas + 1
     if (auxtype == 4)  psplit = psplit + 1
     if (auxtype == 6)  pboundary = pboundary + 1
     if (auxtype == 9)  picm = picm + 1
     if (auxtype == 10) pcdm = pcdm + 1
     if (auxtype >= 1)  ptot = ptot + 1
  end do
! ----------------------------------------------------------------------------

10 rewind(1)

  write(6,*) "Particles   = ", ptot, "   Sinks = ", stot
  write(6,*) "Gas         = ", pgas
  write(6,*) "Boundary    = ", pboundary
  write(6,*) "Intercloud  = ", picm
  write(6,*) "Dark matter = ", pcdm
  write(6,*) "Splitting   = ", psplit
  write(6,*) "Accreted    = ", pdead
  if (psplit /= 0) stop "Fatal error: particle splitting not supported"
  if (ptot /= (pgas + picm + pboundary + pcdm)) &
       &stop "Fatal error: particles do not add up"

! If MPI decomposition read, allocate 'minimal = TRUE'
  call allocate_memory(decomp_read)

! Initialise aux. id variables
  boundaryslot = 0
  icmslot      = pboundary 
  gasslot      = pboundary + picm
  cdmslot      = pboundary + picm + pgas
  sinkslot     = 0

! Now read all particle data into main arrays
! ============================================================================
  do j=1,ptot+stot
     read(1,*,err=20,end=20) auxtype,pdata(2:ndata)

     ! Determine array location of particle depending on type
     if (auxtype == 6) then
        boundaryslot = boundaryslot + 1
        p = boundaryslot
     else if (auxtype == 9) then
        icmslot = icmslot + 1
        p = icmslot
     else if (auxtype == 1) then
        gasslot = gasslot + 1
        p = gasslot
     else if (auxtype == 10) then
        cdmslot = cdmslot + 1
        p = cdmslot
     else if (auxtype == -1) then
        sinkslot = sinkslot + 1
        s = sinkslot
        p = -s
     else
        write(6,*) "Unrecognised particle type in read_data_ascii.F90 :",auxtype
        stop
     end if

     ! If particle is an SPH particle
     ! -----------------------------------------------------------------------
     if (p > 0) then
        do i=2,ndata
           if (data_id(i) == "x") then
              sph(p)%r(1) = pdata(i)
#if NDIM==2 || NDIM==3
           else if (data_id(i) == "y" .and. NDIM > 1) then
              sph(p)%r(2) = pdata(i)
#endif
#if NDIM==3
           else if (data_id(i) == "z" .and. NDIM == 3) then
              sph(p)%r(3) = pdata(i)
#endif
           else if (data_id(i) == "m") then
              sph(p)%m = pdata(i)
           else if (data_id(i) == "h") then
              sph(p)%h = pdata(i)
           else if (data_id(i) == "vx") then
              sph(p)%v(1) = pdata(i)
#if VDIM==2 || VDIM==3
           else if (data_id(i) == "vy" .and. VDIM > 1) then
              sph(p)%v(2) = pdata(i)
#endif
#if VDIM==3
           else if (data_id(i) == "vz" .and. VDIM == 3) then
              sph(p)%v(3) = pdata(i)
#endif
           else if (data_id(i) == "rho") then
              sph(p)%rho = pdata(i)
#if defined(HYDRO)
           else if (data_id(i) == "temp") then
              sph(p)%temp = pdata(i)
#if defined(INTERNAL_ENERGY)
           else if (data_id(i) == "u") then
              sph(p)%u = pdata(i)
#endif
#endif
           end if
        end do

        ! If the temperature is read in but not the internal energy, 
        ! calculate the internal energy from the temperature here.
        uflag = .false.
        tflag = .false.
        do i=2,ndata
           if (data_id(i) == "u") uflag = .true.
           if (data_id(i) == "temp") tflag = .true.
        end do
        if (tflag .and. (.not. uflag)) then
           sph(p)%u = Pconst*sph(p)%temp/(gamma - 1.0_PR)*uscale
        end if

     ! Else, if a sink particle
     ! -----------------------------------------------------------------------
#if defined(SINKS)
     else if (p < 0) then
        do i=2,ndata
           if (data_id(i) == "x") then
              sink(s)%r(1) = pdata(i)
#if NDIM==2 || NDIM==3
           else if (data_id(i) == "y" .and. NDIM > 1) then
              sink(s)%r(2) = pdata(i)
#endif
#if NDIM==3
           else if (data_id(i) == "z" .and. NDIM == 3) then
              sink(s)%r(3) = pdata(i)
#endif
           else if (data_id(i) == "m") then
              sink(s)%m = pdata(i)
           else if (data_id(i) == "h") then
              sink(s)%h = pdata(i)
           else if (data_id(i) == "vx") then
              sink(s)%v(1) = pdata(i)
#if VDIM==2 || VDIM==3
           else if (data_id(i) == "vy" .and. VDIM > 1) then
              sink(s)%v(2) = pdata(i)
#endif
#if VDIM==3
           else if (data_id(i) == "vz" .and. VDIM == 3) then
              sink(s)%v(3) = pdata(i)
#endif
           end if
        end do
        sink(s)%accrete = .true.
        sink(s)%static  = .false.
        sink(s)%radius  = KERNRANGE*sink(s)%h
#endif
     end if
     ! -----------------------------------------------------------------------

  end do
     
#if defined(SINKS)
#if defined(USE_MPI)
! If this is MPI and we have sinks, load them now from the root task
! unless this is a root-only read (for decomposition)
  if (stot > 0) then
     sink(1:stot)%domain = 0 ! Assign all sinks to the root temporarily
     if (.NOT. decomp_read) call sink_share
  end if
#endif
#endif
! ============================================================================

20 close(1)

  return
END SUBROUTINE read_data_ascii
