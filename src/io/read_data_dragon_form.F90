! READ_DATA_DRAGON_FORM.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Reads in initial conditions file in DRAGON ASCII format.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_data_dragon_form(in_file, decomp_read)
  use particle_module
  use hydro_module
  use type_module
  use sink_module
  use time_module
  use scaling_module
#if defined(USE_MPI)
  use mpi
  use mpi_communication_module
#endif
  implicit none

  character(len=*), intent(in) :: in_file    ! formatted DRAGON snapshot
  logical, intent(in)          :: decomp_read ! MPI decomposition read

  integer, allocatable :: auxtype(:)         ! ..
  integer :: boundaryslot                    ! point in list to boundary ptcl
  integer :: cdmslot                         ! point in list to cdm particle
  integer :: gasslot                         ! point in list to insert ptcl
  integer :: icmslot                         ! point in list to icm particle
  integer :: idata(1:20)                     ! dummy integers
  integer :: p                               ! counter to loop over particles
  integer :: psplit                          ! counter for particle splitting
  integer :: pdead                           ! counter for accreted particles
  integer :: sinkslot                        ! point in list sink particle
  integer :: stot_file                       ! sinks in file
  real(kind=SP) :: rdata(1:50)               ! dummy reals
#if defined(SINKS)
  integer :: s                               ! sink counter
#if defined(USE_MPI)
  integer       :: ierr                      ! MPI error value
#endif
#endif
  real(kind=PR) :: raux

  debug1("Reading in formatted data file : "//trim(in_file)//" [read_data_dragon_form.F90]")

! Initialise counters
  pboundary = 0
  psplit    = 0
  pdead     = 0
  picm      = 0
  pgas      = 0
  pcdm      = 0
  pdust     = 0
  pion      = 0
  stot_file = 0

! Open snapshot file
  open(1, file=in_file, status="old", form="formatted")

! First, read in header information 
  do p=1,20
     read(1,*) idata(p)
  end do
  do p=1,50
     read(1,*) rdata(p)
  end do

! Assign variables for important information
  ptot      = idata(1)
  nsteps    = idata(2)
  snapshot  = idata(4)
  pgas_orig = idata(20)
  time      = real(rdata(1),DP)
  lastsnap  = real(rdata(2),DP)
  mgas_orig = real(rdata(50),DP)
  allocate(auxtype(1:ptot))

! First pass to get numbers of particles of each type
! ----------------------------------------------------------------------------
  do p=1,ptot*2
#if NDIM == 1
     read(1,*) rdata(1)
#elif NDIM == 2
     read(1,*) rdata(1), rdata(2)
#else
     read(1,*) rdata(1), rdata(2), rdata(3)
#endif
  end do
  do p=1,ptot*4
     read(1,*) rdata(1)
  end do
  do p=1,ptot
!     read(1,*) auxtype(p)
     read(1,*) raux
     auxtype(p) = int(raux)
#if defined(SINKS)
     ! Only add sinks if we are using them
     if (auxtype(p) == -1) stot_file = stot_file + 1
#endif
     if (auxtype(p) == 0)  pdead = pdead + 1
     if (auxtype(p) == 1)  pgas = pgas + 1
     if (auxtype(p) == 4)  psplit = psplit + 1
     if (auxtype(p) == 6)  pboundary = pboundary + 1
     if (auxtype(p) == 9)  picm = picm + 1
     if (auxtype(p) == 10) pcdm = pcdm + 1
  end do

#if defined(SINKS)
#if defined(USE_MPI)
  ! Sinks will all be on root task, find correct value of stot (for allocating)
  if (rank == 0) stot = stot_file
  if (.NOT. decomp_read) then
     call MPI_BCAST(stot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  end if
#else
  stot = stot_file
#endif
#endif

! Allocate memory now we know 'correct' value of ptot
! If MPI decomposition read, allocate 'minimal = TRUE'
  call allocate_memory(decomp_read)
  if (.NOT. decomp_read) allocate(minimal_sph(1:ptot))

! Reset ptot to only account for gas, boundary and intercloud particles
  ptot = ptot - (stot_file + pdead)
  ! call types ! Why was this called here but not in read_data_dragon_unform?

#if !defined(USE_MPI)
  write(6,*) "Particles   = ", ptot, "   Sinks = ", stot
  write(6,*) "Gas         = ", pgas
  write(6,*) "Boundary    = ", pboundary
  write(6,*) "Intercloud  = ", picm
  write(6,*) "Dark matter = ", pcdm
  write(6,*) "Splitting   = ", psplit
  write(6,*) "Accreted    = ", pdead
#endif
  if (psplit /= 0) stop "Fatal error: particle splitting not supported"
  if (ptot /= (pgas + picm + pboundary + pcdm)) &
       &stop "Fatal error: particles do not add up"


! Second pass to assign data
! ----------------------------------------------------------------------------
  rewind(1)

  do p=1,20
     read(1,*) idata(p)
  end do
  do p=1,50
     read(1,*) rdata(p)
  end do

! Positions
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  cdmslot = pboundary + picm + pgas + 1
  sinkslot = 1
  do p=1,(ptot + (stot_file + pdead))
#if NDIM == 1
     read(1,*) rdata(1)
#elif NDIM == 2
     read(1,*) rdata(1), rdata(2)
#else
     read(1,*) rdata(1), rdata(2), rdata(3)
#endif
     if (auxtype(p) == 6) then
        minimal_sph(boundaryslot)%r(1:NDIM) = rdata(1:NDIM)
        minimal_sph(boundaryslot)%porig = p
        boundaryslot = boundaryslot + 1
     else if (auxtype(p) == 9) then
        minimal_sph(icmslot)%r(1:NDIM) = rdata(1:NDIM)
        minimal_sph(icmslot)%porig = p
        icmslot = icmslot + 1
     else if (auxtype(p) == 1) then
        minimal_sph(gasslot)%r(1:NDIM) = rdata(1:NDIM)
        minimal_sph(gasslot)%porig = p
        gasslot = gasslot + 1
     else if (auxtype(p) == 10) then
        minimal_sph(cdmslot)%r(1:NDIM) = rdata(1:NDIM)
        minimal_sph(cdmslot)%porig = p
        cdmslot = cdmslot + 1
#if defined(SINKS)
     else if (auxtype(p) == -1) then
        sink(sinkslot)%r(1:NDIM) = rdata(1:NDIM)
        sinkslot = sinkslot + 1
#endif
     end if
  end do

! Velocities
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  cdmslot = pboundary + picm + pgas + 1
  sinkslot = 1
  do p=1,(ptot + (stot_file + pdead))
#if NDIM == 1
     read(1,*) rdata(1)
#elif NDIM == 2
     read(1,*) rdata(1), rdata(2)
#else
     read(1,*) rdata(1), rdata(2), rdata(3)
#endif
     if (auxtype(p) == 6) then
        minimal_sph(boundaryslot)%v(1:VDIM) = rdata(1:VDIM)
        boundaryslot = boundaryslot + 1
     else if (auxtype(p) == 9) then
        minimal_sph(icmslot)%v(1:VDIM) = rdata(1:VDIM)
        icmslot = icmslot + 1
     else if (auxtype(p) == 1) then
        minimal_sph(gasslot)%v(1:VDIM) = rdata(1:VDIM)
        gasslot = gasslot + 1
     else if (auxtype(p) == 10) then
        minimal_sph(cdmslot)%v(1:VDIM) = rdata(1:VDIM)
        cdmslot = cdmslot + 1
#if defined(SINKS)
     else if (auxtype(p) == -1) then
        sink(sinkslot)%v(1:VDIM) = rdata(1:VDIM)
        sinkslot = sinkslot + 1
#endif
     end if
  end do

! Temperatures
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  cdmslot = pboundary + picm + pgas + 1
  sinkslot = 1
  do p=1,(ptot + (stot_file + pdead))
     read(1,*) rdata(1)
     if (auxtype(p) == 6) then
        minimal_sph(boundaryslot)%temp = rdata(1)
        boundaryslot = boundaryslot + 1
     else if (auxtype(p) == 9) then
        minimal_sph(icmslot)%temp = rdata(1)
        icmslot = icmslot + 1
     else if (auxtype(p) == 1) then
        minimal_sph(gasslot)%temp = rdata(1)
        gasslot = gasslot + 1
#if defined(SINKS)
     else if (auxtype(p) == 10) then
        minimal_sph(cdmslot)%temp = rdata(1)
        cdmslot = cdmslot + 1
#endif
     end if
  end do

! Smoothing lengths
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  cdmslot = pboundary + picm + pgas + 1
  sinkslot = 1
  do p=1,(ptot + (stot_file + pdead))
     read(1,*) rdata(1)
     if (auxtype(p) == 6) then
        minimal_sph(boundaryslot)%h = rdata(1)
        boundaryslot = boundaryslot + 1
     else if (auxtype(p) == 9) then
        minimal_sph(icmslot)%h = rdata(1)
        icmslot = icmslot + 1
     else if (auxtype(p) == 1) then
        minimal_sph(gasslot)%h = rdata(1)
        gasslot = gasslot + 1
     else if (auxtype(p) == 10) then
        minimal_sph(cdmslot)%h = rdata(1)
        cdmslot = cdmslot + 1
#if defined(SINKS)
     else if (auxtype(p) == -1) then
        sink(sinkslot)%h = rdata(1)
        sinkslot = sinkslot + 1
#endif
     end if
  end do

! Densities
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  cdmslot = pboundary + picm + pgas + 1
  sinkslot = 1
  do p=1,(ptot + (stot_file + pdead))
     read(1,*) rdata(1)
     if (auxtype(p) == 6) then
        minimal_sph(boundaryslot)%rho = rdata(1)
        boundaryslot = boundaryslot + 1
     else if (auxtype(p) == 9) then
        minimal_sph(icmslot)%rho = rdata(1)
        icmslot = icmslot + 1
     else if (auxtype(p) == 1) then
        minimal_sph(gasslot)%rho = rdata(1)
        gasslot = gasslot + 1
     else if (auxtype(p) == 10) then
        minimal_sph(cdmslot)%rho = rdata(1)
        cdmslot = cdmslot + 1
     end if
  end do

! Masses
! ----------------------------------------------------------------------------
  boundaryslot = 1
  icmslot = pboundary + 1
  gasslot = pboundary + picm + 1
  cdmslot = pboundary + picm + pgas + 1
  sinkslot = 1
  do p=1,(ptot + (stot_file + pdead))
     read(1,*) rdata(1)
     if (auxtype(p) == 6) then
        minimal_sph(boundaryslot)%m = rdata(1)
        boundaryslot = boundaryslot + 1
     else if (auxtype(p) == 9) then
        minimal_sph(icmslot)%m = rdata(1)
        icmslot = icmslot + 1
     else if (auxtype(p) == 1) then
        minimal_sph(gasslot)%m = rdata(1)
        gasslot = gasslot + 1
     else if (auxtype(p) == 10) then
        minimal_sph(cdmslot)%m = rdata(1)
        cdmslot = cdmslot + 1
#if defined(SINKS)
     else if (auxtype(p) == -1) then
        sink(sinkslot)%m = rdata(1)
        sinkslot = sinkslot + 1
#endif
     end if
  end do

! Close file once finished
  close(1)
! ----------------------------------------------------------------------------

  deallocate(auxtype)

! Convert temperatures to internal energies if required
! ----------------------------------------------------------------------------
  do p=1,ptot
     if (gamma > 1.0_PR) then
        minimal_sph(p)%u = Pconst*minimal_sph(p)%temp/(gamma - 1.0_PR)*uscale
     end if
  end do

! Initialise some variables not recorded by Dragon format
! ----------------------------------------------------------------------------
#if defined(MHD)
  do p=1,ptot
     minimal_sph(p)%B = 0.0_DP
  end do
#endif

#if defined(SINKS)
  do s=1,stot_file
     sink(s)%accrete = .true.
     sink(s)%static  = .false.
     sink(s)%radius  = KERNRANGE*sink(s)%h
  end do
  
! Copy data from minimal_sph array to sph array if not an MPI decomp read
! ----------------------------------------------------------------------------
  if (.NOT. decomp_read) then
     do p=1,ptot
        sph(p)%porig = minimal_sph(p)%porig
        sph(p)%r = minimal_sph(p)%r
        sph(p)%v = minimal_sph(p)%v
        sph(p)%m = minimal_sph(p)%m
        sph(p)%rho = minimal_sph(p)%rho
        sph(p)%h = minimal_sph(p)%h
        sph(p)%u = minimal_sph(p)%u
        sph(p)%temp = minimal_sph(p)%temp
     end do
     deallocate(minimal_sph)
  end if

#if defined(USE_MPI)
! If this is MPI and we have sinks, load them now from the root task
! unless this is a root-only read (for decomposition)
  if (stot > 0) then
     sink(1:stot)%domain = 0 ! Assign all sinks to the root temporarily
     if (.NOT. decomp_read) call sink_share
  end if
#endif
#endif

! Verification that particles add up correctly
  if (boundaryslot-1 /= pboundary) stop 'boundary particles do not match'
  if (icmslot-1 /= picm + pboundary) stop 'intercloud particles do not match'
  if (gasslot-1 /= pgas + picm + pboundary) stop 'gas particles do not match'
  if (sinkslot-1 /= stot_file) stop 'sinks do not match'

  return
END SUBROUTINE read_data_dragon_form
