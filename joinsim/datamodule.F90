! Program to modify a simulation
! Andrew McLeod

module datamodule

  implicit none

#if NDIM == 1 && defined(IDEAL_MHD)
  integer, parameter :: BDIM=2
  integer, parameter :: VDIM=2
#elif NDIM == 1 && !defined(IDEAL_MHD)
  integer, parameter :: BDIM=1
  integer, parameter :: VDIM=1 
#elif NDIM == 2 && defined(IDEAL_MHD)
  integer, parameter :: BDIM=3
  integer, parameter :: VDIM=3
#elif NDIM == 2 && !defined(IDEAL_MHD)
  integer, parameter :: BDIM=2
  integer, parameter :: VDIM=2 
#elif NDIM == 3
  integer, parameter :: BDIM=3
  integer, parameter :: VDIM=3
#endif

  logical               :: joinmpi            ! Are we joining MPI simulations?
  character (LEN=200)   :: filename           ! Output filename

  integer, parameter :: DP = selected_real_kind(p=15) ! double precision
  integer, parameter :: SP = selected_real_kind(p=6)  ! single precision

!   integer, parameter   :: NDIM = 3

#ifdef DOUBLE_PRECISION
  integer, parameter :: PR = DP                       ! particle precision
#else
  integer, parameter :: PR = SP                       ! default = single
#endif

  integer              :: i, ptot ! Particle counter, particle number
  real (kind=PR), allocatable :: parray(:,:), r(:,:), v(:,:), h(:)
  real (kind=PR), allocatable :: m(:), temp(:), rho(:)
                               ! Particle positions, smoothing length
                               ! mass, temperature, density, B
  integer, allocatable :: ptype(:)

  !integer              :: nsteps, snapshot
  !real (kind=DP)       :: time
  logical              :: formatted_files
  logical              :: seren_format
  logical              :: new_format

  ! MPI stuff

  integer              :: MPInumthreads              ! Number of threads to use
  integer                     :: numtypes(-1:12)     ! Number of each particle type
  integer                     :: slot(-1:12)         ! Current slot numbers for adding particles
  integer                     :: sink_slot           ! Slot for adding SEREN sinks
  integer              :: t_stot, t_pboundary, t_picm, t_pgas, t_pcdm, t_pdust, t_pion
  logical              :: t_is_porig, t_is_parray, t_is_r, t_is_m, t_is_h
  logical              :: t_is_v, t_is_temp, t_is_rho, t_is_u, t_is_sink_v1
  logical              :: t_is_B

end module datamodule

module seren_data_store
  use datamodule
  implicit none
!   integer, parameter :: DP = selected_real_kind(p=15) ! double precision
!   integer, parameter :: SP = selected_real_kind(p=6)  ! single precision
  integer, parameter :: ILP = selected_int_kind(r=15)  ! Integer long precision
  logical           :: doubleprec              ! Is file double precision?
  character(len=20) :: format_id               ! File format (for verification)
  integer           :: nunit                   ! Number of units
  integer           :: ndata                   ! Number of data entries
  integer           :: stot!, ptot              ! Number of particles/sinks
  integer           :: pboundary, picm, pgas   ! Number of each type of particle
  integer           :: pcdm, pdust, pion       ! Number of each type of particle
  integer           :: pgas_orig, pp_gather    ! Simulation variables
  integer (kind=ILP):: nsteps, snapshot, ntempnext      !   ""
  integer (kind=ILP):: ndiagnext, nsnapnext, nsinknext  !   ""
  integer           :: PRtemp, NDIMtemp, VDIMtemp, BDIMtemp ! Important parameters
  integer           :: dmdt_range              ! DMDT_RANGE
  character(len=20) :: data_id(1:50)           ! Char ids of arrays written
  character(len=20) :: unit_data(1:50)         ! Unit data
  integer           :: typedata(1:5,1:50)      ! type data header array
  integer           :: itemp                   ! Particle type
  real (kind=PR)    :: h_fac                   ! Simulation variable
  real (kind=DP)    :: time, lastsnap, mgas_orig ! Simulation variables

  logical           :: is_porig, is_parray, is_r, is_m, is_h
  logical           :: is_v, is_temp, is_rho, is_u, is_B, is_sink_v1

  integer           :: unit_r, unit_m, unit_h, unit_v
  integer           :: unit_temp, unit_rho, unit_u, unit_B

  integer, allocatable :: porig(:)             ! Original particle IDs
  real(kind=PR), allocatable :: u(:)           ! Internal energy
  real(kind=PR), allocatable :: B(:,:)         ! Magnetic fields
  type unknown_type
      integer                    :: width
      integer                    :: type_id
      character(len=20)          :: data_id
      integer, allocatable       :: idata(:)
      real(kind=PR), allocatable :: data(:)
      real(kind=PR), allocatable :: vector(:,:)
  end type unknown_type

  integer           :: nunknown
  type(unknown_type):: unknown_data(1:1000)

  type sink_node
     integer       :: id             ! Id of sink
     integer       :: ncreate        ! nsteps when sink is created
     real(kind=DP) :: tcreate        ! Physical time when sink is created
     real(kind=DP) :: angmom(1:3)    ! Sink internal angular momentum
     real(kind=DP) :: dmdt           ! Accretion rate
     real(kind=DP) :: luminosity     ! Luminosity from unresolved star
     real(kind=DP) :: temperature    ! Surface temperature of unresolved star
     real(kind=DP) :: star_radius    ! Physical radius of unresolved star
     real(kind=DP),allocatable :: tacc(:) ! Times of previous steps
     real(kind=DP),allocatable :: macc(:) ! Masses accreted in previous steps
     real(kind=PR) :: r(1:3)      ! Position vectors (only use NDIM of this)
     real(kind=PR) :: v(1:3)      ! Velocity vectors (only use VDIM of this)
     real(kind=PR) :: m              ! Particle mass
     real(kind=PR) :: h              ! Sink softening length
     real(kind=PR) :: radius         ! Sink accretion radius
     logical       :: accrete        ! Is the sink accreting particles?
     logical       :: static         ! Is the sink static or not?
  end type sink_node

  type(sink_node), allocatable :: sink_array(:)

end module seren_data_store