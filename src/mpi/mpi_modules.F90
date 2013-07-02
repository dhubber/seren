! MPI_MODULES.F90
! A. McLeod - 15/03/2011
! Declares MPI module stuff
! ============================================================================

#include "macros.h"

! ============================================================================

MODULE mpi_particle_types
  use definitions
  use particle_module

! ================== Non-particle data ==================
  integer :: MPI_SYNCHRONISE        ! MPI handle for data type synchronise

  type synchronisetype              ! Type for sharing variables every timestep
    sequence
    real (kind=DP) :: step_min      ! Minimum step size
    integer (kind=ILP) :: nmax      ! Maximum level used
    integer        :: totalptot     ! Number of particles across all domains
    logical        :: loadbalance   ! Whether a load balance is needed
  end type

! ================== Particle data ==================
  integer        :: MPI_PARTICLE       ! MPI handle for data type particletype
  integer (kind=ILP) :: mpi_particle_extent ! Size in bytes of MPI_PARTICLE

  type particleaccess  ! Type to create array with allocatable components
    type(sph_particle), allocatable :: data(:)
  end type particleaccess

! ================ Ghost particles for h_gather and sph ================
  integer :: MPI_GHOST_TYPE         ! MPI handle for data type ghosttype
  type ghosttype                    ! Type for exporting ghost particles
!     sequence
    real (kind=PR) :: r(1:NDIM)     ! Position
    real (kind=PR) :: m             ! Mass
    real (kind=PR) :: v(1:VDIM)     ! Velocity
#ifdef DIV_A
    real (kind=PR) :: a(1:NDIM)     ! Accelerations
#endif
#if defined(SINKS) && defined(SELF_GRAVITY)
    real (kind=PR) :: gpot          ! Gravitational potential
#endif
#ifdef DEBUG_GHOST_TX
    integer        :: porig         ! Original particle number
#endif
    integer        :: ptype         ! Type of particle
  end type ghosttype
  
  type ghostaccess ! Type to create array with allocatable components
    type(ghosttype), allocatable :: data(:)
  end type ghostaccess


! ================== Particles for force calculations ==================
! #ifdef HYDRO
!   integer :: MPI_HYDROTYPE          ! MPI handle for data type exporthydrotype
!   type exporthydrotype              ! Type for exporting particles for hydro
! !     sequence
!     real (kind=PR) :: r(1:NDIM)     ! Position
!     real (kind=PR) :: v(1:VDIM)     ! Velocity
!     real (kind=PR) :: h             ! Smoothing length
!     real (kind=PR) :: rho           ! Density
!     real (kind=PR) :: press         ! Pressure
!     real (kind=PR) :: m             ! Mass
!     real (kind=PR) :: div_v         ! Velocity divergence
! #ifdef ARTIFICIAL_VISCOSITY
!     real (kind=PR) :: sound         ! Sound speed
! #endif
! #ifdef GRAD_H_SPH
!     real (kind=PR) :: omega         ! grad-h correction term
! #endif
! #if defined(VISC_TD) || defined(COND_TD)
!     real (kind=PR) :: talpha        ! Time-dependent alpha
! #endif
! #ifdef VISC_BALSARA
!     real (kind=PR) :: balsara       ! Balsara factor
! #endif
! #ifdef VISC_PATTERN_REC
!     real (kind=PR) :: pattrec       ! Pattern recognition factor
! #endif
! #if defined(INTERNAL_ENERGY) && defined(ARTIFICIAL_CONDUCTIVITY)
!     real (kind=PR) :: u             ! Specific internal energy
! #endif
! #ifdef DIFFUSION
!     real (kind=PR) :: temp          ! Particle temperature
! #endif
!     integer        :: p             ! Particle identifier on originating domain
! #ifdef DEBUG_EXPORT
!     integer        :: dfrom         ! Originating domain
!     integer        :: dto           ! Target domain
! #endif
! ! #ifdef DOUBLE_PRECISION
! !     integer        :: pad           ! Padding to ensure alignment
! ! #endif
!   end type exporthydrotype
! 
! #ifdef USE_F2003
!   type exporthydroaccess  ! Type to create array with allocatable components
!     type(exporthydrotype), allocatable :: data(:)
!   end type exporthydroaccess
! #else
!   type exporthydroaccess  ! Type to create array of pointers to exporthydrotype data
!     type(exporthydrotype), pointer :: data(:) => NULL()
!   end type exporthydroaccess
! #endif
! 
!   integer :: MPI_RETURNHYDRO    ! MPI handle for data type returnhydrotype
!   type returnhydrotype          ! Type for returning particles for hydro
! !     sequence
!     real(kind=PR)  :: acc(1:VDIM) ! Partial acceleration value from export
! #ifdef INTERNAL_ENERGY
!     real(kind=PR)  :: du_dt     ! Partial du_dt from export
! #endif
! #ifdef DIFFUSION
!     real (kind=PR) :: radenergygrad(1:NDIM) ! radiation energy gradient
!     real (kind=PR) :: radenergy ! radiation energy
! #endif
!     integer        :: p         ! Particle identifier on originating domain
! #ifdef DEBUG_EXPORT
!     integer        :: dfrom     ! Originating domain
!     integer        :: calc_in   ! Domain in which calculation occurred
! #endif
! ! #ifdef DOUBLE_PRECISION
! !     integer        :: pad       ! Padding to ensure alignment
! ! #endif
!   end type returnhydrotype
! 
!   type returnhydroaccess  ! Type to create array with allocatable components
!     type(returnhydrotype), allocatable :: data(:)
!   end type returnhydroaccess
! #endif

! ----------------------------------------------------------------------------
#if defined(SELF_GRAVITY)
  integer :: MPI_GRAVITYTYPE   ! MPI handle for data type exportgravitytype
  type exportgravitytype       ! Type for exporting particles for gravity
!     sequence
    real (kind=PR) :: r(1:NDIM)
    real (kind=PR) :: h
#ifdef GRAD_H_SPH
    real (kind=PR) :: zo        ! Zeta over omega
#endif
#if defined(BH_TREE) && !defined(GEOMETRIC_MAC)
    real (kind=PR) :: agravmag  ! agravmag for particle
#endif
    integer        :: p         ! Particle identifier on originating domain
#ifdef DEBUG_EXPORT
    integer        :: d         ! Originating domain
#else
! #ifdef DOUBLE_PRECISION
!     integer        :: pad       ! Padding to ensure alignment
! #endif
#endif
  end type exportgravitytype

  type exportgravityaccess  ! Type to create array with allocatable components
    type(exportgravitytype), allocatable :: data(:)
  end type exportgravityaccess

  integer :: MPI_RETURNGRAVITY   ! MPI handle for data type exportgravitytype
  type returngravitytype         ! Type for returning particles for gravity
!     sequence
    real (kind=DP) :: acc(1:VDIM)  ! Gravitational acceleration
    real (kind=DP) :: pot          ! Gravitational potential
    integer        :: p            ! Particle identifier on originating domain
#ifdef DEBUG_EXPORT
    integer        :: d            ! Originating domain
    integer        :: calc_in      ! Domain in which calculation occurred
#endif
!     integer        :: pad          ! Padding to ensure alignment
  end type returngravitytype

  type returngravityaccess  ! Type to create array with allocatable components
    type(returngravitytype), allocatable :: data(:)
  end type returngravityaccess

#ifdef BH_TREE
  integer                 :: MPI_REMOTE_BH_GRAV   ! MPI handle for BH grav tree
!   type BHgrav_sequence              ! BH gravity tree node
! !     sequence
!     real(kind=PR) :: r(1:NDIM)      ! Centre of mass of cell c
!     real(kind=PR) :: m              ! Total mass of cell c
!     real(kind=PR) :: dminsqd        ! Cell-opening distance squared
! #ifndef GEOMETRIC_MAC
!      real(kind=PR) :: mac
! #endif
! #ifdef QUADRUPOLE
!      real(kind=PR) :: q(1:NQUAD)     ! Quadrupole moment terms
! #endif
! #ifdef OCTUPOLE
!      real(kind=PR) :: s(1:NOCT)      ! Octupole moment terms
! #endif
!     integer :: nextcell             ! Next cell in list if not opening
!     integer :: ifopen               ! First child cell if opened
!   end type BHgrav_sequence
#endif
#endif

#if defined(SINKS)
  integer                 :: MPI_SINK_NODE        ! MPI handle for sink_node

!   integer                 :: MPI_SEND_SINK        ! MPI handle for send_sink
!   type sendsink                     ! Type for sending sink particle data
! !     sequence
!     real(kind=PR) :: r(1:NDIM)      ! position vectors
!     real(kind=PR) :: v(1:VDIM)      ! velocity vectors
!     real(kind=PR) :: a(1:VDIM)      ! Acceleration vectors
!     real(kind=PR) :: rold(1:NDIM)   ! Old position vectors
!     real(kind=PR) :: vold(1:VDIM)   ! Old velocity vectors
!     real(kind=PR) :: gpot           ! Gravitational potential
!     real(kind=PR) :: gpe            ! Gravitational potential energy
! #if defined(RUNGE_KUTTA2)
!     real(kind=PR) :: vhalf(1:VDIM)  ! Half-step velocity arrays
! #endif
! #if defined(LEAPFROG_KDK)
!     real(kind=PR) :: aold(1:VDIM)   ! Old accelerations for PC
! #endif
! #if defined(BH_TREE) && !defined(GEOMETRIC_MAC)
!     real(kind=PR) :: agravmag       ! mag. of grav. acceleration
! #endif
!     integer       :: s              ! Sink number
! ! #ifdef DOUBLE_PRECISION
! !     integer       :: pad            ! Padding to ensure alignment
! ! #endif
!   end type sendsink
#endif

#ifdef CHECK_NEIGHBOUR_TIMESTEPS
integer :: MPI_PARTICLETIMESTEP ! MPI handle for data type particle_timestep
  type particletimesteptype         ! Type for returning particles for gravity
!     sequence
    integer(kind=ILP) :: nstep  ! Minimum neighbour stepsize
    real (kind=PR)    :: parray(1:SMOO) ! Particle position, mass, h
! #ifndef DOUBLE_PRECISION
!     integer        :: pad          ! Padding to ensure alignment
! #endif
  end type particletimesteptype

  type particletimestepaccess  ! Type to create array with allocatable components
    type(particletimesteptype), allocatable :: data(:)
  end type particletimestepaccess
#endif

  type xyz_type
#if NDIM==3
    real (kind=PR) :: xyz(1:3)
#elif NDIM==2
    real (kind=PR) :: xyz(1:2)
#else
    real (kind=PR) :: xyz(1:1)
#endif
  end type xyz_type

  type integer_access ! Type to create integer array with allocatable components
    integer, allocatable :: data(:)
  end type integer_access
  type dp_access ! Type to create double precision array with allocatable components
    double precision, allocatable :: data(:)
  end type dp_access
  type tree_xyz ! Type to create real array with x,y,z allocatable components
    type (xyz_type), allocatable :: data(:)
  end type tree_xyz
  type real_access ! Type to create real array with allocatable components
    real (kind=PR), allocatable :: data(:)
  end type real_access
  type reals_access ! Type to create 2D real array with allocatable components
    real (kind=PR), allocatable :: data(:,:)
  end type reals_access


END MODULE mpi_particle_types


! ============================================================================


MODULE mpi_communication_module
  use mpi_particle_types
  use tree_module
  use definitions

  integer, allocatable :: exportlist(:,:)! Identifiers of particles to be exported
                                         ! (for hydro, either gather or
                                         ! gather/scatter) and domain into which
                                         ! the particle overlaps/has crossed into

! Arrays for sending and receiving entire particles when they cross boundaries
! ----------------------------------------------------------------------------
  type(sph_particle), allocatable    :: sendparticles(:)     ! Particles to be
                                         ! entirely sent to other domains
  type(particleaccess), allocatable  :: receiveparticles(:)  ! Particles to be
                                         ! entirely received from other domains

! #ifdef HYDRO
! ! Arrays for sending and receiving, and storing for calculation, particles that
! ! need hydro calculation on remote domains
! ! ----------------------------------------------------------------------------
!   type(exporthydrotype), allocatable   :: sendhydro(:)   ! Particles to be
!                                          ! exported for hydro force calculation
!   type(exporthydroaccess), allocatable :: receivehydro(:) ! Receive buffer for
!                                          ! exported particles
! 
! ! Arrays for the sending and receiving the results of hydro calculation after
! ! exportation to other domains
! ! ----------------------------------------------------------------------------
!   type(returnhydrotype), allocatable   :: returnhydro(:) ! Exported particles with
!                                          ! calculated hydro forces to be
!                                          ! returned to their host domain - also
!                                          ! a send buffer
!   type(returnhydroaccess), allocatable :: receivereturnhydro(:) ! Exported
!                                          ! particles that have had hydro forces
!                                          ! calculated for them on remote domains
! ! ----------------------------------------------------------------------------
! #endif

#if defined(SELF_GRAVITY)
! Arrays for sending and receiving, and storing for calculation, particles that
! need gravity calculation on remote domains
! ----------------------------------------------------------------------------
  type(exportgravitytype), allocatable   :: sendgravity(:)   ! Particles to be
                                         ! exported for gravity force calculation
  type(exportgravityaccess), allocatable :: receivegravity(:) ! Receive buffer for
                                         ! exported particles

! Arrays for the sending and receiving the results of gravity calculation after
! exportation to other domains
! ----------------------------------------------------------------------------
  type(returngravitytype), allocatable   :: returngravity(:) ! Exported particles with
                                         ! calculated gravity forces to be
                                         ! returned to their host domain - also
                                         ! a send buffer
  type(returngravityaccess), allocatable :: receivereturngravity(:) ! Exported
                                         ! particles that have had gravity forces
                                         ! calculated for them on remote domains
#endif

  character(len=6)   :: MPI_ext     ! MPI task number
  integer :: lastrank               ! id of final proc. (numtasks - 1)
  integer :: lastrank2              ! max(1,lastrank)
  integer :: endranklist            ! for list from 0:lastrank-1
  integer :: numtasks               ! Total number of processors
  integer :: totalptot              ! Total number of particles

  integer :: pexport                ! No. of received exported hydro particles
  integer :: pexportgrav            ! No. of received exported grav particles

  integer, allocatable :: grav_fromeach(:,:)  ! Which domain and how many
                                              ! particles from it

  real(kind=PR), allocatable :: bbmin(:,:)        ! Min extent of neib search bb
  real(kind=PR), allocatable :: bbmax(:,:)        ! Max extent of neib search bb
  real(kind=PR), allocatable :: activemin(:,:)    ! Min extent of active particles
  real(kind=PR), allocatable :: activemax(:,:)    ! Max extent of active particles
  real(kind=PR), allocatable :: domain_bbmin(:,:) ! Min extent of domain bb
  real(kind=PR), allocatable :: domain_bbmax(:,:) ! Max extent of domain bb
  real(kind=PR) :: local_min(1:NDIM) ! Min extent of non-overlapped domain
  real(kind=PR) :: local_max(1:NDIM) ! Max extent of non-overlapped

  integer, allocatable       :: ghost_dom(:,:)    ! Ghost particle id & task
  type(integer_access), allocatable :: ghostfrom(:) ! Ghost particle original p
  integer              :: sent_to_tasks      ! How many tasks we sent ghosts to
  logical, allocatable :: ghosts_from_here(:) ! Did we get ghosts from task

#if defined(SELF_GRAVITY)
#ifdef BH_TREE
  type(BHgrav_node), allocatable :: BHlocal_grav(:)
                                        ! Pruned version of local gravity tree
  type(BHgrav_node), allocatable :: BHremote_grav(:,:)
                                        ! Imported BH grav tree

  integer, parameter      :: remotetreedepth=3
                                        ! Depth of imported trees
  integer                 :: cmax_remotegrav
                                        ! Maximum number of cells in pruned trees
  integer                 :: ctot_localgrav
                                        ! Number of cells in local pruned gravity tree
  integer, allocatable    :: treesend(:) ! MPI handle for pruned tree sends
  integer, allocatable    :: treerecv(:) ! MPI handle for pruned tree receives
#endif
#endif

! Storage for periodic ghost particles not needed in this domain
! ----------------------------------------------------------------------------

#if defined(PERIODIC)
  type ghost_sph_reference_type
     integer :: porig                  ! ID of main particle
     real(kind=PR) :: r(1:NDIM)        ! Ghost position
  end type ghost_sph_reference_type
  
  type ghost_sph_reference_access  ! Type to create array with allocatable components
    type(ghost_sph_reference_type), allocatable :: data(:)
  end type ghost_sph_reference_access
  
  type(ghost_sph_reference_access), allocatable :: ghost_sph_reference(:)
                              ! Store for ghost particles to be sent externally
#endif


! Message size limits for transferring particles
! ----------------------------------------------------------------------------
  integer, parameter :: MAX_CHUNKS=1000    ! Maximum number of chunks to send
  integer, parameter :: safety_factor=10   ! Safety factor for chunks
  integer            :: max_message_size   ! Maximum message size (PARTICLES)

! Load balancing timing variables
! ----------------------------------------------------------------------------

  double precision :: walltime, non_start_time, calcstart, calctime ! Time variables
  double precision :: exportstart, exporttime ! Time in foreign particle calc
  double precision :: waitstart, waittime ! Time spent MPI waiting
  double precision :: lastwrite ! Time of last output to screen
  double precision :: start_calctime ! calctime at start of integration step
  double precision :: diag_calctime ! Time of this integration step
  double precision :: since_loadbalance_time ! All time elapsed since loadbalance
  double precision :: last_loadbalance_time ! Last time between loadbalances
  logical          :: do_load_balance  ! Do a load balance this step
  logical          :: demand_load_balance ! Require load balance next step

! Load balancing parameters
! ----------------------------------------------------------------------------
  double precision, parameter :: aim_time_load_balance=60._DP
                ! Aimed for CALC time between loadbalancing
  integer (kind=ILP) :: loadbalance_nsteps
  integer (kind=ILP), parameter :: min_steps_load_balance=8_ILP
                                 ! Best if an even number?
                                 ! Minimum nsteps between loadbalancing
  integer (kind=ILP), parameter :: max_steps_load_balance=64_ILP
                                 ! Best if an even number?
                                 ! Maximum nsteps between loadbalancing

! Particle work tracking parameters
! ----------------------------------------------------------------------------
  integer              :: sum_acctot    ! Sum of acceleration calculations
  real(kind=PR)        :: predict_acctot ! Predicted acctot
  integer              :: sum_second_h  ! Sum of second h_gather calculations
  real(kind=PR)        :: sum_costs     ! Sum of cost calculations

END MODULE mpi_communication_module


! ============================================================================


MODULE mpi_decomposition_module
  use mpi_particle_types
  use definitions

  integer              :: MPItreedepth           ! Depth of binary tree required

  ! MPI Geometry
  integer, allocatable :: MPIgeometry(:)
                           ! Mapping of each task number to tree number
  integer, allocatable :: MPIrevgeometry(:)
                           ! MPI geometry reversed
  type(integer_access), allocatable :: MPItreeoccupation(:)
                           ! Number of occupied leaves through the tree
!   type(integer_access), allocatable :: MPItreeptot(:)
!                            ! Number of particles through the tree
  type(tree_xyz), allocatable       :: MPItreemin(:)
                           ! Minimum x,y,z of boxes within tree
  type(tree_xyz), allocatable       :: MPItreemax(:)
                           ! Minimum x,y,z of boxes within tree
!   type(tree_xyz), allocatable    :: minbox(:)
!                            ! (1:NDIM, task) minimum x,y,z corner of box
!   type(tree_xyz), allocatable    :: maxbox(:)
!                            ! (1:NDIM, task) maximum x,y,z corner of box
  type(dp_access), allocatable :: taskbias_tree(:)
                             ! Tree of summed taskbias values

  ! Particle density grid
  integer, parameter    :: gridsize=128   ! Size of particle density grid
#if NDIM==3
  real (kind=PR), allocatable :: MPIpdensitygrid(:,:,:)
                           ! Particle density grid
#elif NDIM==2
  real (kind=PR), allocatable :: MPIpdensitygrid(:,:)
                           ! Particle density grid
#else
  real (kind=PR), allocatable :: MPIpdensitygrid(:)
                           ! Particle density grid
#endif


END MODULE mpi_decomposition_module



! ============================================================================

! #ifdef LEAPFROG_KDK
!   integer, parameter :: SINK_ACCRETE_DATA=4+NDIM+5*VDIM
! #else
!   integer, parameter :: SINK_ACCRETE_DATA=4+NDIM+4*VDIM
! #endif
!                                      ! Width of data field for accretion data
!                                      ! to be merged at root task
! 
!   type(sink_node), allocatable :: sink(:)  ! Main sink data structure
!   type(sink_node), allocatable :: newsinks(:) ! Data structure for new sinks
!   integer                      :: new_sink_number ! How many new sinks
!   integer                      :: next_sink_id ! Next unique id for a sink
! 
!   real(kind=DP) :: store_data(1:SINK_ACCRETE_DATA,1:SMAX)
!                                       ! Accretion data to be merged at root
!   real(kind=DP), allocatable :: task_accrete_data(:,:,:)
!                                       ! Receive buffer for accretion data
! 
!   integer                 :: accretesend    ! MPI handle for accretion send
!   integer, allocatable    :: accreterecv(:) ! MPI handle for accretion receives
! 
!   integer        :: createdead                 ! create_sink dead ghost number
!   integer, allocatable :: createdeadlist(:) ! create_sink dead ghost particles

! ============================================================================
