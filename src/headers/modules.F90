! MODULES.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Contains all global variable modules
! - constant_module
! - Eos_module; Tprof_module; Eos_functions
! - ewald_module
! - filename_module
! - HIIregion_module
! - hydro_module
! - kernel_module
! - Nbody_module
! - neighbour_modules
! - particle_module
! - periodic_module
! - scaling_module
! - sink_module
! - sink_correction_module
! - time_module
! - timing_module
! - turbulence_module
! - tree_module
! - type_module
! ============================================================================

#include "macros.h"

! ============================================================================
MODULE constant_module
  use definitions

  ! Physical constants in SI units (unless stated otherwise)
  real(kind=DP),parameter :: r_pc    = 3.08568E16_DP     ! Parsec
  real(kind=DP),parameter :: r_au    = 1.49597870E11_DP  ! Astronomical unit
  real(kind=DP),parameter :: r_sun   = 6.96E8_DP         ! Solar radius
  real(kind=DP),parameter :: r_earth = 6.371E6_DP        ! Earth radius
  real(kind=DP),parameter :: pc_au   = 206265.0_DP       ! Parsec (in AU)
  real(kind=DP),parameter :: km_cm   = 1.E5_DP           ! km (in cm)
  real(kind=DP),parameter :: m_sun   = 1.98892E30_DP     ! Solar mass
  real(kind=DP),parameter :: m_jup   = 1.8986E27_DP      ! Jupiter mass
  real(kind=DP),parameter :: m_earth = 5.9736E24_DP      ! Earth mass
  real(kind=DP),parameter :: myr = 3.1556952E13_DP       ! Megayear
  real(kind=DP),parameter :: yr  = 3.1556952E7_DP        ! Year
  real(kind=DP),parameter :: day = 8.64E4_DP             ! Day
  real(kind=DP),parameter :: amu = 1.660538782E-27_DP    ! Atomic mass unit
  real(kind=DP),parameter :: m_hydrogen = 1.66054E-27_DP ! Hydrogen mass
  real(kind=DP),parameter :: G_const = 6.67428E-11_DP    ! Grav. constant
  real(kind=DP),parameter :: k_boltzmann = 1.3807E-23_DP ! Boltzmann constant 
  real(kind=DP),parameter :: stefboltz = 5.6704E-8_DP    ! Stefan-boltzmann
  real(kind=DP),parameter :: e_charge = 1.6021765E-19_DP ! Electron charge
  real(kind=DP),parameter :: mu_0 = 1.25663706144E-6_DP  ! Permeability of
                                                         ! free space
  real(kind=DP),parameter :: kappa_const = 2.09E-4_DP    ! Opacity kappa (cgs)
  real(kind=DP),parameter :: L_sun = 3.839E26_DP         ! Solar luminosity

  ! Numerical constants
  real(kind=PR),parameter :: pi = 3.1415926536_PR
  real(kind=PR),parameter :: twopi = 6.283185307_PR
  real(kind=PR),parameter :: invpi = 0.318309886_PR
  real(kind=PR),parameter :: invlogetwo = 1.442695041_PR
  real(kind=PR),parameter :: invlog10two = 3.321928095_PR
  real(kind=PR),parameter :: invsqrttwo = 0.707106781_PR
  real(kind=PR),parameter :: onethird = 0.333333333_PR
  real(kind=PR),parameter :: onesixth = 0.166666666_PR
  real(kind=PR),parameter :: twothirds = 0.666666666_PR
  real(kind=PR),parameter :: big_number = 9.9e20_PR
  real(kind=PR),parameter :: small_number = 1.0e-20_PR

  real(kind=DP),parameter :: pi_dp = 3.1415926536_DP
  real(kind=DP),parameter :: twopi_dp = 6.283185307_DP
  real(kind=DP),parameter :: invpi_dp = 0.318309886_DP
  real(kind=DP),parameter :: invlogetwo_dp = 1.442695041_DP
  real(kind=DP),parameter :: invlog10two_dp = 3.321928095_DP
  real(kind=DP),parameter :: invsqrttwo_dp = 0.707106781_DP
  real(kind=DP),parameter :: onethird_dp = 0.333333333_DP
  real(kind=DP),parameter :: onesixth_dp = 0.166666666_DP
  real(kind=DP),parameter :: twothirds_dp = 0.666666666_DP
  real(kind=DP),parameter :: big_number_dp = 9.9e20_DP
  real(kind=DP),parameter :: small_number_dp = 1.0e-20_DP  

END MODULE constant_module


! ============================================================================
MODULE diagnostics_module
  use definitions

  real(kind=DP) :: etot                       ! Total energy
  real(kind=DP) :: etot0                      ! Total initial energy
  real(kind=DP) :: mtot                       ! Total mass
  real(kind=DP) :: mtot0                      ! Total initial amount of mass
  real(kind=DP) :: rcom(1:NDIM)               ! Position of centre of mass
  real(kind=DP) :: rcom0(1:NDIM)              ! Initial position of COM
  real(kind=DP) :: vcom(1:VDIM)               ! Velocity of centre of mass
  real(kind=DP) :: vcom0(1:VDIM)              ! Initial velocity of COM

END MODULE diagnostics_module


! ============================================================================
MODULE Eos_module
  use definitions

  integer :: dim_temp                           ! Dimension of temp table
  integer :: dim_dens                           ! Dimension of dens table

  real(kind=PR) :: bdens                        ! log10 density factor
  real(kind=PR) :: btemp                        ! log10 temp. factor
  real(kind=PR) :: densmax                      ! Max. density in eos table
  real(kind=PR) :: densmin                      ! Min. density in eos table
  real(kind=PR) :: fcolumn                      ! Column density
                                                ! Polytrope correction
  real(kind=PR) :: rad_const                    ! Scaling constant
  real(kind=PR) :: tempmin                      ! Max. temp. in eos table
  real(kind=PR) :: tempmax                      ! Min. temp. in eos table

  real(kind=PR), allocatable :: eos_dens(:)     ! Density of EOS table
  real(kind=PR), allocatable :: eos_temp(:)     ! Temperature of EOS table
  real(kind=PR), allocatable :: eos_energy(:,:) ! Energy from EOS table
  real(kind=PR), allocatable :: eos_mu(:,:)     ! mu from EOS table
  real(kind=PR), allocatable :: kappa(:,:)      ! Mean opacity EOS table
  real(kind=PR), allocatable :: kappap(:,:)     ! Planck opacity EOS table
  real(kind=PR), allocatable :: kappar(:,:)     ! Rosseland opacity EOS table
#if defined(DEBUG_RAD)
  real(kind=PR), allocatable :: rad_info(:,:)   ! Debug info
#endif

END MODULE Eos_module


! ============================================================================
MODULE Tprof_module
  use definitions

  ! Temperature profile away from a luminosity source
  real(kind=PR) :: ptemp0        ! Temperature at r=1AU from the star
  real(kind=PR) :: ptemp_r0      ! Temperature softening radius (<<1AU)
  real(kind=PR) :: ptemp_q       ! Temperature power law index
  real(kind=PR) :: temp_inf      ! Temperature at infinity

END MODULE Tprof_module


! ============================================================================
MODULE ewald_module
  use definitions
  implicit none

  ! Ewald correction grid size
#if NDIM==2
  integer, parameter :: ewsize(1:NDIM)=(/32,32/)
#elif NDIM==3
  integer, parameter :: ewsize(1:NDIM)=(/32,32,32/)
#endif
  real(kind=PR) :: eforce(1:NDIM)              ! Ewald correction acceleration
  real(kind=PR) :: ewsizeil(1:NDIM)            ! (ewsize - 1)/L
  real(kind=PR) :: L(1:NDIM)                   ! Periodic box size
  real(kind=PR), allocatable :: fcorr(:,:,:,:) ! Ewald correction grids

END MODULE ewald_module


! ============================================================================
MODULE filename_module
  use definitions

  logical :: incomplete_params        ! Accept an incomplete params file?
  logical :: inifile                  ! Has .ini file been written yet?
  logical :: restart                  ! Is this simulation a restart?
  integer :: nparams                  ! No. of params available
  integer :: ntemp                    ! Next temp file to be written
  integer :: ptrack                   ! (Original) id of particle to track
  character(len=256) :: eos_opa_file  ! Filename of eos/opacity table
  character(len=256) :: error_file    ! Error report file
  character(len=256) :: in_file       ! Initial conditions file
  character(len=256) :: run_id        ! Simulation run identifier (name)
  character(len=256) :: run_dir       ! Name of run directory
  character(len=256) :: fileform_ext  ! File format extension for filenames
  character(len=256) :: in_file_form  ! Format of initial conditions file
  character(len=256) :: out_init      ! Initial snapshot output
  character(len=256) :: out_final     ! Final snapshot output
  character(len=256) :: out_temp      ! Temporary snapshot output
  character(len=256) :: out_file_form ! Format of snapshot files
  character(len=256) :: out_temp1     ! Temporary snapshot 1
  character(len=256) :: out_temp2     ! Temporary snapshot 2
  character(len=256) :: restart_log   ! Restart file that records last snapshot
  character(len=256) :: param_file    ! Parameter file
  real(kind=PR) :: rzero(1:NDIM)      ! Position of origin for debug output

  type seren_param                             ! Seren parameter structure
     character(len=256)          :: var_name   ! Variable name
     character(len=1)            :: var_type   ! Variable type
     logical                     :: done       ! Has variable been read?
     character(len=256), pointer :: var_c      ! 256-character string
     character(len=20), pointer  :: var_u      ! 20-character unit
     logical, pointer            :: var_l      ! Logical
     integer, pointer            :: var_i      ! Integer
     integer(kind=ILP), pointer  :: var_j      ! Long integer
     real(kind=PR), pointer      :: var_r      ! PR real
     real(kind=DP), pointer      :: var_d      ! DP real
  end type seren_param
  type(seren_param) :: params(1:256)           ! Main parameter array

END MODULE filename_module


! ============================================================================
MODULE HP_module
  use definitions
  use healpix_types

  logical :: HP_ionize                    ! Ionizing radiation?
  integer :: imax                         ! Maximum no. of rays
  integer :: lmax_hp                      ! Max number of HEALPix levels
  integer :: ltot_hp                      ! Total number of HEALPix levels
  integer :: HPtot                        ! No. of UV sources
  integer :: nnewwind                     ! No. of new wind particles
  integer(kind=i4b) :: x2pix(1:128)       ! Internal HEALPix array
  integer(kind=i4b) :: y2pix(1:128)       ! ""
  integer(kind=i4b) :: pix2x(0:1023)      ! ""
  integer(kind=i4b) :: pix2y(0:1023)      ! ""
  real(kind=DP) :: a_star                 ! Recombination coefficient
  real(kind=PR) :: f1                     ! Integration step factor
  real(kind=PR) :: f2                     ! Opening criterion factor
  real(kind=PR) :: f3                     ! Temperature smoothing factor
  real(kind=PR) :: f4                     ! Density interpolation factor
  real(kind=PR) :: HPmaxres               ! Maximum HEALPix resolution
  real(kind=DP) :: intmax                 ! Value of integral at IF
  real(kind=DP) :: N_LyC                  ! Flux of UV photons
  real(kind=PR) :: rstatic(1:3)           ! Position of single static source
  real(kind=PR) :: Tion                   ! Temp. of the ionized gas
  real(kind=PR) :: Tneut                  ! Temp. of the neutral gas
  real(kind=PR) :: Xfrac                  ! Fraction by mass of hydrogen
  real(kind=PR) :: Yfrac                  ! Fraction by mass of helium
  real(kind=PR) :: M_loss                 ! Stellar wind mass loss rate (dmdt)
  real(kind=PR) :: v_wind                 ! Stellar wind velocity
#if defined(DEBUG_HP_WALK_ALL_RAYS)
  integer, allocatable :: whichHPlevel(:) ! HP level of particle p
#endif
#if defined(STELLAR_WIND)
  integer, allocatable :: windsurface(:)  ! Wind surface info
#endif

#if defined(HEALPIX)
  type HPlevel_node                           ! HEALPix level
     integer :: ifirst                        ! First ray on current level
     integer :: ilast                         ! Last ray on current level
  end type HPlevel_node
  type(HPlevel_node) :: HPlevel(0:HP_LEVELS)  ! HEALPix level info

  type HPray_node                             ! HEALPix ray 
     logical :: done                          ! Is this ray finished?
     integer(kind=I4B) :: ipix                ! HEALPix pixel id
     integer :: first                         ! 1st particle in linked list
     integer :: last                          ! Last particle in lnked list
     integer :: rayend                        ! Final particle in ray
     real(kind=PR) :: hep                     ! Smoothing length at e.p.
     real(kind=PR) :: rhoep                   ! Density at e.p.
     real(kind=PR) :: rep(1:NDIM)             ! Position of e.p.
#if defined(IONIZING_UV_RADIATION)
     logical :: UVdone                        ! Finished UV radiation
     real(kind=PR) :: integral                ! Current ionization integral
#endif
#if defined(STELLAR_WIND)
     logical :: winddone                      ! Finished stellar winds
#endif
  end type HPray_node
  type(HPray_node), allocatable :: HPray(:)   ! HEALPix rays

  type HPsource_node                          ! HEALPix source
     integer :: sinkid                        ! Sink id of source
     integer :: Nlist                         ! Size of list order
     real(kind=DP) :: arot(1:3,1:3)           ! Rotation matrix
     real(kind=DP) :: Aangle                  ! Rotation angle
     real(kind=DP) :: Bangle                  !    "       "
     real(kind=DP) :: Cangle                  !    "       "
     real(kind=PR) :: r(1:NDIM)               ! Position of HP source
     integer, allocatable :: distorder(:)     ! Distance order array
#if defined(IONIZING_UV_RADIATION)
     real(kind=PR) :: intmax                  ! Stromgren integral value
     real(kind=DP) :: N_LyC                   ! Ionizing photons per sec.
     real(kind=PR) :: HIIangle                ! Solid angle of HII region
     real(kind=PR) :: HIIvolume               ! Volume of HII region
#endif
#if defined(STELLAR_WIND)
     real(kind=PR) :: M_loss                  ! Stellar mass loss rate
     real(kind=PR) :: v_wind                  ! Radial wind speed
     real(kind=PR) :: windangle               ! Solid angle of HII region
     real(kind=PR) :: windvolume              ! Volume of HII region
     real(kind=DP) :: tlastwind               ! Time of last particle injection
#endif
#if defined(STELLAR_LUMINOSITY)
     real(kind=PR) :: L_star                  ! Stellar luminosity
#endif
#if defined(PARTICLE_INJECTION_WINDS)
     integer :: numb_inj                      ! Marker for no. of injections  
#endif
  end type HPsource_node
  type(HPsource_node) :: HPsource(1:SMAX)     ! HP source array

#endif


END MODULE HP_module


! ============================================================================
MODULE hydro_module
  use definitions

  character(len=256) :: cooling_law           ! Cooling law
  character(len=256) :: energy_integration    ! Energy equation integration 
  real(kind=PR) :: alpha                      ! alpha viscosity parameter
  real(kind=PR) :: alpha_min                  ! alpha_min
  real(kind=PR) :: beta                       ! beta viscosity parameter
  real(kind=PR) :: gamma                      ! Ratio of specific heats
  real(kind=PR) :: gammaone                   ! (gamma - 1)
  real(kind=PR) :: isotemp                    ! Isothermal temperature
  real(kind=PR) :: Kpoly                      ! Polytropic constant
  real(kind=PR) :: mu_bar                     ! Mean mol. weight
  real(kind=PR) :: mu_ion                     ! mu of ionised gas
  real(kind=PR) :: newsound_const             ! Sound speed scale
  real(kind=PR) :: Pconst                     ! Ideal gas law scaling const.
  real(kind=PR) :: Pconst2                    !  "" for RAD_WS EoS.
  real(kind=PR) :: Pext                       ! External pressure
  real(kind=PR) :: rhobary                    ! Barotropic density
  real(kind=PR) :: sound_const                ! Sound speed scale
#if defined(ARTIFICIAL_CONDUCTIVITY)
  real(kind=PR) :: alpha_cond                 ! Conductivity alpha
#endif

END MODULE hydro_module


! ============================================================================
MODULE kernel_module
  use definitions

#if defined(KERNEL_TABLES)
  real(kind=PR), allocatable :: w0table(:)     ! Kernel function
  real(kind=PR), allocatable :: w1table(:)     ! Kernel derivative (TC)
  real(kind=PR), allocatable :: w2table(:)     ! Kernel derivaitive (orig)
  real(kind=PR), allocatable :: wgravtable(:)  ! Grav. force kernel
  real(kind=PR), allocatable :: wpottable(:)   ! Grav. potential kernel
  real(kind=PR), allocatable :: womegatable(:) ! 'grad-h' correction kernel
  real(kind=PR), allocatable :: wzetatable(:)  ! 'grad-h' gravity kernel
#endif

END MODULE kernel_module


! ============================================================================
MODULE Nbody_module
  use definitions

  integer :: nbin                    ! No. of mutiple systems in simulation
  integer :: npec                    ! No. of P(EC)^n iterations
  real(kind=DP) :: gammapertmax      ! Max. allowed perturbation for binary
  real(kind=DP) :: nbody_endtime     ! Final runtime of N-body simulation
  real(kind=DP) :: nbody_frac        ! Frac. of gas removed before N-body sim
  real(kind=DP) :: nbody_lastsnap    ! Time of last nbody snapshot
  real(kind=DP) :: nbody_timemult    ! N-body timestep multiplier

  type star_node
     logical :: accdo                ! Acceleration step?
     integer(kind=ILP) :: ncreate    ! nsteps when sink is created
     integer(kind=ILP) :: nlast      ! n of beginning of current timestep
     integer(kind=ILP) :: nlevel     ! Timestep level of star
     real(kind=DP) :: tcreate        ! Physical time when sink is created
     real(kind=DP) :: r(1:NDIM)      ! Position vectors
     real(kind=DP) :: v(1:VDIM)      ! Velocity vectors
     real(kind=DP) :: m              ! Particle mass
     real(kind=DP) :: h              ! Star softening length
     real(kind=DP) :: radius         ! Sink accretion radius
     real(kind=DP) :: rold(1:NDIM)   ! Old position vectors
     real(kind=DP) :: vold(1:VDIM)   ! Old velocity vectors
     real(kind=DP) :: angmom(1:3)    ! Star internal angular momentum
     real(kind=DP) :: gpe            ! Gravitational potential energy
     real(kind=DP) :: gpot           ! Gravitational potential
     real(kind=DP) :: agravmag       ! Magnitude of grav. acceleration
     real(kind=DP) :: dmdt           ! Accretion rate
     real(kind=DP) :: laststep       ! Last timestep
     real(kind=DP) :: a(1:VDIM)      ! Acceleration vectors
     real(kind=DP) :: a0(1:VDIM)     ! Acceleration at start of timestep
     real(kind=DP) :: adot(1:VDIM)   ! Jerk vector
     real(kind=DP) :: adot0(1:VDIM)  ! Jerk at start of timestep
     real(kind=DP) :: a2dot(1:VDIM)  ! 2nd deriv accel
     real(kind=DP) :: a2dot0(1:VDIM) ! 2nd deriv accel at start of timestep
     real(kind=DP) :: a3dot(1:VDIM)  ! 3rd deriv accel
     real(kind=DP) :: luminosity     ! Luminosity from unresolved star
     real(kind=DP) :: temperature    ! Surface temperature of unresolved star
     real(kind=DP) :: star_radius    ! Physical radius of unresolved star
     real(kind=DP) :: macc(1:DMDT_RANGE) ! Masses accreted in previous steps
     real(kind=DP) :: tacc(1:DMDT_RANGE) ! Times of previous steps
  end type star_node
  type(star_node), allocatable   :: star(:)   ! Star array (N-body)

  type binary_node                   ! Binary star data structure
     integer :: id                   ! System id
     integer :: s1                   ! id of sink/binary 1
     integer :: s2                   ! id of sink/binary 2
     real(kind=DP) :: r(1:NDIM)      ! Position of COM of binary
     real(kind=DP) :: v(1:VDIM)      ! Velocity of COM of binary
     real(kind=DP) :: m              ! Total mass
     real(kind=DP) :: angmom(1:3)    ! Magnitude of angular momentum
     real(kind=DP) :: binen          ! Binding energy
     real(kind=DP) :: ecc            ! Eccentricity
     real(kind=DP) :: period         ! Period
     real(kind=DP) :: q              ! Mass-ratio (=m2/m1)
     real(kind=DP) :: sma            ! Semi-major axis
     real(kind=DP) :: drmag          ! Instantaneous dist. between components
  end type binary_node
  type(binary_node), allocatable :: binary(:)  ! Main binary array

END MODULE Nbody_module


! ============================================================================
MODULE neighbour_module
  use definitions

  integer :: listtot                      ! Max. no. of lists available
  integer :: pp_gather                    ! neighbours wanted
  integer :: pp_limit                     ! max. neighbours limit
  real(kind=PR) :: hmin                   ! Minimum allowed smoothing length
  real(kind=PR) :: h_fac                  ! grad-h density-h factor

  integer, allocatable :: pptot(:)        ! number of neighbours
  integer, allocatable :: pplist(:,:)     ! list of neighbours

END MODULE neighbour_module


! ============================================================================
MODULE particle_module
  use definitions

  logical :: com_frame                        ! Flag to change to COM frame
  logical :: rho_remove                       ! Remove particle below min. rho?
  logical :: energy_remove                    ! Remove escaping particles?
  logical :: rad_remove                       ! Remove distant particles?
  integer :: pghostmax                        ! Max. no. of ghost particles
  integer :: pmax                             ! Length of particle arrays
  integer :: ptot                             ! Total number of particles
  integer :: rseed                            ! Random number seed
  real(kind=PR) :: rholost                    ! Particle removal density
  real(kind=PR) :: rad_lost                   ! Particle removal radius
  real(kind=PR) :: rmax(1:NDIM)               ! Max. extent of particles
  real(kind=PR) :: rmin(1:NDIM)               ! Min. extent of particles
#if defined(H_RHO)
  real(kind=PR) :: rextent                    ! Max extent in any 1 dimension
#endif
#if defined(REMOVE_OUTLIERS)
  real(kind=DP) :: mlost                      ! Mass removed from simulation
  real(kind=DP) :: rlost(1:NDIM)              ! COM of removed mass
  real(kind=DP) :: momlost(1:NDIM)            ! Momentum of removed mass
  real(kind=DP) :: angmomlost(1:3)            ! Ang mom. of removed mass
#endif

  real(kind=PR), allocatable :: parray(:,:)   ! r, m and h particle data
                                              ! grouped into one array

  ! Main SPH particle data structure
  ! --------------------------------------------------------------------------
  type sph_particle
     logical :: accdo                  ! Update particle properties
     integer :: porig                  ! Original id
     integer :: ptype                  ! Particle type
     integer(kind=ILP) :: nlast        ! Integer time of last update
     integer(kind=ILP) :: nlevel       ! Integer time step level
#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
     integer(kind=ILP) :: nminneib     ! Minimum neighbour stepsize
#endif
     real(kind=PR) :: r(1:NDIM)        ! Position
     real(kind=PR) :: v(1:VDIM)        ! Velocity vector
     real(kind=PR) :: a(1:VDIM)        ! Acceleration vector
     real(kind=PR) :: m                ! Mass
     real(kind=PR) :: div_v            ! Velocity divergence
     real(kind=PR) :: rho              ! Density
#if !defined(LOW_MEM)
     real(kind=PR) :: invrho           ! 1 / rho
#endif
     real(kind=PR) :: h                ! Smoothing length
#if !defined(LOW_MEM)
     real(kind=PR) :: invh             ! 1 / h
     real(kind=PR) :: hfactor          ! (1 / h)^(NDIM + 1)
#endif
     real(kind=PR) :: u                ! Specific internal energy
     real(kind=PR) :: Aent             ! Entropic function
     real(kind=PR) :: press            ! Pressure
     real(kind=PR) :: sound            ! Sound speed
     real(kind=PR) :: temp             ! Temperature
#if defined(TWO_FLUIDS)
     real(kind=PR) :: h2               ! 2nd fluid smoothing length
     real(kind=PR) :: invh2            ! 1 / h2
#endif
#if defined(SM2012_SPH)
     real(kind=PR) :: q                ! Internal energy density
#endif
#if defined(GRAD_H_SPH)
     real(kind=PR) :: omega            ! 'grad-h' correction
#endif
#if defined(GRAD_H_SPH) && defined(SELF_GRAVITY)
     real(kind=PR) :: zo               ! Gravity 'grad-h' correction
#endif
#if defined(GRAVITY)
     real(kind=PR) :: gpot             ! Gravitational potential
#endif
#if defined(SELF_GRAVITY)
     real(kind=PR) :: agravmag         ! Mag. of grav. acceleration
#endif
     real(kind=DP) :: laststep         ! Previous step size
     real(kind=PR) :: r_old(1:NDIM)    ! Old particle positions
     real(kind=PR) :: v_old(1:VDIM)    ! Old particle velocities
#if defined(RUNGE_KUTTA2)
     real(kind=PR) :: v_half(1:VDIM)   ! Half-step velocities
#endif
#if defined(LEAPFROG_KDK)
     real(kind=PR) :: a_old(1:VDIM)    ! Old particle accelerations
#endif
#if defined(SMOOTHED_VELOCITY)
     real(kind=PR) :: v_smooth(1:VDIM) ! Smoothed velocity
#endif
#if defined(RAD_WS) && defined(SELF_GRAVITY)
     integer :: idens                  ! Indices for dens from table
     integer :: itemp                  ! Indices for temp from table
     real(kind=PR) :: sphgpot          ! Grav. pot just due to SPH ptcls
     real(kind=PR) :: column2          ! Squared column density
     real(kind=PR) :: ueq              ! Equilibrium u
     real(kind=PR) :: dt_therm         ! Thermal cooling timescale
#endif
#if defined(SINKS) && defined(SELF_GRAVITY)
     logical :: ispotmin               ! Does this particle have the min
                                       ! potential of its neibs?
#endif
#if defined(STELLAR_WIND)
     logical :: windflag               ! Is particle ejected from wind?
     real(kind=PR) :: a_wind(1:NDIM)   ! Momentum driven wind acceleration
#endif
#if defined(HYDRO)
#if defined(ENERGY_EQN)
     real(kind=PR) :: u_old            ! Old internal energy
     real(kind=PR) :: dudt             ! Compressional heating rate
#if defined(LEAPFROG_KDK)
     real(kind=PR) :: dudt_old         ! Old dudt
#endif
#if defined(DIFFUSION)
     real(kind=PR) :: du_dt_diff       ! Diffused energy rate
     real(kind=PR) :: k_cond           ! Diffusion conductivity
#endif
#if defined(COOLING_HEATING)
     real(kind=PR) :: dudt_cool        ! Cooling rate
#endif
#endif
#if defined(ENTROPIC_FUNCTION)
     real(kind=PR) :: Aold             ! Old entropic function
     real(kind=PR) :: dAdt             ! Rate of change of entropic fn.
#if defined(LEAPFROG_KDK)
     real(kind=PR) :: Ahalf            ! Half-step entropic function
#endif
#endif
#if defined(VISC_TD)
     real(kind=PR) :: talpha           ! Time-dependent alpha parameter
     real(kind=PR) :: talpha_old       ! alpha at start of timestep
     real(kind=PR) :: dalpha_dt        ! Rate of change of alpha
#endif
#if defined(VISC_BALSARA)
     real(kind=PR) :: balsara          ! Balsara factor
#endif
#endif
#if defined(IONIZING_UV_RADIATION)
     integer :: newtemp                ! New temperature calculated
     real(kind=PR) :: tempmin          ! Minimum temperature
     real(kind=PR) :: tempaux          ! Aux. temp array
#endif
#if defined(DIV_A)
     real(kind=PR) :: div_a            ! Acceleration divergence
#endif
#if defined(DEBUG_FORCES)
     real(kind=PR) :: a_hydro(1:VDIM)  ! Hydro acceleration
     real(kind=PR) :: a_grav(1:VDIM)   ! Gravitational acceleration
     real(kind=PR) :: a_mag(1:VDIM)    ! Magnetic acceleration
     real(kind=PR) :: a_visc(1:VDIM)   ! Viscous acceleration
#endif
#if defined(DEBUG_DUDTRAD)
     real(kind=PR) :: dudt_rad         ! Radiative cooling rate
#endif
#if defined(USE_MPI) && !defined(LOW_MEM)
#if defined(HYDRO)
     real(kind=PR) :: hydro_calcs      ! No. of hydro calcs for particle
#endif
#if defined(GRAVITY)
     real(kind=PR) :: gravity_calcs    ! No. of grav. calcs for particle
#endif
#endif
#if defined(CHEMCOOL)
     real(kind=PR) :: abundances(1:NCHEM)       ! Chemical abundances
#endif
  end type sph_particle
  
  type minimal_sph_particle
     integer :: porig                  ! Original id
     real(kind=PR) :: r(1:NDIM)        ! Position
     real(kind=PR) :: v(1:VDIM)        ! Velocity vector
     real(kind=PR) :: m                ! Mass
     real(kind=PR) :: rho              ! Density
     real(kind=PR) :: h                ! Smoothing length
     real(kind=PR) :: u                ! Specific internal energy
     real(kind=PR) :: temp             ! Temperature
#if defined(ENTROPIC_FUNCTION)
     real(kind=PR) :: Aent             ! Entropic function
#endif
#if defined(CHEMCOOL)
     real(kind=PR) :: abundances(1:NCHEM)       ! Chemical abundances
#endif
  end type minimal_sph_particle
  
  type(sph_particle), allocatable :: sph(:)  ! Main sph particle array
  type(minimal_sph_particle), allocatable :: minimal_sph(:) ! Minimal SPH array

END MODULE particle_module


! ============================================================================
MODULE periodic_module
  use definitions

  integer :: psphere                     ! Id of sphere particle
  real(kind=PR) :: periodic_min(1:3)     ! Minimum extent of periodic box
  real(kind=PR) :: periodic_max(1:3)     ! Maximum extent of periodic box
  real(kind=PR) :: periodic_size(1:3)    ! Size of periodic box
  real(kind=PR) :: periodic_half(1:3)    ! Half-size of periodic box
  real(kind=PR) :: rspheremax            ! Radius of spherical mirror
#if defined(USE_MPI)
  real(kind=PR) :: periodic_mpi_min(1:3) ! Point at which to wrap particles
  real(kind=PR) :: periodic_mpi_max(1:3) ! Point at which to wrap particles
#endif
  real(kind=PR) :: ghost_search_min(1:3) ! Search box min for ghost particles
  real(kind=PR) :: ghost_search_max(1:3) ! Search box max for ghost particles
  logical       :: leftwall(1:3)         ! Is left wall a boundary?
  logical       :: rightwall(1:3)        ! Is right wall a boundary?
  real(kind=PR) :: periodic_half_minval  ! Smallest 1/2 of box

END MODULE periodic_module


! ============================================================================
MODULE scaling_module
  use definitions

  logical :: dimensionless           ! Flag if using dimensionless units

! Scaling unit strings
  character(len=20) :: runit         ! length unit
  character(len=20) :: munit         ! mass unit
  character(len=20) :: tunit         ! time unit
  character(len=20) :: vunit         ! velocity unit
  character(len=20) :: aunit         ! acceleration unit
  character(len=20) :: rhounit       ! density unit
  character(len=20) :: sigmaunit     ! column density unit
  character(len=20) :: Punit         ! pressure unit
  character(len=20) :: funit         ! force unit
  character(len=20) :: Eunit         ! energy unit
  character(len=20) :: momunit       ! momentum unit
  character(len=20) :: angmomunit    ! angular momentum unit
  character(len=20) :: angvelunit    ! angular vel unit
  character(len=20) :: dmdtunit      ! accretion rate unit
  character(len=20) :: Lunit         ! luminosity unit
  character(len=20) :: kappaunit     ! opacity unit
  character(len=20) :: Bunit         ! magnetic field unit
  character(len=20) :: Qunit         ! electronic charge unit
  character(len=20) :: Junit         ! current density unit
  character(len=20) :: uunit         ! specific internal energy unit
  character(len=20) :: dudtunit      ! rate of change of u unit
  character(len=20) :: tempunit      ! temperature unit

! Scaling conversion factors
  real(kind=DP) :: rscale            ! length scaling factor
  real(kind=DP) :: mscale            ! mass scaling factor
  real(kind=DP) :: tscale            ! time scaling factor
  real(kind=DP) :: vscale            ! velocity scaling factor
  real(kind=DP) :: ascale            ! acceleration scaling factor
  real(kind=DP) :: rhoscale          ! density scaling factor
  real(kind=DP) :: sigmascale        ! column density scaling factor
  real(kind=DP) :: Pscale            ! pressure scaling factor
  real(kind=DP) :: fscale            ! force scaling factor
  real(kind=DP) :: Escale            ! energy scaling factor
  real(kind=DP) :: momscale          ! momentum scaling factor
  real(kind=DP) :: angmomscale       ! angular momentum scaling factor
  real(kind=DP) :: angvelscale       ! angular velocity scaling factor
  real(kind=DP) :: dmdtscale         ! accretion rate scaling factor
  real(kind=DP) :: Lscale            ! luminosity scaling factor
  real(kind=DP) :: kappascale        ! opacity scaling factor
  real(kind=DP) :: Bscale            ! magnetic field scaling factor
  real(kind=DP) :: Qscale            ! electronic charge scaling factor
  real(kind=DP) :: Jscale            ! current density scaling factor
  real(kind=DP) :: uscale            ! specific internal energy scaling factor
  real(kind=DP) :: dudtscale         ! dudt scaling factor
  real(kind=DP) :: tempscale         ! temperature scaling factor

! Scaling factors for cgs units
  real(kind=DP) :: rcgs              ! length unit in cgs
  real(kind=DP) :: mcgs              ! mass unit in cgs
  real(kind=DP) :: tcgs              ! time unit in cgs
  real(kind=DP) :: vcgs              ! velocity unit in cgs
  real(kind=DP) :: acgs              ! acceleration unit in cgs
  real(kind=DP) :: rhocgs            ! density unit in cgs
  real(kind=DP) :: sigmacgs          ! column density unit in cgs
  real(kind=DP) :: Pcgs              ! pressure unit in cgs
  real(kind=DP) :: fcgs              ! force unit in cgs
  real(kind=DP) :: Ecgs              ! energy unit in cgs
  real(kind=DP) :: momcgs            ! momentum unit in cgs
  real(kind=DP) :: angmomcgs         ! angular momentum unit in cgs
  real(kind=DP) :: angvelcgs         ! angular velocity unit in cgs
  real(kind=DP) :: dmdtcgs           ! accretion rate unit in cgs
  real(kind=DP) :: Lcgs              ! luminosity unit in cgs
  real(kind=DP) :: kappacgs          ! opacity unit in cgs
  real(kind=DP) :: Bcgs              ! magnetic field unit in cgs
  real(kind=DP) :: Qcgs              ! electronic charge unit in cgs
  real(kind=DP) :: Jcgs              ! current density unit in cgs
  real(kind=DP) :: ucgs              ! specific internal energy unit in cgs
  real(kind=DP) :: dudtcgs           ! dudt unit in cgs
  real(kind=DP) :: tempcgs           ! temperature unit in cgs

! Scaling factors for S.I. units
  real(kind=DP) :: r_SI              ! length unit in S.I.
  real(kind=DP) :: m_SI              ! mass unit in S.I.
  real(kind=DP) :: t_SI              ! time unit in S.I.
  real(kind=DP) :: v_SI              ! velocity unit in S.I.
  real(kind=DP) :: a_SI              ! acceleration unit in S.I.
  real(kind=DP) :: rho_SI            ! density unit in S.I.
  real(kind=DP) :: sigma_SI          ! column density unit in S.I.
  real(kind=DP) :: P_SI              ! pressure unit in S.I.
  real(kind=DP) :: f_SI              ! force unit in S.I.
  real(kind=DP) :: E_SI              ! energy unit in S.I.
  real(kind=DP) :: mom_SI            ! momentum unit in S.I.
  real(kind=DP) :: angmom_SI         ! angular momentum unit in S.I.
  real(kind=DP) :: angvel_SI         ! angular velocity unit in S.I.
  real(kind=DP) :: dmdt_SI           ! accretion rate unit in S.I.
  real(kind=DP) :: L_SI              ! luminosity unit in S.I.
  real(kind=DP) :: kappa_SI          ! opacity unit in S.I.
  real(kind=DP) :: B_SI              ! magnetic field unit in S.I.
  real(kind=DP) :: Q_SI              ! electronic charge unit in S.I.
  real(kind=DP) :: J_SI              ! current density unit in S.I.
  real(kind=DP) :: u_SI              ! specific internal energy unit in S.I.
  real(kind=DP) :: dudt_SI           ! dudt unit in S.I.
  real(kind=DP) :: temp_SI           ! temperature unit in S.I.

END MODULE scaling_module


! ============================================================================
MODULE seren_sim_module
  use definitions

  logical :: nbody_sim               ! N-body simulation flag
  logical :: nbody_sph_sim           ! Hybrid SPH/N-body simulation flags
  logical :: sph_sim                 ! SPH simulation flag

END MODULE seren_sim_module


! ============================================================================
MODULE sink_module
  use definitions

  logical :: rho_search                 ! Use density to identify sinks?
  logical :: potmin_search              ! Should sink be bottom of pot. well?
  logical :: hill_sphere_search         ! Use hill sphere for selecting sinks
  logical :: energy_search              ! Should new sink be grav. bound?
  logical :: thermal_search             ! Should thermal energy < grav. energy?
  logical :: div_v_search               ! Should new sink have div_v < 0?
  logical :: div_a_search               ! Should new sink have div_a < 0?
  logical :: timescale_search           ! Compare timescales for selecting sinks
  logical :: energy_accrete             ! Should accreted particle be bound?
  logical :: first                      ! ???
  logical :: accdo_sinks                ! accdo variable for all sinks
  integer :: stot                       ! Total number of sinks
  integer(kind=ILP) :: nlast_sinks      ! Beginning of current sink timestep
  integer(kind=ILP) :: nlevel_sinks     ! Step size of current sink timestep
  real(kind=PR) :: alpha_ss             ! Sunyaev-Shakura viscosity parameter
  real(kind=PR) :: f_accretion          ! Frac. of accretion lum. radiated away
  real(kind=PR) :: feedback_tdelay      ! Time delay to switch on 
                                        ! stellar radiative feeback
  real(kind=PR) :: feedback_minmass     ! Min. sink mass for switching on 
                                        ! stellar radiative feeback
  real(kind=PR) :: star_radius          ! (unresolved) star radius
  real(kind=PR) :: alpha_EA             ! S-S viscosity for episodic accretion
  real(kind=DP) :: dmdt_regular         ! dmdt when no EA
  real(kind=PR) :: z_factor             ! ??? (Ask Dimitri!!)
  real(kind=DP) :: laststep_sinks       ! Previous sink timestep
  real(kind=PR) :: rhosink              ! Sink creation density
  real(kind=PR) :: sinkrad              ! Sink creation radius
  real(kind=PR) :: sink_frac            ! Frac. of gas mass in sinks
  real(kind=DP) :: smooth_accrete_frac  ! Critical mass frac. before accretion
  real(kind=DP) :: smooth_accrete_dt    ! Critical timestep before accretion
  real(kind=DP) :: smooth_accrete_rho   ! Critical dens. frac before accretion

  type sink_node
     logical :: accrete                 ! Is the sink accreting particles?
     logical :: static                  ! Is the sink static or not?
#if defined(USE_MPI)
     integer :: slot                    ! Lazy way to find where sink should go
     integer :: domain                  ! MPI task that 'owns' sink
#endif
     integer :: id                      ! Sink id
     integer(kind=ILP) :: ncreate       ! nsteps when sink is created
     real(kind=DP) :: tcreate           ! Physical time when sink is created
     real(kind=PR) :: r(1:NDIM)         ! Position vectors
     real(kind=PR) :: v(1:VDIM)         ! Velocity vectors
     real(kind=PR) :: m                 ! Particle mass
     real(kind=PR) :: h                 ! Sink softening length
     real(kind=PR) :: radius            ! Sink accretion radius
     real(kind=PR) :: a(1:VDIM)         ! Acceleration vectors
     real(kind=PR) :: rold(1:NDIM)      ! Old position vectors
     real(kind=PR) :: vold(1:VDIM)      ! Old velocity vectors
     real(kind=DP) :: angmom(1:3)       ! Sink internal angular momentum
     real(kind=DP) :: angmomnet(1:3)    ! Total sink ang mom accreted
     real(kind=PR) :: gpot              ! Gravitational potential energy
     real(kind=PR) :: gpe               ! Gravitational potential energy
     real(kind=DP) :: dmdt              ! Accretion rate
     real(kind=DP) :: luminosity        ! Luminosity from unresolved star
     real(kind=DP) :: temperature       ! Surface temp. of unresolved star
     real(kind=DP) :: star_radius       ! Physical radius of unresolved star
     real(kind=DP) :: macc(1:DMDT_RANGE) ! Masses accreted in previous steps
     real(kind=DP) :: tacc(1:DMDT_RANGE) ! Times of previous steps
     real(kind=DP) :: menc              ! Mass enclosed by sink (sink + gas)
     real(kind=DP) :: taccrete          ! Mass accretion rate
     real(kind=PR) :: mmax              ! Accretion mass limit for sink
#if defined(SMOOTH_ACCRETION)
     real(kind=DP) :: trot              ! Rotational period at edge of sink
     real(kind=DP) :: racc              ! Accretion radius
     real(kind=DP) :: rho               ! SPH gas density at location of sink
#endif
#if defined(SMOOTH_ACCRETION) || defined(SINK_REMOVE_ANGMOM)
     real(kind=PR) :: cmean             ! Mean sound speed of particles in sink
     real(kind=DP) :: tvisc             ! Viscous timescale for sink
#endif
#if defined(HEALPIX)
     integer :: HPid                    ! id of corresponding UV source
#endif
#if defined(RUNGE_KUTTA2)
     real(kind=PR) :: vhalf(1:VDIM)     ! Half-step velocity arrays
#endif
#if defined(LEAPFROG_KDK)
     real(kind=PR) :: aold(1:VDIM)      ! Old accelerations for PC
#endif
!!#if !defined(GEOMETRIC_MAC)
     real(kind=PR) :: agravmag          ! Mag. of grav. acceleration
!!#endif
#if defined(DEBUG_FORCES)
     real(kind=PR) :: agrav(1:VDIM)     ! Gravitational acceleration
     real(kind=PR) :: ahydro(1:VDIM)    ! Hydro force (of accreted particles)
#endif
#if defined(USE_MPI)
#if defined(SELF_GRAVITY)
     real(kind=DP) :: remote_agrav(1:VDIM) ! Grav. accel. from remote tasks
     real(kind=DP) :: remote_gpot          ! Grav. pot. from remote tasks
#endif
#endif
  end type sink_node
  type(sink_node), allocatable :: sink(:)  ! Main sink data structure

END MODULE sink_module


! ============================================================================
MODULE stellar_module
  use definitions

  integer :: Ntable                            ! No. of table elements
  type table_node                              ! Look-up table node
     real(kind=PR) :: mass                     ! Star mass (M_sun)
     real(kind=PR) :: log_L                    ! log10 (L/L_sun)
     real(kind=PR) :: log_N_LyC                ! log10 (N_LyC/s^{-1})
     real(kind=PR) :: radius                   ! Stellar radius
     real(kind=PR) :: Teff                     ! Effective temperature (K)
     real(kind=PR) :: M_loss                   ! Mass loss rate
     real(kind=PR) :: v_wind                   ! Wind speed
  end type table_node
  type(table_node), allocatable :: stellar_table(:)  ! Table of stellar props.
  
END MODULE stellar_module


! ============================================================================
MODULE time_module
  use definitions

  logical :: synchronise_all         ! Synchronise all particle timesteps

  integer(kind=ILP) :: level_max     ! Maximum timestep level occupied
  integer(kind=ILP) :: level_step    ! Level of smallest (half) step
  integer(kind=ILP) :: n             ! Current integer time
  integer(kind=ILP) :: nbuild        ! Int. time for next tree build
  integer(kind=ILP) :: nbuildstep    ! Int. time interval between tree builds
  integer(kind=ILP) :: ndiagnext     ! Integer time for next diagnostic
                                     ! (measured in numbers of max. steps)
  integer(kind=ILP) :: ndiagstep     ! Integer steps between diagnostic
                                     ! (measured in numbers of max. steps)
  integer(kind=ILP) :: nionall       ! Int. time for next 'all' ionization
  integer(kind=ILP) :: nionallstep   ! Int. time between 'all' ionization calc.
  integer(kind=ILP) :: nionize       ! Int. time for next HEALPix walk
  integer(kind=ILP) :: nlevels       ! Number of quantized timestep levels
  integer(kind=ILP) :: nmaxsteps     ! Number of maximum timesteps
  integer(kind=ILP) :: noutput       ! Int. time for next screen output
  integer(kind=ILP) :: noutputstep   ! Int. time interval of screen outputs
  integer(kind=ILP) :: nresync       ! Integer time to resynchronise
  integer(kind=ILP) :: nsearchnext   ! Integer time for next sink search
  integer(kind=ILP) :: nsearchstep   ! Int. time interval between sink searches
  integer(kind=ILP) :: nsinknext     ! Integer time for next sink output
  integer(kind=ILP) :: nsinkstep     ! Int. time interval between sink output
                                     ! (measured in numbers of steps)
  integer(kind=ILP) :: nsnapnext     ! Integer time for next snapshot
  integer(kind=ILP) :: nsnapstep     ! Int. time interval between snapshots
                                     ! (measured in numbers of steps)
  integer(kind=ILP) :: nsteps        ! Number of steps (Euler) or
                                     ! half-steps (RK, LF..)
  integer(kind=ILP) :: nstepsize     ! Current integer step size
  integer(kind=ILP) :: nstock        ! Integer time for next tree stock
  integer(kind=ILP) :: ntempnext     ! Integer time for next temp snapshot
                                     ! (measured in numbers of steps)
  integer(kind=ILP) :: ntempstep     ! Integer steps between temp files
                                     ! (measured in numbers of steps)
#if defined(USE_MPI)
  integer(kind=ILP) :: nloadbalance  ! Integer time of last loadbalance
#endif
  integer(kind=ILP) :: snapshot      ! Snapshot number

  real(kind=DP) :: accel_mult        ! Acceleration timestep multiplier
  real(kind=DP) :: courant_mult      ! Courant timestep multiplier
  real(kind=DP) :: dt_fixed          ! Reference timestep for creating levels
  real(kind=DP) :: dt_max            ! Maximum stepsize
  real(kind=DP) :: endtime           ! End time
  real(kind=DP) :: firstsnap         ! Time for first snapshot file
  real(kind=DP) :: lastsnap          ! Time of last snapshot
  real(kind=DP) :: nextsnap          ! Time for next snapshot
  real(kind=DP) :: sink_mult         ! Sink timestep multiplier
  real(kind=DP) :: snaptime          ! Snapshot interval (time)
  real(kind=DP) :: time              ! Physical time
  real(kind=DP) :: timestep          ! Real timestep of smallest integer step

  real(kind=DP) :: sph_endtime       ! End time of SPH simulation
  real(kind=DP) :: nbody_sph_endtime ! End time of N-body/SPH simulation

END MODULE time_module


! ============================================================================
MODULE timing_module
  use definitions
  implicit none

  character(len=50) :: marker_id(1:NBLOCKS)  ! Array of marker strings
  integer           :: last_id               ! Number of last marker
  integer(kind=ILP) :: last_itime            ! Last integer time mark
  integer(kind=ILP) :: itime                 ! Total integer time taken
  integer           :: mark_tot              ! Toal number of markers
  integer(kind=ILP) :: iblock(1:NBLOCKS)     ! Integer time taken by each block
  integer(kind=ILP) :: ngravcomp             ! No. of grav. computations
  integer(kind=ILP) :: nhydrocomp            ! No. of hydro computations
  integer(kind=ILP) :: nsphcomp              ! No. of SPH computations
  real(kind=DP) :: last_rtime                ! Last real time mark
  real(kind=DP) :: rtime                     ! Total real time taken
  real(kind=DP) :: rblock(1:NBLOCKS)         ! Real time taken by each block

END MODULE timing_module


! ============================================================================
MODULE turbulence_module
  use definitions
  implicit none
  
#if defined(TURBULENT_FORCING)
  complex(kind=PR), allocatable :: turb_last(:,:,:,:)
                                       ! Turbulent spectrum at time = t
  complex(kind=PR), allocatable :: turb_next(:,:,:,:)
                                       ! Turbulent spectrum at time = t + dt
  real(kind=PR), allocatable    :: afield_last(:,:,:,:)
                                       ! Forcing field at time = t
  real(kind=PR), allocatable    :: afield_next(:,:,:,:)
                                       ! Forcing field at time = t + dt
  real(kind=PR), allocatable    :: afield_now(:,:,:,:)
                                       ! Forcing field now
  real(kind=PR), allocatable    :: power_spec(:,:,:)
                                       ! Power spectrum of turbulence
  real(kind=PR)    :: sol_frac         ! Solenoidal fraction
  real(kind=DP)    :: turb_last_time   ! Time of old turbulent field
  real(kind=DP)    :: turb_next_time   ! Time of next turbulent field
  real(kind=DP)    :: turb_dt          ! Turbulent velocity evolution timestep
  real(kind=PR)    :: turb_decay_frac  ! Decay fraction per dt
  real(kind=PR)    :: turb_space(1:3)  ! Grid spacing
  real(kind=PR)    :: turb_gs_real     ! real(TURB_GS, PR)
  real(kind=PR)    :: turb_norm        ! Normalizing constant from combination
                                       ! of Ornstein-Uhlenbeck process, initial
                                       ! power spectrum and projection
  integer          :: nturbtemp        ! Either 1 or 2 (redundant file output)
  logical          :: turb_changed     ! Have we updated fields since output?
  character(len=256) :: turb_file_last    ! filename for 'last' field
  character(len=256) :: turb_file_next    ! filename for 'last' field
  character(len=256) :: turb_file_header  ! filename for header file
#endif
  real(kind=PR)    :: comp_frac        ! Compressive fractions
  real(kind=DP)    :: turb_T           ! Turbulent velocity autocorrelation time
  integer          :: turb_Ndt         ! Number of timesteps per autocorr. time
  real(kind=PR)    :: turb_rms         ! rms turbulent forcing acceleration
  real(kind=PR)    :: turb_min(1:3)    ! Minimum extent of forcing box
  real(kind=PR)    :: turb_max(1:3)    ! Maximum extent of forcing box
  
END MODULE turbulence_module


! ============================================================================
MODULE tree_module
  use definitions

  real(kind=PR) :: thetamaxsqd       ! Opening angle criterion squared
  real(kind=PR) :: abserror          ! Absolute error parameter (Gadget MAC)


! BH tree variables
! ----------------------------------------------------------------------------
#if defined(BH_TREE)
  integer :: cmax_grav                   ! Maximum allowed no. of cells
  integer :: ctot_grav                   ! Total number of cells
  integer :: ltot_grav                   ! Bottom level of gravity tree
  integer :: first_cell_grav(0:LMAX)     ! id of first cell on level
  integer :: last_cell_grav(0:LMAX)      ! id of last cell on level

  integer :: cmax_hydro                  ! Maximum possible number of cells
  integer :: ctot_hydro                  ! Total number of cells
  integer :: ltot_hydro                  ! Bottom level of hydro tree
  integer :: first_cell_hydro(0:LMAX)    ! id of first cell on level
  integer :: last_cell_hydro(0:LMAX)     ! id of last cell on level

  integer :: cmax_skeleton               ! Maximum possible number of cells
  integer :: ctot_skeleton               ! Total number of cells
  integer :: ltot_skeleton               ! Bottom level of hydro tree
  integer :: first_cell_skeleton(0:LMAX) ! id of first cell on level
  integer :: last_cell_skeleton(0:LMAX)  ! id of last cell on level

#if defined(GHOST_PARTICLES)
  integer :: cmax_ghost                  ! Maximum possible number of cells
  integer :: ctot_ghost                  ! Total number of cells
  integer :: ltot_ghost                  ! Bottom level of hydro tree
  integer :: first_cell_ghost(0:LMAX)    ! id of first cell on level
  integer :: last_cell_ghost(0:LMAX)     ! id of last cell on level
#if defined(PERIODIC)
  integer :: ctot_periodic               ! Total number of cells
  integer :: ltot_periodic               ! Bottom level of hydro tree
  integer :: first_cell_periodic(0:LMAX) ! id of first cell on level
  integer :: last_cell_periodic(0:LMAX)  ! id of last cell on level
#endif
#endif

#if defined(REORDER_TREE)
  integer, allocatable :: gravorder(:)   ! Order of grav tree cells in memory
  integer, allocatable :: hydroorder(:)  ! Order of hydro tree cells in memory
#endif

  type BHgrav_node                      ! BH gravity tree node
     integer :: leaf                    ! Number of particles if leaf
     integer :: plist(1:LEAFMAX)        ! List of particle ids
     integer :: nextcell                ! Next cell in list if not opening
     integer :: ifopen                  ! First child cell if opened
     real(kind=PR) :: r(1:NDIM)         ! Centre of mass of cell c
     real(kind=PR) :: m                 ! Total mass of cell c
     real(kind=PR) :: dminsqd           ! Cell-opening distance squared
#ifndef GEOMETRIC_MAC
     real(kind=PR) :: mac               ! Multipole-acceptance criterion
#endif
#if defined(CELL_VELOCITIES)
     real(kind=PR) :: v(1:NDIM)         ! Centre of mass velocity of cell c
#endif
#if defined(QUADRUPOLE)
     real(kind=PR) :: q(1:NQUAD)        ! Quadrupole moment terms
#endif
#if defined(OCTUPOLE)
     real(kind=PR) :: s(1:NOCT)         ! Octupole moment terms
#endif
  end type BHgrav_node
  type(BHgrav_node), allocatable :: BHgrav(:)  ! BH gravity tree array

  type BHhydro_node                     ! BH hydro tree node
     integer :: leaf                    ! Number of particles if leaf
     integer :: nextcell                ! Next cell in list if not opening
     integer :: ifopen                  ! First child cell if opened
     integer :: plist(1:LEAFMAX)        ! List of particle ids
     real(kind=PR) :: r(1:NDIM)         ! Centre of position of cell c
     real(kind=PR) :: rmax              ! Distance of furthest particle in c
     real(kind=PR) :: hrangemax         ! Maximum extent of kernel
  end type BHhydro_node
  type(BHhydro_node), allocatable :: BHhydro(:) ! BH hydro tree array
#if defined(GHOST_PARTICLES)
  type(BHhydro_node), allocatable :: BHghost(:) ! BH ghost tree array
#if defined(PERIODIC)
  type(BHhydro_node), allocatable :: BHperiodic(:) ! BH periodic ghost tree
#endif
#endif

  type BH_skeleton_node                 ! Skeleton tree node
     integer :: pfirst                  ! First particle in linked list
     integer :: plast                   ! Last particle in linked list
     integer :: leaf                    ! Number of particles if leaf
     integer :: nextcell                ! Next cell in list if not opening
     integer :: ifopen                  ! First child cell if opened
     integer :: plist(1:LEAFMAX)        ! List of particle ids
     integer :: childof(1:NCHILD)       ! List if child cell ids
     real(kind=PR) :: r(1:NDIM)         ! Centre of position of cell c
  end type BH_skeleton_node
  type(BH_skeleton_node), allocatable :: BHskeleton(:) ! Skeleton tree array

  integer, allocatable :: cellof(:)     ! Cell particle p is in
  integer, allocatable :: BHnextptcl(:) ! Next particle in linked list
  integer, allocatable :: whichchild(:) ! Child cell particle p is in

  type auxilary_node                    ! Auxilary cell node
     real(kind=PR) :: hmax              ! Max value of h in a cell
     real(kind=PR) :: bbmin(1:NDIM)     ! rmin
     real(kind=PR) :: bbmax(1:NDIM)     ! rmax
  end type auxilary_node
  type(auxilary_node), allocatable :: BHstock(:)  ! bounding box array

#endif

END MODULE tree_module


! ============================================================================
MODULE type_module
  use definitions

  integer, parameter :: boundaryid = 1     ! Boundary particle type id
  integer, parameter :: icmid = 2          ! ICM particle type id
  integer, parameter :: gasid = 3          ! Gas particle type id
  integer, parameter :: cdmid = 4          ! CDM particle type id
  integer, parameter :: dustid = 5         ! Dust particle type id
  integer, parameter :: ionid = 6          ! Ion particle type id
  integer, parameter :: nmaxtypes = 6      ! Max. no of particle types

  integer :: ntypes                        ! No. of types
  integer :: typeorder(1:nmaxtypes)        ! Order of particle types in memory

  integer :: pboundary                     ! No. of boundary ptcls
  integer :: pcdm                          ! No. of cold dark matter ptcls
  integer :: pdust                         ! No. of dust ptcls
  integer :: pgas                          ! No. of gas ptcls
  integer :: picm                          ! No. of ICM ptcls
  integer :: pion                          ! No. of ion ptcls
  integer :: pghost                        ! No. of ghost ptcls
  integer :: pperiodic                     ! No. of periodic ghost ptcls
                                           ! (= pghost for non-MPI)
#if defined(USE_MPI)
  integer :: last_pghost                   ! Previous number of ghosts
#endif

  integer :: pboundarystart                ! first boundary ptcl
  integer :: pboundaryend                  ! last boundary ptcl
  integer :: pcdmstart                     ! first cdm ptcl
  integer :: pcdmend                       ! last cdm ptcl
  integer :: pduststart                    ! first dust ptcls
  integer :: pdustend                      ! last dust ptcls
  integer :: pgas_orig                     ! Original no. of gas ptcls
  integer :: pgasstart                     ! first gas ptcl
  integer :: pgasend                       ! last gas ptcl
  integer :: picmstart                     ! first ICM ptcl
  integer :: picmend                       ! last ICM ptcl
  integer :: pionstart                     ! first ion ptcl
  integer :: pionend                       ! last ion ptcl
  integer :: pionized                      ! No. of ionized ptcls in HII sim
  integer :: pgravityend                   ! last ptcl in gravity loop
  integer :: pgravitystart                 ! first ptcl in gravity loop
  integer :: phydroend                     ! last ptcl in hydro loop
  integer :: phydrostart                   ! first ptcl in hydro loop
  integer :: pionizestart                  ! first ptcl in ionization routine

  integer :: phtreestart                   ! first ptcl in neighbour tree
  integer :: phtreeend                     ! last ptcl in neighbour tree
  integer :: pgtreestart                   ! first ptcl in gravity tree
  integer :: pgtreeend                     ! last ptcl in gravity tree
  integer :: psftreestart                  ! first ptcl in 2nd fluid tree
  integer :: psftreeend                    ! last ptcl in 2nd fluid tree

  real(kind=DP) :: mgas                    ! Total gas mass
  real(kind=DP) :: mgas_orig               ! Initial Total gas mass
  real(kind=DP) :: micm                    ! Total mass of icm ptcls
  real(kind=DP) :: mboundary               ! Total mass of boundary ptcls
  real(kind=DP) :: mcdm                    ! Total mass of cdm ptcls
  real(kind=DP) :: mmean                   ! Mean SPH gas ptcl mass 
                                           ! (N.B. not mu_bar)

  logical :: sphmask(1:nmaxtypes)          ! SPH mask
  logical :: hydromask(1:nmaxtypes)        ! Hydro force mask
  logical :: gravmask(1:nmaxtypes)         ! Gravity force mask
  logical :: dragmask(1:nmaxtypes)         ! Gas-dust drag term mask

  type particle_type                       ! Particle type node
     character(len=256) :: name            ! Particle type name
     character(len=256) :: eos             ! EOS for particle type
     integer :: dragonid                   ! Dragon id for particle id
     integer :: N                          ! No. of particles of type
     integer :: pfirst                     ! first particle
     integer :: plast                      ! last particle
     logical :: static                     ! Do particles move?
     logical :: h(1:nmaxtypes)             ! Mask to determine h
     logical :: h2(1:nmaxtypes)            ! Mask to determine second h
     logical :: hydro(1:nmaxtypes)         ! Mask to compute hydro forces
     logical :: grav(1:nmaxtypes)          ! Mask to compute grav forces
     logical :: drag(1:nmaxtypes)          ! Mask to compute drag terms
  end type particle_type
  type(particle_type) :: typeinfo(1:nmaxtypes)  ! Array of particle type info


END MODULE type_module


! ============================================================================
