# ----------------------------------------------------------------------------
# Seren Makefile version 1.5.1
# Date : 01/07/2013
# ----------------------------------------------------------------------------
F90                       = ifort
MPIF90                    = mpif90
VERSION_NO                = 1.5.1
SRCDIR                    = $(PWD)/src
EXEDIR                    = $(PWD)
OPENMP                    = 1
MPI                       = 0
MPI_LIBRARY               = mpich2
COMPILER_MODE             = FASTDEBUG
OUTPUT_LEVEL              = 2
DIAGNOSTIC_OUTPUT         = 1
NDIM                      = 3
PRECISION                 = SINGLE
PERIODIC                  = 0
X_BOUNDARY                = 0
Y_BOUNDARY                = 0
Z_BOUNDARY                = 0
GHOST_PARTICLES           = 0

# ----------------------------------------------------------------------------
# Simulation selection options
# ----------------------------------------------------------------------------
SPH_SIMULATION            = 1
NBODY_SIMULATION          = 0

# ----------------------------------------------------------------------------
# SPH simulation options
# ----------------------------------------------------------------------------
SPH                       = GRAD_H_SPH
MHD                       = 1
SPH_INTEGRATION           = LFDKD
KERNEL                    = M4
HFIND                     = H_RHO
MINIMUM_H                 = 0
HYDRO                     = 1
ENERGY_EQN                = 0
ENTROPY_EQN               = 0
ARTIFICIAL_VISCOSITY      = MON97
VISC_TD                   = 1
BALSARA                   = 0
ARTIFICIAL_CONDUCTIVITY   = 0
EXTERNAL_PRESSURE         = 0

# ----------------------------------------------------------------------------
# Heating, cooling and feedback options
# ----------------------------------------------------------------------------
RAD_WS                    = 0
FLUX_LIMITED_DIFFUSION    = 0
SINK_POTENTIAL_WS         = 0
AMBIENT_HEATING_WS        = 0
SINK_HEATING_WS           = 0
COOLING_HEATING           = 0
IONIZING_RADIATION        = 0
STELLAR_WIND              = 0

# ----------------------------------------------------------------------------
# SPH gravity options
# ----------------------------------------------------------------------------
EXTERNAL_FORCE            = 0
SELF_GRAVITY              = KS
MEANH_GRAVITY             = 0
EWALD                     = 0

# ----------------------------------------------------------------------------
# Sink and N-body options
# ----------------------------------------------------------------------------
SINKS                     = SIMPLE
SINK_RADIUS               = HMULT
SINK_REMOVE_ANGMOM        = 0
SINK_GRAVITY_ONLY         = 0
NBODY_INTEGRATION         = HERMITE4
BINARY_STATS              = 0

# ----------------------------------------------------------------------------
# Tree options
# ----------------------------------------------------------------------------
TREE                      = BH
MULTIPOLE                 = QUADRUPOLE
MAC                       = GADGET
REORDER                   = 0

# ----------------------------------------------------------------------------
# Misc. options
# ----------------------------------------------------------------------------
SORT                      = INSERTION
TIMESTEP                  = ADAPTIVE
CHECK_NEIB_TIMESTEP       = 1
NEIGHBOURLISTS            = 0
KERNEL_TABLES             = 0
REMOVE_OUTLIERS           = 0
TURBULENT_FORCING         = 0
TIMING_CODE               = 1
TEST                      = 0

# ----------------------------------------------------------------------------
# Debugging options
# ----------------------------------------------------------------------------
ifneq ($(OUTPUT_LEVEL),0)
#DFLAGS += -DIEEE_EXCEPTION_HANDLING
#DFLAGS += -DDIV_A
#DFLAGS += -DDEBUG_ACCRETE
#DFLAGS += -DDEBUG_ALLOCATE_MEMORY
#DFLAGS += -DDEBUG_BHTREEBUILD
#DFLAGS += -DDEBUG_BHTREESTOCK
#DFLAGS += -DDEBUG_BHTREEWALK
#DFLAGS += -DDEBUG_BHTREEGRAVITY
#DFLAGS += -DDEBUG_BINARY_PROPERTIES
#DFLAGS += -DDEBUG_BINARY_SEARCH
#DFLAGS += -DDEBUG_BLOCK_TIMESTEPS
#DFLAGS += -DDEBUG_COPY_PARTICLE_DATA
#DFLAGS += -DDEBUG_CREATE_SINK
#DFLAGS += -DDEBUG_CREATE_HP_SOURCE
#DFLAGS += -DDEBUG_DENSITY
#DFLAGS += -DDEBUG_DIV_V
#DFLAGS += -DDEBUG_DUDTRAD
#DFLAGS += -DDEBUG_ENERGY_EQN
#DFLAGS += -DDEBUG_FOLIATE
DFLAGS += -DDEBUG_FORCES
#DFLAGS += -DDEBUG_FREEFALL
#DFLAGS += -DDEBUG_GATHER_NEIB
#DFLAGS += -DDEBUG_GET_NEIB
#DFLAGS += -DDEBUG_GHOST_PARTICLES
#DFLAGS += -DDEBUG_GRAD_H_SPH
#DFLAGS += -DDEBUG_GRID_RESULTS
#DFLAGS += -DDEBUG_HEAPSORT
#DFLAGS += -DDEBUG_HERMITE4
#DFLAGS += -DDEBUG_H_GATHER
#DFLAGS += -DDEBUG_H_GATHER_DENSITY
#DFLAGS += -DDEBUG_H_GUESS
#DFLAGS += -DDEBUG_H_RHO_ITERATION
#DFLAGS += -DDEBUG_HP_IF
#DFLAGS += -DDEBUG_HP_OUTPUT
#DFLAGS += -DDEBUG_HP_SPLIT_ACTIVE_RAYS
#DFLAGS += -DDEBUG_HP_WALK_ALL_RAYS
#DFLAGS += -DDEBUG_HP_WALK_RAY
#DFLAGS += -DDEBUG_INTEGRATE
#DFLAGS += -DDEBUG_KERNEL
#DFLAGS += -DDEBUG_NBODYSETUP
#DFLAGS += -DDEBUG_OUTPUT_STAR_DATA
#DFLAGS += -DDEBUG_PARAMETERS
#DFLAGS += -DDEBUG_PLOT_DATA
#DFLAGS += -DDEBUG_RAD
#DFLAGS += -DDEBUG_REDUCE_TIMESTEP
#DFLAGS += -DDEBUG_REMOVE_OUTLIERS
#DFLAGS += -DDEBUG_RSPH_OUTPUT
#DFLAGS += -DDEBUG_SINK_REMOVE_ANGMOM
#DFLAGS += -DDEBUG_SINK_SEARCH
#DFLAGS += -DDEBUG_SINK_TIMESTEP
#DFLAGS += -DDEBUG_SKELETON
#DFLAGS += -DDEBUG_SMOOTH_ACCRETE_PARTICLES
#DFLAGS += -DSMOOTHED_VELOCITY
#DFLAGS += -DDEBUG_SPH_UPDATE
#DFLAGS += -DDEBUG_SWAP_PARTICLE_DATA
DFLAGS += -DDEBUG_TIMESTEP_SIZE
#DFLAGS += -DDEBUG_TRACK_PARTICLE
#DFLAGS += -DDEBUG_TREE_BUILD
#DFLAGS += -DDEBUG_TREE_GRAVITY
#DFLAGS += -DDEBUG_TREESTOCK
#DFLAGS += -DDEBUG_TREEWALK
#DFLAGS += -DDEBUG_TYPES
#DFLAGS += -DDEBUG_VISC_BALSARA
#DFLAGS += -DDEBUG_VISC_PATTERN_REC
#DFLAGS += -DDEBUG_WRITE_MPI_TASK
#DFLAGS += -C=all
#DFLAGS += -C=undefined
#DFLAGS += -Wall -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace
CFLAGS += $(DFLAGS)
endif


# Now include the Makefile tail which contains all the processed options 
# and generates list of object files for compilation
# ----------------------------------------------------------------------------
include makefiletail.mk


# List of possible flags for making the code
# ----------------------------------------------------------------------------
# F90                     : FORTRAN compiler
#		            f95      = NAG f95 compiler
#                           g95      = free f95 compiler (not gnu)
#                           gfortran = gnu f95 compiler
#                           pgf90    = Portland group compiler (coma)
#                           pgf95    = Portland group compiler (iceberg)
#                           ifort    = Intel Fortran compiler
#                           If using MPI, leave this set to the non-mpi version
#                           of the compiler from the list above

# MPIF90                  : MPI FORTRAN compiler 
#                           Set this to the MPI wrapper compiler, normally
#                           mpif90 or mpiifort

# SRCDIR                  : Absolute path to main Seren directory

# EXEDIR                  : Absolute path to location of Seren executable

# OPENMP                  : OpenMP parallelisation (0 or 1)

# MPI                     : MPI parallelisation (0 or 1; still in development)

# MPI_LIBRARY             : MPI library for compiling code
#                           Set to 'intel' if using Intel MPI

# OPTIMISE                : Compiler optimisation level (0, 1, 2 or 3)

# COMPILER_MODE           : Select flags for compiler
#                           STANDARD = Standard 'O3' compilation flags
#                           FAST     = O3 + fast (potentially unsafe) flags
#                           DEBUG    = Debugging flags to use gdb

# OUTPUT_LEVEL            : 0 = No standard output
#                           1 = Simple output to screen (time, snapshots, etc)
#                           2 = Code section location information
#                           3 = Output on particle loops

# NDIM                    : Number of dimensions (1, 2 or 3)

# PRECISION               : Precision of floating point variables
#                           SINGLE    = 32-bit floating point precision
#                           DOUBLE    = 64-bit floating point precision

# INFILE_FORMAT           : Available input file format options in SEREN
#                           ALL    = All possible file formats enabled
#                           DRAGON = DRAGON file format
#                           SEREN  = SEREN file format
#                           ASCII  = Simple ASCII file format

# OUTFILE_FORMAT          : Available output file format options in SEREN
#                           ALL    = All possible file formats enabled
#                           SEREN  = SEREN file format
#                           ASCII  = Simple ASCII file format

# PERIODIC                : Periodic boundary conditions (0 or 1)

# X_BOUNDARY              : Boundary-type in x-direction
#                           0        = No boundaries
#                           PERIODIC = Periodic boundaries
#                           WALL     = LHS and RHS walls in x-dimension
#                           WALL_LHS = LHS wall in x-dimension
#                           WALL_RHS = RHS wall in x-dimension

# Y_BOUNDARY              : as X_BOUNDARY
# Z_BOUNDARY              : as X_BOUNDARY

# SPHERICAL_WALL          : (0 or 1)

# CYLINDRICAL_WALL        : (0 or 1)

# GHOST_PARTICLES         : Use ghost particles for periodic/domain boundaries
#                           (Experimental)
#                           0 = No ghosts; use standard periodic approach
#                           1 = Use ghost particles (Not working fully yet)

# SPH_SIMULATION          : (0 or 1)

# NBODY_SPH_SIMULATION    : (0 or 1)

# NBODY_SIMULATION        : (0 or 1)

# SPH                     : SPH mode
#                           STANDARD   = Traditional SPH (e.g. Monaghan 1992)
#                           GRAD_H_SPH = 'grad-h' conservative SPH 
#                                        (e.g. Price & Monaghan 2005)
#                           SM2012_SPH = Saitoh & Makino (2012) SPH.
#                           RTSPH      = Ritchie & Thomas (2001) SPH.  
#                                        Same as STANDARD with modified density
#                           RPSPH      = 'Relative pressure' SPH (Abel 2010)

# SPH_INTEGRATION         : Integration scheme used
#                           EULER = 1st order Euler scheme
#                           RK2   = 2nd order Runge-Kutta
#                           LFKDK = 2nd order Kick-drift-kick Leapfrog
#                           LFDKD = 2nd order Drift-kick-drift Leapfrog

# KERNEL                  : Kernel used to calculate SPH quantities
#                           M4          = M4 kernel (Monaghan & Lattanzio 1985)
#                           M4TC        = Modified M4 kernel gradient 
#                                         (Thomas & Couchman 1992)
#                           QUINTIC     = Quintic polynomial kernel 
#                                         (Morris 1996)
#                           QUINTICTC   = As QUINTIC, but with modified grad.
#                           GAUSSIAN_3H = Gaussian kernel truncated at 3h
#                           LINEAR      = Simple linear function designed 
#                                         for simple testing purposes.

# HFIND                   : Method to calculate h in standard SPH formulation
#                           NUMBER   = h contains specified no. of neibs
#                           MASS     = h contains specified mass
#                           CONSTANT = Use a constant value of h
#                           H_RHO    = Iterate h and rho to become consistent 
#                                      with relation : h = const*(m/rho)^{1/D)

# MINIMUM_H               : Enforce a minimum smoothing length (0 or 1)

# HYDRO                   : Hydrodynamical forces (0 or 1)

# ENERGY_EQN              : Activate Energy equation in compliation (0 or 1)

# ENTROPY_EQN             : Activate Entropy equation in compliation (0 or 1)

# RAD_WS                  : Activate Radiative cooling method of 
#                           Stamatellos et al. (2007) in compilation (0 or 1)

# FLUX_LIMITED_DIFFUSION  : Hybrid flux-limited diffusion method (0 or 1)

# SINK_POTENTIAL_WS       : Use sink potential in RAD_WS method (0 or 1)

# AMBIENT_HEATING_WS      : External heating source (0 or 1)

# SINK_HEATING_WS         : Options for heating gas from sinks
#                           0                   = No heating from sinks
#                           STAR_HEATING        = ..
#                           STAR_SIMPLE_HEATING = ..
#                           HDISC_HEATING = heating by central star in disc 
#                                (assumes one star/disc with disc in x-y plane)
#                                T = sqrt{To^2 [(R^2+Ro^2)/AU^2]^(-q)+Tinf^2}
#                                i.e. T ~ To*(R/Ro)^(-q), 
#                                where R is the distance in disc midplane 
#                                (To, Ro (in AU), Tinf set in params.dat)

# COOLING_HEATING         : Method for integrating cooling/heating terms 
#                           (Experimental)
#                           0           = No cooling/heating
#                           EXPLICIT    = Include cooling/heating in energy 
#                                         equation for explicit integration.
#                           EXPONENTIAL = Exponential integration for 
#                                         simple functions

# IONIZING_RADIATION      : Include sources of ionizing radiation or not
#                           0                     = No sources of 
#                                                   ionizing radiation
#                           SINGLE_STATIC_SOURCE  = Single static source 
#                                                   located at position given  
#                                                   in params file.
#                           SINGLE_SINK_SOURCE    = Sink particle 1 becomes 
#                                                   source of ionizing rad.

# STELLAR_WIND            : Include sources of momentum winds using HEALPix
#                           0                     = No winds
#                           SINGLE_STATIC_SOURCE  = Single static wind source 
#                                                   located at position given 
#                                                   in params file
#                           SINGLE_SINK_SOURCE    = Sink particle 1 becomes 
#                                                   source of momentum wind.

# ARTIFICIAL_VISCOSITY    : Artificial viscosity formulation
#                           0      = No artificial viscosity used
#                           AB     = Alpha-Beta (standard)
#                           MON97  = Monaghan Riemann-viscosity

# BALSARA                 : Balsara switch for artificial viscosity (0 or 1)

# VISC_TD                 : Time-dependent artificial viscosity (0 or 1)

# PATTERN_REC             : Keplerian pattern recognition visc. switch (0 or 1)

# ARTIFICIAL_CONDUCTIVITY : Artificial conductivity formulation
#                           PRICE2008   = Price (2008) conductivity
#                           WADSLEY2008 = Wadsley et al. (2008) conductivity

# EXTERNAL_PRESSURE       : Simple external pressure formulation (0 or 1)

# EXTERNAL_FORCE          : External gravitational force options
#                           0        = No external gravitational force
#                           PLUMMER  = Plummer potential
#                           UDS      = Uniform-density sphere
#                           NFW1996  = Navarro, Frenk & White (1996) potential

# SELF_GRAVITY            : Method of computing gravitational force
#		            0     = No gravitational forces computed
#		            KS    = Kernel-softened gravity for 2-body forces
#                           NBODY = Newton's grav. law for all 2-body forces

# MEANH_GRAVITY           : Use mean-h method of symmeterising gravitational 
#                           forces (0 or 1)

# EWALD                   : Ewald periodic gravity forces (0 or 1)

# SINKS                   : Sink particle options
#                           0           = No sinks used
#                           SIMPLE      = Simple sink particles; accretion only
#                           NO_ACC      = Simple sinks with no accretion
#                           SMOOTH_ACC  = Sinks with smooth accretion 

# SINK_RADIUS             : Method of selecting the sink radius for new sinks
#                           FIXED_ABSOLUTE = Absolute (constant) sink radius
#                           FIXED_HMULT    = Multiple of h 
#                                            (same sinkrad for all)
#                           HMULT          = Multiple of h 
#                                            (indiv. sinkrad values)

# SINK_REMOVE_ANGMOM      : Redistribute angular momentum from sinks (0 or 1)

# SINK_GRAVITY_ONLY       : Sink gravity options
#                           0 = Self-gravity from all SPH gas/sink particles
#                           1 = No self-gravity; sink-gravity only 
#                              (experimental)

# NBODY_INTEGRATION       : Integration scheme for N-body sim
#                           HERMITE4    = 4th order Hermite scheme
#                           HERMITE4_TS = 4th order Hermite time-symmetric
#                           LFKDK       = 2nd order Leapfrog kick-drift-kick

# BINARY_STATS            : Automatically calculate binary statistics for 
#                           stars during N-body integrator (0 or 1)

# TREE                    : Tree options
#                           0         = No tree (direct-sum only)
#                           BH        = Barnes-Hut Octal-spatial tree
#                           BINARY    = Binary-number tree (not fully working)

# MULTIPOLE               : Gravity tree multipole expansion options
#                           0          = Monopole moments only
#                           QUADRUPOLE = Expand to quadrupole order
#                           OCTUPOLE   = Expand to octupole order

# MAC                     : Multipole-acceptance criterion used
#                           GEOMETRIC = Original Barnes-Hut geometric criterion
#                           GADGET    = GADGET-style MAC
#                           GADGET2   = GADGET2-style MAC
#                           EIGEN     = Eigenvalue MAC

# REORDER                 : Options for re-ordering particles/tree
#                           0         = No reordering of particles
#                           PARTICLES = Only reorder particle array
#                           TREE      = Only reorder tree cells
#                           ALL       = Reorder particle and tree arrays 

# CELL_WALK               : 0 or 1 (Experimental)

# SORT                    : Sorting algorithm 
#                           INSERTION = Insertion sort
#                           HEAP      = Heapsort

# TIMESTEP                : Options for multiple-particle timestepping
#                           ADAPTIVE   = Block timestep levels adjusted at 
#                                        resync
#                           RESTRICTED = Timestep levels only take certain 
#                                        values (dtmax parameter times 
#                                        integer power of 2)
#                           FIXED      = Fixed block timestep levels

# CHECK_NEIB_TIMESTEP     : Options for reducing particles timesteps
#                           0 = No checking of neighbour timesteps
#                           1 = Only reduce timesteps at end of current 
#                               timestep
#                           2 = Reduce timestep in middle of current timestep 
#                               if required

# NEIGHBOURLISTS          : Store neighbour lists in memory
#                           0         = Neighbour lists not stored
#                           PARTICLES = Particle neighbour lists stored

# KERNEL_TABLES           : Tabulate kernels for fast referencing (0 or 1)

# REMOVE_OUTLIERS         : 0 or 1 (Experimental)

# TURBULENT_FORCING       : Impose time-dependent turbulent forcing (0 or 1)

# TIMING_CODE             : Use internal timing routines (0 or 1)

# TEST                    : Specific test options
#                           0        = No test options
#                           SPIEGEL  = Spiegel (ref??) test
#                           FREEFALL = Freefall collapse test
#                           BINARY   = Orbitting binary stars test
#                           PLUMMER  = Plummer sphere stability test
#                           ENTROPY  = Entropy-core test
#                           WIND     = Momentum winds test

