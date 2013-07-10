! MACROS.H
! C. P. Batty & D. A. Hubber - 19/11/2007
! Definition of macros used for 
! 1. Array-element ids
! 2. Common constants
! 3. Concise debugging statements (Function-like macros)
! ============================================================================

! Seren version string
! ----------------------------------------------------------------------------
#define SEREN_VERSION "1.6.0"

! Numerical constants
! ----------------------------------------------------------------------------
#define PI 3.141592653589793_PR
#define TWOPI 6.283185307179586_PR
#define INVPI 0.318309886183791_PR
#define INVLOGETWO 1.442695040888963_PR
#define INVLOG10TWO 3.321928094887362_PR
#define INVSQRTTWO 0.707106781186548_PR
#define ONETHIRD 0.333333333333333_PR
#define ONESIXTH 0.166666666666666_PR
#define TWOTHIRDS 0.666666666666666_PR

#define PI_DP 3.141592653589793_DP
#define TWOPI_DP 6.283185307179586_DP
#define ONETHIRD_DP 0.333333333333333333_DP
#define ONESIXTH_DP 0.166666666666666666_DP
#define TWOTHIRDS_DP 0.666666666666666666_DP

! Macros for referencing grouped arrays
! ----------------------------------------------------------------------------
#if NDIM == 1
#define MASS 2
#define SMOO 3
#define ZETA 4
#define NDIMPR 1.0_PR
#define INVNDIM 1.0_PR
#define NDIMPLUS1 2
#elif NDIM == 2
#define MASS 3
#define SMOO 4
#define ZETA 5
#define NDIMPR 2.0_PR
#define INVNDIM 0.5_PR
#define NDIMPLUS1 3
#elif NDIM == 3
#define MASS 4
#define SMOO 5
#define ZETA 6
#define NDIMPR 3.0_PR
#define INVNDIM ONETHIRD
#define NDIMPLUS1 4
#endif
#if defined(GRAVITY)
#define DATATOT ZETA
#else
#define DATATOT SMOO
#endif

! No. of velocity and B-field dimensions depending on compiler options
! ----------------------------------------------------------------------------
#if NDIM == 1 && defined(IDEAL_MHD)
#define BDIM 2
#define VDIM 2
#elif NDIM == 1 && !defined(IDEAL_MHD)
#define BDIM 1
#define VDIM 1 
#elif NDIM == 2 && defined(IDEAL_MHD)
#define BDIM 3
#define VDIM 3
#elif NDIM == 2 && !defined(IDEAL_MHD)
#define BDIM 2
#define VDIM 2 
#elif NDIM == 3
#define BDIM 3
#define VDIM 3
#endif

! Standard big and small numbers for max/min determinations 
! ----------------------------------------------------------------------------
#define BIG_NUMBER 9.9E20_PR
#define BIG_NUMBER_DP 9.9E20_DP
#define SMALL_NUMBER 1.0E-20_PR
#define SMALL_NUMBER_DP 1.0E-20_DP

! SPH particle types
! ----------------------------------------------------------------------------
#define SINKID    -1
#define DEADID     0
#define GASID      1
#define SPLITID    4
#define BOUNDARYID 6
#define ICMID      9
#define CDMID      10

! SPH variables
! ----------------------------------------------------------------------------
#define ETA_SQD 0.01_PR
#define HMULT 1.1_PR

! Kernel variables
! ----------------------------------------------------------------------------
#if defined(M4_KERNEL)
#define KERNTOT         1000
#define HALFKERNTOT     500.0_PR
#define KERN_H          500.0_PR
#define KERNRANGE       2.0_PR
#define KERNRANGESQD    4.0_PR
#define INVKERNRANGE    0.5_PR
#define KERNRANGE_DP    2.0_DP
#define INVKERNRANGE_DP 0.5_DP
#if NDIM==1
#define KERNNORM 0.666666666666666666_PR
#elif NDIM==2
#define KERNNORM (INVPI*10.0_PR/7.0_PR)
#elif NDIM==3
#define KERNNORM INVPI
#endif
#elif defined(QUINTIC_KERNEL)
#define KERNTOT         1000
#define HALFKERNTOT     333.3_PR
#define KERN_H          333.3_PR
#define KERNRANGE       3.0_PR
#define KERNRANGESQD    9.0_PR
#define INVKERNRANGE    0.333333333333333333_PR
#define KERNRANGE_DP    3.0_DP
#define INVKERNRANGE_DP 0.333333333333333333_DP
#if NDIM==1
#define KERNNORM (1.0_PR/120.0_PR)
#elif NDIM==2
#define KERNNORM (INVPI*7.0_PR/478.0_PR)
#elif NDIM==3
#define KERNNORM (INVPI*3.0_PR/359.0_PR)
#endif
#elif defined(GAUSSIAN_KERNEL) && defined(GAUSSIAN_3H)
#define KERNTOT         1000
#define HALFKERNTOT     333.3_PR
#define KERN_H          333.3_PR
#define KERNRANGE       3.0_PR
#define KERNRANGESQD    9.0_PR
#define INVKERNRANGE    0.333333333333333333_PR
#define KERNRANGE_DP    3.0_DP
#define INVKERNRANGE_DP 0.333333333333333333_DP
#elif defined(LINEAR_KERNEL)
#define KERNTOT         1000
#define HALFKERNTOT     500.0_PR
#define KERN_H          500.0_PR
#define KERNRANGE       2.0_PR
#define KERNRANGESQD    4.0_PR
#define INVKERNRANGE    0.5_PR
#define KERNRANGE_DP    2.0_DP
#define INVKERNRANGE_DP 0.5_DP
#define KERNNORM        (3.0_PR*INVPI/(KERNRANGE**4))
#endif

! Hydro variables
! ----------------------------------------------------------------------------
#define C_1 0.1_PR
#define TVISC_FAC 0.6_DP
#define BAL_DENOM 1.0E-3_PR

! BH tree variables
! ----------------------------------------------------------------------------
#define LEAFMAX 8
#define LMAX 30
#if NDIM==1
#define NCHILD 2
#define NQUAD 1
#define NOCT 1
#elif NDIM==2
#define NCHILD 4
#define NQUAD 3
#define NOCT 4
#elif NDIM==3
#define NCHILD 8
#define NQUAD 5
#define NOCT 10
#endif

! Sink variables
! ----------------------------------------------------------------------------
#define SMAX 512
#define DMDT_RANGE 8
#define NEW_SINK_RMAX 1.0_PR
#define ANGMOMRAD 2.0_PR

#if defined(SINK_REMOVE_ANGMOM)
#define REXTENTSQD (ANGMOMRAD*ANGMOMRAD)
#else
#define REXTENTSQD 1.0_PR
#endif

! grad-h SPH macros
! ----------------------------------------------------------------------------
#if defined(NBODY_SPH_SIMULATION)
#define H_CON 0.0000000000001_PR
#define GRADH_ITERATION_MAX 150
#else
#define H_CON 0.02_PR
#define GRADH_ITERATION_MAX 30
#endif

! Ionization macros
! ----------------------------------------------------------------------------
#define RMIN_MULT 12.0_PR
#define HP_LEVELS 13

! Turbulent forcing macros
! ----------------------------------------------------------------------------
#define TURB_GS 64

! Timing macros
! ----------------------------------------------------------------------------
#define NBLOCKS 100

! Timestep macros
! ----------------------------------------------------------------------------
#define TIMESTEP_LEVEL_DIFF_MAX 2_ILP

! Memory macros
! ----------------------------------------------------------------------------
#if defined(USE_MPI) 
#define PMAXMULT 1.5_PR
#elif defined(GHOST_PARTICLES)
#define PMAXMULT 2.0_PR
#elif defined(PARTICLE_INJECTION_WINDS)
#define PMAXMULT 1.5_PR
#else
#define PMAXMULT 1.0_PR
#endif
#define LISTSIZE 512
#define GLISTSIZE 1024

! Parallelisation macros
! ----------------------------------------------------------------------------
#define CHUNKFRAC 0.002_PR

! Debug macros to reduce line numbers for debug statements 
! ----------------------------------------------------------------------------
#if defined(USE_MPI) && defined(DEBUG1)
#define debug1(x)   if (rank == 0) write (6,*) x
#elif defined(DEBUG1)
#define debug1(x)   write (6,*) x
#else
#define debug1(x)
#endif

#ifdef DEBUG2
#define debug2(x)   write (6,*) x
#else
#define debug2(x)
#endif

#ifdef DEBUG3
#define debug3(x, y)   write (6,*) x, y
#else
#define debug3(x, y)
#endif

#ifdef DEBUG4
#define debug4(x, y)   write (6,*) x, y
#else
#define debug4(x, y)
#endif

! Timing call macro
! ----------------------------------------------------------------------------
#ifdef TIMING
#define debug_timing(x)  call timing(x)
#else
#define debug_timing(x)
#endif

! MPI macros
! ----------------------------------------------------------------------------
#ifdef USE_MPI
#define MPI_ROOT if (rank == 0) then
#define MPI_NOT_ROOT if (rank /= 0) then
#define MPI_END end if
#else
#define MPI_ROOT
#define MPI_END
#define MPI_NOT_ROOT
#endif

! Macros for MPI tags
! ----------------------------------------------------------------------------
#define TRANSFER_NUM_TAG 40
#define TRANSFER_TAG 50
#define GHOST_NUM_TAG 90
#define GHOST_TAG 100
#define SCATTER_GHOST_TAG 110
#define SINK_GHOST_TAG 120
#define DIFFUSION_GHOST_TAG 130
#define GHOST_RHO_TAG 140
#define HYDRO_EXPORT_TAG 300
#define HYDRO_RETURN_TAG 310
#define GRAV_EXPORT_TAG 400
#define GRAV_RETURN_TAG 410
#define GRAVTREES_TAG 420
#define SINK_TAG 900
#define PARTICLETIMESTEP_TAG 700

! MPI timings macro
! ----------------------------------------------------------------------------
#define WAIT_TIME_MACRO waitstart = MPI_WTIME(); calctime = calctime + waitstart - calcstart
#define CALC_TIME_MACRO calcstart = MPI_WTIME(); waittime = waittime + calcstart - waitstart
#define EXPORT_START_TIME_MACRO exportstart=MPI_WTIME()
#define EXPORT_STOP_TIME_MACRO exporttime=exporttime+MPI_WTIME()-exportstart

! MPI periodic boundary macros
! ----------------------------------------------------------------------------
#define MPI_PERIODIC_OVERLAP 0.25_PR

! Cost definitions
! ----------------------------------------------------------------------------
#define COST_BASE 0.0_PR

#define COST_MONOPOLE 1.0_PR
#define COST_QUADRUPOLE 2.5_PR
#define COST_OCTUPOLE 3.0_PR

#define COST_HYDRO 0.01_PR
#define COST_GRAVITY 0.00001_PR

#define COST_SINK 1000.0_PR
