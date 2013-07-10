! WRITE_MAKEFILE_OPTIONS
! D. A. Hubber - 24/06/2010
! Write all Makefile options to the selected output (e.g. 6 for screen).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_makefile_options(unitno)
  use definitions
  implicit none

  integer, intent(in) :: unitno        ! Unit no. to write information to

! Compiler flag output 
! ----------------------------------------------------------------------------
  write(unitno,'(A)') "======================="
  write(unitno,'(A)') "MAKEFILE COMPILER FLAGS"
  write(unitno,'(A)') "======================="

#if defined(COMPILER_GFORTRAN)
  write(unitno,'(A)') "F90                         = gfortran"
#elif defined(COMPILER_IFORT)
  write(unitno,'(A)') "F90                         = ifort"
#elif defined(COMPILER_F95)
  write(unitno,'(A)') "F90                         = f95"
#elif defined(COMPILER_G95)
  write(unitno,'(A)') "F90                         = g95"
#elif defined(COMPILER_PGF90)
  write(unitno,'(A)') "F90                         = pgf90"
#elif defined(COMPILER_PGF95)
  write(unitno,'(A)') "F90                         = pgf95"
#elif defined(COMPILER_MPIF90)
  write(unitno,'(A)') "F90                         = mpif90"
#elif defined(COMPILER_SUNF95)
  write(unitno,'(A)') "F90                         = sunf95"
#endif

  write(unitno,'(A)') "VERSION_NO                  = " // SEREN_VERSION

#if defined(OPENMP)
  write(unitno,'(A)') "OPENMP                      = 1"
#else
  write(unitno,'(A)') "OPENMP                      = 0"
#endif

#if defined(COMPILER_MODE_STANDARD)
  write(unitno,'(A,A)') "COMPILER_MODE               = STANDARD"
#elif defined(COMPILER_MODE_FAST)
  write(unitno,'(A,A)') "COMPILER_MODE               = FAST"
#elif defined(COMPILER_MODE_DEBUG)
  write(unitno,'(A,A)') "COMPILER_MODE               = DEBUG"
#elif defined(COMPILER_MODE_FASTDEBUG)
  write(unitno,'(A,A)') "COMPILER_MODE               = FASTDEBUG"
#elif defined(COMPILER_MODE_0)
  write(unitno,'(A,A)') "COMPILER_MODE               = 0"
#endif

  write(unitno,'(A,I1)') "OUTPUT_LEVEL                = ",OUTPUT_LEVEL

  write(unitno,'(A,I1)') "DIAGNOSTIC_OUTPUT           = ",DIAGNOSTIC_OUTPUT 

  write(unitno,'(A,I1)') "NDIM                        = ",NDIM

#if defined(QUADRUPLE_PRECISION)
  write(unitno,'(A)') "PRECISION                   = QUADRUPLE"
#elif defined(DOUBLE_PRECISION)
  write(unitno,'(A)') "PRECISION                   = DOUBLE"
#else
  write(unitno,'(A)') "PRECISION                   = SINGLE"
#endif

#if defined(DRAGON_INPUT) && defined(SEREN_INPUT) && defined(ASCII_INPUT)
  write(unitno,'(A)') "INFILE_FORMAT               = ALL"
#elif defined(DRAGON_INPUT)
  write(unitno,'(A)') "INFILE_FORMAT               = DRAGON"
#elif defined(SEREN_INPUT)
  write(unitno,'(A)') "INFILE_FORMAT               = SEREN"
#elif defined(ASCII_INPUT)
  write(unitno,'(A)') "INFILE_FORMAT               = ASCII"
#endif

#if defined(DRAGON_OUTPUT) && defined(SEREN_OUTPUT)
  write(unitno,'(A)') "OUTFILE_FORMAT              = ALL"
#elif defined(DRAGON_OUTPUT)
  write(unitno,'(A)') "OUTFILE_FORMAT              = DRAGON"
#elif defined(SEREN_OUTPUT)
  write(unitno,'(A)') "OUTFILE_FORMAT              = SEREN"
#endif

#if defined(PERIODIC) && defined(BOUNDARY_CONDITIONS)
  write(unitno,'(A)') "PERIODIC                    = 1"
#else
  write(unitno,'(A)') "PERIODIC                    = 0"
#endif
#if defined(PERIODIC_X)
  write(unitno,'(A)') "X_BOUNDARY                  = PERIODIC"
#else
  write(unitno,'(A)') "X_BOUNDARY                  = 0"
#endif
#if defined(PERIODIC_Y)
  write(unitno,'(A)') "Y_BOUNDARY                  = PERIODIC"
#else
  write(unitno,'(A)') "Y_BOUNDARY                  = 0"
#endif
#if defined(PERIODIC_Z)
  write(unitno,'(A)') "Z_BOUNDARY                  = PERIODIC"
#else
  write(unitno,'(A)') "Z_BOUNDARY                  = 0"
#endif

#if defined(SPH_SIMULATION)
  write(unitno,'(A)') "SPH_SIMULATION              = 1"
#else
  write(unitno,'(A)') "SPH_SIMULATION              = 0"
#endif
#if defined(NBODY_SIMULATION)
  write(unitno,'(A)') "NBODY_SIMULATION            = 1"
#else
  write(unitno,'(A)') "NBODY_SIMULATION            = 0"
#endif

#if defined(GRAD_H_SPH)
  write(unitno,'(A)') "SPH                         = GRAD_H_SPH"
#else
  write(unitno,'(A)') "SPH                         = STANDARD"
#endif

#if defined(EULER)
  write(unitno,'(A)') "SPH_INTEGRATION             = EULER"
#elif defined(RUNGE_KUTTA2)
  write(unitno,'(A)') "SPH_INTEGRATION             = RK2"
#elif defined(LEAPFROG_KDK)
  write(unitno,'(A)') "SPH_INTEGRATION             = LFKDK"
#elif defined(LEAPFROG_DKD)
  write(unitno,'(A)') "SPH_INTEGRATION             = LFDKD"
#endif

#if defined(M4_KERNEL) && defined(TC_KERNEL)
  write(unitno,'(A)') "KERNEL                      = M4TC"
#elif defined(M4_KERNEL)
  write(unitno,'(A)') "KERNEL                      = M4"
#elif defined(QUINTIC_KERNEL) && defined(TC_KERNEL)
  write(unitno,'(A)') "KERNEL                      = QUINTICTC"
#elif defined(QUINTIC_KERNEL)
  write(unitno,'(A)') "KERNEL                      = QUINTIC"
#elif defined(GAUSSIAN_3H_KERNEL)
  write(unitno,'(A)') "KERNEL                      = QUINTIC"
#endif

#if defined(HGATHER)
  write(unitno,'(A)') "HFIND                       = NUMBER"
#elif defined(HMASS)
  write(unitno,'(A)') "HFIND                       = MASS"
#elif defined(CONSTANT_H)
  write(unitno,'(A)') "HFIND                       = CONSTANT"
#endif

#if defined(MINIMUM_H)
  write(unitno,'(A)') "MINIMUM_H                   = 1"
#else
  write(unitno,'(A)') "MINIMUM_H                   = 0"
#endif

#if defined(HYDRO)
  write(unitno,'(A)') "HYDRO                       = 1"
#else
  write(unitno,'(A)') "HYDRO                       = 0"
#endif

#if defined(ENERGY_EQN)
  write(unitno,'(A)') "ENERGY_EQN                  = 1"
#else
  write(unitno,'(A)') "ENERGY_EQN                  = 0"
#endif

#if defined(ENTROPY_EQN)
  write(unitno,'(A)') "ENTROPY_EQN                 = 1"
#else
  write(unitno,'(A)') "ENTROPY_EQN                 = 0"
#endif

#if defined(RAD_WS)
  write(unitno,'(A)') "RAD_WS                      = 1"
#else
  write(unitno,'(A)') "RAD_WS                      = 0"
#endif

#if defined(RAD_WS) && defined(RAD_WS_SINK_POT)
  write(unitno,'(A)') "SINK_POTENTIAL_WS           = 1"
#else
  write(unitno,'(A)') "SINK_POTENTIAL_WS           = 0"
#endif

#if defined(RAD_WS) && defined(AMBIENT_HEATING) && defined(CONST_HEATING)
  write(unitno,'(A)') "AMBIENT_HEATING_WS          = 1"
#else
  write(unitno,'(A)') "AMBIENT_HEATING_WS          = 0"
#endif

#if defined(RAD_WS) && defined(HDISC_HEATING)
  write(unitno,'(A)') "SINK_HEATING_WS             = HDISC_HEATING"
#elif defined(RAD_WS) && defined(HDISC_HEATING_3D_SINGLE)
  write(unitno,'(A)') "SINK_HEATING_WS             = HDISC_HEATING_3D_SINGLE"
#elif defined(RAD_WS) && defined(STAR_HEATING)
  write(unitno,'(A)') "SINK_HEATING_WS             = STAR_HEATING"
#elif defined(RAD_WS) && defined(STAR_SIMPLE_HEATING)
  write(unitno,'(A)') "SINK_HEATING_WS             = STAR_SIMPLE_HEATING"
#else
  write(unitno,'(A)') "SINK_HEATING_WS             = 0"
#endif

#if defined(DIFFUSION)
  write(unitno,'(A)') "FLUX_LIMITED_DIFFUSION      = 1"
#else
  write(unitno,'(A)') "FLUX_LIMITED_DIFFUSION      = 0"
#endif

#if defined(IONIZING_UV_RADIATION) && defined(SINGLE_STATIC_SOURCE)
  write(unitno,'(A)') "IONIZING_UV_RADIATION       = SINGLE_STATIC_SINGLE"
#elif defined(IONIZING_UV_RADIATION) && defined(MULTIPLE_SINK_SOURCE)
  write(unitno,'(A)') "IONIZING_UV_RADIATION       = MULTIPLE_SINK_SOURCES"
#else
  write(unitno,'(A)') "IONIZING_UV_RADIATION       = 0"
#endif

#if defined(RIEMANN_SOLVER)
  write(unitno,'(A)') "RIEMANN_SOLVER              = 1"
#else
  write(unitno,'(A)') "RIEMANN_SOLVER              = 0"
#endif

#if defined(ARTIFICIAL_VISCOSITY) && defined(VISC_AB)
  write(unitno,'(A)') "ARTIFICIAL_VISCOSITY        = AB"
#elif defined(ARTIFICIAL_VISCOSITY) && defined(VISC_MON97)
  write(unitno,'(A)') "ARTIFICIAL_VISCOSITY        = MON97"
#else
  write(unitno,'(A)') "ARTIFICIAL_VISCOSITY        = 0"
#endif

#if defined(ARTIFICIAL_VISCOSITY) && defined(VISC_TD)
  write(unitno,'(A)') "VISC_TD                     = 1"
#else
  write(unitno,'(A)') "VISC_TD                     = 0"
#endif

#if defined(ARTIFICIAL_VISCOSITY) && defined(VISC_BALSARA)
  write(unitno,'(A)') "BALSARA                     = 1"
#else
  write(unitno,'(A)') "BALSARA                     = 0"
#endif

#if defined(ARTIFICIAL_CONDUCTIVITY) && defined(COND_PRICE2008)
  write(unitno,'(A)') "ARTIFICIAL_CONDUCTIVITY     = PRICE2008"
#elif defined(ARTIFICIAL_CONDUCTIVITY) && defined(COND_WADSLEY2008)
  write(unitno,'(A)') "ARTIFICIAL_CONDUCTIVITY     = WADSLEY2008"
#else
  write(unitno,'(A)') "ARTIFICIAL_CONDUCTIVITY     = 0"
#endif

#if defined(EXTERNAL_PRESSURE)
  write(unitno,'(A)') "EXTERNAL_PRESSURE           = 1"
#else
  write(unitno,'(A)') "EXTERNAL_PRESSURE           = 0"
#endif

#if defined(EXTERNAL_FORCE) && defined(PLUMMER_POTENTIAL)
  write(unitno,'(A)') "EXTERNAL_FORCE              = PLUMMER"
#elif defined(EXTERNAL_FORCE) && defined(UDS_POTENTIAL)
  write(unitno,'(A)') "EXTERNAL_FORCE              = UDS"
#elif defined(EXTERNAL_FORCE) && defined(NFW1996_POTENTIAL)
  write(unitno,'(A)') "EXTERNAL_FORCE              = NFW1996"
#else
  write(unitno,'(A)') "EXTERNAL_FORCE              = 0"
#endif

#if defined(SELF_GRAVITY) && defined(N_BODY)
  write(unitno,'(A)') "SELF_GRAVITY                = NBODY"
#elif defined(SELF_GRAVITY)
  write(unitno,'(A)') "SELF_GRAVITY                = KS"
#else
  write(unitno,'(A)') "SELF_GRAVITY                = 0"
#endif

#if defined(MEANH_GRAVITY)
  write(unitno,'(A)') "MEANH_GRAVITY               = 1"
#else
  write(unitno,'(A)') "MEANH_GRAVITY               = 0"
#endif

#if defined(GRAVITY) && defined(EWALD)
  write(unitno,'(A)') "EWALD                       = 1"
#else
  write(unitno,'(A)') "EWALD                       = 0"
#endif

#if defined(GRAVITY) && defined(SINKS) && defined(SMOOTH_ACCRETION)
  write(unitno,'(A)') "SINKS                       = SMOOTH_ACC"
#elif defined(GRAVITY) && defined(SINKS) && defined(NO_ACCRETION)
  write(unitno,'(A)') "SINKS                       = NO_ACC"
#elif defined(GRAVITY) && defined(SINKS)
  write(unitno,'(A)') "SINKS                       = SIMPLE"
#else
  write(unitno,'(A)') "SINKS                       = 0"
#endif

#if defined(GRAVITY) && defined(SINKS) && defined(FIXED_ABSOLUTE_SINKRAD)
  write(unitno,'(A)') "SINK_RADIUS                 = FIXED_ABSOLUTE"
#elif defined(GRAVITY) && defined(SINKS) && defined(FIXED_HMULT_SINKRAD)
  write(unitno,'(A)') "SINK_RADIUS                 = FIXED_HMULT"
#elif defined(GRAVITY) && defined(SINKS) && defined(HMULT_SINKRAD)
  write(unitno,'(A)') "SINK_RADIUS                 = HMULT"
#endif

#if defined(GRAVITY) && defined(SINKS) && defined(SINK_REMOVE_ANGMOM)
  write(unitno,'(A)') "SINK_REMOVE_ANGMOM          = 1"
#else
  write(unitno,'(A)') "SINK_REMOVE_ANGMOM          = 0"
#endif

#if defined(SINKS) && defined(SINK_GRAVITY_ONLY)
  write(unitno,'(A)') "SINK_GRAVITY_ONLY           = 1"
#else
  write(unitno,'(A)') "SINK_GRAVITY_ONLY           = 0"
#endif

#if defined(NBODY_SIMULATION) || defined(NBODY_SPH_SIMULATION)
  write(unitno,'(A)') "NBODY_INTEGRATION           = HERMITE4"
#else
  write(unitno,'(A)') "NBODY_INTEGRATION           = 0"
#endif

#if defined(BINARY_STATS)
  write(unitno,'(A)') "BINARY_STATS                = 1"
#else
  write(unitno,'(A)') "BINARY_STATS                = 0"
#endif

#if defined(BH_TREE)
  write(unitno,'(A)') "TREE                        = BH"
#else
  write(unitno,'(A)') "TREE                        = 0"
#endif

#if defined(SELF_GRAVITY) && defined(BH_TREE) && defined(QUADRUPOLE)
  write(unitno,'(A)') "MULTIPOLE                   = QUADRUPOLE"
#elif defined(SELF_GRAVITY) && defined(BH_TREE) && defined(OCTUPOLE)
  write(unitno,'(A)') "MULTIPOLE                   = OCTUPOLE"
#else
  write(unitno,'(A)') "MULTIPOLE                   = 0"
#endif

#if defined(SELF_GRAVITY) && defined(BH_TREE) && defined(GEOMETRIC_MAC)
  write(unitno,'(A)') "MAC                         = GEOMETRIC"
#elif defined(SELF_GRAVITY) && defined(BH_TREE) && defined(GADGET_MAC)
  write(unitno,'(A)') "MAC                         = GADGET"
#elif defined(SELF_GRAVITY) && defined(BH_TREE) && defined(GADGET2_MAC)
  write(unitno,'(A)') "MAC                         = GADGET2"
#elif defined(SELF_GRAVITY) && defined(BH_TREE) && defined(EIGEN_MAC)
  write(unitno,'(A)') "MAC                         = EIGEN"
#else
  write(unitno,'(A)') "MAC                         = 0"
#endif

#if defined(BH_TREE) && defined(REORDER_PARTICLES)
  write(unitno,'(A)') "REORDER                     = PARTICLES"
#else
  write(unitno,'(A)') "REORDER                     = 0"
#endif

#if defined(HEAPSORT)
  write(unitno,'(A)') "SORT                        = HEAP"
#else
  write(unitno,'(A)') "SORT                        = INSERTION"
#endif

#if defined(ADAPTIVE_TIMESTEP_LEVELS)
  write(unitno,'(A)') "TIMESTEP                    = ADAPTIVE"
#elif defined(RESTRICTED_TIMESTEP_LEVELS)
  write(unitno,'(A)') "TIMESTEP                    = RESTRICTED"
#elif defined(FIXED_TIMESTEP_LEVELS)
  write(unitno,'(A)') "TIMESTEP                    = FIXED"
#endif

#if defined(CHECK_NEIGHBOUR_TIMESTEPS) && defined(IMMEDIATE_TIMESTEP_REDUCTION)
  write(unitno,'(A)') "CHECK_NEIB_TIMESTEPS        = 2"
#elif defined(CHECK_NEIGHBOUR_TIMESTEPS)
  write(unitno,'(A)') "CHECK_NEIB_TIMESTEPS        = 1"
#else
  write(unitno,'(A)') "CHECK_NEIB_TIMESTEPS        = 0"
#endif

#if defined(NEIGHBOUR_LISTS)
  write(unitno,'(A)') "NEIGHBOURLISTS              = 1"
#else
  write(unitno,'(A)') "NEIGHBOURLISTS              = 0"
#endif

#if defined(KERNEL_TABLES)
  write(unitno,'(A)') "KERNEL_TABLES               = 1"
#else
  write(unitno,'(A)') "KERNEL_TABLES               = 0"
#endif

#if defined(REMOVE_OUTLIERS)
  write(unitno,'(A)') "REMOVE_OUTLIERS             = 1"
#else
  write(unitno,'(A)') "REMOVE_OUTLIERS             = 0"
#endif

#if defined(TURBULENT_FORCING)
  write(unitno,'(A)') "TURBULENT_FORCING           = 1"
#else
  write(unitno,'(A)') "TURBULENT_FORCING           = 0"
#endif

#if defined(TIMING)
  write(unitno,'(A)') "TIMING_CODE                 = 1"
#else
  write(unitno,'(A)') "TIMING_CODE                 = 0"
#endif

#if defined(SPIEGEL_TEST)
  write(unitno,'(A)') "TEST                        = SPIEGEL"
#elif defined(FREEFALL_TEST)
  write(unitno,'(A)') "TEST                        = FREEFALL"
#elif defined(BINARY_TEST)
  write(unitno,'(A)') "TEST                        = BINARY"
#elif defined(PLUMMER_TEST)
  write(unitno,'(A)') "TEST                        = PLUMMER"
#else
  write(unitno,'(A)') "TEST                        = 0"
#endif


  return
END SUBROUTINE write_makefile_options
