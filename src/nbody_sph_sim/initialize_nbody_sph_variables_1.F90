! INITIALIZE_NBODY_SPH_VARIABLES_1.F90
! D. A. Hubber - 1/10/2007
! Sets values for particular variables that need to be initialized 
! BEFORE the first force calculations in setup.  Other variables 
! (after the first force calculation) are initialized in 
! initialize_variables_2.F90
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_nbody_sph_variables_1
  use particle_module
  use hydro_module
  use periodic_module
  use time_module
  use scaling_module
  use filename_module
  use sink_module
  use neighbour_module
  use type_module
  use Nbody_module
  use timing_module
  implicit none

  integer :: p                ! Particle counter
  integer :: s                ! Star counter

  debug2("Initializing variables [initialize_nbody_sph_variables_1.F90]")


! Set accdo variables to ensure all particles are included in 
! sph and force subroutines
! ----------------------------------------------------------------------------
  do p=1,ptot
     sph(p)%accdo    = .true.
     sph(p)%nlevel   = 0
     sph(p)%nlast    = n
     sph(p)%laststep = 0.0_DP
  end do
  do s=1,stot
     star(s)%accdo    = .true.
     star(s)%nlevel   = 0
     star(s)%nlast    = n
     star(s)%laststep = 0.0_DP
#if defined(FORCE_SPLITTING)
     star(s)%accsph   = .true.
     star(s)%lastsph  = time
#endif
  end do
  timestep = 0.0_DP


! Artificial viscosity and conductivity variables
! ----------------------------------------------------------------------------
#if defined(VISC_TD)
  sph(1:ptot)%talpha    = alpha_min
  sph(1:ptot)%talpha_old = alpha_min
  sph(1:ptot)%dalpha_dt = 0.0_PR
#endif
#if defined(VISC_BALSARA)
  sph(1:ptot)%balsara = 0.0_PR
#endif
#if defined(VISC_PATTERN_REC)
  sph(1:ptot)%pattrec = 0.0_PR
#endif


! Grad-h SPH correction terms
! ----------------------------------------------------------------------------
#if defined(GRAD_H_SPH)
  do p=1,ptot
     if (sph(p)%rho == 0.0_PR) sph(p)%rho = 1.0_PR
     sph(p)%omega = 1.0_PR
#if defined(GRAVITY) && defined(SELF_GRAVITY)
     sph(p)%zo = 0.0_PR
#endif
  end do
#endif


! Zero acceleration arrays
! ----------------------------------------------------------------------------
  do p=1,ptot
     sph(p)%a(1:VDIM) = 0.0_PR
#if defined(SELF_GRAVITY) && !defined(GEOMETRIC_MAC)
     sph(p)%agravmag = 0.0_PR
#endif
#if defined(DEBUG_FORCES) && defined(GRAVITY)
     sph(p)%a_grav(1:VDIM) = 0.0_PR
#endif
#if defined(DEBUG_FORCES) && defined(HYDRO)
     sph(p)%a_hydro(1:VDIM) = 0.0_PR
#endif
  end do

  do s=1,stot
     star(s)%agravmag = 0.0_PR
     star(s)%a = 0.0_DP
     star(s)%adot = 0.0_DP
  end do


! Misc variables
! ----------------------------------------------------------------------------
  sph(1:ptot)%div_v = 0.0_PR
#if defined(DIV_A)
  sph(1:ptot)%div_a = 0.0_PR
#endif
#if defined(GRAVITY)
  sph(1:ptot)%gpot = 0.0_PR
#endif
#if defined(SINKS) && defined(SELF_GRAVITY)
  sph(1:ptot)%ispotmin = .false.
#endif
#if defined(IONIZING_UV_RADIATION)
  sph(1:ptot)%tempmin = 0.0_PR
#endif
#if defined(TIMING)
  ngravcomp = 0_ILP
  nhydrocomp = 0_ILP
  nsphcomp = 0_ILP
#endif

! Tree-building, stocking and ionization time variables
  nbuild = nsteps
  nstock = nsteps
  nionize = nsteps
  nionall = nsteps
#if defined(BINARY_TREE)
  nskeleton = nsteps
#endif

  pghost = 0

  return
END SUBROUTINE initialize_nbody_sph_variables_1
