! INITIALIZE_SPH_VARIABLES_1.F90
! D. A. Hubber - 1/10/2007
! Sets values for particular variables that need to be initialized BEFORE the 
! first force calculations in sph_setup.F90.  Other variables (after the first 
! force calculation) are initialized in initialize_sph_variables_2.F90.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_sph_variables_1
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

  integer :: p             ! Particle counter

  debug2("Initializing variables [initialize_sph_variables_1.F90]")


! Set accdo variables to ensure all particles are included in 
! sph and force subroutines
! ----------------------------------------------------------------------------
  do p=1,ptot
     sph(p)%accdo    = .true.
     sph(p)%nlevel   = 0
     sph(p)%nlast    = n
     sph(p)%laststep = 0.0_DP
  end do
  timestep = 0.0_DP
#if defined(SINKS)
  accdo_sinks    = .true.
  nlevel_sinks   = 0
  laststep_sinks = 0.0_DP
#endif


! Artificial viscosity and conductivity variables
! ----------------------------------------------------------------------------
#if defined(VISC_TD)
  sph(1:ptot)%talpha     = alpha_min
  sph(1:ptot)%talpha_old = alpha_min
  sph(1:ptot)%dalpha_dt  = 0.0_PR
#endif
#if defined(VISC_BALSARA)
  sph(1:ptot)%balsara = 0.0_PR
#endif


! Grad-h SPH correction terms
! ----------------------------------------------------------------------------
#if defined(GRAD_H_SPH)
  do p=1,ptot
     if (sph(p)%rho == 0.0_PR) sph(p)%rho = 1.0_PR
     sph(p)%omega = 1.0_PR
#if defined(SELF_GRAVITY)
     sph(p)%zo = 0.0_PR
#endif
  end do
#endif


! Zero acceleration arrays
! ----------------------------------------------------------------------------
  do p=1,ptot
     sph(p)%a(1:VDIM) = 0.0_PR
#if defined(STELLAR_WIND)
     sph(p)%a_wind(1:NDIM) = 0.0_PR
#endif
#if defined(DEBUG_FORCES) && defined(GRAVITY)
     sph(p)%a_grav(1:VDIM) = 0.0_PR
#endif
#if defined(DEBUG_FORCES) && defined(HYDRO)
     sph(p)%a_hydro(1:VDIM) = 0.0_PR
#endif
#if defined(DEBUG_FORCES) && defined(ARTIFICIAL_VISCOSITY)
     sph(p)%a_visc(1:VDIM) = 0.0_PR
#endif
  end do
#if defined(SELF_GRAVITY) && !defined(GEOMETRIC_MAC)
  sph(1:ptot)%agravmag = 0.0_PR
#if defined(SINKS)
  sink(1:stot)%agravmag = 0.0_PR
#endif
#endif


! MHD arrays
! ----------------------------------------------------------------------------
#if defined(MHD)
  do p=1,ptot
     sph(p)%B_old = sph(p)%B
     sph(p)%dBdt = 0.0_PR
     sph(p)%phi = 0.0_PR
     sph(p)%dphi_dt = 0.0_PR
     sph(p)%phi_old = 0.0_PR
  end do
#endif


! Misc variables
! ----------------------------------------------------------------------------
  do p=1,ptot
     sph(p)%div_v = 0.0_PR
#if defined(GRAVITY)
     sph(p)%gpot = 0.0_PR
#if defined(RAD_WS)
     sph(p)%sphgpot = 0.0_PR
#endif
#endif
#if defined(DIV_A)
     sph(p)%div_a = 0.0_PR
#endif
#if defined(OSPH) && defined(ENTROPIC_FUNCTION) && !defined(ENTROPY_EQN)
     sph(p)%Aent = 1.0_PR
#endif
#if defined(SINKS) && defined(SELF_GRAVITY)
     sph(p)%ispotmin = .false.
#endif
#if defined(IONIZING_UV_RADIATION)
     sph(p)%tempmin = 0.0_PR
#endif
#if defined(STELLAR_WIND)
     sph(p)%windflag = .false.
#endif
  end do

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

  pghost = 0
  pperiodic = 0

  return
END SUBROUTINE initialize_sph_variables_1
