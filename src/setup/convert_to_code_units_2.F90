! CONVERT_TO_CODE_UNITS_2.F90
! C. P. Batty & D. A. Hubber - 11/12/2006
! Converts all simulation variables from physical units to dimensionless code 
! units using scaling factors calculated in units.F90.  Variables are 
! converted to dimensionless form by dividing by Xscale 
! (where X is r,t,m etc..). 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE convert_to_code_units_2
  use constant_module
  use scaling_module
  use particle_module
  use periodic_module
  use sink_module
  use hydro_module
  use time_module
  use type_module
  use neighbour_module
  use Nbody_module
  use HP_module
  use turbulence_module
  implicit none

  integer :: p        ! counter to loop over particles
#if defined(SINKS)
  integer :: s        ! sink counter
#endif

  debug1("Converting to dimensionless code units [convert_to_code_units_2.F90]")

! Scale periodic variables
  periodic_min(1:3)  = periodic_min(1:3) / real(rscale,PR)
  periodic_max(1:3)  = periodic_max(1:3) / real(rscale,PR)
  periodic_half(1:3) = periodic_half(1:3) / real(rscale,PR)
  periodic_size(1:3) = periodic_size(1:3) / real(rscale,PR)

! Scale SPH simulation time variables
  dt_fixed    = dt_fixed  / tscale
  firstsnap   = firstsnap / tscale
  lastsnap    = lastsnap / tscale
  snaptime    = snaptime  / tscale
  sph_endtime = sph_endtime / tscale
  time        = time / tscale

! Scale N-body/SPH simulation time variables
  nbody_sph_endtime = nbody_sph_endtime / tscale

! Scale N-body simulation time variables
  nbody_endtime  = nbody_endtime / tscale

! Scale density variables 
  rhobary = rhobary / real(rhoscale*rhocgs,PR)
  rhosink = rhosink / real(rhoscale*rhocgs,PR)

! Other misc. variables
  Pext       = Pext / real(Pscale,PR)
  mgas_orig  = mgas_orig / mscale
  hmin       = hmin / real(rscale,PR)
  rspheremax = rspheremax / real(rscale,PR)
  rholost    = rholost / real(rhoscale*rhocgs,PR)
  rad_lost   = rad_lost / real(rscale,PR)

! HEALPix variables
  M_loss = M_loss / real(dmdtscale,PR)
  v_wind = v_wind / real(vscale,PR)
  
! Turbulent forcing parameters
#if defined(TURBULENT_FORCING)
  turb_T = turb_T / tscale
  turb_rms = turb_rms / ascale
#endif

! If using fixed absolute sink radius, scale variable from au to code units
#if defined(FIXED_ABSOLUTE_SINKRAD) && !defined(DIMENSIONLESS)
  sinkrad = sinkrad * real(r_au/(rscale*r_SI),PR)
#endif

  return
END SUBROUTINE convert_to_code_units_2
