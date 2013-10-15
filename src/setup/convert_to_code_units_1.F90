! CONVERT_TO_CODE_UNITS_1.F90
! C. P. Batty & D. A. Hubber - 11/12/2006
! Converts all particle and sink variables from physical units to dimensionless 
! code units using scaling factors calculated in units.F90.  Variables are 
! converted to dimensionless form by dividing by Xscale 
! (where X is r,t,m etc..). 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE convert_to_code_units_1(minimal)
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
  implicit none

  logical, intent(in) :: minimal

  integer :: p        ! counter to loop over particles
#if defined(SINKS)
  integer :: s        ! sink counter
#endif

  debug1("Converting to dimensionless code units [convert_to_code_units_1.F90]")

! Scale main particle data to code units 
! ----------------------------------------------------------------------------
  if (minimal) then
     do p=1,ptot
        minimal_sph(p)%r(1:NDIM) = minimal_sph(p)%r(1:NDIM) / real(rscale,PR)
        minimal_sph(p)%m         = minimal_sph(p)%m / real(mscale,PR)
        minimal_sph(p)%h         = minimal_sph(p)%h / real(rscale,PR)
        minimal_sph(p)%v(1:VDIM) = minimal_sph(p)%v(1:VDIM) / real(vscale,PR)
#if defined(MHD)
        minimal_sph(p)%B(1:BDIM) = minimal_sph(p)%B(1:BDIM) / real(Bscale,PR)
#endif
        minimal_sph(p)%rho       = minimal_sph(p)%rho / real(rhoscale,PR)
#if defined(HYDRO)
        minimal_sph(p)%u         = minimal_sph(p)%u / real(uscale,PR)
#endif
     end do
  else
     do p=1,ptot
        sph(p)%r(1:NDIM) = sph(p)%r(1:NDIM) / real(rscale,PR)
        sph(p)%m         = sph(p)%m / real(mscale,PR)
        sph(p)%h         = sph(p)%h / real(rscale,PR)
        sph(p)%v(1:VDIM) = sph(p)%v(1:VDIM) / real(vscale,PR)
#if defined(MHD)
        sph(p)%B(1:BDIM) = sph(p)%B(1:BDIM) / real(Bscale,PR)
#endif
        sph(p)%rho       = sph(p)%rho / real(rhoscale,PR)
#if defined(HYDRO)
        sph(p)%u         = sph(p)%u / real(uscale,PR)
#endif
     end do
  end if

! Scale sink variables
! ----------------------------------------------------------------------------
#if defined(SINKS)
  do s=1,stot
     sink(s)%tcreate     = sink(s)%tcreate / tscale
     sink(s)%r(1:NDIM)   = sink(s)%r(1:NDIM) / real(rscale,PR)
     sink(s)%v(1:VDIM)   = sink(s)%v(1:VDIM) / real(vscale,PR)
     sink(s)%m           = sink(s)%m / real(mscale,PR)
     sink(s)%h           = sink(s)%h / real(rscale,PR)
     sink(s)%radius      = sink(s)%radius / real(rscale,PR)
     sink(s)%angmom(1:3) = sink(s)%angmom(1:3) / angmomscale
     sink(s)%dmdt        = sink(s)%dmdt / dmdtscale
     sink(s)%star_radius = sink(s)%star_radius / rscale
     sink(s)%luminosity  = sink(s)%luminosity / Lscale
     sink(s)%mmax        = sink(s)%mmax / mscale
     sink(s)%macc(1:DMDT_RANGE) = sink(s)%macc(1:DMDT_RANGE) / mscale
     sink(s)%tacc(1:DMDT_RANGE) = sink(s)%tacc(1:DMDT_RANGE) / tscale
  end do
#endif

  return
END SUBROUTINE convert_to_code_units_1
