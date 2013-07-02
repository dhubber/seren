! SPH_TURBULENCE_FORCES.F90
! 
! 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_turbulence_forces
  use interface_module, only : active_particle_list,hydro,hydro_gradh
  use type_module
  use particle_module
  use timing_module
  use type_module
#if defined(OPENMP)
  use omp_lib
#endif
  implicit none

#if defined(DRIVEN_TURBULENCE)
  integer :: acctot                   ! No. of particles on acc. step
  integer :: i                        ! Aux. particle counter
  integer :: p                        ! particle counter
  integer, allocatable :: acclist(:)  ! List of particles on acc. step
#if defined(OPENMP)
  integer :: chunksize                ! Data packet size for dynamic OpenMP
#endif

  debug2("[sph_turbulence_forces.F90]")

! For multiple particle timesteps, first make a list of all hydro SPH 
! particles on an acceleration step, and then parallelize over that list.
  allocate(acclist(1:ptot))
  call active_particle_list(acctot,acclist,hydromask)
#if defined(OPENMP)
  chunksize = int(CHUNKFRAC*real(acctot,PR)/real(omp_get_max_threads(),PR)) + 1
#endif

! ----------------------------------------------------------------------------
  if (acctot > 0) then
     debug_timing("HYDRO_FORCES")

     ! Loop over all hydro particles (excluding boundary particles) and 
     ! calculate hydro forces.
     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) &
     !$OMP PRIVATE(aturb,p,rp)
     do i=1,acctot
        p = acclist(i)

        rp(1:NDIM) = sph(p)%r(1:NDIM)

        aturb(1:NDIM) = 0.0_PR

        
        sph(p)%a(1:NDIM) = sph(p)%a(1:NDIM) + aturb(1:NDIM)

     end do
     !$OMP END PARALLEL DO
     ! ------------------------------------------------------------------------

  end if
! -----------------------------------------------------------------------------

#endif

  return
END SUBROUTINE sph_turbulence_forces
