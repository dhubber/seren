! SPH_HYDRO_FORCES.F90
! C. P. Batty & D. A. Hubber - 11/1/2007
! Calculates hydrodynamical accelerations for all particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_hydro_forces
  use interface_module, only : active_particle_list,hydro,hydro_gradh
  use type_module
  use particle_module
  use time_module
  use timing_module
#if defined(OPENMP)
  use omp_lib
#endif
#if defined(CELL_WALK)
  use tree_module
#endif
#if defined(USE_MPI)
  use mpi_communication_module
#endif
  implicit none

#if defined(HYDRO)
  integer :: acctot                   ! No. of particles on acc. step
  integer :: i                        ! Aux. particle counter
  integer :: p                        ! particle counter
  integer, allocatable :: acclist(:)  ! List of particles on acc. step
#if defined(OPENMP)
  integer :: chunksize                ! Data packet size for dynamic OpenMP
#endif
#if defined(DEBUG_VISC_BALSARA) || defined(DEBUG_VISC_PATTERN_REC)
  real(kind=DP) :: mean               ! Mean Balsara/pattern recog. factor
#endif

  debug2("Calculating hydro forces for all particles [sph_hydro_forces.F90]")

! For multiple particle timesteps, first make a list of all hydro SPH 
! particles on an acceleration step, and then parallelize over that list.
  allocate(acclist(1:ptot))
  call active_particle_list(acctot,acclist,hydromask)
#if defined(OPENMP)
  chunksize = int(CHUNKFRAC*real(acctot,PR)/real(omp_get_max_threads(),PR)) + 1
#endif
#if defined(TIMING)
  nhydrocomp = nhydrocomp + int(acctot,ILP)
#endif

! ----------------------------------------------------------------------------
  if (acctot > 0) then
     debug_timing("HYDRO_FORCES")

     ! Loop over all hydro particles (excluding boundary particles) and 
     ! calculate hydro forces.
     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
#if defined(USE_MPI) && !defined(LOW_MEM)
        sph(p)%hydro_calcs = 0.0_PR
#endif
#if defined(GRAD_H_SPH)
        call hydro_gradh(p)
#elif defined(OSPH)
        call hydro_osph(p)
#else
        call hydro(p)
#endif
#if defined(MHD)
        call mhd(p)
#endif
     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------

  end if
! ----------------------------------------------------------------------------

! Calculate and write out mean pattern recognition factor
#if defined(VISC_PATTERN_REC) && defined(DEBUG_VISC_PATTERN_REC)
  mean = 0.0_PR
  do p=1,ptot
     mean = mean + sph(p)%pattrec
  end do
  open(1, file="vg.dat", position="append")
  write(1,*) nsteps,mean/real(ptot,PR)
  close(1)
#endif

! Calculate and write out mean Balsara factor
#if defined(VISC_BALSARA) && defined(DEBUG_VISC_BALSARA)
  mean = 0.0_PR
  do p=1,ptot
     mean = sph(p)%balsara + mean
  end do
  open(1, file="bal.dat", position="append")
  write(1,*) nsteps,mean/real(ptot,PR)
  close(1)
#endif

  deallocate(acclist)

#endif

! ----------------------------------------------------------------------------
  ! Add costs for MPI loadbalancing
#if defined(USE_MPI) && defined(HYDRO) && !defined(LOW_MEM)
  sum_costs = sum_costs + COST_HYDRO * &
     & sum(sph(phydrostart:ptot)%hydro_calcs,mask=sph(phydrostart:ptot)%accdo)
#endif

  return
END SUBROUTINE sph_hydro_forces
