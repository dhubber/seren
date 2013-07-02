! TURB_FORCES.F90
! A.McLeod - 26/11/2012
! Loop over particles and apply turbulent forcing
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE turb_forces
   use interface_module, only : active_particle_list
   use definitions
   use type_module
   use particle_module
   use time_module
   use turbulence_module
#if defined(OPENMP)
   use omp_lib
#endif
   implicit none
   
   integer :: acctot                   ! No. of particles on acc. step
   integer :: i                        ! Aux. particle counter
   integer :: p                        ! particle counter
   integer, allocatable :: acclist(:)  ! List of particles on acc. step
#if defined(OPENMP)
   integer :: chunksize                ! Data packet size for dynamic OpenMP
#endif
   real(kind=PR) :: turb_last_tfrac    ! Time fraction since last
   real(kind=PR) :: turb_next_tfrac    ! Time fraction until next
   real(kind=PR) :: a_turb(1:NDIM)     ! Turbulent forcing acceleration

   debug2("Calculating turbulent forcing for all particles [turb_forces.F90]")

! For multiple particle timesteps, first make a list of all SPH 
! particles on an acceleration step, and then parallelize over that list.
   allocate(acclist(1:ptot))
   call active_particle_list(acctot,acclist)
#if defined(OPENMP)
   chunksize = int(CHUNKFRAC*real(acctot,PR)/real(omp_get_max_threads(),PR)) + 1
#endif
   
! ----------------------------------------------------------------------------
   if (acctot > 0) then
      debug_timing("TURB_FORCES")
      
      turb_last_tfrac = real((time - turb_last_time) / turb_dt, PR)
      turb_next_tfrac = 1.0_PR - turb_last_tfrac
      
      afield_now = turb_last_tfrac*afield_last + turb_next_tfrac*afield_next

      ! Loop over all hydro particles (excluding boundary particles) and 
      ! calculate hydro forces.
      ! -----------------------------------------------------------------------
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p, a_turb)
      do i=1,acctot
         p = acclist(i)
         call turb_force_apply(p, a_turb)
         sph(p)%a(1:NDIM) = sph(p)%a(1:NDIM) + a_turb
      end do
!$OMP END PARALLEL DO
      ! -----------------------------------------------------------------------

  end if
! ----------------------------------------------------------------------------


   deallocate(acclist)
   
   return
  
END SUBROUTINE turb_forces