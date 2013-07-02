! UPDATE_THERMAL_PROPERTIES.F90
! D. A. Hubber - 17/9/2009
! Update the thermal properties of all SPH particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE update_thermal_properties
  use interface_module, only : active_particle_list,HP_walk_all_rays,thermal
  use particle_module, only : ptot
  use hydro_module
  use type_module
  implicit none

  integer :: acctot                    ! No. of particles on acc. step
  integer :: i                         ! Aux. loop counter
  integer :: p                         ! Particle counter
  integer, allocatable :: acclist(:)   ! List of active particles
  
  debug2("Calculating thermal properties of all particles [update_thermal_properties.F90]")


! For multiple particle timesteps, first make a list of all hydro SPH 
! particles on an acceleration step, and then parallelize over that list.
! ----------------------------------------------------------------------------
  debug_timing("THERMAL")
  allocate(acclist(1:ptot))
  call active_particle_list(acctot,acclist,hydromask)


! Loop over all active particles and calculate the thermal properties.
! ----------------------------------------------------------------------------
  if (acctot > 0) then     
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
        call thermal(p)
     end do
     !$OMP END PARALLEL DO
  end if

  deallocate(acclist)


  return
END SUBROUTINE update_thermal_properties
