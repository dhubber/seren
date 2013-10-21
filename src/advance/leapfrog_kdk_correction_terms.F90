! LEAPFROG_KDK_CORRECTION_TERMS.F90
! D. A. Hubber - 25/07/2011
! Compute all correction terms for SPH particles using the Leapfrog-KDK
! integration scheme.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE leapfrog_kdk_correction_terms
  use interface_module, only : active_particle_list
  use particle_module
  use type_module
  use sink_module
  use hydro_module
  use time_module
  implicit none

  integer :: acctot                   ! No. of particles on acc. step
  integer :: i                        ! Aux. particle counter
  integer :: p                        ! particle counter
  integer, allocatable :: acclist(:)  ! List of particles on acc. step
#if defined(SINKS)
  integer :: s                        ! Sink counter
#endif

  debug2("Compute end-step corrections [leapfrog_kdk_correction_terms.F90]")

  ! Compute list of all active SPH particles and compute corrections
  allocate(acclist(1:ptot))
  call active_particle_list(acctot,acclist,sphmask)

  ! Calculate 'correction' step for velocities in Leapfrog-KDK scheme
  ! --------------------------------------------------------------------------
  do i=1,acctot
     p = acclist(i)
     sph(p)%v_old(1:VDIM) = sph(p)%v_old(1:VDIM) + 0.5_PR*&
          &(sph(p)%a(1:NDIM) - sph(p)%a_old(1:VDIM))*real(sph(p)%laststep,PR)
#if defined(BOUNDARY_CONDITIONS)
     call check_boundary_conditions(sph(p)%r(1:NDIM),sph(p)%v_old(1:VDIM))
#endif
     sph(p)%v(1:VDIM) = sph(p)%v_old(1:VDIM)
     sph(p)%a_old(1:VDIM) = sph(p)%a(1:VDIM)
  end do

  ! Calculate corrections for sink velocities
  ! --------------------------------------------------------------------------
#if defined(SINKS)
  if (accdo_sinks) then
     do s=1,stot
        sink(s)%vold(1:VDIM) = sink(s)%vold(1:VDIM) + &
             & 0.5_PR*(sink(s)%a(1:VDIM) - sink(s)%aold(1:VDIM))*&
             & real(laststep_sinks,PR)
#if defined(BOUNDARY_CONDITIONS)
        call check_boundary_conditions(sink(s)%r(1:NDIM),sink(s)%vold(1:VDIM))
#endif
        sink(s)%v(1:NDIM) = sink(s)%vold(1:NDIM)
        sink(s)%aold(1:NDIM) = sink(s)%a(1:NDIM)
     end do
  end if
#endif
  
  ! Compute correction terms for all active hydro quantities
  ! --------------------------------------------------------------------------
  call active_particle_list(acctot,acclist,hydromask)
  do i=1,acctot
     p = acclist(i)
#if defined(ENTROPIC_FUNCTION) && defined(ENTROPY_EQN)
     if (typeinfo(sph(p)%ptype)%eos == "entropy_eqn" .and. &
          &energy_integration == "explicit") then
        sph(p)%Aold = &
             &sph(p)%Ahalf + 0.5_PR*sph(p)%dAdt*real(sph(p)%laststep,PR)
        sph(p)%Aent = sph(p)%Aold
        call thermal(p)
     end if
#elif defined(HYDRO) && defined(ENERGY_EQN)
     if (typeinfo(sph(p)%ptype)%eos == "energy_eqn") then
        sph(p)%u_old = sph(p)%u_old + &
             & 0.5_PR*(sph(p)%dudt - sph(p)%dudt_old)*real(sph(p)%laststep,PR)
        sph(p)%u = sph(p)%u_old
        call thermal(p)
     end if
#endif
#if defined(MHD)
     sph(p)%B_old = sph(p)%B_old + &
          & 0.5_PR*(sph(p)%dBdt - sph(p)%dBdt_old)*real(sph(p)%laststep,PR)
     sph(p)%B = sph(p)%B_old
     sph(p)%phi_old = sph(p)%phi_old + &
          & 0.5_PR*(sph(p)%dphi_dt - &
          & sph(p)%dphi_dt_old)*real(sph(p)%laststep,PR)
     sph(p)%phi = sph(p)%phi_old
#endif
  end do

  deallocate(acclist)

  return
END SUBROUTINE leapfrog_kdk_correction_terms
