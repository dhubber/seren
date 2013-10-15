! CREATE_NEW_SPH_PARTICLE.F90
! D. A. Hubber - 12/10/2011
! Creates new SPH particle in memory
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE create_new_sph_particle(ptype,pnew,part)
  use particle_module
  use type_module
  use time_module
  implicit none

  integer, intent(in) :: ptype              ! New particle type
  integer, intent(in) :: pnew               ! New particle array id
  type(sph_particle), intent(in) :: part    ! New particle data

! Stop simulation if we have run out of memory
  if (ptot == pmax) &
       &stop 'Run out of memory!  Cannot add new SPH particle to arrays!'

! Check the particle has a valid type
  if (ptype <= 0 .and. ptype > nmaxtypes) &
       &stop 'Invalid particle type when adding new particle'


! If we have enough memory, add particle to end of sph arrays
  sph(pnew) = part
  sph(pnew)%accdo     = .true.
  sph(pnew)%porig     = pnew
  sph(pnew)%ptype     = ptype
  sph(pnew)%nlast     = n
  sph(pnew)%nlevel    = level_max
#if defined(CHECK_NEIGHBOUR_STEPS)
  sph(pnew)%nminneib  = level_max
#endif
  sph(pnew)%laststep  = 0.0_DP
  sph(pnew)%a(1:VDIM) = 0.0_DP
#if defined(MHD)
  sph(pnew)%B(1:BDIM) = 0.0_DP
#endif
  sph(pnew)%rho       = part%m*(1.2_PR/part%h)**(NDIM)
#if !defined(LOW_MEM)
  sph(pnew)%invrho    = 1.0_PR/sph(pnew)%rho
#endif
  sph(pnew)%h         = part%h
#if !defined(LOW_MEM)
  sph(pnew)%invh      = 1.0_PR/part%h
#endif
  sph(pnew)%div_v     = 0.0_PR
#if defined(GRAVITY)
  sph(pnew)%gpot      = 0.0_PR
#endif
#if defined(GRAD_H_SPH)
  sph(pnew)%omega     = 1.0_PR
#endif
#if defined(GRAD_H_SPH) && defined(SELF_GRAVITY)
  sph(pnew)%zo        = 0.0_PR
#endif
  sph(pnew)%r_old(1:NDIM) = part%r(1:NDIM)
  sph(pnew)%v_old(1:VDIM) = part%v(1:NDIM)
#if defined(LEAPFROG_KDK)
  sph(pnew)%a_old(1:VDIM) = 0.0_PR
#endif
#if defined(HYDRO)
  sph(pnew)%temp          = part%temp
#endif
#if defined(HYDRO) && defined(INTERNAL_ENERGY)
  sph(pnew)%u_old         = part%u
#endif


  return
END SUBROUTINE create_new_sph_particle
