! AMBIENTTEMP.F90
! D. Stamatellos - 3/1/2008
! Gets ambient (i.e. heating) temperature of particle. 
! TODO : Scaling should be done 'properly' in this rouinte (Check with Dimitri)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE ambient_temp(p,atemp)
  use particle_module
  use sink_module
  use Tprof_module
  use scaling_module
  implicit none

  integer, intent(in) :: p             ! particle id
  real(kind=PR), intent(out) :: atemp  ! Ambient temperature

#if defined(HDISC_HEATING) || defined(HDISC_HEATING_3D_SINGLE) || defined(STAR_SIMPLE_HEATING)
  integer :: s                         ! sink id
  real(kind=PR) :: drsqd               ! radius squared
  real(kind=PR) :: dr(1:NDIM)          ! relative displacement vector
#endif


! Assumes one star/disc with disc in the x-y plane 
! ----------------------------------------------------------------------------
#if defined(HDISC_HEATING)  
  do s=1,stot
     if (sink(s)%m > 0.3) then
        dr(1:NDIM) = sink(s)%r(1:NDIM) - sph(p)%r(1:NDIM)
        drsqd= (dr(1)*dr(1) + dr(2)*dr(2))*(6.684e-14*rscale*rcgs)**2 
        ! x-y distance in AU
        atemp = sqrt(ptemp0*ptemp0*&
             &(drsqd + ptemp_r0*ptemp_r0)**(-ptemp_q) + &
             &temp_inf*temp_inf)
     end if
  end do


! Assumes one star/disc with disc in the x-y plane 
! ----------------------------------------------------------------------------
#elif defined(HDISC_HEATING_3D_SINGLE)
  s = 1
  dr(1:NDIM) = sink(s)%r(1:NDIM) - sph(p)%r(1:NDIM)
  drsqd= (dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))*&
       &(6.684e-14*rscale*rcgs)**2 
  atemp = sqrt(ptemp0*ptemp0*&
       &(drsqd + ptemp_r0*ptemp_r0)**(-ptemp_q) + temp_inf*temp_inf)


! ..
! ----------------------------------------------------------------------------
#elif defined(STAR_SIMPLE_HEATING)
  atemp = 0.0_PR
  do s=1,stot       
     dr(1:NDIM) = sink(s)%r(1:NDIM) - sph(p)%r(1:NDIM)
     drsqd= dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
     atemp = (atemp**4 + 0.25_PR*(sink(s)%star_radius**2/drsqd)&
          &*sink(s)%temperature**4)**(0.25_PR)
  end do

#if defined(AMBIENT_HEATING)
  atemp = (temp_inf**4 + atemp**4)**0.25_PR
#endif

! ...
! ----------------------------------------------------------------------------
#elif defined(STAR_HEATING)
  write(*,*) 'Not done yet'
  stop

! Constant background temperature
! ----------------------------------------------------------------------------
#elif defined(AMBIENT_HEATING)
  atemp = temp_inf

#endif
! ----------------------------------------------------------------------------

  return
END SUBROUTINE ambient_temp

