! DIRECT_GRAVITY.F90
! C. P. Batty & D. A. Hubber - 24/8/2007
! Calculates gravitational accelerations for particle p (or sink s = -p) 
! using direct summation.  Also calculates gravitational potential.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE direct_sph_gravity(p,invhp,zo_p,rp,agravp,potp)
  use interface_module, only : gravity_gradh,gravity_gradh_meanh,&
       &gravity_nbody,gravity_sph,wpot
  use definitions
  use particle_module
  use type_module, only : pgravitystart, pgravityend
#if defined(SINKS) 
  use sink_module
#endif
#if defined(USE_MPI)
  use mpi_communication_module
#endif
  implicit none

  integer, intent(in) :: p                      ! Id of current particle
  real(kind=PR), intent(in) :: invhp            ! Smoothing length of p
  real(kind=PR), intent(in) :: zo_p             ! zeta/omega for p
  real(kind=PR), intent(in) :: rp(1:NDIM)       ! Position of particle p
  real(kind=DP), intent(out) :: agravp(1:NDIM)  ! Gravitational accelertation
  real(kind=DP), intent(out) :: potp            ! Gravitational potential

  integer       :: pp                ! second particle identifier
  real(kind=PR) :: atemp(1:NDIM)     ! temp. grav. accel. variable
  real(kind=PR) :: dpotp             ! grav potential of p due to pp
#if defined(MEANH_GRAVITY)
  real(kind=PR) :: hp                ! Smoothing length of particle
#endif

  debug3("Calculating gravitational force [direct_sph_gravity.F90] for particle",p)

! Initialise gravitational acceleration and potential to self-contribution
  agravp(1:NDIM) = 0.0_PR
  potp = 0.0_PR
  if (p > 0) potp = sph(p)%m*invhp*wpot(0.0_PR)
#if defined(MEANH_GRAVITY)
  hp = 1.0_PR / invhp
#endif


! Loop over all other gravitating particles for SPH particles
! ----------------------------------------------------------------------------
  do pp=pgravitystart,pgravityend
     if (pp == p) cycle
#if defined(N_BODY)
     call gravity_nbody(sph(pp)%m,rp(1:NDIM),&
          &sph(pp)%r(1:NDIM),atemp(1:NDIM),dpotp)
#elif defined(GRAD_H_SPH) && defined(MEANH_GRAVITY)
     if (p < 0) then
        call gravity_meanh(0.5_PR*(hp + sph(pp)%h),sph(pp)%m,&
             &rp(1:NDIM),sph(pp)%r(1:NDIM),atemp(1:NDIM),dpotp)
     else
        call gravity_gradh_meanh(hp,sph(pp)%h,sph(pp)%m,&
             &rp(1:NDIM),sph(pp)%r(1:NDIM),zo_p,sph(pp)%zo,&
             &atemp(1:NDIM),dpotp)
     end if
#elif defined(GRAD_H_SPH)
     if (p < 0) then
        call gravity_sph(invhp,sph(pp)%h,sph(pp)%m,&
             &rp(1:NDIM),sph(pp)%r(1:NDIM),atemp(1:NDIM),dpotp)
     else
        call gravity_gradh(invhp,sph(pp)%h,sph(pp)%m,rp(1:NDIM),&
             &sph(pp)%r(1:NDIM),zo_p,sph(pp)%zo,atemp(1:NDIM),dpotp)
     end if
#elif defined(MEANH_GRAVITY)
     call gravity_meanh(0.5_PR*(hp + sph(pp)%h),sph(pp)%m,&
          &rp(1:NDIM),sph(pp)%r(1:NDIM),atemp(1:NDIM),dpotp)
#else
     call gravity_sph(invhp,sph(pp)%h,sph(pp)%m,&
          &rp(1:NDIM),sph(pp)%r(1:NDIM),atemp(1:NDIM),dpotp)
#endif
     agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),DP)
     potp = potp + real(dpotp,DP)
     
  end do
! ----------------------------------------------------------------------------

! Add to particle gravity_calcs for MPI loadbalancing instrumentation
#if defined(USE_MPI) && !defined(LOW_MEM)
  if (p > 0) sph(p)%gravity_calcs = &
       & sph(p)%gravity_calcs + COST_MONOPOLE*(ptot-1)
#endif

  return
END SUBROUTINE direct_sph_gravity
