! NBODY_ACCRETE_BOUND_PARTICLES.F90
! D. A. Hubber - 23/7/2010
! Accrete all SPH particles that are bound to any sink particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_accrete_bound_particles
  use particle_module
  use sink_module
  use Nbody_module
  use type_module
  implicit none

#if defined(GRAVITY)
  integer :: ndead                    ! Number of newly accreted particles
  integer :: p                        ! Particle counter
  integer :: s                        ! Sink counter
  integer :: ss                       ! Secondary sink counter
  integer :: sinkid_p(1:SMAX)         ! Auxilary sorting arrays
  integer, allocatable :: deadlist(:) ! List of dead (accreted) particles
  integer, allocatable :: sinkid(:)   ! Sink that particle p is most bound to
  real(kind=DP) :: atemp(1:NDIM)      ! Auxilary accel variable 
  real(kind=PR) :: atemp_pr(1:NDIM)   ! ..
  real(kind=PR) :: dpotp              ! Aux. sink pot variable
  real(kind=DP) :: dr(1:NDIM)         ! Relative position vector (DP)
  real(kind=DP) :: drsqd              ! Distance squared (DP)
  real(kind=PR) :: dv(1:NDIM)         ! Relative velocity vector
  real(kind=DP) :: dv2(1:NDIM)        ! Relative velocity vector (DP)
  real(kind=PR) :: energy_p(1:SMAX)   ! Energy of p-s system
  real(kind=PR) :: invhs              ! Smoothing length of sink s
  real(kind=PR) :: invhp              ! 1 / hp
  real(kind=DP) :: invdrmag           ! 1 / drmag
  real(kind=DP) :: invdrsqd           ! 1 / drsqd
  real(kind=PR) :: mp                 ! Mass of particle p
  real(kind=PR) :: ms                 ! Mass of sink s
  real(kind=DP) :: potp               ! Potential of sink s
  real(kind=PR) :: reduced_mass       ! Reduced mass of sink-particle system
  real(kind=PR) :: rp(1:NDIM)         ! Position of particle p
  real(kind=PR) :: rs(1:NDIM)         ! Position of sink s
  real(kind=DP) :: rss(1:NDIM)        ! Position of sink ss
  real(kind=PR) :: vp(1:NDIM)         ! Velocity of particle p
  real(kind=PR) :: vs(1:NDIM)         ! Velocity of sink s
  real(kind=DP) :: vsqd               ! Velocity squared (DP)


! Remove gas if there is any present in the simulation
! ============================================================================
  if (pgas > 0) then    

     debug2("Adding bound SPH particles to stars [nbody_setup.F90]")

     allocate(deadlist(1:ptot))
     allocate(sinkid(1:ptot))
     sinkid(1:ptot) = -1
     ndead = 0

     ! Loop over all particles and find sink each particle is most bound to.
     ! -----------------------------------------------------------------------
     do p=pgravitystart,pgravityend
        rp(1:NDIM) = sph(p)%r(1:NDIM)
        mp         = sph(p)%m
        invhp      = 1.0_PR / sph(p)%h
        vp(1:NDIM) = sph(p)%v(1:NDIM)
        energy_p(1:stot) = BIG_NUMBER
        
        ! Loop over all sinks and calculate the 2-body energies
        do s=1,stot
           sinkid_p(s) = s
           dv(1:NDIM) = sink(s)%v(1:NDIM) - vp(1:NDIM)
           reduced_mass = mp*sink(s)%m / (mp + sink(s)%m)
#if defined(N_BODY)
           call gravity_nbody(sink(s)%m,rp(1:NDIM),&
                &sink(s)%r(1:NDIM),atemp_pr(1:NDIM),dpotp)
#elif defined(GRAD_H_SPH)
           call gravity_gradh(invhp,sink(s)%h,sink(s)%m,rp(1:NDIM),&
                &sink(s)%r(1:NDIM),0.0_PR,0.0_PR,atemp_pr(1:NDIM),dpotp)
#else
           call gravity_sph(invhp,sink(s)%h,sink(s)%m,rp(1:NDIM),&
                &sink(s)%r(1:NDIM),atemp_pr(1:NDIM),dpotp)
#endif
           energy_p(s) = 0.5_PR*reduced_mass*&
                &dot_product(dv(1:NDIM),dv(1:NDIM)) - dpotp*mp
        end do
        
        ! Sort energies into ascending order.  Only accrete particle to 
        ! sink if it is bound to only one sink.
        call insertion_sort_real(stot,sinkid_p(1:stot),energy_p(1:stot))
        if (energy_p(1) < 0.0_PR .and. energy_p(2) > 0.0_PR) &
             & sinkid(p) = sinkid_p(1)
 
     end do
     ! -----------------------------------------------------------------------


     ! Now loop back over sinks and add mass and momentum from particles
     ! -----------------------------------------------------------------------
     do s=1,stot
        ms = sink(s)%m
        rs(1:NDIM) = ms*sink(s)%r(1:NDIM)
        vs(1:NDIM) = ms*sink(s)%v(1:NDIM)

        do p=pgravitystart,pgravityend
           if (sinkid(p) /= s) cycle
           ms = ms + sph(p)%m
           rs(1:NDIM) = rs(1:NDIM) + sph(p)%m*sph(p)%r(1:NDIM)
           vs(1:NDIM) = vs(1:NDIM) + sph(p)%m*sph(p)%v(1:NDIM)
           ndead = ndead + 1
           deadlist(ndead) = p
        end do

        ! Record new sink info in main arrays
        sink(s)%m = ms
        sink(s)%r(1:NDIM) = rs(1:NDIM) / ms
        sink(s)%v(1:NDIM) = vs(1:NDIM) / ms
     end do
     ! -----------------------------------------------------------------------


     ! If particles have been accreted, remove from main arrays
     if (ndead > 0) then
        call insertion_sort_int(ndead,deadlist)
        call remove_from_list(ndead,deadlist)
     end if

     if (allocated(sinkid)) deallocate(sinkid)
     if (allocated(deadlist)) deallocate(deadlist)

  end if
! ============================================================================
#endif

  return
END SUBROUTINE nbody_accrete_bound_particles
