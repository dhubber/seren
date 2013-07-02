! MASS_LOADING_WINDS.F90
! J. E. Dale, D. A. Hubber, J. Ngoumou - 15/08/2011
! Pressure-driven winds algorithm that loads mass as well as momentum 
! onto wind-surface particles.  Uses HEALPix structures to determine which 
! particles are wind surface particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE mass_loading_winds(isource)
  use particle_module
  use HP_module
  use type_module
  use hydro_module
  use scaling_module
  use time_module
  implicit none

  integer, intent(in) :: isource          ! HP source id
  
  integer :: i                            ! Aux. counter
  integer :: iwind                        ! loop index
  integer :: Nrays                        ! No. of Healpix rays
  integer :: nwind                        ! no. injected particles 
  integer :: p                            ! Particle counter
  real(kind=PR) :: dput                   ! distance fr source to put new part
  real(kind=PR) :: dr(1:NDIM)             ! vect. between source & target part
  real(kind=PR) :: dsource                ! Distance to source
  real(kind=PR) :: mdotwind               ! mass loss rate
  real(kind=PR) :: mwind                  ! mass of each wind particle
  real(kind=PR) :: rwsource(1:NDIM)       ! location of wind source
  real(kind=PR) :: vwind                  ! v_infinity
  real(kind=PR) :: Twind                  ! shocked wind temp.    
  real(kind=PR) :: uwind                  ! shocked wind internal energy
  real(kind=PR) :: wtstep                 ! wind timestep
  type(sph_particle) :: windpart          ! struct for new wind particle

  debug2("Computing mass loading for pressure-driven wind [mass_loading_winds.F90]")

! Obtain location of wind source
  rwsource(1:NDIM) = HPsource(isource)%r(1:NDIM)
  write(*,*) rwsource(1:NDIM)

! Obtain mass flux and vinfinity - NB these must agree with mass
! and vinfinity being used by the momentum winds code
  mdotwind = HPsource(isource)%M_loss
  vwind = HPsource(isource)%v_wind

! Determine number of wind particles
  nwind = 0

! Obtain wind timestep
  if (HPsource(isource)%tlastwind .lt. 0.0_DP)then
     HPsource(isource)%tlastwind = time
     return
  else
     wtstep = time - HPsource(isource)%tlastwind
     HPsource(isource)%tlastwind = time
  end if

! Compute all properties of the wind, including the temperature of shocked 
! wind (see Lamers and Cassinelli, eqn 12.3)
  mwind = mdotwind*wtstep
  Twind = 1.4e5_PR*(vwind/(100.0_PR*1000.0_PR/vscale/v_SI))**2
  uwind = Pconst*Twind/(gamma - 1.0_PR)

!#if defined(DEBUG_MASS_LOADING)
  write(*,*) 'mdotwind : ',mdotwind*dmdtscale,trim(dmdtunit)
  write(*,*) 'vwind    : ',vwind*vscale,trim(vunit)
  write(*,*) 'mwind : ',mwind*mscale,trim(munit),mwind/mmean
  write(*,*) 'Twind : ',Twind
  write(*,*) 'wind int. energy: ',uwind*uscale,trim(uunit)

!#endif


! Loop over all particles and identify total wind surface
! =============================================================================
  do p=1,ptot

     ! Ignore if not a surface particle
     if (windsurface(p) == -1) cycle

     ! find unit vector between source and target ISM particle
     dr(1:NDIM) = sph(p)%r(1:NDIM) - rwsource(1:NDIM)
     dsource = sqrt(dot_product(dr(1:NDIM),dr(1:NDIM))) + SMALL_NUMBER
     dr(1:NDIM) = dr(1:NDIM)/dsource

     ! Compute position of new wind particle relative to source
     dput = max(0.5_PR*dsource,dsource - sph(p)%h)
     Nrays = 12*4**(windsurface(p))

     
     ! If the windsurface particle is not a wind particle, or is a wind 
     ! particle that has exceeded its mass threshold, then create a new 
     ! particle.  Otherwise, add the mass, energy and momentum to the 
     ! wind surface particle.
     ! ------------------------------------------------------------------------
     if (sph(p)%windflag .and. &
         &(sph(p)%m < mmean .or. dsource < 4.0_PR*KERNRANGE*sph(p)%h)) then

#if defined(DEBUG_MASS_LOADING)
        write(6,*) "Adding mass to particle ",p,&
            &mwind/Nrays/mmean,sph(p)%m/mmean
#endif

        sph(p)%u = sph(p)%u*sph(p)%m
        sph(p)%m = sph(p)%m + mwind/Nrays
        sph(p)%u = (sph(p)%u + mwind*uwind)/sph(p)%m


     ! ------------------------------------------------------------------------
     else

        ! Create new wind particle with the same h as the target particle.
        windpart%r(1:NDIM) = rwsource(1:NDIM) + dput*dr(1:NDIM)
        windpart%v(1:NDIM) = sph(p)%v(1:NDIM)
        windpart%m         = mwind/Nrays
        windpart%h         = sph(p)%h
        windpart%temp      = Twind
        windpart%u         = uwind
        windpart%windflag  = .true.

        ! Changes counters for new particle
        nwind = nwind + 1
        iwind = ptot + nwind

#if defined(DEBUG_MASS_LOADING)
        write(6,*) "Created new wind particle : ",nwind,sph(p)%windflag,&
           &dput,dsource,4.0_PR*KERNRANGE*sph(p)%h,windpart%m/mmean
#endif

        ! insert new wind particle into main arrays
        call create_new_sph_particle(gasid,iwind,windpart)
     
     end if
     ! ------------------------------------------------------------------------

  end do
! =============================================================================
  
!#if defined(DEBUG_MASS_LOADING)
  write(*,*) 'Injected ',nwind,' new wind particles.'
!#endif

  nnewwind = nnewwind + nwind

#if defined(DEBUG_MASS_LOADING)
  nwind = 0
  do p=1,ptot
     if (windsurface(p) /= -1) nwind = nwind + 1
  end do
  write(6,*) "No. of wind surface particles : ",nwind
  write(6,*) "No. of wind partcles : ",ptot - 100000
#endif


  return
END SUBROUTINE mass_loading_winds
