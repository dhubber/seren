! HP_STELLAR_FEEDBACK.F90
! J. Ngoumou & D. A. Hubber - 23/03/2011
! Calculate all stellar feedback terms (i.e. momentum-driven winds, 
! stellar irradiation, etc..) using HEALPix sources.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE HP_stellar_feedback(i,isource,level,rfirst,rsource,nextptcl)
  use interface_module, only : distance2,distance3
  use particle_module
  use hydro_module
  use HP_module
  implicit none
  
  integer, intent(in) :: i                     ! HP ray id
  integer, intent(in) :: isource               ! HP source id
  integer, intent(in) :: level                 ! Current HEALPix level
  integer, intent(in) :: nextptcl(1:ptot)      ! forward linked list
  real(kind=PR), intent(out) :: rfirst         ! dist. of 1st ptcl from source
  real(kind=PR), intent(in) :: rsource(1:NDIM) ! position of source

  integer :: clevel                            ! Aux. HP level variable
  integer :: j                                 ! Aux. loop counter
  integer :: nlist                             ! No. of particles in ray list
  integer :: Nrays                             ! Number of Healpix rays
  integer :: p                                 ! Particle id
  integer :: pfirst                            ! id of first particle in ray
  integer, allocatable :: windlist(:)          ! list of particles in ray
  real(kind=PR) :: dr(1:NDIM)                  ! Relative displacment vector
  real(kind=PR) :: drsqd                       ! distance from the source
  real(kind=PR) :: hfirst                      ! Smoothing length of 1st ptcl
  real(kind=PR) :: sumelt                      ! Wind normalisation variable

  real(kind=PR) :: mdotwind                    ! mass loss rate
  real(kind=PR) :: vwind                       ! v_infinity
  real(kind=PR) :: V_each(1:NDIM)              ! ..
  real(kind=PR) :: vshell                      ! ..
  real(kind=PR) :: FACT                        ! ..

! If the ray has already been done, exit routine.
  if (HPray(i)%winddone) return

! Set wind properties
  mdotwind = HPsource(isource)%M_loss
  vwind = HPsource(isource)%v_wind

! Find the required ray resolution from the first particle's position and 
! smoothing length
  pfirst = HPray(i)%first
  call distance2(rsource(1:NDIM),pfirst,dr(1:NDIM),drsqd)
  rfirst = sqrt(drsqd) 
  hfirst = sph(pfirst)%h
  clevel = int(log10(rfirst/(f2*hfirst))*INVLOG10TWO + 0.5_PR)

! If the ray resolution is too low, return to main healpix routine and wait 
! for next level for increased resolution.
  if (clevel > level) return

  allocate(windlist(1:pmax))
  windlist(1:pmax) = -1
  Nrays = 12*4**(level)

! Start with first particle in list.
  nlist = 1
  windlist(nlist) = pfirst
  sumelt = sph(pfirst)%m/drsqd
  p = pfirst

! Record level if particle is on surface of wind
#if defined(STELLAR_WIND) && defined(PARTICLE_INJECTION_WINDS)
  windsurface(pfirst) = level
#endif

! Follow linked list to find all other particles in list within 2*hfirst 
! of the first particle.  Determine the normalisation factor.
  do 
     if (p == HPray(i)%last) exit
     p =  nextptcl(p)  
     call distance2(rsource(1:NDIM),p,dr(1:NDIM),drsqd)
     if (drsqd > (rfirst + KERNRANGE*hfirst)**2) exit
     nlist = nlist + 1
     windlist(nlist) = p
     sumelt = sumelt + sph(p)%m/drsqd
  end do
  
! Now we have calculated normalisation factor, loop over particles again and 
! compute wind acceleration on each particle.
  do j=1,nlist
     p = windlist(j)
     call distance2(rsource(1:NDIM),p,dr(1:NDIM),drsqd)
#if defined(MOMENTUM_WIND) && !defined(PARTICLE_INJECTION_WINDS)
     sph(p)%a_wind(1:NDIM) = sph(p)%a_wind(1:NDIM) + &
          & HPsource(isource)%M_loss*HPsource(isource)%v_wind/&
          & (Nrays*drsqd*sumelt)*dr(1:NDIM)/sqrt(drsqd)
#endif
#if defined(STELLAR_LUMINOSITY)
     sph(p)%dudt_rad = sph(p)%dudt_rad + HPsource(isource)%L_star/&
          & (Nrays*drsqd*sumelt)
#endif
#if defined(PARTICLE_INJECTION_WINDS)
     if (HPsource(isource)%numb_inj .lt. 1)then
        FACT = 0.0_PR  ! max(1.0_PR,2.4_PR*time**(0.1_PR))
     else 
        FACT = 0.0_PR
     end if
     sph(p)%a_wind(1:NDIM) = sph(p)%a_wind(1:NDIM) + &
          & (FACT)*HPsource(isource)%M_loss*HPsource(isource)%v_wind/&
          & (Nrays*drsqd*sumelt)*dr(1:NDIM)/sqrt(drsqd)
#endif
  end do
  
! Flag ray as done
  HPray(i)%winddone = .true.
#if !defined (IONIZING_UV_RADIATION)
  HPray(i)%done = .true.
#endif
  deallocate(windlist)
  
  return
END SUBROUTINE HP_stellar_feedback




