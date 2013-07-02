! BINARY_PROPERTIES.F90
! D. A. Hubber - 12/3/2008
! Calculates binary properties from mechanical information (m,r,v).
! Also stores binary properties in binary array.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE binary_properties(s1,s2,binen,m1,m2,r1,r2,v1,v2)
  use interface_module, only : distance3_dp
  use definitions
  use sink_module
  use Nbody_module
  use scaling_module
  use constant_module
  implicit none

  integer, intent(in) :: s1                ! id of system 1
  integer, intent(in) :: s2                ! id of system 2
  real(kind=DP),intent(in) :: binen        ! Two-body energy
  real(kind=DP),intent(in) :: m1           ! Mass of system 1
  real(kind=DP),intent(in) :: m2           ! Mass of system 2
  real(kind=DP),intent(in) :: r1(1:NDIM)   ! Position of system 1
  real(kind=DP),intent(in) :: r2(1:NDIM)   ! Position of system 2
  real(kind=DP),intent(in) :: v1(1:NDIM)   ! Velocity of system 1
  real(kind=DP),intent(in) :: v2(1:NDIM)   ! Velocity of system 2

  real(kind=DP) :: dr(1:NDIM)              ! Relative displacement vector
  real(kind=DP) :: drmag                   ! Distance
  real(kind=DP) :: drsqd                   ! Distance squared
  real(kind=DP) :: dv(1:NDIM)              ! Relative velocity
  real(kind=DP) :: eccent                  ! Eccentricity
  real(kind=DP) :: L(1:3)                  ! Angular momentum vector
  real(kind=DP) :: Lsqd                    ! Angular momentum squared
  real(kind=DP) :: period                  ! Period
  real(kind=DP) :: q                       ! Mass-ratio
  real(kind=DP) :: rbin(1:NDIM)            ! Position of centre of mass
  real(kind=DP) :: reduced_mass            ! Reduced mass
  real(kind=DP) :: sma                     ! Semi-major axis
  real(kind=DP) :: vbin(1:NDIM)            ! Velocity of centre of mass

! Calculate position and velocity and centre of mass and the relative velocity
  rbin(1:NDIM) = (m1*r1(1:NDIM) + m2*r2(1:NDIM))/(m1 + m2)
  vbin(1:NDIM) = (m1*v1(1:NDIM) + m2*v2(1:NDIM))/(m1 + m2)
  dv(1:NDIM) = v2(1:NDIM) - v1(1:NDIM)
  call distance3_dp(r1(1:NDIM),r2(1:NDIM),dr(1:NDIM),drsqd)
  drmag = sqrt(drsqd)
           
! Calculate angular momentum squared of binary
  reduced_mass = m1*m2 / (m1 + m2)
#if NDIM==3
  L(1) = (dr(2)*dv(3) - dr(3)*dv(2))*reduced_mass
  L(2) = (dr(3)*dv(1) - dr(1)*dv(3))*reduced_mass
  L(3) = (dr(1)*dv(2) - dr(2)*dv(1))*reduced_mass
  Lsqd = L(1)*L(1) + L(2)*L(2) + L(3)*L(3)
#elif NDIM == 2
  L(3) = (dr(1)*dv(2) - dr(2)*dv(1))*reduced_mass
  Lsqd = L(3)*L(3)
#endif
        
! Now calculate all binary parameters
  sma    = -0.5_DP*m1*m2/binen
  eccent = 1.0_DP - Lsqd/(m1 + m2)/sma/(reduced_mass**2)
  eccent = max(0.0_DP,eccent)
  eccent = sqrt(eccent)
  period = 2.0_DP*PI*sqrt(sma**3/(m1 + m2))
  if (m1 > m2) then
     q = m2/m1
  else
     q = m1/m2
  end if

! Record all binary information in main arrays
  binary(nbin)%id          = stot + nbin
  binary(nbin)%s1          = s1
  binary(nbin)%s2          = s2
  binary(nbin)%r(1:NDIM)   = rbin(1:NDIM)
  binary(nbin)%v(1:NDIM)   = vbin(1:NDIM)
  binary(nbin)%m           = m1 + m2
  binary(nbin)%angmom(1:3) = L(1:3)
  binary(nbin)%binen       = binen
  binary(nbin)%ecc         = eccent
  binary(nbin)%period      = period
  binary(nbin)%q           = q
  binary(nbin)%sma         = sma
  binary(nbin)%drmag       = drmag


! Output information to screen
! ----------------------------------------------------------------------------
#if defined(DEBUG_BINARY_PROPERTIES)
  if (s1 <= stot .and. s2 <= stot) then
     write(6,*) "Binary ",stot+nbin," from stars ",s1," and ",s2
  else if (s1 > stot .and. s2 > stot) then
     write(6,*) "Quadrupole ",stot+nbin," formed from binaries ",s1," and ",s2
  else if (s1 <= stot .and. s2 > stot) then
     write(6,*) "Triple ",stot+nbin," formed from star ",s1," and binary ",s2
  else if (s2 <= stot .and. s1 > stot) then
     write(6,*) "Triple ",stot+nbin," formed from star ",s2," and binary ",s1
  end if

  write(6,*) "Masses :",m1*mscale,m2*mscale,trim(munit)
  write(6,*) "binen  :",binen*Escale,trim(Eunit)
  write(6,*) "sma    :",sma*rscale,trim(runit)
  write(6,*) "period :",period*tscale,trim(tunit)
  write(6,*) "eccent :",eccent
  write(6,*) "q      :",q
  write(6,*)
#endif


  return
END SUBROUTINE binary_properties
