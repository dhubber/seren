! IC_SUBROUTINES.F90
! D. A. Hubber - 08/02/2010
! Various subroutines used for generating initial conditions.
! - add_azimuthal_density_perturbation
! - add_radial_density_gradient
! - add_rotational_velocity_field
! - add_turbulent_velocity_field
! - sort_radial_particle_types
! - solve_isothermal_eqn
! - solve_lane_emden_eqn
! ============================================================================

#include "macros.h"

! ============================================================================
! ADD_AZIMUTHAL_DENSITY_PERTURBATION
! Add an mpert azimutal density perturbation of amplitude amp to the 
! density field centred on r0.
! ============================================================================
SUBROUTINE add_azimuthal_density_perturbation(pertmode,mpert,&
     &pstart,pend,amp,rinner,r0)
  use particle_module
  use type_module
  implicit none

  character(len=50), intent(in) :: pertmode  ! Name of IC file
  integer, intent(in) :: mpert               ! Mode of perturbation
  integer, intent(in) :: pend                ! Last particle
  integer, intent(in) :: pstart              ! First particle
  real(kind=PR), intent(in) :: amp           ! Amplitude of perturbation
  real(kind=PR), intent(in) :: r0(1:NDIM)    ! Centre of sphere to perturb
  real(kind=PR), intent(in) :: rinner        ! 'Core' radius for Machida core

  logical :: flag                            ! Aux. flag
  integer :: i                               ! Aux. loop counter
  integer :: p                               ! Particle counter
  integer :: tabtot                          ! Total number of entries in table
  real(kind=PR) :: ampaux                    ! Local azimuthal amplitude
  real(kind=PR) :: phi                       ! Angle
  real(kind=PR) :: phi1                      ! Aux. table angle 1
  real(kind=PR) :: phi2                      ! Aux. table abgle 2
  real(kind=PR) :: phiprime                  ! Perturbed angle
  real(kind=PR) :: rmag                      ! Distance 
  real(kind=PR) :: rp(1:NDIM)                ! Position of particle p
  real(kind=PR) :: rsqd                      ! Distance squared
  real(kind=PR) :: spacing                   ! Table spacing

  write(6,*) "Adding azimuthal density perturbation"
  write(6,*) "mpert : ",mpert,"     amp : ",amp

! Set up search table
  tabtot  = 10000
  spacing = 2.0_PR*PI / real(tabtot - 1,PR)

! Perturb positions of particles in cloud
! ----------------------------------------------------------------------------
  do p=pstart,pend
     rp(1:NDIM) = sph(p)%r(1:NDIM) - r0(1:NDIM)

     ! Perturbation mode
     if (pertmode == "machida") then
        rsqd = dot_product(rp(1:NDIM),rp(1:NDIM))
        if (rsqd < rinner*rinner) then
           ampaux = amp*rsqd/rinner/rinner
        else
           ampaux = amp
        end if
     else
        ampaux = amp
     end if
     if (ampaux > 1.0_PR) ampaux = 0.999999999_PR

     ! Make local copy of position and obtain axial distance from z-axis
     rsqd = rp(1)*rp(1) + rp(2)*rp(2)
     rmag = sqrt(rsqd)

     ! Find azimuthal angle around z-axis
     if (rmag > SMALL_NUMBER) then
        phi = asin(abs(rp(2))/rmag)
     else
        phi = 0.0_PR
     end if

     if (rp(1) < 0.0_PR .and. rp(2) > 0.0_PR) then 
        phi = PI - phi
     else if (rp(1) < 0.0_PR .and. rp(2) < 0.0_PR) then
        phi = PI + phi
     else if (rp(1) > 0.0_PR .and. rp(2) < 0.0_PR) then
        phi = 2.0_PR*PI - phi
     end if

     ! Wrap around if table begins after phi
     if (phi < ampaux/real(mpert,PR)) phi = phi + 2.0_PR*PI
     phiprime = phi

     ! Numerically find new phi angle for perturbation.  Search through 
     ! grid of values, find upper and lower bounds, then use linear 
     ! interpolation to find new value of phi
     flag = .false.

     do i=1,tabtot-1
        phi1 = spacing*real(i - 1,PR)
        phi2 = spacing*real(i,PR)
        phi1 = phi1 + ampaux*cos(real(mpert,PR)*phi1)/real(mpert,PR)
        phi2 = phi2 + ampaux*cos(real(mpert,PR)*phi2)/real(mpert,PR)

        if (phi2 >= phi .and. phi1 < phi) then
           phiprime = spacing*real(i - 1,PR) + &
                &spacing*(phi - phi1) / (phi2 - phi1)
           flag = .true.
        end if
        if (flag) exit
     end do

     ! Reposition particle with new angle
     sph(p)%r(1) = r0(1) + rmag*cos(phiprime)
     sph(p)%r(2) = r0(2) + rmag*sin(phiprime)

  end do
! ----------------------------------------------------------------------------

  return
END SUBROUTINE add_azimuthal_density_perturbation



! ============================================================================
! ADD_RADIAL_DENSITY_GRADIENT
! Add a radial density profile of the form 'rho = rho0 + drhodr*r' 
! for a circular (NDIM=2) or spherical (NDIM=3) distribution of unity radius.
! ============================================================================
SUBROUTINE add_radial_density_gradient(pstart,pend,rho0,drhodr,rcloud,r0)
  use particle_module
  implicit none

  integer, intent(in) :: pstart            ! First particle
  integer, intent(in) :: pend              ! Last particle
  real(kind=PR), intent(in) :: rho0        ! Central density
  real(kind=PR), intent(in) :: drhodr      ! Density gradient
  real(kind=PR), intent(in) :: rcloud        ! Maximum radius of cloud
  real(kind=PR), intent(in) :: r0(1:NDIM)  ! Centre of radial distribution

  integer :: p                             ! Particle counter
  real(kind=PR) :: dr(1:NDIM)              ! Relative displacement
  real(kind=PR) :: drmag                   ! Distance
  real(kind=PR) :: drsqd                   ! Distance squared
  real(kind=PR) :: rhigh                   ! Upper limit of binary search
  real(kind=PR) :: rlow                    ! Lower limit of binary search
  real(kind=PR) :: rguess                  ! 'Guess' of radial distance
  real(kind=PR) :: m0                      ! Total mass in density field
  real(kind=PR) :: menc                    ! Enclosed mass
  real(kind=PR) :: mp                      ! Mean particle mass

! Calculate total mass depending on dimensionality
#if NDIM==2
  m0 = PI*rho0*rcloud*rcloud + TWOTHIRDS*PI*drhodr*(rcloud**3)
#elif NDIM==3
  m0 = 4.0_PR*ONETHIRD*PI*rho0*(rcloud**3) + PI*drhodr*(rcloud**4)
#endif
  mp = m0 / real(pend - pstart + 1,PR)

  write(6,*) "m0 : ",m0,"    mp : ",mp,"    plimits : ",pstart,pend
  write(6,*) "r0 : ",r0(1:NDIM),"    rho : ",rho0,"   drhodr : ",drhodr

! Reposition all particles to reproduce density gradient
! ----------------------------------------------------------------------------
  do p=pstart,pend

     dr(1:NDIM) = sph(p)%r(1:NDIM)
     drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
     drmag = sqrt(drsqd) + SMALL_NUMBER
     rlow = 0.0_PR
     rhigh = rcloud

     ! Binary chop to find correct radial distance of particle
     ! -----------------------------------------------------------------------
     do
        rguess = 0.5_PR*(rlow + rhigh)
#if NDIM==2
        menc = PI*rho0*rguess*rguess + TWOTHIRDS*PI*drhodr*(rguess**3)
#elif NDIM==3
        menc = 4.0_PR*ONETHIRD*PI*rho0*(rguess**3) + PI*drhodr*(rguess**4)
#endif
        
        ! Alter limits of binary chop
        if (menc > m0*drmag**(NDIMPR)) then
           rhigh = rguess
        else
           rlow = rguess
        end if

        ! Check if convergence has been achieved
        if ((rhigh - rlow)/rcloud < 1.0E-6) exit

     end do
     ! -----------------------------------------------------------------------

     sph(p)%r(1:NDIM) = r0(1:NDIM) + rguess*dr(1:NDIM)/drmag
     sph(p)%m = mp

  end do
! ----------------------------------------------------------------------------

  return
END SUBROUTINE add_radial_density_gradient




! ============================================================================
! ADD_RADIAL_VELOCITY_FIELD
! Add a radial velocity field to the selected particle distribution.
! ============================================================================
SUBROUTINE add_radial_velocity_field(velmode,beta0,dvdr,pstart,pend,r0)
  use interface_module, only : distance2
  use particle_module
  implicit none

  character(len=50), intent(in) :: velmode  ! ..
  integer, intent(in) :: pstart             ! First particle
  integer, intent(in) :: pend               ! Last particle
  real(kind=PR), intent(in) :: beta0        ! Desired rotational k.e.
  real(kind=PR), intent(in) :: dvdr         ! ..
  real(kind=PR), intent(in) :: r0(1:NDIM)   ! Centre of sphere to perturb

  integer :: p
  real(kind=DP) :: kerad                    ! ..
  real(kind=PR) :: drsqd                    ! Distance from origin
  real(kind=PR) :: dr(1:NDIM)               ! Position of particle p
  real(kind=PR) :: vsqd                     ! ..
  real(kind=PR), allocatable :: vrad(:,:)   ! ..

  write(6,*) "Adding radial velocity field : ",velmode
  write(6,*) "beta0 : ",beta0,"     dvdr : ",dvdr

  allocate(vrad(1:VDIM,1:ptot))
  vrad(1:VDIM,1:ptot) = 0.0_PR


! Calculate radial velocity components
! ----------------------------------------------------------------------------
  do p=pstart,pend
     call distance2(r0(1:NDIM),p,dr(1:NDIM),drsqd)
     vrad(1:NDIM,p) = dvdr*dr(1:NDIM)
  end do
! ----------------------------------------------------------------------------


! Scale velocities to give required rotational kinetic energy
! ----------------------------------------------------------------------------
  if (velmode == "energy") then
     kerad = 0.0_DP
     do p=pstart,pend
        vsqd = dot_product(vrad(1:NDIM,p),vrad(1:NDIM,p))
        kerad = kerad + real(sph(p)%m*vsqd,DP)
     end do
     kerad = 0.5_DP*kerad
     write(6,*) "kerad : ",kerad,"    scaling factor : ",sqrt(beta0/kerad)
     
     ! Now rescale velocities and add to particle velocity
     do p=pstart,pend
        sph(p)%v(1:NDIM) = sph(p)%v(1:NDIM) + vrad(1:NDIM,p)*sqrt(beta0/kerad)
     end do

! Scale velocities to give required angular momentum
! ----------------------------------------------------------------------------
  else if (velmode == "dvdr") then
     stop 'dvdr option not written yet'

     do p=pstart,pend
        sph(p)%v(1:NDIM) = sph(p)%v(1:NDIM) + vrad(1:NDIM,p)
     end do
     
! If no valid option is selected, end program here.
! ----------------------------------------------------------------------------
  else if (velmode /= "none") then
     stop 'Invalid rotational velocity option selected'

  end if
! ----------------------------------------------------------------------------

  deallocate(vrad)

  return
END SUBROUTINE add_radial_velocity_field




! ============================================================================
! ADD_ROTATIONAL_VELOCITY_FIELD
! Add a rotational velocity field to the selected particle distribution.
! ============================================================================
SUBROUTINE add_rotational_velocity_field(velmode,beta0,&
     &angmom,angpower,angvel,pstart,pend,r0)
  use particle_module
  implicit none

  character(len=50), intent(in) :: velmode  ! ..
  integer, intent(in) :: pstart             ! First particle
  integer, intent(in) :: pend               ! Last particle
  real(kind=PR), intent(in) :: beta0        ! Desired rotational k.e.
  real(kind=PR), intent(in) :: angmom       ! ..
  real(kind=PR), intent(in) :: angpower     ! ..
  real(kind=PR), intent(in) :: angvel       ! ..
  real(kind=PR), intent(in) :: r0(1:NDIM)   ! Centre of sphere to perturb

  integer :: p                              ! Particle counter
  real(kind=DP) :: kerot                    ! Rotational kinetic energy
  real(kind=PR) :: phi                      ! Azimuthal angle
  real(kind=PR) :: rmag                     ! Distance from origin
  real(kind=PR) :: rp(1:NDIM)               ! Position of particle p
  real(kind=PR) :: rsqd                     ! Distance squared
  real(kind=PR) :: vmag                     ! Speed
  real(kind=PR) :: vsqd                     ! Speed squared
  real(kind=PR), allocatable :: vrot(:,:)   ! Rotational velocities

  write(6,*) "Adding rotational velocity field : ",velmode
  write(6,*) "beta0 : ",beta0,"    angmom : ",angmom
  write(6,*) "angpower : ",angpower,"     angvel : ",angvel

  allocate(vrot(1:VDIM,1:ptot))
  vrot(1:VDIM,1:ptot) = 0.0_PR

! Loop over all selected particles
! ----------------------------------------------------------------------------
  do p=pstart,pend

     ! Make local copy of position and obtain axial distance from z-axis
     rp(1:NDIM) = sph(p)%r(1:NDIM) - r0(1:NDIM)
     rsqd = rp(1)*rp(1) + rp(2)*rp(2) 
     rmag = sqrt(rsqd)

     ! Find azimuthal angle around z-axis
     if (rmag > SMALL_NUMBER) then
        phi = asin(abs(rp(2))/rmag)
     else
        phi = 0.0_PR
     end if

     ! Make sure angle is in correct quadrant.
     if (rp(1) < 0.0_PR .and. rp(2) > 0.0_PR) then 
        phi = PI - phi
     else if (rp(1) < 0.0_PR .and. rp(2) < 0.0_PR) then
        phi = PI + phi
     else if (rp(1) > 0.0_PR .and. rp(2) < 0.0_PR) then
        phi = 2.0_PR*PI - phi
     end if

     ! Set-up rotational field depending on mode chosen
     if (velmode == "angvel") then
        vmag = angvel*rmag
     else
        vmag = rmag**(angpower + 1.0_PR)
     end if

     ! Set velocities in temporary array before scaling
     vrot(1,p) = -vmag*sin(phi)
     vrot(2,p) = vmag*cos(phi)
     vrot(3,p) = 0.0_PR

  end do
! ----------------------------------------------------------------------------


! Scale velocities to give required rotational kinetic energy
! ----------------------------------------------------------------------------
  if (velmode == "energy") then
     kerot = 0.0_DP
     do p=pstart,pend
        vsqd = dot_product(vrot(1:NDIM,p),vrot(1:NDIM,p))
        kerot = kerot + real(sph(p)%m*vsqd,DP)
     end do
     kerot = 0.5_DP*kerot
     write(6,*) "kerot : ",kerot,"    scaling factor : ",sqrt(beta0/kerot)
     
     ! Now rescale velocities and add to particle velocity
     do p=pstart,pend
        sph(p)%v(1:NDIM) = sph(p)%v(1:NDIM) + vrot(1:NDIM,p)*sqrt(beta0/kerot)
     end do

! Scale velocities to give required angular momentum
! ----------------------------------------------------------------------------
  else if (velmode == "angmom") then
     stop 'angmom option not written yet'

     do p=pstart,pend
        sph(p)%v(1:NDIM) = sph(p)%v(1:NDIM) + vrot(1:NDIM,p)*angmom
     end do
     
! Add velocities for angvel case
! ----------------------------------------------------------------------------
  else if (velmode == "angvel") then

     ! Now rescale velocities and add to particle velocity
     do p=pstart,pend
        sph(p)%v(1:NDIM) = sph(p)%v(1:NDIM) + vrot(1:NDIM,p)
     end do

! If no valid option is selected, end program here.
! ----------------------------------------------------------------------------
  else if (velmode /= "none") then
     stop 'Invalid rotational velocity option selected'

  end if
! ----------------------------------------------------------------------------

  deallocate(vrot)

  return
END SUBROUTINE add_rotational_velocity_field




! ============================================================================
! ADD_TURBULENT_VELOCITY_FIELD
! Add a turbulent velocity field to the selected particle distribution.
! ============================================================================
SUBROUTINE add_turbulent_velocity_field(eturb,vpower,&
     pstart,pend,ngrid,iseed1,iseed2)
  use particle_module
  use filename_module
  use time_module
  implicit none

  integer, intent(in) :: iseed1         ! Random number seed 1
  integer, intent(in) :: iseed2         !   "      "     "   2
  integer, intent(in) :: ngrid          ! No. of grid points in each dimension
  integer, intent(in) :: pstart         ! First particle
  integer, intent(in) :: pend           ! Last particle
  real(kind=PR), intent(in) :: eturb    ! Total turbulent kinetic energy
  real(kind=PR), intent(in) :: vpower   ! Power-spectrum index

  integer :: i                          ! x-grid coordinate
  integer :: j                          ! y-grid coordinate
  integer :: k                          ! z-grid coordinate
  integer :: ngridtemp                  ! ..
  integer :: p                          ! Particle counter
  real(kind=PR) :: dx(1:NDIM)           ! ..
  real(kind=PR) :: dxgrid               ! Grid spacing
  real(kind=DP) :: keturb               ! Total kinetic energy
  !real(kind=PR) :: rmax(1:NDIM)         ! Maximum bounding box extent
  !real(kind=PR) :: rmin(1:NDIM)         ! Minimum bounding box extent
  real(kind=PR) :: vint(1:8)            ! ..
  real(kind=PR) :: vsqd                 ! Speed squared
  real(kind=PR), allocatable :: vtable(:,:,:,:)  ! Velocity field grid
  real(kind=PR), allocatable :: vturb(:,:)       ! Turbulent velocity array
  real(kind=PR), allocatable :: r(:,:)

  write(6,*) "Adding turbulent velocity field"
  write(6,*) "eturb : ",eturb,"     vpower : ",vpower,"    ngrid : ",ngrid

! Compute turbulent velocity field onto a grid
! ----------------------------------------------------------------------------
  if (eturb > SMALL_NUMBER) then
     allocate(vturb(1:NDIM,1:ptot))
#if defined(FFTW)
     allocate(vtable(1:NDIM,1:ngrid,1:ngrid,1:ngrid))
     call turbsub(2,1.0_PR,vpower/2.0_PR,ngrid,ngrid,vtable)
#else
     stop 'FFTW flag not activated.  Cannot generate turbulent velocity field'
#endif
     allocate(r(1:NDIM,1:ptot))
     do p=pstart,pend
        r(1:NDIM,p+1-pstart) = sph(p)%r(1:NDIM)
     end do
     call bounding_box(pend-pstart+1,r,rmax(1:NDIM),rmin(1:NDIM))
     deallocate(r)
     dxgrid = maxval(rmax(1:NDIM) - rmin(1:NDIM))/real(ngrid - 1,PR)

     ! Now interpolate gridded velocity field onto particle positions
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dx,i,j,k,vint)
     do p=pstart,pend
        dx(1) = (sph(p)%r(1) - maxval(rmin))/dxgrid
        dx(2) = (sph(p)%r(2) - maxval(rmin))/dxgrid
        dx(3) = (sph(p)%r(3) - maxval(rmin))/dxgrid
        i = int(dx(1)) + 1
        j = int(dx(2)) + 1
        k = int(dx(3)) + 1
        if (i >= ngrid .or. j >= ngrid .or. k >= ngrid) then
           write(6,*) "Grid too big! : ",&
                &i,j,k,ngrid,sph(p)%r(1:NDIM),rmin,rmax,dx(1:NDIM)
           stop
        end if
        dx(1) = dx(1) - int(dx(1))
        dx(2) = dx(2) - int(dx(2))
        dx(3) = dx(3) - int(dx(3))

        ! Interpolate to get more accurate velocities
#if NDIM==3
        vint(1) = (1.0_PR - dx(1))*(1.0_PR - dx(2))*(1.0_PR - dx(3))
        vint(2) = (1.0_PR - dx(1))*(1.0_PR - dx(2))*(dx(3))
        vint(3) = (1.0_PR - dx(1))*(dx(2))*(1.0_PR - dx(3))
        vint(4) = (1.0_PR - dx(1))*(dx(2))*(dx(3))
        vint(5) = (dx(1))*(1.0_PR - dx(2))*(1.0_PR - dx(3))
        vint(6) = (dx(1))*(1.0_PR - dx(2))*(dx(3))
        vint(7) = (dx(1))*(dx(2))*(1.0_PR - dx(3))
        vint(8) = (dx(1))*(dx(2))*(dx(3))

        vturb(1,p)=(vint(1)*vtable(1,i,j,k) + vint(2)*vtable(1,i,j,k+1)+&
             & vint(3)*vtable(1,i,j+1,k) + vint(4)*vtable(1,i,j+1,k+1)+&
             & vint(5)*vtable(1,i+1,j,k) + vint(6)*vtable(1,i+1,j,k+1)+&
             & vint(7)*vtable(1,i+1,j+1,k) + vint(8)*vtable(1,i+1,j+1,k+1))
        vturb(2,p)=(vint(1)*vtable(2,i,j,k) + vint(2)*vtable(2,i,j,k+1)+&
             & vint(3)*vtable(2,i,j+1,k) + vint(4)*vtable(2,i,j+1,k+1)+&
             & vint(5)*vtable(2,i+1,j,k) + vint(6)*vtable(2,i+1,j,k+1)+&
             & vint(7)*vtable(2,i+1,j+1,k) + vint(8)*vtable(2,i+1,j+1,k+1))
        vturb(3,p)=(vint(1)*vtable(3,i,j,k) + vint(2)*vtable(3,i,j,k+1)+&
             & vint(3)*vtable(3,i,j+1,k) + vint(4)*vtable(3,i,j+1,k+1)+&
             & vint(5)*vtable(3,i+1,j,k) + vint(6)*vtable(3,i+1,j,k+1)+&
             & vint(7)*vtable(3,i+1,j+1,k) + vint(8)*vtable(3,i+1,j+1,k+1))
#endif
     end do
     !$OMP END PARALLEL DO

     ! Calculate total kinetic energy of turbulent velocity field and 
     ! the total gravitational potential energy
     keturb = 0.0_DP
     do p=pstart,pend
        vsqd = dot_product(vturb(1:NDIM,p),vturb(1:NDIM,p))
        keturb = keturb + real(sph(p)%m*vsqd,DP)
     end do
     keturb = 0.5_DP*keturb

     write(6,*) "keturb : ",keturb
     write(6,*) "scaling factor : ",sqrt(eturb/keturb)

     ! Now rescale velocities and add to particle velocity
     do p=pstart,pend
        sph(p)%v(1:NDIM) = sph(p)%v(1:NDIM) + vturb(1:NDIM,p)*sqrt(eturb/keturb)
     end do

     deallocate(vtable)
     deallocate(vturb)

  end if
! ----------------------------------------------------------------------------

  return
END SUBROUTINE add_turbulent_velocity_field




! ============================================================================
! SORT_RADIAL_PARTICLE_TYPES
! Create spherically symmetric distribution with different particle types.
! (i.e. gas, icm and boundary particles).
! ============================================================================
SUBROUTINE sort_radial_particle_types(radgas,radicm,radboundary,r0)
  use interface_module, only : distance2
  use particle_module
  use hydro_module
  use type_module
  implicit none

  real(kind=PR), intent(in) :: radgas      ! Radial extent of gas particles
  real(kind=PR), intent(in) :: radicm      ! Radial extent of icm particles
  real(kind=PR), intent(in) :: radboundary ! Radial extent of boundary ptcls
  real(kind=PR), intent(in) :: r0(1:NDIM)  ! Position of origin

  integer :: p                                 ! Particle counter
  real(kind=PR) :: dr(1:NDIM)                  ! Relative position
  real(kind=PR) :: drsqd                       ! Distance squared
  real(kind=PR), allocatable :: rboundary(:,:) ! Positions of boundary ptcls
  real(kind=PR), allocatable :: rgas(:,:)      ! Positions of gas particles
  real(kind=PR), allocatable :: ricm(:,:)      ! Positions of icm particles

  write(6,*) "Sorting particles radially [sort_radial_particle_types]"

! Calculate average value of smoothing length and gas cloud radius  
! to know where to select boundary particles.
  pboundary = 0
  picm      = 0
  pgas      = 0
  pcdm      = 0
  pdust     = 0
  pion      = 0
  allocate(rgas(1:NDIM,1:ptot))
  allocate(ricm(1:NDIM,1:ptot))
  allocate(rboundary(1:NDIM,1:ptot))
  
! Now sort particles into gas particles and boundary particles
! ----------------------------------------------------------------------------
  do p=1,ptot
     call distance2(r0(1:NDIM),p,dr(1:NDIM),drsqd)
     if (drsqd <= radgas*radgas) then
        pgas = pgas + 1
        rgas(1:NDIM,pgas) = sph(p)%r(1:NDIM)
     else if (drsqd <= radicm*radicm) then
        picm = picm + 1
        ricm(1:NDIM,picm) = sph(p)%r(1:NDIM)
     else if (drsqd <= radboundary*radboundary) then
        pboundary = pboundary + 1
        rboundary(1:NDIM,pboundary) = sph(p)%r(1:NDIM)
     else
        write(6,*) "Particle outside range : ",sqrt(drsqd)
        stop
     end if
  end do
! ----------------------------------------------------------------------------
  
  write(6,*) "pboundary        :",pboundary
  write(6,*) "picm             :",picm
  write(6,*) "pgas             :",pgas
  write(6,*) "ptot             :",ptot
  
  ! Now reorder particle positions for boundary and gas particles
  if (pboundary > 0) then
     do p=1,pboundary
        sph(p)%r(1:NDIM) = rboundary(1:NDIM,p)
     end do
  end if
  if (picm > 0) then
     do p=1,picm
        sph(p+pboundary)%r(1:NDIM) = ricm(1:NDIM,p)
     end do
  end if
  if (pgas > 0) then
     do p=1,pgas
        sph(p+picm+pboundary)%r(1:NDIM) = rgas(1:NDIM,p)
     end do
  end if
  
  deallocate(rboundary)
  deallocate(ricm)
  deallocate(rgas)

  return
END SUBROUTINE sort_radial_particle_types




! ============================================================================
! SORT_PARTICLE_TYPES
! Sort all particles into the correct 'order' in the main arrays according 
! to their types given by the ptype array.
! ============================================================================
!SUBROUTINE sort_particle_types(ptype)
!  use particle_module, only : ptot
!  use hydro_module
!  use type_module
!  implicit none
!
!  integer, intent(in) :: ptype(1:ptot)         ! Particle types
!
!  integer :: boundaryslot                      ! Aux boundary ptcl counter
!  integer :: gasslot                           ! Aux gas ptcl counter
!  integer :: icmslot                           ! Aux icm ptcl counter
!  integer :: p                                 ! Particle counter
!  integer, allocatable :: porder(:)            ! New particle order
!
!  write(6,*) "Sorting particles types [sort_particle_types]"
!
!! Calculate average value of smoothing length and gas cloud radius  
!! to know where to select boundary particles.
!  pboundary = 0
!  picm      = 0
!  pgas      = 0
!  pcdm      = 0
!  pdust     = 0
!  pion      = 0
!  allocate(porder(1:ptot))
!  
!
!! Find how many particle there are of each type
!! ----------------------------------------------------------------------------
!  do p=1,ptot
!     if (ptype(p) == BOUNDARYID) then
!        pboundary = pboundary + 1
!     else if (ptype(p) == ICMID) then
!        picm = picm + 1
!     else if (ptype(p) == GASID) then
!        pgas = pgas + 1
!     else
!        write(6,*) "Unidentified particle type : ",p,ptype(p)
!        stop
!     end if
!  end do
!
!  write(6,*) "pboundary        :",pboundary
!  write(6,*) "picm             :",picm
!  write(6,*) "pgas             :",pgas
!  write(6,*) "ptot             :",ptot
!
!
!! Now make new ordered list of particles
!! ----------------------------------------------------------------------------
!  boundaryslot = 0
!  icmslot = pboundary
!  gasslot = pboundary + picm
!
!  do p=1,ptot
!     if (ptype(p) == BOUNDARYID) then
!        boundaryslot = boundaryslot + 1
!        porder(boundaryslot) = p
!     else if (ptype(p) == ICMID) then
!        icmslot = icmslot + 1
!        porder(icmslot) = p
!     else if (ptype(p) == GASID) then
!        gasslot = gasslot + 1
!        porder(gasslot) = p
!     end if
!  end do
!
!! Now reorder all arrays
!  call reorder_particle_arrays(1,ptot,porder(1:ptot))
!  
!  deallocate(porder)
!
!  return
!END SUBROUTINE sort_particle_types




! ============================================================================
! SOLVE_ISOTHERMAL_EQN
! Solve the isothermal Lane-Emden equation using a 4th-order Runge-Kutta 
! integrator.
! ============================================================================
SUBROUTINE solve_isothermal_eqn(nmax,delta_xi,xi_array,&
    &psi_array,phi_array,mu_array)
  use definitions
  implicit none

  integer, intent(in) :: nmax                         ! ..
  real(kind=PR), intent(in) :: delta_xi               ! ..
  real(kind=PR), intent(inout) :: mu_array(1:nmax)    ! ..
  real(kind=PR), intent(inout) :: phi_array(1:nmax)   ! ..
  real(kind=PR), intent(inout) :: psi_array(1:nmax)   ! ..
  real(kind=PR), intent(inout) :: xi_array(1:nmax)    ! ..

  integer :: i                       ! ..
  real(kind=PR) :: phi               ! ..
  real(kind=PR) :: psi               ! ..
  real(kind=PR) :: xi                ! ..
  real(kind=PR) :: k1_phi, k1_psi    ! 4th order Runge-Kutta values
  real(kind=PR) :: k2_phi, k2_psi    ! for phi, psi and theta
  real(kind=PR) :: k3_phi, k3_psi    ! 
  real(kind=PR) :: k4_phi, k4_psi    ! 

  debug2("Solving isothermal Lane-Emden equation [solve_isothermal_eqn]")

  ! Tabulate central values using boundary conditions
  xi_array(1)  = 0.0_PR
  psi_array(1) = 0.0_PR
  phi_array(1) = 0.0_PR
  mu_array(1)  = 0.0_PR
  
  ! Use first few terms of series solution for first step 
  ! (due to singularity in differential equation at xi = 0)
  xi  = delta_xi
  psi = (ONESIXTH)*(xi**2) - (1.0_PR/120.0_PR)*(xi**4) + &
       &(1.0_PR/1890.0_PR)*(xi**6)
  phi = (ONETHIRD)*(xi) - (1.0_PR/30.0_PR)*(xi**3) + (1.0_PR/315.0_PR)*(xi**5)
  xi_array(2)  = xi
  psi_array(2) = psi
  phi_array(2) = phi
  mu_array(2)  = phi*(xi**2)
  
  
  ! Now loop over all over integration points
  ! --------------------------------------------------------------------------
  do i=3,nmax

     ! Solve using 4th order Runge-Kutta method
     k1_phi = delta_xi*(exp(-psi) - 2*phi/xi)
     k1_psi = delta_xi*phi
     
     k2_phi = delta_xi*(exp(-psi - 0.5_PR*k1_psi) - &
          &2.0_PR*(phi + 0.5_PR*k1_phi)/(xi + delta_xi/2))
     k2_psi = delta_xi*(phi + k1_phi/2)
     
     k3_phi = delta_xi*(exp(-psi - 0.5_PR*k2_psi) - &
          &2.0_PR*(phi + 0.5_PR*k2_phi)/(xi + delta_xi/2))
     k3_psi = delta_xi*(phi + k2_phi/2)
     
     k4_phi = delta_xi*(exp(-psi - k3_psi) - &
          &2.0_PR*(phi + k3_phi)/(xi + delta_xi))
     k4_psi = delta_xi*(phi + k3_phi)
     
     phi = phi + ONESIXTH*(k1_phi + k4_phi) + ONETHIRD*(k2_phi + k3_phi)
     psi = psi + ONESIXTH*(k1_psi + k4_psi) + ONETHIRD*(k2_psi + k3_psi)
     
     xi = real(i - 1,PR)*delta_xi
     
     ! Tabulate values
     xi_array(i)  = xi
     psi_array(i) = psi
     phi_array(i) = phi
     mu_array(i)  = phi*(xi**2)
  end do
  ! --------------------------------------------------------------------------

  return
END SUBROUTINE solve_isothermal_eqn




! ============================================================================
! SOLVE_LANE_EMDEN_EQN
! Solve the general Lane-Emden equation using a 4th-order Runge-Kutta 
! integrator.
! ============================================================================
SUBROUTINE solve_lane_emden_eqn(ibound,nmax,delta_xi,&
     &npoly,xi_array,theta_array,phi_array,mu_array)
  use definitions
  implicit none

  integer, intent(inout) :: ibound                    ! ..
  integer, intent(in) :: nmax                         ! ..
  real(kind=PR), intent(in) :: delta_xi               ! ..
  real(kind=PR), intent(in) :: npoly                  ! ..
  real(kind=PR), intent(inout) :: mu_array(1:nmax)    ! ..
  real(kind=PR), intent(inout) :: phi_array(1:nmax)   ! ..
  real(kind=PR), intent(inout) :: theta_array(1:nmax) ! ..
  real(kind=PR), intent(inout) :: xi_array(1:nmax)    ! ..

  integer :: i                          ! ..
  real(kind=PR) :: phi                  ! ..
  real(kind=PR) :: theta                ! ..
  real(kind=PR) :: xi                   ! ..
  real(kind=PR) :: k1_phi, k1_theta     ! 4th order Runge-Kutta values
  real(kind=PR) :: k2_phi, k2_theta     ! for phi, psi and theta
  real(kind=PR) :: k3_phi, k3_theta     ! 
  real(kind=PR) :: k4_phi, k4_theta     ! 

  debug2("Solving Lane-Emden equation for polytrope [solve_lane_emden_eqn]")

! Tabulate central values using boundary conditions
  mu_array(1)    = 0.0_PR
  phi_array(1)   = 0.0_PR
  theta_array(1) = 1.0_PR
  xi_array(1)    = 0.0_PR
  
! Use first few terms of series solution for first step 
! (due to singularity in differential equation at xi = 0)
  xi             = delta_xi
  theta          = 1.0_PR - (ONESIXTH)*(xi**2) + (npoly/120.0_PR)*(xi**4)
  phi            = (-ONETHIRD)*(xi) + (npoly/30.0_PR)*(xi**3)
  phi_array(2)   = phi
  mu_array(2)    = -phi*(xi**2)
  theta_array(2) = theta
  xi_array(2)    = xi
  
  
! Now loop over all over integration points
! ----------------------------------------------------------------------------
  do i=3,nmax
     ibound = i - 1
     
     ! Solve using 4th order Runge-Kutta method
     k1_phi   = delta_xi*(-theta**(npoly) - 2.0_PR*phi/xi)
     k1_theta = delta_xi*phi
     
     k2_phi   = delta_xi*(-(theta + 0.5_PR*k1_theta)**(npoly) - &
          &2.0_PR*(phi + 0.5_PR*k1_phi)/(xi + 0.5_PR*delta_xi))
     k2_theta = delta_xi*(phi + 0.5_PR*k1_phi)
     
     k3_phi   = delta_xi*(-(theta + 0.5_PR*k2_theta)**(npoly) - &
          &2.0_PR*(phi + 0.5_PR*k2_phi)/(xi + 0.5_PR*delta_xi))
     k3_theta = delta_xi*(phi + 0.5_PR*k2_phi)
     
     k4_phi   = delta_xi*(-(theta + k3_theta)**(npoly) - &
          &2.0_PR*(phi + k3_phi)/(xi + delta_xi))
     k4_theta = delta_xi*(phi + k3_phi)
     
     phi   = phi + ONESIXTH*(k1_phi + k4_phi) + ONETHIRD*(k2_phi + k3_phi)
     theta = theta + ONESIXTH*(k1_theta + k4_theta) + &
          &ONETHIRD*(k2_theta + k3_theta)
     
     ! Tabulate values
     xi             = real(i - 1,PR)*delta_xi
     xi_array(i)    = xi
     theta_array(i) = theta
     phi_array(i)   = phi
     mu_array(i)    = -phi*(xi**2)
     if (.not.(theta > 0.0_PR .and. theta < BIG_NUMBER)) exit
     
  end do
! ----------------------------------------------------------------------------
  
  return
END SUBROUTINE solve_lane_emden_eqn




! ============================================================================
! ADD_UNIFORM_3D_LATTICE
! Add a uniform 2-D grid of nx*ny particles to specified domain 
! ([xmin,ymin] to [xmax,ymax]).  Particles are placed in data arrays 
! starting at p=pfirst.
! ============================================================================
SUBROUTINE add_uniform_2d_lattice(xmin,xmax,ymin,ymax,raux,plat,nx,ny)
  use definitions
  implicit none

  integer, intent(in) :: nx              ! No. of ptcl layers in x-dimension
  integer, intent(in) :: ny              ! No. of ptcl layers in x-dimension
  integer, intent(in) :: plat            ! ..
  real(kind=PR), intent(in) :: xmin      ! LHS extent of grid
  real(kind=PR), intent(in) :: xmax      ! ..
  real(kind=PR), intent(in) :: ymin      ! ..
  real(kind=PR), intent(in) :: ymax      ! ..
  real(kind=PR), intent(out) :: raux(1:NDIM,1:plat)  ! ..

  integer :: i                           ! ..
  integer :: j                           ! ..
  integer :: p                           ! ..
  
  do i=0,nx-1
     do j=0,ny-1
        p = j*nx + i + 1
        if (p > plat) then
           write(6,*) 'p > plat : ',p,plat,nx,ny,i,j,nx*ny
           stop
        end if
        raux(1,p) = xmin + &
             &(real(i,PR) + 0.5_PR)*(xmax - xmin)/real(nx,PR)
        raux(2,p) = ymin + &
             &(real(j,PR) + 0.5_PR)*(ymax - ymin)/real(ny,PR)
     end do
  end do

  return
END SUBROUTINE add_uniform_2d_lattice




! ============================================================================
! ADD_UNIFORM_3D_LATTICE
! Add a uniform 2-D grid of nx*ny particles to specified domain 
! ([xmin,ymin] to [xmax,ymax]).  Particles are placed in data arrays 
! starting at p=pfirst.
! ============================================================================
SUBROUTINE add_uniform_3d_lattice(xmin,xmax,ymin,ymax,zmin,zmax,&
     &raux,plat,nx,ny,nz)
  use definitions
  implicit none

  integer, intent(in) :: nx              ! No. of ptcl layers in x-dimension
  integer, intent(in) :: ny              ! No. of ptcl layers in x-dimension
  integer, intent(in) :: nz              ! No. of ptcl layers in x-dimension
  integer, intent(in) :: plat            ! ..
  real(kind=PR), intent(in) :: xmin      ! LHS extent of grid
  real(kind=PR), intent(in) :: xmax      ! ..
  real(kind=PR), intent(in) :: ymin      ! ..
  real(kind=PR), intent(in) :: ymax      ! ..
  real(kind=PR), intent(in) :: zmin      ! ..
  real(kind=PR), intent(in) :: zmax      ! ..
  real(kind=PR), intent(out) :: raux(1:NDIM,1:plat)  ! ..

  integer :: i                           ! ..
  integer :: j                           ! ..
  integer :: k                           ! ..
  integer :: p                           ! ..
  
  do i=0,nx-1
     do j=0,ny-1
        do k=0,nz-1
           p = k*nx*ny + j*nx + i + 1
           if (p > plat) then
              write(6,*) 'p > plat : ',p,plat,nx,ny,nz,i,j,k,nx*ny*nz
              stop
           end if
           raux(1,p) = xmin + &
                &(real(i,PR) + 0.5_PR)*(xmax - xmin)/real(nx,PR)
           raux(2,p) = ymin + &
                &(real(j,PR) + 0.5_PR)*(ymax - ymin)/real(ny,PR)
           raux(3,p) = zmin + &
                &(real(k,PR) + 0.5_PR)*(zmax - zmin)/real(nz,PR)
        end do
     end do
  end do

  return
END SUBROUTINE add_uniform_3d_lattice




! ============================================================================
! ADD_UNIFORM_2D_CIRCLE_LATTICE
! Add a uniform 2-D grid of nx*ny particles to specified domain 
! ([xmin,ymin] to [xmax,ymax]).  Particles are placed in data arrays 
! starting at p=pfirst.
! ============================================================================
SUBROUTINE add_uniform_2d_circle_lattice(x0,y0,rcirc,raux,pcirc,plat,nx,ny)
  use definitions
  implicit none

  integer, intent(in) :: nx              ! No. of ptcl layers in x-dimension
  integer, intent(in) :: ny              ! No. of ptcl layers in x-dimension
  integer, intent(in) :: plat            ! ..
  integer, intent(out) :: pcirc          ! ..
  real(kind=PR), intent(in) :: x0        ! ..
  real(kind=PR), intent(in) :: y0        ! ..
  real(kind=PR), intent(in) :: rcirc     ! ..
  real(kind=PR), intent(out) :: raux(1:NDIM,1:plat)  ! ..

  integer :: i                           ! ..
  integer :: j                           ! ..
  integer :: p                           ! ..
  real(kind=PR) :: dr(1:NDIM)            ! ..
  real(kind=PR) :: xmin                  ! LHS extent of grid
  real(kind=PR) :: xmax                  ! ..
  real(kind=PR) :: ymin                  ! ..
  real(kind=PR) :: ymax                  ! ..

  pcirc = 0
  xmin = x0 - rcirc
  xmax = x0 + rcirc
  ymin = y0 - rcirc
  ymax = y0 + rcirc
  
  do i=0,nx-1
     do j=0,ny-1
        pcirc = pcirc + 1
        raux(1,p) = xmin + &
             &(real(i,PR) + 0.5_PR)*(xmax - xmin)/real(nx,PR)
        raux(2,p) = ymin + &
             &(real(j,PR) + 0.5_PR)*(ymax - ymin)/real(ny,PR)
        dr(1) = raux(1,p) - x0; dr(2) = raux(2,p) - y0
        if (dr(1)*dr(1) + dr(2)*dr(2) >= rcirc*rcirc) pcirc = pcirc - 1
     end do
  end do

  return
END SUBROUTINE add_uniform_2d_circle_lattice



! ============================================================================
! ADD_UNIFORM_3D_SPHERE_LATTICE
! Add a uniform 2-D grid of nx*ny particles to specified domain 
! ([xmin,ymin] to [xmax,ymax]).  Particles are placed in data arrays 
! starting at p=pfirst.
! ============================================================================
SUBROUTINE add_uniform_3d_sphere_lattice(x0,y0,z0,rcirc,raux,&
     &pcirc,plat,nx,ny,nz)
  use definitions
  implicit none

  integer, intent(in) :: nx              ! No. of ptcl layers in x-dimension
  integer, intent(in) :: ny              ! No. of ptcl layers in x-dimension
  integer, intent(in) :: nz              ! No. of ptcl layers in x-dimension
  integer, intent(in) :: plat            ! ..
  integer, intent(out) :: pcirc          ! ..
  real(kind=PR), intent(in) :: x0        ! ..
  real(kind=PR), intent(in) :: y0        ! ..
  real(kind=PR), intent(in) :: z0        ! ..
  real(kind=PR), intent(in) :: rcirc     ! ..
  real(kind=PR), intent(out) :: raux(1:NDIM,1:plat)  ! ..

  integer :: i                           ! ..
  integer :: j                           ! ..
  integer :: k                           ! ..
  integer :: p                           ! ..
  real(kind=PR) :: dr(1:NDIM)            ! ..
  real(kind=PR) :: xmin                  ! LHS extent of grid
  real(kind=PR) :: xmax                  ! ..
  real(kind=PR) :: ymin                  ! ..
  real(kind=PR) :: ymax                  ! ..
  real(kind=PR) :: zmin                  ! ..
  real(kind=PR) :: zmax                  ! ..

  pcirc = 0
  xmin = x0 - rcirc
  xmax = x0 + rcirc
  ymin = y0 - rcirc
  ymax = y0 + rcirc
  zmin = z0 - rcirc
  zmax = z0 + rcirc
  
  do i=0,nx-1
     do j=0,ny-1
        do k=0,nz-1
           pcirc = pcirc + 1
           raux(1,pcirc) = xmin + &
                &(real(i,PR) + 0.5_PR)*(xmax - xmin)/real(nx,PR)
           raux(2,pcirc) = ymin + &
                &(real(j,PR) + 0.5_PR)*(ymax - ymin)/real(ny,PR)
           raux(3,pcirc) = zmin + &
                &(real(k,PR) + 0.5_PR)*(zmax - zmin)/real(nz,PR)
           dr(1) = raux(1,pcirc) - x0
           dr(2) = raux(2,pcirc) - y0
           dr(3) = raux(3,pcirc) - z0
           if (dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) >= rcirc*rcirc) &
                &pcirc = pcirc - 1
        end do
     end do
  end do

  return
END SUBROUTINE add_uniform_3d_sphere_lattice



! ============================================================================
! SOLVE_BONDI_ACCRETION_EQN
! '''DOESN'T WORK!!'''
! ============================================================================
SUBROUTINE solve_bondi_accretion_eqn(nmax,gamma,fraclambda,xmax,&
     &mu_array,u_array,x_array,y_array,z_array)
  use definitions
  implicit none

  integer, intent(in) :: nmax
  real(kind=PR), intent(in) :: gamma
  real(kind=PR), intent(in) :: fraclambda
  real(kind=PR), intent(in) :: xmax
  real(kind=PR), intent(inout) :: mu_array(1:nmax)
  real(kind=PR), intent(inout) :: u_array(1:nmax)
  real(kind=PR), intent(inout) :: x_array(1:nmax)
  real(kind=PR), intent(inout) :: y_array(1:nmax)
  real(kind=PR), intent(inout) :: z_array(1:nmax)

  integer :: i
  real(kind=PR) :: fu
  real(kind=PR) :: lambda
  real(kind=PR) :: lambdaprime
  real(kind=PR) :: gfactor
  real(kind=PR) :: gx
  real(kind=PR) :: u
  real(kind=PR) :: uold
  real(kind=PR) :: x

  write(6,*) "Solving Bondi problem : ",nmax,xmax

  i = 1
  x = xmax/real(nmax,PR)

  ! Calculate value for Bondi accretion factor depending on value of gamma
  if (abs(1.6666666666666_PR - gamma) < 0.000001_PR) then
     lambda = 0.25_PR*fraclambda
     gfactor = 0.25_PR
     lambdaprime = lambda**(-2.0_PR*gfactor)

     ! Calculate properties in first table region
     u = 1.0_PR
     y_array(i) = 0.5_PR*sqrt(2.0_PR/x)
     z_array(i) = 0.125_PR*(0.5_PR*x)**(1.5_PR)
     mu_array(i) = ONETHIRD*(2.0_PR*x)**(1.5_PR)

  else if (abs(1.0_PR - gamma) < 0.000001_PR) then
     lambda = 0.25_PR*exp(1.5_PR)*fraclambda
     gfactor = 0.0_PR
     lambdaprime = 1.0_PR

     ! Calculate properties in first table region
     u = (4.0_PR/lambda)**(0.5_PR*(gamma - 1.0_PR))*&
          &(0.5_PR*x)**(0.75_PR*(gamma - 1.0_PR) - 0.5_PR)
     y_array(i) = sqrt(2.0_PR/x)
     z_array(i) = 0.25_PR*lambda*(0.5_PR*x)**(-1.5_PR)
     mu_array(i) = TWOTHIRDS*lambda*(2.0_PR*x)**(1.5_PR)

  else
     lambda = 0.25_PR*fraclambda*(2.0_PR/(5.0_PR - 3.0_PR*gamma))**&
          &((5.0_PR - 3.0_PR*gamma)/(2.0_PR*gamma - 2.0_PR))
     gfactor = (gamma - 1.0_PR)/(gamma + 1.0_PR)
     lambdaprime = lambda**(-2.0_PR*gfactor)

     ! Calculate properties in first table region
     u = (4.0_PR/lambda)**(0.5_PR*(gamma - 1.0_PR))*&
          &(0.5_PR*x)**(0.75_PR*(gamma - 1.0_PR) - 0.5_PR)
     y_array(i) = sqrt(2.0_PR/x)
     z_array(i) = 0.25_PR*lambda*(0.5_PR*x)**(-1.5_PR)
     mu_array(i) = TWOTHIRDS*lambda*(2.0_PR*x)**(1.5_PR)
  end if

  x_array(i) = x
  u_array(i) = u


  write(6,*) "STUFF : ",lambda,gfactor,lambdaprime
  write(6,*) "BONDI : ",i,x,u,y_array(i),z_array(i),mu_array(i)
  
  
  ! Calculate solutions for all other table elements
  ! --------------------------------------------------------------------------
  do i=2,nmax

     x = real(i,PR)*xmax/real(nmax,PR)
     !fu = lambdaprime*((x**(4.0*gfactor))/(gamma - 1.0_PR) + &
     !     &x**((3.0_PR*gamma - 5.0_PR)/(gamma + 1.0_PR)))

     fu = lambdaprime*((x**(4.0_PR*gfactor))*(1.0_PR/x + 1.0_PR/(gamma - 1.0_PR)))

     u = 0.5_PR*u

     ! Iterate u towards consistent solution
     do 
        uold = u
        !u = (2.0_PR*(fu - (u**(-2.0_PR*gfactor))/(gamma - 1.0_PR)))**&
        !     &(0.25_PR*(gamma + 1.0_PR))
        !u = ((gamma - 1.0_PR)*(fu - 0.5_PR*u**(4.0_PR/(gamma + 1.0_PR))))&
        !     &**(-0.5_PR/gfactor)
        u = (fu/(0.5_PR + 1.0_PR/(gamma - 1.0_PR)/(u**2)))**(0.25_PR*(gamma + 1.0_PR))
        write(6,*) "Iterating .... ",u,uold,abs((u - uold)/uold)
        !read(5,*)
        if (abs((u - uold)/uold) < 1.0E-10) exit
     end do

     write(6,*) "Complete? : ",0.5_PR*u**(4.0_PR/(gamma + 1.0_PR)) + &
          &u**(-2.0_PR*gfactor)/(gamma - 1.0_PR),fu

     ! Once consistent solution is obtained, record values in tables
     x_array(i) = x
     u_array(i) = u
     y_array(i) = (u**(2.0_PR/(gamma - 1.0_PR)))*(lambda/x**2)**(gfactor)
     z_array(i) = (lambda/(u*x*x))**(2.0_PR/(gamma - 1.0_PR))

     mu_array(i) = mu_array(i-1) + 0.5_PR*(x_array(i-1) - x_array(i))*&
          &(z_array(i) + z_array(i-1))

     write(6,*) "BONDI    : ",i,x,u,y_array(i),z_array(i),mu_array(i)
     if (x < 0.2_PR) write(6,*) "COMPARE? : ",i,x,&
          &(4.0_PR/lambda)**(0.5_PR*(gamma - 1.0_PR))*&
          &(0.5_PR*x)**(0.75_PR*(gamma - 1.0_PR) - 0.5_PR),&
          &sqrt(2.0_PR/x),0.25_PR*lambda*(0.5_PR*x)**(-1.5_PR),&
          &TWOTHIRDS*lambda*(2.0_PR*x)**(1.5_PR)
     read(5,*)


  end do
  ! --------------------------------------------------------------------------


  return
END SUBROUTINE solve_bondi_accretion_eqn
