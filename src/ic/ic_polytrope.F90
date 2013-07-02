! IC_POLYTROPE.F90
! D. A. Hubber - 30/09/2008
! Solves Lane-Emden equation (or isothermal equation) to obtain the
! Lane-Emden functions for various Polytropes.  Next, reads in a uniform 
! density sphere and stretches particles radially and rescales masses to 
! reproduce required polytropic density distribution.
! Reads in a parameters file 'polytrope.dat' with the following data : 
! in_file       : Filename of relaxed unit-sphere
! in_file_form  : Format of input sphere file 
!                 (e.g. dragon_form, dragon_unform, seren_form, seren_unform)
! out_file      : Output filename
! out_file_form : Format of filename 
!                 (e.g. dragon_form, dragon_unform, seren_form, seren_unform)
! isocloud      : Isothermal eqn (.TRUE.) or Lane-Emden eqn (.FALSE.)
! etapoly       : Polytropic index
! xi_bound      : Dimensionless cloud boundary
! delta_xi      : Integration step for integrating Lane-Emden Eqn.
! mpoly         : Mass of polytrope
! munit         : Mass unit
! rho0          : Central density of polytrope
! rhounit       : Density unit
! mflag         : 'rho0' or 'mass' : Set either mass or central rho
! Kpoly         : Kpoly (a_0^2 for isothermal option)
! vunit         : Velocity unit
! menvelope     : Mass of gas envelope around polytrope
! micmenv       : Mass of IcM envelope around polytrope/gas envelope.
! hboundary     : Size of static boundary zone (in units of mean h)
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_polytrope
  use interface_module, only : distance2,paramerror,&
       &paramstore,read_data,write_data
  use particle_module
  use filename_module
  use scaling_module
  use time_module
  use type_module
  use neighbour_module
  use hydro_module
  implicit none
  
  character(len=256) :: out_file    ! Name of output file
  character(len=50) :: mflag        ! Set 'mass' of 'rho0' in options files
  character(len=256) :: store_file  ! parameters output
  logical :: isocloud               ! Flag if isothermal cloud or not
  integer :: i                      ! Integration counter
  integer :: ibound                 ! Table element of boundary
  integer :: ip                     ! ??
  integer :: nmax                   ! No. of integration points
  integer :: p                      ! Particle counter
  real(kind=PR) :: a0sqd            ! Isothermal sound speed squared
  real(kind=PR) :: delta_xi         ! Step-size of xi
  real(kind=PR) :: dr(1:NDIM)       ! Relative position
  real(kind=PR) :: drmag            ! Distance
  real(kind=PR) :: drsqd            ! Distance squared
  real(kind=PR) :: etapoly          ! Polytropic exponent/index (??)
  real(kind=PR) :: hboundary        ! Thickness of boundary region (in h)
  real(kind=PR) :: hmean            ! Mean smoothing length for uniform sphere
  real(kind=PR) :: mpoly            ! Mass of polytrope
  real(kind=PR) :: menvelope        ! Mass of envelope
  real(kind=PR) :: micmenv          ! Mass of IcM envelope
  real(kind=PR) :: mint             ! Mass interior to radius
  real(kind=PR) :: mp               ! Mass of particle p
  real(kind=PR) :: mu               ! Dimensionless mass
  real(kind=PR) :: npoly            ! Polytropic index
  real(kind=PR) :: radp             ! Aux. radius variable
  real(kind=PR) :: rboundary        ! Radius containing boundary particles
  real(kind=PR) :: rcentre(1:NDIM)  ! Centre of polytrope
  real(kind=PR) :: rcloud           ! Radius containing polytrope particles
  real(kind=PR) :: remainder        ! Remainder of interpolation
  real(kind=PR) :: renvelope        ! Radius containing gas envelope
  real(kind=PR) :: rgas             ! Radius (of initial sphere) containing gas
  real(kind=PR) :: rho0             ! Central density
  real(kind=PR) :: rhofrac          ! Density contrast of external medium
  real(kind=PR) :: ricm             ! Radius containing IcM particles
  real(kind=PR) :: rpoly            ! Physical radius of polytrope
  real(kind=PR) :: R0               ! Scaling radius of polytrope
  real(kind=PR) :: xi               ! Dimensionless distance
  real(kind=PR) :: xi_bound         ! Dimensionless radius/boundary
  real(kind=PR), allocatable :: mass_array(:)  ! Interior mass table
  real(kind=PR), allocatable :: density(:)     ! Polytrope density table
  real(kind=PR), allocatable :: mu_array(:)    ! Dimensionless mass table
  real(kind=PR), allocatable :: phi_array(:)   ! ?? table
  real(kind=PR), allocatable :: pressure(:)    ! Pressure table
  real(kind=PR), allocatable :: psi_array(:)   ! ?? table 
  real(kind=PR), allocatable :: theta_array(:) ! ?? table
  real(kind=PR), allocatable :: xi_array(:)    ! Dimensionless radius table
  
  debug1("Create polytrope density distribution [ic_polytrope.F90]")
  
! Read parameters file
  call default_parameters
  
! Read polytrope set-up file  
  open(unit=1,file='polytrope.dat',status='old')
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) in_file      
  read(1,*) in_file_form 
  read(1,*) out_file     
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) isocloud
  read(1,*) etapoly        
  read(1,*) xi_bound 
  read(1,*) delta_xi
  read(1,*) mpoly
  read(1,*) munit
  read(1,*) rho0
  read(1,*) rhounit
  read(1,*) mflag
  read(1,*) Kpoly
  read(1,*) vunit
  read(1,*) menvelope
  read(1,*) micmenv
  read(1,*) hboundary
  read(1,*) rhofrac
  close(1)

! Check parameter file is consistent with chosen EOS.
  if (xi_bound < 0.0_PR) call paramerror("xi_bound < 0")
  if (mflag /= "mass" .and. mflag /= "rho0") &
       &call paramerror("Unrecognised option for mflag")

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types

  in_file       = trim(adjustl(in_file))
  in_file_form  = trim(adjustl(in_file_form))
  out_file      = trim(adjustl(out_file))
  out_file_form = trim(adjustl(out_file_form))
  restart       = .false.
  run_id        = ''
  rcentre(1:3)  = 0.0_PR

! Calculate scaling variables
  call units

! Add 2 for round-off error and end points
  nmax   = int(xi_bound/delta_xi) + 2
  ibound = 1
  write(6,*) "nmax : ",nmax
  write(6,*) "Mass : ",mpoly

! Allocate arrays now number of array points is known
  allocate(xi_array(1:nmax))
  allocate(psi_array(1:nmax))
  allocate(phi_array(1:nmax))
  allocate(mu_array(1:nmax))
  allocate(pressure(1:nmax))
  allocate(density(1:nmax))
  allocate(theta_array(1:nmax))
  allocate(mass_array(1:nmax))
  ibound = nmax

  menvelope = menvelope / mscale
  mpoly     = mpoly / mscale
  micmenv   = micmenv / mscale

  !gamma = etapoly
  if (etapoly /= 0.0_PR) npoly = 1.0_PR / (etapoly - 1.0_PR)


! Isothermal sphere solution
! ----------------------------------------------------------------------------
  if (isocloud) then     
    call solve_isothermal_eqn(nmax,delta_xi,&
          &xi_array,psi_array,phi_array,mu_array)

!    mu = mu_array(i) + (mu_array(i+1) - mu_array(i))*remainder
    mu    = mu_array(nmax)
    a0sqd = Kpoly / (vscale*vscale)
    if (mflag == 'mass') then
       mpoly = mpoly / mscale
       rho0  = a0sqd*a0sqd*a0sqd*mu*mu / (4.0_PR*PI*mpoly*mpoly)
    else
       rho0  = rho0 / rhoscale
       mpoly = a0sqd*mu*sqrt(a0sqd/(4.0_PR*PI*rho0))
    end if
    R0 = sqrt(a0sqd / (4.0_PR*PI*rho0))
    rpoly = xi_bound*R0

    ! Scale dimensionless solution to physical units
    do i=1,nmax
       density(i)    = rho0*exp(-psi_array(i))
       mass_array(i) = mpoly*mu_array(i)/mu
       pressure(i)   = a0sqd*density(i)
       if (density(i) < 0.0_PR .and. i > ibound) ibound = i
    end do

! General polytrope solution
! ----------------------------------------------------------------------------
  else
     call solve_lane_emden_eqn(ibound,nmax,delta_xi,&
          &npoly,xi_array,theta_array,phi_array,mu_array)

     etapoly = 1.0_PR + 1.0_PR/npoly
     xi      = xi_array(ibound)
     mu      = mu_array(ibound)

     if (mflag == 'mass') then
        mpoly = mpoly / mscale
        rho0  = (4.0_PR*PI*(Kpoly*(npoly + 1.0_PR)/(4.0_PR*PI))**&
             &(1.5_PR)*mu/mpoly)**(2.0_PR*npoly/(npoly - 3.0_PR))
     else
        rho0  = rho0 / rhoscale
        !mpoly = 1.0_PR/(4.0_PR*PI*(Kpoly*(npoly + 1.0_PR)/(4.0_PR*PI))**&
        !     &(1.5_PR)/mu)*rho0**((npoly - 3.0_PR)/2.0_PR*npoly)
        mpoly = 4.0_PR*PI*(Kpoly*(npoly + 1.0_PR)/(4.0_PR*PI))**&
             &(1.5_PR)*rho0**((3.0_PR - npoly)/2.0_PR*npoly)*mu
     end if

     R0 = sqrt((Kpoly*(npoly + 1)*rho0**(1.0_PR/npoly - 1.0_PR)) / &
          &(4.0_PR*PI))
     rpoly  = xi*R0
     write(6,*) "R0 : ",R0,"   xi : ",xi,"    rpoly : ",rpoly
     write(6,*) "mu : ",mpoly,sqrt(TWOPI)
  
     ! Scale dimensionless solution to physical units
     do i=1,ibound
        density(i)    = rho0*(theta_array(i)**(npoly))
        mass_array(i) = mpoly*mu_array(i)/mu
        pressure(i)   = Kpoly*(density(i)**(etapoly))
     end do
     do i=ibound+1,nmax
        density(i)    = 0.0_PR
        mass_array(i) = mpoly
        pressure(i)   = 0.0_PR
     end do

  end if
! ----------------------------------------------------------------------------

  write(6,*) "R0                :",R0*rscale,runit
  write(6,*) "Radius of cloud   :",rpoly*rscale,runit
  write(6,*) "Mass of cloud     :",mpoly*mscale,munit
  write(6,*) "Central density   :",rho0*rhoscale,rhounit
  write(6,*) "External pressure :",pressure(nmax)*Pscale,Punit

! Load in uniform density sphere file
  write(6,*) "in_file : ",trim(in_file)
  call read_data(in_file,in_file_form)

! Writing compiler flags and parameters to file
  store_file = trim("poly.params")
  call paramstore(store_file)

! Sort particles into gas, icm and boundary zones
  hmean     = (pp_gather/8.0_PR/real(ptot,PR))**(ONETHIRD)
  rboundary = 1.0000001_PR
  ricm      = rboundary - hboundary*hmean
  rgas      = ricm*((mpoly + menvelope)/&
       &(mpoly + menvelope + micmenv))**(ONETHIRD)
  rcloud    = ricm*(mpoly/(mpoly + menvelope + micmenv))**(ONETHIRD)

  write(6,*) "hmean     : ",hmean
  write(6,*) "rgas      : ",rgas
  write(6,*) "ricm      : ",ricm
  write(6,*) "rboundary : ",rboundary
  write(6,*) "rcloud    : ",rcloud

  if (ricm < 0.0_PR .or. rgas < 0.0_PR) then
     stop 'Negative ricm or rgas'
  else if (rboundary - rgas < SMALL_NUMBER) then
     write(6,*) "All gas particles : No radial ordering needed"
     pgas = ptot
     picm = 0
     pcdm = 0
     pboundary = 0
     pion = 0
     pdust = 0
  else
     call sort_radial_particle_types(rgas,ricm,rboundary,rcentre(1:NDIM))
  end if

! Setting up for different particle types
  call types

#if defined(INTERNAL_ENERGY)
  typeinfo(gasid)%eos = "energy_eqn"
  typeinfo(boundaryid)%eos = "energy_eqn"
  typeinfo(icmid)%eos = "energy_eqn"
#endif


! Now stretch all particles and rescale masses to reproduce the density 
! distribution of the polytrope (and the boundary envelope).
! ----------------------------------------------------------------------------
  mp = (mpoly + menvelope) / real(pgas,PR)

  do p=1,ptot
     call distance2(rcentre(1:NDIM),p,dr(1:NDIM),drsqd)
     drmag = sqrt(drsqd) + SMALL_NUMBER
     if (drmag >= rcloud) then
        radp = (3.0_PR*mpoly)/(4.0_PR*PI*density(nmax))&
             &*((drmag/rcloud)**3 - 1.0_PR) + rpoly**3
        radp = radp**(ONETHIRD)
        sph(p)%r(1:NDIM) = sph(p)%r(1:NDIM)*radp/drmag
        sph(p)%m = rhofrac*mp
#if defined(INTERNAL_ENERGY)
        sph(p)%u = 1.0_PR/rhofrac
#endif
     else
        mint = mpoly*(drmag/rcloud)**3
        ip = 1
        do i=1,ibound
           ip = i
           if (mass_array(i) >= mint) exit
        end do
        remainder = (mint - mass_array(ip))/(mass_array(ip-1) - mass_array(ip))
        xi = xi_array(ip - 1) + remainder*delta_xi
        radp = R0*xi
        sph(p)%r(1:NDIM) = sph(p)%r(1:NDIM)*radp/drmag
        sph(p)%m = mp
#if defined(INTERNAL_ENERGY)
        sph(p)%u = 1.0_PR
#endif
     end if
  end do

! Now set other variables
  do p=1,ptot
     sph(p)%v(1:VDIM) = 0.0_PR
  end do

! Initialising kernel tables
#if defined(KERNEL_TABLES)
  call tabulate_kernel_functions
#endif

! Initialize variables
  call initialize_seren_variables_2
  call initialize_sph_variables_1

! Build and stock tree
  call tree_update(nbuild,nstock)

! Make initial guesses of h either using tree or a global estimate
#if defined(BH_TREE) && !defined(CONSTANT_H)
  if (.not. restart) call BHhydro_hguess
#elif !defined(CONSTANT_H)
  if (.not. restart) call h_guess
#endif

! Calculating initial SPH quantities
  call sph_update

! Calculate thermal properties
  do p=1,ptot
     sph(p)%press = Kpoly*sph(p)%rho**gamma
     sph(p)%temp = sph(p)%press/sph(p)%rho/Pconst
#if defined(INTERNAL_ENERGY)
     if (isocloud) then
        sph(p)%u = 1.5_PR*sph(p)%u*sph(p)%press/sph(p)%rho
     else
        sph(p)%u = sph(p)%u*sph(p)%press/sph(p)%rho/(gamma - 1.0_PR)
     end if
#endif
  end do
  call initialize_thermal_properties

! Calculate hydro forces on all SPH particles
#if defined(HYDRO)
  call sph_hydro_forces
#endif

! Calculate gravitational forces on all SPH particles
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  call sph_grav_forces
#endif

! Initialize other variables
  call initialize_sph_variables_2

! Write tables to file
  open(unit=1,file='polysolution.dat',status='unknown')
  do i=1,ibound
     write (1,'(9E15.7)') xi_array(i),psi_array(i),phi_array(i),&
          &mu_array(i),theta_array(i),xi_array(i)*R0*rscale,&
          &density(i)*rhoscale*rhocgs,pressure(i)*Pscale,&
          &mu_array(i)*mpoly*mscale/mu
  end do
  close(1)

! Write particle data to file
  out_file = trim(adjustl(out_file))
  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  call write_data_debug("ICPOLY.debug.dat",rcentre(1:NDIM))
#endif

  renvelope = ((3.0_PR*mpoly)/(4.0_PR*PI*density(nmax))&
       &*((rgas/rcloud)**3 - 1.0_PR) + rpoly**3)**(ONETHIRD)
  ricm      = ((3.0_PR*mpoly)/(4.0_PR*PI*density(nmax))&
       &*((ricm/rcloud)**3 - 1.0_PR) + rpoly**3)**(ONETHIRD)
  rboundary = ((3.0_PR*mpoly)/(4.0_PR*PI*density(nmax))&
       &*((rboundary/rcloud)**3 - 1.0_PR) + rpoly**3)**(ONETHIRD)

  write(6,*) "rpoly     : ",rpoly*rscale,runit
  write(6,*) "renvelope : ",renvelope*rscale,runit
  write(6,*) "ricm      : ",ricm*rscale,runit
  write(6,*) "rboundary : ",rboundary*rscale,runit

! Clean-up all memory
  call clean_up
  deallocate(mass_array)
  deallocate(theta_array)
  deallocate(density)
  deallocate(pressure)
  deallocate(mu_array)
  deallocate(phi_array)
  deallocate(psi_array)
  deallocate(xi_array)
  
  stop
END PROGRAM ic_polytrope
