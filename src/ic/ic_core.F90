! IC_CORE.F90
! D. A. Hubber - 30/09/2008
! ..
! Reads in a parameters file 'polytrope.dat' with the following data : 
! in_file       : Filename of relaxed unit-sphere
! in_file_form  : Format of input sphere file 
! out_file      : Output filename
! out_file_form : Format of output file
! runit         : Length unit
! munit         : Mass unit
! rhounit       : Density unit
! vunit         : Velocity unit
! rho0          : 'core' density (for 'powerlaw', rho0 = density at r=rcore)
! rcore         : Radius of 'core' (for 'powerlaw', this is a reference 
!                                   radius to scale the density with rho0)
! radmax        : Maximum extent of cloud
! pofr          : Radius powerlaw exponent
! coretype      : Type of density profile added to core
!                 plummer : Plummer-like sphere
!                 powerlaw : radial-power law density profile
! menvelope     : Mass of gas envelope around polytrope
! micm          : Mass of IcM region around gas
! hboundary     : Size of static boundary zone (in units of mean h)
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_core
  use interface_module, only : distance2,paramerror,read_data,write_data
  use hydro_module, only : isotemp
  use particle_module
  use filename_module
  use scaling_module
  use time_module
  use type_module
  use neighbour_module
  use constant_module
  implicit none
  
  character(len=256) :: out_file    ! Name of output file
  character(len=50) :: coretype     ! Set 'mass' of 'rho0' in options files
  integer :: i                      ! Integration counter
  integer :: ip                     ! ..
  integer :: nmax                   ! No. of integration points
  integer :: p                      ! Particle counter
  real(kind=PR) :: dr(1:NDIM)       ! Relative position
  real(kind=PR) :: drmag            ! Distance
  real(kind=PR) :: drsqd            ! Distance squared
  real(kind=PR) :: hboundary        ! Thickness of boundary region (in h)
  real(kind=PR) :: hicm             ! Thickness of ICM region (in h)
  real(kind=PR) :: hmean            ! Mean smoothing length for uniform sphere
  real(kind=PR) :: menvelope        ! Mass of envelope
  real(kind=PR) :: micmenv          ! Mass of IcM region
  real(kind=PR) :: mint             ! Mass interior to radius
  real(kind=PR) :: mp               ! Mass of particle p
  real(kind=PR) :: mtot_temp        ! Mass
  real(kind=PR) :: old_mtot_temp    ! ..
  real(kind=PR) :: pofr             ! ..
  real(kind=PR) :: r1               ! ..
  real(kind=PR) :: r2               ! ..
  real(kind=PR) :: radmax           ! ..
  real(kind=PR) :: radp             ! Aux. radius variable
  real(kind=PR) :: rboundary        ! Radius containing boundary particles
  real(kind=PR) :: rcentre(1:NDIM)  ! Centre of polytrope
  real(kind=PR) :: rcloud           ! Radius containing polytrope particles
  real(kind=PR) :: rcore            ! ..
  real(kind=PR) :: remainder        ! Remainder of interpolation
  real(kind=PR) :: renvelope        ! Radius containing gas envelope
  real(kind=PR) :: rgas             ! Radius (of initial sphere) containing gas
  real(kind=PR) :: rho0             ! Central density
  real(kind=PR) :: ricm             ! Radius containing IcM particles
  real(kind=PR), allocatable :: mass_array(:)  ! Interior mass table
  real(kind=PR), allocatable :: density(:)     ! Polytrope density table
  real(kind=PR), allocatable :: radius(:)      ! Dimensionless radius table
  
  debug1("Create polytrope density distribution [ic_polytrope.F90]")
  
! Read parameters file
  call default_parameters
  
! Read polytrope set-up file  
  open(unit=1,file='core.dat',status='old')
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) in_file      
  read(1,*) in_file_form 
  read(1,*) out_file     
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) runit
  read(1,*) munit
  read(1,*) rhounit
  read(1,*) vunit
  read(1,*) coretype
  read(1,*) rho0
  read(1,*) rcore
  read(1,*) radmax
  read(1,*) pofr
  read(1,*) mtot_temp
  read(1,*) menvelope
  read(1,*) micmenv
  read(1,*) hboundary
  close(1)

! Check that parameters are sensible
  if (coretype /= "plummer" .and. coretype /= "powerlaw" .and. &
       &coretype /= "uniform") &
       &call paramerror("Unrecognised option for mflag")
  if (coretype == "powerlaw" .and. pofr <= -3.0_PR) &
       &call paramerror("Invalid power law index")
  if (rcore > radmax) call paramerror("rcore > radmax")

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

! Scale input parameters
  rcore     = rcore * real(r_au/(rscale*r_SI),PR)
  radmax    = radmax * real(r_au/(rscale*r_SI),PR)
  rho0      = rho0 / rhoscale !real(rhoscale*rhocgs,PR)
  menvelope = menvelope / mscale
  micmenv   = micmenv / mscale
  mtot_temp = mtot_temp / mscale

! Initialise variables
  nmax = 2


! Iterate towards converged solution of mass distribution if density function
! is not integratable.
! ============================================================================
  if (coretype == "plummer") then

     ! -----------------------------------------------------------------------
     do        
        old_mtot_temp = mtot_temp
        mtot_temp = 0.0_PR
        
        ! Allocate arrays now number of array points is known
        allocate(density(1:nmax))
        allocate(mass_array(1:nmax))
        allocate(radius(1:nmax))
        
        ! Find mass-table using Simpson's rule
        ! --------------------------------------------------------------------
        do i=1,nmax
           r1 = radmax*real(i - 1,PR)/real(nmax,PR)
           r2 = radmax*real(i,PR)/real(nmax,PR)
           if (coretype == "plummer") then
              mtot_temp = mtot_temp + &
                   & 0.5_PR*(r2 - r1)*(fmplummer(r1,rcore,rho0,pofr) &
                   & + fmplummer(r2,rcore,rho0,pofr))
              density(i) = fplummer(r2,rcore,rho0,pofr)
           end if
           radius(i) = r2
           mass_array(i) = mtot_temp
        end do
        ! --------------------------------------------------------------------
        
        write(6,*) "Iteration : ",nmax,mtot_temp,&
             &abs((mtot_temp - old_mtot_temp)/mtot_temp)
        if (abs((mtot_temp - old_mtot_temp)/mtot_temp) < 1.0E-8_PR) exit
        if (nmax > 1000000) then
           write(6,*) "Solution not converging!  Exitting now"
           stop
        end if
        
        deallocate(radius)
        deallocate(mass_array)
        deallocate(density)
        nmax = 2*nmax
        
     end do
     ! -----------------------------------------------------------------------


! If we have an exact, integratable density/mass function, then simply 
! prepare a dense table of values without iteration.
! ============================================================================
  else if (coretype == "powerlaw") then

     nmax = 10000

     ! Allocate arrays now number of array points is known
     allocate(density(1:nmax))
     allocate(mass_array(1:nmax))
     allocate(radius(1:nmax))
     
     ! Find mass-table using Simpson's rule
     ! -----------------------------------------------------------------------
     do i=1,nmax
        r1            = radmax*real(i - 1,PR)/real(nmax,PR)
        r2            = radmax*real(i,PR)/real(nmax,PR)
        radius(i)     = r2
        mtot_temp     = mpowerlaw(r2,rcore,rho0,pofr)
        density(i)    = fpowerlaw(r2,rcore,rho0,pofr)
        mass_array(i) = mtot_temp
     end do
     ! -----------------------------------------------------------------------
     

! If we have an exact, integratable density/mass function, then simply 
! prepare a dense table of values without iteration.
! ============================================================================
  else if (coretype == "uniform") then

     nmax = 10000
     rho0 = 3.0_PR*mtot_temp/(4.0_PR*PI*(radmax**3))

     ! Allocate arrays now number of array points is known
     allocate(density(1:nmax))
     allocate(mass_array(1:nmax))
     allocate(radius(1:nmax))
     
     ! ..
     ! -----------------------------------------------------------------------
     do i=1,nmax
        r1            = radmax*real(i - 1,PR)/real(nmax,PR)
        r2            = radmax*real(i,PR)/real(nmax,PR)
        radius(i)     = r2
        mtot_temp     = 4.0_PR*PI*rho0*(r2**3)/3.0_PR
        density(i)    = rho0
        mass_array(i) = mtot_temp
     end do
     ! -----------------------------------------------------------------------

  end if
! ============================================================================


  write(6,*) "Radius of cloud   : ",radmax*rscale,runit
  write(6,*) "Radius of core    : ",rcore*rscale,runit
  write(6,*) "Mass of cloud     : ",mtot_temp*mscale,munit
  write(6,*) "Central density   : ",rho0*rhoscale,rhounit

! Load in uniform density sphere file
  write(6,*) "in_file : ",trim(in_file)
  call read_data(in_file,in_file_form)

! Sort particles into gas, icm and boundary zones
  hmean     = (pp_gather/8.0_PR/real(ptot,PR))**(ONETHIRD)
  rboundary = 1.000000001_PR
  ricm      = rboundary - hboundary*hmean
  rgas      = ricm*((mtot_temp + menvelope)/&
       &(mtot_temp + menvelope + micmenv))**(ONETHIRD)
  rcloud    = ricm*(mtot_temp/(mtot_temp + menvelope + micmenv))**(ONETHIRD)

  write(6,*) "hmean     : ",hmean
  write(6,*) "rgas      : ",rgas
  write(6,*) "ricm      : ",ricm
  write(6,*) "rboundary : ",rboundary

  if (ricm < 0.0_PR .or. rgas < 0.0_PR) then
     stop 'Negative ricm or rgas'
  else if (rboundary - rgas < SMALL_NUMBER) then
     write(6,*) "All gas particles : No radial ordering needed"
  else
     call sort_radial_particle_types(rgas,ricm,rboundary,rcentre(1:NDIM))
  end if

! Setting up for different particle types
  call types


! Now stretch all particles and rescale masses to reproduce the density 
! distribution of the polytrope (and the boundary envelope).
! ----------------------------------------------------------------------------
  mp = (mtot_temp + menvelope) / real(pgas,PR)

  do p=1,ptot
     call distance2(rcentre(1:NDIM),p,dr(1:NDIM),drsqd)
     drmag = sqrt(drsqd) + SMALL_NUMBER
     if (drmag >= rcloud) then
        radp = (3.0_PR*mtot_temp)/(4.0_PR*PI*density(nmax))&
             &*((drmag/rcloud)**3 - 1.0_PR) + radmax**3
        radp = radp**(ONETHIRD)
        sph(p)%r(1:NDIM) = sph(p)%r(1:NDIM)*radp/drmag
     else
        mint = mtot_temp*(drmag/rcloud)**3
        ip = 1
        do i=1,nmax
           ip = i
           if (mass_array(i) >= mint) exit
        end do
        remainder = (mint - mass_array(ip))/(mass_array(ip-1) - mass_array(ip))
        radp = radius(ip - 1) + remainder*(radius(ip) - radius(ip - 1))
        sph(p)%r(1:NDIM) = sph(p)%r(1:NDIM)*radp/drmag
     end if
     sph(p)%m = mp
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
  call initialize_thermal_properties

! Calculate hydro forces on all SPH particles
#if defined(HYDRO)
  call sph_hydro_forces
#endif

! Calculate gravitational forces on all SPH particles
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  call sph_grav_forces
#endif

! Calculate gravitational forces on all sink particles
#if defined(SINKS)
  call sph_sink_forces
#endif

! Initialize other variables
  call initialize_seren_variables_2

! Write particle data to file
  out_file = trim(adjustl(out_file))
  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  call write_data_debug("ICCORE.debug.dat",rcentre(1:NDIM))
#endif

  renvelope = ((3.0_PR*mtot_temp)/(4.0_PR*PI*density(nmax))&
       &*((rgas/rcloud)**3 - 1.0_PR) + radmax**3)**(ONETHIRD)
  ricm      = ((3.0_PR*mtot_temp)/(4.0_PR*PI*density(nmax))&
       &*((ricm/rcloud)**3 - 1.0_PR) + radmax**3)**(ONETHIRD)
  rboundary = ((3.0_PR*mtot_temp)/(4.0_PR*PI*density(nmax))&
       &*((rboundary/rcloud)**3 - 1.0_PR) + radmax**3)**(ONETHIRD)

  write(6,*) "rcore     : ",rcore*rscale,runit
  write(6,*) "radmax    : ",radmax*rscale,runit
  write(6,*) "renvelope : ",renvelope*rscale,runit
  write(6,*) "ricm      : ",ricm*rscale,runit
  write(6,*) "rboundary : ",rboundary*rscale,runit


! Clean-up all memory
  call clean_up
  deallocate(mass_array)
  deallocate(density)
  deallocate(radius)
  
  stop
  

contains


! ============================================================================
! FPLUMMER
! User defined density distribution
! ============================================================================
FUNCTION fplummer(rad,rcore,rho0,pofr)
  use definitions
  implicit none

  real(kind=PR), intent(in) :: rad,rcore,rho0,pofr
  real(kind=PR) :: fplummer

  fplummer = rho0/(( 1.0_PR + (rad/rcore)**2 )**(pofr))

  return
END FUNCTION fplummer



! ============================================================================
! FMPLUMMER
! Function 4\pi r^2 \rho(r) to be integrated to get mass within r
! ============================================================================
FUNCTION fmplummer(rad,rcore,rho0,pofr)
  use definitions
  implicit none

  real(kind=PR), intent(in) :: rad,rcore,rho0,pofr
  real(kind=PR) :: fmplummer

  fmplummer = 16.0_PR*atan(1.0_PR)*rad*rad*fplummer(rad,rcore,rho0,pofr)

  return
END FUNCTION fmplummer



! ============================================================================
! FPOWERLAW
! User defined density distribution
! ============================================================================
FUNCTION fpowerlaw(rad,rcore,rho0,pofr)
  use definitions
  implicit none

  real(kind=PR), intent(in) :: rad,rcore,rho0,pofr
  real(kind=PR) :: fpowerlaw

  fpowerlaw = rho0*(rad/rcore)**pofr

  return
END FUNCTION fpowerlaw



! ============================================================================
! MPOWERLAW
! Total mass contained within radius rad of power law density distribution
! ============================================================================
FUNCTION mpowerlaw(rad,rcore,rho0,pofr)
  use definitions
  implicit none

  real(kind=PR), intent(in) :: rad,rcore,rho0,pofr
  real(kind=PR) :: mpowerlaw

  mpowerlaw = 4.0*PI*rho0*(rad**(3.0_PR + pofr))/(3.0_PR + pofr)/rcore**pofr

  return
END FUNCTION mpowerlaw



END PROGRAM ic_core
