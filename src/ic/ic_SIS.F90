! IC_SIS.F90
! D. A. Hubber - 21/7/2010
! Create a singular isothermal sphere density profile (i.e. sph(p)%rho = k*rho^-2) 
! from a unit-radius uniform density sphere.
! Reads in a parameters file 'SISparams.dat' with the following data : 
! in_file       : Filename of relaxed unit-sphere
! in_file_form  : Format of input sphere file 
! out_file      : Output filename
! out_file_form : Format of output file
! runit         : Length unit
! munit         : Mass unit
! rhounit       : Density unit
! vunit         : Velocity unit
! radmax        : Maximum extent of cloud
! a0            : Sound speed
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_SIS
  use interface_module, only : comperror,create_sink,distance2,&
       &paramerror,read_data,write_data
  use particle_module
  use filename_module
  use scaling_module
  use time_module
  use type_module
  use neighbour_module
  use constant_module
  use sink_module
  use hydro_module, only : isotemp, temp
  implicit none
  
  character(len=256) :: out_file    ! Name of output file
  integer :: p                      ! Particle counter
  integer :: pmin                   ! ..
  real(kind=PR) :: a0               ! ..
  real(kind=PR) :: dr(1:NDIM)       ! Relative position
  real(kind=PR) :: dr_unit(1:NDIM)  ! ..
  real(kind=PR) :: drmag            ! Distance
  real(kind=PR) :: drsqd            ! Distance squared
  real(kind=PR) :: hboundary        ! Thickness of boundary region (in h)
  real(kind=PR) :: hmean            ! Mean smoothing length for uniform sphere
  real(kind=PR) :: menvelope        ! Mass of envelope
  real(kind=PR) :: micmenv          ! Mass of IcM region
  real(kind=PR) :: mint             ! Mass interior to radius
  real(kind=PR) :: mp               ! Mass of particle p
  real(kind=PR) :: mcloud           ! Mass
  real(kind=PR) :: radmax           ! ..
  real(kind=PR) :: radp             ! Aux. radius variable
  real(kind=PR) :: rboundary        ! Radius containing boundary particles
  real(kind=PR) :: rcentre(1:NDIM)  ! Centre of polytrope
  real(kind=PR) :: rcloud           ! Radius containing polytrope particles
  real(kind=PR) :: rgas             ! Radius (of initial sphere) containing gas
  real(kind=PR) :: rhoenv           ! ..
  real(kind=PR) :: ricm             ! Radius containing IcM particles
  real(kind=PR) :: rsqdmin          ! ..

  debug1("Create polytrope density distribution [ic_SIS.F90]")
  
! Read parameters file
  call default_parameters
  
! Read polytrope set-up file  
  open(unit=1,file='SISparams.dat',status='old')
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
  read(1,*) radmax
  read(1,*) a0
  read(1,*) sinkrad
  read(1,*) menvelope
  read(1,*) micmenv
  read(1,*) hboundary
  close(1)

! Check that parameters are sensible
#if !defined(ISOTHERMAL)
  call comperror("ISOTHERMAL flag not switched on")
#endif
#if !defined(SINKS)
  call comperror("SINKS not switched on")
#endif
  if (a0 <= 0.0_PR) call paramerror("a0 <= 0")
  if (radmax <= 0.0_PR) call paramerror("radmax <= 0")

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
  radmax    = radmax / rscale
  a0        = a0 / vscale
  menvelope = menvelope / mscale
  micmenv   = micmenv / mscale
  mcloud      = 2.0_PR*a0*a0*radmax
  rhoenv    = a0*a0/(TWOPI*radmax*radmax)

  write(6,*) "Radius of cloud   : ",radmax*rscale,runit
  write(6,*) "Mass of cloud     : ",mcloud*mscale,munit

! Load in uniform density sphere file
  write(6,*) "in_file : ",trim(in_file)
  call read_data(in_file,in_file_form)

! Sort particles into gas, icm and boundary zones
  hmean     = (pp_gather/8.0_PR/real(ptot,PR))**(ONETHIRD)
  rboundary = 1.000000001_PR
  ricm      = rboundary - hboundary*hmean
  rgas      = ricm*((mcloud + menvelope)/&
       &(mcloud + menvelope + micmenv))**(ONETHIRD)
  rcloud    = ricm*(mcloud/(mcloud + menvelope + micmenv))**(ONETHIRD)

  write(6,*) "hmean     : ",hmean
  write(6,*) "rgas      : ",rgas
  write(6,*) "ricm      : ",ricm
  write(6,*) "rboundary : ",rboundary
  write(6,*) "rcloud    : ",rcloud

  if (ricm < 0.0_PR .or. rgas < 0.0_PR) then
     stop 'Negative ricm or rgas'
  else if (rboundary - rgas < SMALL_NUMBER) then
     write(6,*) "All gas particles : No radial ordering needed"
  else
     call sort_radial_particle_types(rgas,ricm,rboundary,rcentre(1:NDIM))
  end if

! Setting up for different particle types
  call types

!  com_frame = .true.
!  call COM


! Now stretch all particles and rescale masses to reproduce the density 
! distribution of the polytrope (and the boundary envelope).
! ----------------------------------------------------------------------------
!  mp = mcloud / real(pgas,PR)
  mp = (mcloud + menvelope) / real(pgas,PR)

  do p=1,ptot
     call distance2(rcentre(1:NDIM),p,dr(1:NDIM),drsqd)
     drmag = sqrt(drsqd) + SMALL_NUMBER
     dr_unit(1:NDIM) = dr(1:NDIM) / drmag

     if (drmag >= rcloud) then
        radp = (3.0_PR*mcloud)/(4.0_PR*PI*rhoenv)&
             &*((drmag/rcloud)**3 - 1.0_PR) + radmax**3
        radp = radp**(ONETHIRD)
        sph(p)%r(1:NDIM) = sph(p)%r(1:NDIM)*radp/drmag
     else
        mint = mcloud*(drmag/rcloud)**3
        sph(p)%r(1:NDIM) = radmax*dr_unit(1:NDIM)*(drmag/rcloud)**3
     end if

     sph(p)%m = mp
     sph(p)%v(1:NDIM) = 0.0_PR
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
  call sph_hydro_forces

! Calculate gravitational forces on all SPH particles
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  call sph_grav_forces
#endif

! Initialize other variables
  call initialize_sph_variables_2

! Place sink particle near centre of cloud
  pmin = 1
  rsqdmin = BIG_NUMBER
  do p=1,ptot
     call distance2(rcentre(1:NDIM),p,dr(1:NDIM),drsqd)
     if (drsqd < rsqdmin) then
        pmin = p
        rsqdmin = drsqd
     end if
  end do
  write(6,*) "rmin : ",rsqdmin,pmin

  call create_sink(pmin)
  write(6,*) "rs : ",sink(1)%r(1:NDIM),sink(1)%radius,sinkrad
  call diagnostics

! Give all particles near or inside sink a small inwards velocity perturbation.
  do p=1,ptot
     call distance2(rcentre(1:NDIM),p,dr(1:NDIM),drsqd)
     if (drsqd < 9.0_PR*sink(1)%radius**2) then
        sph(p)%v(1:VDIM) = 0.1_PR*a0*sph(p)%r(1:NDIM)*&
             &(ONETHIRD*sqrt(drsqd) - sink(1)%radius)/&
             &(sink(1)%radius)/sqrt(drsqd)
     end if
  end do

! Write particle data to file
  out_file = trim(adjustl(out_file))
  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  call write_data_debug("ICSIS.debug.dat",rcentre(1:NDIM))
#endif

! Clean-up all memory
  call clean_up

  stop
END PROGRAM ic_SIS
