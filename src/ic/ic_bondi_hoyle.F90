! IC_BONDI_HOYLE.F90
! D. A. Hubber - 25/02/2011
! ..
! Reads in a parameters file 'BHparams.dat' with the following data : 
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
PROGRAM ic_bondi_hoyle
  use interface_module, only : comperror,distance2,&
       &paramerror,read_data,write_data
  use particle_module
  use filename_module
  use scaling_module
  use time_module
  use type_module
  use neighbour_module
  use constant_module
  use periodic_module
  use sink_module
  use hydro_module
  implicit none
  
  character(len=256) :: out_file      ! Name of output file
  integer :: p                        ! Particle counter
  integer :: pmin                     ! ..
  real(kind=PR) :: agas               ! ..
  real(kind=PR) :: boxsize            ! ..
  real(kind=PR) :: dmdt               ! ..
  real(kind=PR) :: dr(1:NDIM)         ! Relative position
  real(kind=PR) :: dr_unit(1:NDIM)    ! ..
  real(kind=PR) :: drmag              ! Distance
  real(kind=PR) :: drsqd              ! Distance squared
  real(kind=PR) :: mp                 ! Mass of particle p
  real(kind=PR) :: msink              ! Mass
  real(kind=PR) :: qs                 ! ..
  real(kind=PR) :: racc               ! Accretion radius
  real(kind=PR) :: rcentre(1:NDIM)    ! ..
  real(kind=PR) :: rhogas             ! ..
  real(kind=PR) :: rsink              ! ..
  real(kind=PR) :: rsonic             ! ..
  real(kind=PR) :: vwind              ! ..

  debug1("Create polytrope density distribution [ic_bondi_hoyle.F90]")
  
! Read parameters file
  call default_parameters
  
! Read polytrope set-up file  
  open(unit=1,file='BHparams.dat',status='old')
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
  read(1,*) gamma
  read(1,*) rhogas
  read(1,*) mu_bar
  read(1,*) agas
  read(1,*) vwind
  read(1,*) boxsize
  read(1,*) rsink
  read(1,*) msink
  close(1)

! Check that parameters are sensible
#if !defined(SINKS)
  call comperror("SINKS not switched on")
#endif
  if (gamma > 1.6666666666666_PR) call paramerror("gamma > 5/3")
  if (agas <= 0.0_PR) call paramerror("agas <= 0")
  if (rhogas <= 0.0_PR) call paramerror("rhogas <= 0")

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
  rhogas    = rhogas  / rhoscale
  rsink     = rsink   / rscale
  agas      = agas    / vscale
  vwind     = vwind   / vscale
  msink     = msink   / mscale
  boxsize   = boxsize / rscale

  periodic_max(1:NDIM) = 0.5_PR*boxsize
  periodic_min(1:NDIM) = -0.5_PR*boxsize
  periodic_size(1:NDIM) = boxsize
  periodic_half(1:NDIM) = 0.5_PR*boxsize
#if defined(MINIMUM_H)
  hmin = INVKERNRANGE*rsink
#endif

! Load in uniform density sphere file
  write(6,*) "in_file : ",trim(in_file)
  call read_data(in_file,in_file_form)

! Setting up for different particle types
  pgas = ptot
  picm = 0
  pboundary = 0
  pcdm = 0
  call types

! Calculate important quantities
  racc = 2.0_PR*msink/agas/agas
  rsonic = racc*(5.0_PR - 3.0_PR*gamma)/8.0_PR
  if (abs(1.6666666666666_PR - gamma) < 0.0001_PR) then
     qs = 0.25_PR
  else if (abs(1.0_PR - gamma) < 0.0001_PR) then
     qs = 0.25_PR*exp(1.5_PR)
  else
     qs = 0.25_PR*(2.0_PR/(5.0_PR - 3.0_PR*gamma))**&
          &((5.0_PR - 3.0_PR*gamma)/(2.0_PR*gamma - 2.0_PR))
  end if
  dmdt = 4.0_PR*PI*qs*(msink*msink*rhogas)/(agas*agas + vwind*vwind)**(1.5_PR)
  mp = rhogas*boxsize**(NDIM) / real(pgas,PR)

  write(6,*) "Gas density         : ",rhogas*rhoscale,rhounit
  write(6,*) "Total gas mass      : ",rhogas*boxsize**(NDIM)*mscale,munit
  write(6,*) "SPH particle mass   : ",mp*mscale,munit
  write(6,*) "Gas sound speed     : ",agas*vscale,vunit
  write(6,*) "Sink radius         : ",rsink*rscale,runit
  write(6,*) "Sink mass           : ",msink*mscale,munit
  write(6,*) "Accretion radius    : ",racc*rscale,runit
  write(6,*) "Sonic radius        : ",rsonic*rscale,runit
  write(6,*) "boxsize/volume      : ",boxsize*rscale,boxsize**(NDIM)
  write(6,*) "qs                  : ",qs
  write(6,*) "Accretion rate      : ",dmdt*dmdtscale,dmdtunit
  write(6,*) "Accretion timescale : ",rhogas*boxsize**(NDIM)/dmdt*tscale,tunit
  write(6,*) "Freefall time       : ",(1.0_PR/sqrt(rhogas))*tscale,tunit
  write(6,*) "Sound crossing time : ",boxsize/agas*tscale,tunit
  write(6,*) "Temperature         : ",agas**2/gamma/Pconst

! Set SPH particle properties
! ----------------------------------------------------------------------------
  do p=1,ptot
     sph(p)%r(1:NDIM) = boxsize*(sph(p)%r(1:NDIM) - 0.5_PR)
     sph(p)%m = mp
     sph(p)%v(1:NDIM) = 0.0_PR
     sph(p)%v(1) = vwind
     sph(p)%press = rhogas*agas**2/gamma
    sph(p)%temp = sph(p)%press / (Pconst*rhogas)
#if defined(ENERGY_EQN)
     sph(p)%u = sph(p)%press/(gamma - 1.0_PR)/rhogas
#endif
  end do

! Place sink particle at centre of cloud
! ----------------------------------------------------------------------------
  stot = 1
  sink(stot)%r(1:NDIM) = 0.0_PR
  sink(stot)%v(1:VDIM) = 0.0_PR
  sink(stot)%m = msink
  sink(stot)%radius = rsink
  sink(stot)%h = INVKERNRANGE*rsink
  sink(stot)%accrete = .true.
  sink(stot)%static = .false.
  sink(stot)%tcreate = 0.0_DP
  sink(stot)%ncreate = 0
  sink(stot)%macc(1:DMDT_RANGE) = 0.0_DP
  sink(stot)%tacc(1:DMDT_RANGE) = 0.0_DP
  sink(stot)%angmom(1:3) = 0.0_DP
  sink(stot)%angmomnet(1:3) = 0.0_DP

! Initialising kernel tables
#if defined(KERNEL_TABLES)
  call tabulate_kernel_functions
#endif

! Create Ewald correction table
#if defined(GRAVITY) && defined(EWALD)
  call ewald_init
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
  call diagnostics

! Write particle data to file
  out_file = trim(adjustl(out_file))
  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  call write_data_debug("ICBH.debug.dat",rcentre(1:NDIM))
#endif

! Clean-up all memory
  call clean_up

  stop
END PROGRAM ic_bondi_hoyle
