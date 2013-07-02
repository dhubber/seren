! IC_BONDI.F90
! D. A. Hubber - 29/04/2011
! ..
! Reads in a parameters file 'bondiparams.dat' with the following data : 
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
  logical :: isothermal               ! ..
  integer :: i                        ! ..
  integer :: ip                       ! ..
  integer :: nmax                     ! ..
  integer :: p                        ! Particle counter
  integer :: pmin                     ! ..
  real(kind=PR) :: agas               ! ..
  real(kind=PR) :: dmdt               ! ..
  real(kind=PR) :: dr(1:NDIM)         ! Relative position
  real(kind=PR) :: dr_unit(1:NDIM)    ! ..
  real(kind=PR) :: drmag              ! Distance
  real(kind=PR) :: drsqd              ! Distance squared
  real(kind=PR) :: mcloud             ! ..
  real(kind=PR) :: mint               ! ..
  real(kind=PR) :: mp                 ! Mass of particle p
  real(kind=PR) :: msink              ! Mass
  real(kind=PR) :: qs                 ! ..
  real(kind=PR) :: racc               ! Accretion radius
  real(kind=PR) :: radp               ! ..
  real(kind=PR) :: rbondi             ! ..
  real(kind=PR) :: rcentre(1:NDIM)    ! ..
  real(kind=PR) :: rcloud             ! ..
  real(kind=PR) :: remainder          ! ..
  real(kind=PR) :: rhogas             ! ..
  real(kind=PR) :: rsink              ! ..
  real(kind=PR) :: rsonic             ! ..
  real(kind=PR) :: vwind              ! ..
  real(kind=PR), allocatable :: mass_array(:)  ! Interior mass table
  real(kind=PR), allocatable :: radius(:)      ! Dimensionless radius table
  real(kind=PR), allocatable :: v_array(:)     ! ..
  real(kind=PR), allocatable :: a_array(:)     ! ..
  real(kind=PR), allocatable :: rho_array(:)   ! ..
  real(kind=PR), allocatable :: mu_array(:)    ! ..
  real(kind=PR), allocatable :: u_array(:)     ! ..
  real(kind=PR), allocatable :: w_array(:)     ! ..
  real(kind=PR), allocatable :: x_array(:)     ! ..
  real(kind=PR), allocatable :: y_array(:)     ! ..
  real(kind=PR), allocatable :: z_array(:)     ! ..

  debug1("Create polytrope density distribution [ic_bondi.F90]")
  
! Read parameters file
  call default_parameters
  
! Read Bondi model set-up file  
  open(unit=1,file='bondiparams.dat',status='old')
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) in_file      
  read(1,*) in_file_form 
  read(1,*) out_file     
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) isothermal
  read(1,*) runit
  read(1,*) munit
  read(1,*) rhounit
  read(1,*) vunit
  read(1,*) gamma
  read(1,*) rhogas
  read(1,*) mu_bar
  read(1,*) agas
  read(1,*) rsink
  read(1,*) msink
  read(1,*) mcloud
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
  rhogas = rhogas / rhoscale
  agas   = agas   / vscale
  msink  = msink  / mscale

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
  if (isothermal) then
     gamma = 1.0_PR
     racc = 2.0_PR*msink/agas/agas
     rsonic = 0.5_PR*msink/agas/agas
     dmdt = exp(1.5_PR)*PI*msink*msink*rhogas/(agas**3)
     !rsink = rsink*rsonic
  else
     racc = 2.0_PR*msink/agas/agas
     rsonic = racc*(5.0_PR - 3.0_PR*gamma)/8.0_PR
     if (abs(1.6666666666666_PR - gamma) < 0.0001_PR) then
        qs = 0.25_PR
     else if (abs(1.0_PR - gamma) < 0.0001_PR) then
        stop 'gamma = 1 without using isothermal flag'
     else
        qs = 0.25_PR*(2.0_PR/(5.0_PR - 3.0_PR*gamma))**&
             &((5.0_PR - 3.0_PR*gamma)/(2.0_PR*gamma - 2.0_PR))
     end if
     dmdt = 4.0_PR*PI*qs*(msink*msink*rhogas)/(agas**3)
     !rsink = rsink*racc
     rcloud = rcloud*racc
  end if

  ! Allocate arrays now number of array points is known
  nmax = 1000
  allocate(mass_array(0:nmax))
  allocate(radius(0:nmax))
  allocate(v_array(0:nmax))
  allocate(a_array(0:nmax))
  allocate(mu_array(0:nmax))
  allocate(rho_array(0:nmax))
  allocate(u_array(0:nmax))
  allocate(w_array(0:nmax))
  allocate(x_array(0:nmax))
  allocate(y_array(0:nmax))
  allocate(z_array(0:nmax))


  ! If imposing an isothermal EOS
  ! --------------------------------------------------------------------------
  if (isothermal) then
     call bondi_isothermal_solution(nmax,w_array,x_array,y_array,z_array)

     do i=0,nmax
        drmag = rsonic*x_array(i)
        radius(i) = drmag
        mass_array(i) = z_array(i)
        v_array(i) = w_array(i)
        a_array(i) = agas
     end do

  end if
  ! --------------------------------------------------------------------------


  ! Check required mass of cloud is not greater than that provided by table
  write(6,*) "mcloud : ",mcloud,"   z_array : ",10.0_DP**(z_array(nmax))
  if (log10(mcloud) > z_array(nmax)) stop 'Cloud mass too big'


  ! Find radius of cloud
  ip = 1
  do i=1,nmax
     ip = i
     if (z_array(i) >= log10(mcloud)) exit
  end do
  remainder = (log10(mcloud) - z_array(ip - 1))/(z_array(ip) - z_array(ip - 1))
  rcloud = x_array(ip - 1) + remainder*(x_array(ip) - x_array(ip - 1))
  rcloud = 10.0_DP**(rcloud)
  write(6,*) "rcloud : ",rcloud,ip,nmax,remainder,&
       &10.0_PR**(z_array(ip)),10.0_PR**(z_array(ip - 1)),mcloud


  ! Now stretch all particles and rescale masses to reproduce the density 
  ! distribution of the polytrope (and the boundary envelope).
  ! ---------------------------------------------------------------------------
  mp = 4.0_PR*PI*(rsonic**3)*rhogas*mcloud / real(pgas,PR)
  !mp = 0.5_PR*PI*rhogas*mcloud*msink**3 / (agas**6) / real(pgas,PR)
  !mp = mcloud/real(pgas,PR)
  write(6,*) "mp : ",mp,&
       &0.5_PR*PI*rhogas*mcloud*msink**3 / (agas**6) / real(pgas,PR)

  do p=1,ptot
     call distance2(rcentre(1:NDIM),p,dr(1:NDIM),drsqd)
     drmag = sqrt(drsqd) + SMALL_NUMBER
     mint = log10(mcloud*(drmag)**3)
     ip = 1
     do i=1,nmax
        ip = i
        if (z_array(i) >= mint) exit
     end do
     remainder = (mint - z_array(ip - 1))/(z_array(ip) - z_array(ip - 1))
     radp = x_array(ip - 1) + remainder*(x_array(ip) - x_array(ip - 1))
     sph(p)%r(1:NDIM) = rsonic*sph(p)%r(1:NDIM)*10.0_PR**(radp)/drmag
     sph(p)%m = mp
     sph(p)%v(1:NDIM) = w_array(ip - 1) + &
          &remainder*(w_array(ip) - w_array(ip - 1))
     sph(p)%v(1:NDIM) = -agas*10.0_PR**(sph(p)%v(1:NDIM))*dr(1:NDIM)/drmag
  end do
  

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

! Find minimum smoothing length
  hmin = BIG_NUMBER
  do p=1,ptot
     hmin = min(hmin,sph(p)%h)
  end do
  !rsink = rsink*KERNRANGE*hmin
  rsink = rsink !*racc

! Place sink particle at centre of cloud
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
  sink(stot)%mmax = 0.0_DP

! Now find total mass inside sinkand set to mmax
  do p=1,ptot
     call distance2(rzero(1:NDIM),p,dr(1:NDIM),drsqd)
     if (drsqd > rsink*rsink) cycle
     sink(stot)%mmax = sink(stot)%mmax + sph(p)%m
  end do

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

  ! Output important info to screen
  ! ---------------------------------------------------------------------------
  if (isothermal) then
     write(6,*) "Isothermal?              : ",isothermal
     write(6,*) "Gas density              : ",rhogas*rhoscale,rhounit
     write(6,*) "Total gas mass           : ",&
          &0.5_PR*PI*rhogas*mcloud*msink**3*mscale/(agas**6),munit
     write(6,*) "SPH particle mass        : ",mp*mscale,munit
     write(6,*) "Gas sound speed          : ",agas*vscale,vunit
     write(6,*) "Minimum smoothing length : ",hmin*rscale,runit
     write(6,*) "Sink radius              : ",rsink*rscale,runit
     write(6,*) "Sink mass                : ",msink*mscale,munit
     write(6,*) "Accretion radius         : ",racc*rscale,runit
     write(6,*) "Sonic radius             : ",rsonic*rscale,runit
     write(6,*) "Cloud radius             : ",rcloud*rscale,runit
     write(6,*) "qs                       : ",qs
     write(6,*) "Accretion rate           : ",dmdt*dmdtscale,dmdtunit
     write(6,*) "Accretion timescale      : ",&
          &rhogas*4.0_PR*PI*ONETHIRD*rcloud**3/dmdt*tscale,tunit
     write(6,*) "Freefall time            : ",&
          &(1.0_PR/sqrt(rhogas))*tscale,tunit
     write(6,*) "Sound crossing time      : ",rcloud/agas*tscale,tunit
     write(6,*) "Temperature              : ",agas**2/gamma/Pconst
     write(6,*) "Mass test                : ",4.0_PR*PI*rsonic**3,&
          &0.5_PR*PI*msink**3*rhogas/(agas**6)
  end if

! Write particle data to file
  out_file = trim(adjustl(out_file))
  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  call write_data_debug("ICBH.debug.dat",rcentre(1:NDIM))
#endif

! Clean-up all memory
  call clean_up

  deallocate(radius)
  deallocate(mass_array)

  stop
END PROGRAM ic_bondi_hoyle




