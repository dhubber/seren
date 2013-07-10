! GRAVTEST.F90
! D. A. Hubber - 7/4/2008
! Calculates fractional error of gravity tree compared to direct sum 
! gravity for a range of opening angles.  
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM gravtest
  use interface_module
  use particle_module
  use tree_module
  use time_module
  use filename_module
  implicit none

  character(len=256) :: filename               ! ..
  character(len=256) :: filename2              ! ..
  character(len=10) :: file_ext                ! filename ext. for data output
  character(len=2) :: file_numb                ! filename ext. for data output
  character(len=256) :: store_file             ! parameters output
  integer :: i                                 ! aux. counter
  integer :: j                                 ! aux. counter
  integer :: p                                 ! particle counter
  integer :: ttot                              ! no. of
  integer, allocatable :: forder(:)            ! order of max force error
  integer, allocatable :: porder(:)            ! order of max pot error
  real(kind=DP) :: absavaccel                  ! ..
  real(kind=PR) :: adirectmag                  ! ..
  real(kind=DP) :: adirectsqd                  ! ..
  real(kind=DP) :: agravp(1:NDIM)              ! grav. accel. of particle p
  real(kind=DP) :: afactor                     ! ..
  real(kind=PR) :: da(1:NDIM)                  ! ..
  real(kind=DP) :: damean(1:NDIM)              ! ..
  real(kind=DP) :: dasqd                       ! ..
  real(kind=DP) :: dpot                        ! ..
  real(kind=PR) :: error                       ! ..
  real(kind=PR) :: maxerror                    ! max. accel error
  real(kind=PR) :: maxpoterror                 ! max. error in grav. pot.
  real(kind=DP) :: meanabsdev                  ! ..
  real(kind=DP) :: potp                        ! grav. pot. of particle p
  real(kind=PR) :: power                       ! index of abserror value
  real(kind=PR) :: powermin                    ! min index
  real(kind=PR) :: powermax                    ! max index
  real(kind=DP) :: rmsferror                   ! root-mean squared frac. error
  real(kind=DP) :: rmspoterror                 ! ..
  real(kind=PR) :: theta                       ! opening angle
  real(kind=DP) :: t0                          ! start time
  real(kind=DP) :: t1                          ! end time
  real(kind=DP) :: tdiff                       ! time difference
  real(kind=DP) :: tdirect                     ! direct-sum time
  real(kind=PR), allocatable :: adirect(:,:)   ! direct-sum grav. accel
  real(kind=PR), allocatable :: errorarray(:)  ! grav. accel. errors
  real(kind=PR), allocatable :: potdirect(:)   ! direct-sum grav. pot.
  real(kind=PR), allocatable :: perrorarray(:) ! grav. pot errors

! Read in an process command line arguments
  call read_arguments

! Set default parameters
  call default_parameters

  write(6,*) "------------------------------"
  write(6,*) "          gravtest            "
  write(6,*) "------------------------------"
  open(unit=1,file="gravtest.dat",status="old")
  read(1,*) in_file
  read(1,*) in_file_form
  read(1,*) run_id
  read(1,*) ttot
  read(1,*) powermin
  read(1,*) powermax

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types

! Checking compiler flags and parameter values
  call sanitycheck

! Setting up scaling units for simulation
  call units

! Reading in formatted data file (DRAGON format)
  call read_data(in_file,in_file_form)

! Setting up for different particle types
  call types

! Writing compiler flags and parameters to file
  store_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".params"
  call paramstore(store_file)

! Converting to dimensionless code units
  call convert_to_code_units

! Initialising SPH kernel tables
#if defined(KERNEL_TABLES)
  call tabulate_kernel_functions
#endif

! Initialize certain variables before first force calculation
  call initialize_seren_variables_2
  call initialize_sph_variables_1

! Build and stock trees for first time
  call tree_update(nbuild,nstock)

! Make initial guesses of h either using tree or a global estimate
#if defined(BH_TREE) && !defined(CONSTANT_H)
  if (.NOT. restart) call BHhydro_hguess
#else
  if (.NOT. restart) call h_guess
#endif

! Calculating initial SPH quantities (h, neibs, rho, press etc.)
  call sph_update

! Updating tree properties
#if defined(BH_TREE) && defined(SELF_GRAVITY)
  call BHgrav_stock
#endif

! Allocate temporary array for direct sum case
  allocate(adirect(1:NDIM,1:ptot))
  allocate(errorarray(1:ptot))
  allocate(perrorarray(1:ptot))
  allocate(potdirect(1:ptot))
  allocate(forder(1:ptot))
  allocate(porder(1:ptot))

! Calculate direct sum gravity and record in array
  write(6,*) "Calculating direct sum gravity ..."

  ! Record start time before tree stocking and calculating grav. force
  call cpu_time(t0)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(agravp,potp)
  do p=1,ptot
     call direct_sph_gravity(p,1.0_PR/sph(p)%h,sph(p)%r(1:NDIM),agravp,potp)
     adirect(1:NDIM,p) = real(agravp(1:NDIM),PR)
     potdirect(p) = real(potp,PR)
  end do
!$OMP END PARALLEL DO

  ! Calculate time taken to calculate gravitational forces
  call cpu_time(t1)
  tdirect = t1 - t0

#if defined(DEBUG_FORCES)
  a_grav(1:NDIM,1:ptot) = adirect(1:NDIM,1:ptot)
#endif
#if !defined(GEOMETRIC_MAC)
  do p=1,ptot
     sph(p)%agravmag = sqrt(dot_product(adirect(1:NDIM,p),adirect(1:NDIM,p)))
  end do
#endif

! Open main statistics file here for writing
  filename = trim(run_id)//".dat"
  open(1,file=filename,status="unknown")

#if defined(BH_TREE)
  call BHgrav_build
#endif


! Loop over various values of theta/abserror and calculate statistics
! ============================================================================
  do i=1,ttot
     power = powermin + real(i - 1,PR)*(powermax - powermin)/real(ttot - 1,PR)
#if defined(GEOMETRIC_MAC)
     theta = 10.0_PR**(power)
     thetamaxsqd = theta*theta
     write(6,*) "theta :",theta,"   thetamaxsqd :",thetamaxsqd
#else
     abserror = 10.0_PR**(power)
     theta = abserror
     write(6,*) "abserror :",abserror
     do p=1,ptot
        sph(p)%a(1:NDIM) = adirect(1:NDIM,p)
     end do
#endif
     sph(1:ptot)%gpot = potdirect(1:ptot)

     ! Record start time before tree stocking and calculating grav. force
     call cpu_time(t0)

#if defined(BH_TREE)
     call BHgrav_stock
#endif

     ! Calculate gravitational force of all particles via tree
     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(agravp,potp)
     do p=1,ptot
#if defined(BH_TREE)
        call BHgrav_accel(p,1.0_PR/sph(p)%h,&
             &sph(p)%r(1:NDIM),agravp,potp)
#endif
        sph(p)%a(1:NDIM) = real(agravp(1:NDIM),PR)
        sph(p)%gpot = real(potp,PR)
     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------

     ! Calculate time taken to calculate gravitational forces
     call cpu_time(t1)
     tdiff = t1 - t0

     ! Calculate distribution of fractional errors
     ! First zero all relevant variables
     absavaccel         = 0.0_DP
     damean(1:NDIM)     = 0.0_DP
     maxerror           = 0.0_DP
     maxpoterror        = 0.0_DP
     meanabsdev         = 0.0_DP
     rmsferror          = 0.0_DP
     rmspoterror        = 0.0_DP
     errorarray(1:ptot) = 0.0_DP

     ! Now calculate error contributions from all particles 
     ! -----------------------------------------------------------------------
     do p=1,ptot
        da(1:NDIM)     = real(sph(p)%a(1:NDIM) - adirect(1:NDIM,p),DP)
        damean(1:NDIM) = damean(1:NDIM) + da(1:NDIM)
        dasqd          = dot_product(da(1:NDIM),da(1:NDIM))
        adirectsqd     = dot_product(adirect(1:NDIM,p),adirect(1:NDIM,p))
        rmsferror      = rmsferror + dasqd/adirectsqd
        absavaccel     = absavaccel + sqrt(adirectsqd)
        error          = sqrt(dasqd/adirectsqd)
        maxerror       = max(maxerror,error)
        errorarray(p)  = max(errorarray(p),error)
        dpot           = sph(p)%gpot - potdirect(p)
        error          = abs(dpot)/potdirect(p)
        rmspoterror    = rmspoterror + error*error
        perrorarray(p) = error
        maxpoterror    = max(maxpoterror,error)
        forder(p)      = p
     end do
     ! -----------------------------------------------------------------------

     rmsferror = sqrt(rmsferror / real(ptot,DP))
     rmspoterror = sqrt(rmspoterror / real(ptot,DP))

     damean(1:NDIM) = damean(1:NDIM) / real(ptot,DP)
     absavaccel = absavaccel / real(ptot,DP)
     do p=1,ptot
        da(1:NDIM) = real(sph(p)%a(1:NDIM) - adirect(1:NDIM,p),DP) - damean(1:NDIM)
        meanabsdev = meanabsdev + &
             &real(sqrt(dot_product(da(1:NDIM),da(1:NDIM))),DP)
     end do
     meanabsdev = meanabsdev / real(ptot,DP)

     ! Sort errors for easy analysis and visualisation
     call heapsort_real(ptot,errorarray(1:ptot),forder(1:ptot))
     call heapsort_real(ptot,perrorarray(1:ptot),porder(1:ptot))

     ! Write all error information to file
!     write(1,'(14E15.7)') theta,tdiff/tdirect,rmsferror,maxerror,&
!          &errorarray(forder(ptot/2)),meanabsdev/absavaccel,&
!          &rmspoterror,maxpoterror,perrorarray(porder(ptot/2))
     write(1,'(14E15.7)') theta,tdiff/tdirect,rmsferror,maxerror,&
          &errorarray(ptot/2),meanabsdev/absavaccel,&
          &rmspoterror,maxpoterror,perrorarray(ptot/2)

     ! Data output needs extension
     if (i < 10) then
        write(file_numb,"(I1)") i
        file_ext="ferror.0"//file_numb
     else if (snapshot < 100) then
        write(file_numb,"(I2)") snapshot
        file_ext="ferror."//file_numb
     end if
     filename2 = trim(adjustl(run_id))//"."//trim(adjustl(file_ext))

     ! Write error information for each accuracy parameter and all particles
     open(unit=2,file=filename2,form='formatted')
     do j=1,ptot
        p = forder(j)
        write(2,'(17E15.7)') real(j,PR)/real(ptot,PR),errorarray(p),&
             &perrorarray(j),sph(p)%r(1:NDIM),&
             &adirect(1:NDIM,p),sqrt(dot_product(adirect(1:NDIM,p),&
             &adirect(1:NDIM,p))),sph(p)%a(1:NDIM),sqrt(dot_product(sph(p)%a(1:NDIM),&
             &sph(p)%a(1:NDIM))),potdirect(p),sph(p)%gpot
     end do
     close(2)

     write(6,*) rmsferror,maxerror,errorarray(ptot/2),&
          &rmspoterror,maxpoterror,perrorarray(ptot/2)

  end do
! ============================================================================

  close(1)

! Release all allocated memory
  deallocate(porder)
  deallocate(forder)
  deallocate(potdirect)
  deallocate(errorarray)
  deallocate(adirect)
  call clean_up

  stop
END PROGRAM gravtest

