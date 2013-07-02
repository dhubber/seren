! ERROR_NORM.F90
! D. A. Hubber - 30/09/2010
! Calculates the L1 and L2 error norms for some property of a snapshot 
! compared to a given analytical solution (from a file).
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM error_norm
  use interface_module
  use particle_module
  use filename_module
  use scaling_module
  use time_module
  use type_module
  use neighbour_module
  use constant_module
  use sink_module
  use hydro_module
  implicit none

  character(len=256) :: file1                 ! File to be analysed
  character(len=256) :: file1_form            ! File format
  character(len=256) :: analytic_file         ! Analytic solution filename
  character(len=256) :: params_file           ! Parameters file name
  character(len=256) :: testid                ! Test identity
  integer :: i                                ! Table element counter
  integer :: iaux                             ! Aux integer variable
  integer :: idebug                           ! ..
  integer :: isol                             ! ..
  integer :: jdebug                           ! ..
  integer :: jsol                             ! ..
  integer :: ndebug                           ! ..
  integer :: nelements                        ! ..
  integer :: normorder                        ! ..
  integer :: nsol                             ! ..
  integer :: p                                ! Particle counter
  real(kind=PR) :: alldata(1:100)             ! ..
  real(kind=DP) :: L1denom                    ! L1 error denominator term
  real(kind=DP) :: L1norm                     ! L1 error norm
  real(kind=DP) :: L1num                      ! L1 error numerator term
  real(kind=DP) :: L2denom                    ! L2 error denominator term
  real(kind=DP) :: L2norm                     ! L2 error norm
  real(kind=DP) :: L2num                      ! L2 error numerator term
  real(kind=PR) :: raux                       ! ..
  real(kind=PR) :: rorigin(1:3)               ! Position of origin
  real(kind=PR), allocatable :: solution(:,:) ! Analytic solution array

! Set default parameters
  call default_parameters

! Read polytrope set-up file  
  open(unit=1,file='errornorm.dat',status='old')
  read(1,*) file1
  read(1,*) file1_form
  read(1,*) params_file
  read(1,*) testid
  read(1,*) analytic_file
  read(1,*) normorder

! Read-in parameters file
  call read_parameters(params_file)
  restart = .false.

! Initialise some variables using parameters
  call initialize_seren_variables_1

! Checking compiler flags and parameter values
  call sanitycheck

! Setting up scaling units for simulation
  call units

! Reading in formatted data file
  call read_data(file1,file1_form)

! Setting up for different particle types
  call types

! Converting to dimensionless code units
  call convert_to_code_units

! Initialising kernel tables
  call kernel

! Read in opacity tables for radiation transport
#if defined(HYDRO) && defined(GRAVITY) && defined(RAD_WS)
  call read_eos
#endif

! Create Ewald correction table
#if defined(GRAVITY) && defined(EWALD)
  call ewald_init
#endif

! Calculate COM of system, and if required, change to COM
  com_frame = .true.
  call COM

! Initialise all other particles arrays
  call initialize_seren_variables_2

! Perform all necessary SPH calculations
  call sph_setup


! Now read in analytic solution files
! ----------------------------------------------------------------------------
  if (testid == "polyrho") then

     open(unit=1,file=analytic_file,status="old")
     nelements = 9
     isol = 6
     jsol = 7
     idebug = 17
     jdebug = 11
     nsol = 0
     rorigin(1:NDIM) = 0.0_PR

     ! Find out how many elements are in solutions file
     do
        read(1,*,err=5,end=5) alldata(1:nelements)
        nsol = nsol + 1
     end do

5    allocate(solution(1:nelements,1:nsol))
     rewind(1)
     do i=1,nsol
        read(1,*,err=10,end=10) solution(1:nelements,i)
     end do

10   continue

! If no recognised test selected, quit program
! ----------------------------------------------------------------------------
  else
     stop 'No valid test selected'

  end if
! ----------------------------------------------------------------------------


  L1num   = 0.0_DP
  L1denom = 0.0_DP
  L2num   = 0.0_DP
  L2denom = 0.0_DP

  
! Sum error terms over all particles
! ----------------------------------------------------------------------------
  do p=1,ptot
     call record_particle_data(p,ndebug,alldata,rorigin)

     ! Find table limits with solution
     do i=1,nsol-1
        iaux = i
        if (solution(isol,i+1) > alldata(idebug)) exit
     end do

     ! Interpolate table elements for more accurate solution
     raux = solution(jsol,iaux) + &
          &(solution(jsol,iaux+1) - solution(jsol,iaux))*&
          &(alldata(idebug) - solution(isol,iaux))/&
          &(solution(isol,iaux+1) - solution(isol,iaux))
     
     L1num   = L1num   + abs(alldata(jdebug) - raux)
     L1denom = L1denom + abs(raux)
     L2num   = L2num   + abs(alldata(jdebug) - raux)**2
     L2denom = L2denom + abs(raux)**2

  end do
! ----------------------------------------------------------------------------

  L1norm = L1num/real(ptot,PR)
  L2norm = sqrt(L2num/real(ptot,PR))

! Output all results to screen
  write(6,*) "L1norm : ",L1norm
  write(6,*) "L2norm : ",L2norm

  deallocate(solution)

  stop
END PROGRAM error_norm
