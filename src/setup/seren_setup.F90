! SEREN_SETUP.F90
! C. P. Batty, D. A. Hubber, A.McLeod & A. P. Whitworth - 8/12/2006
! Run all main Seren set-up routines.  Used for all simulation modes 
! (i.e. SPH, N-body or hybrid SPH/N-body mode).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE seren_setup
  use interface_module, only : paramstore,read_data,read_parameters,&
       &w0,w1,w2,wgrav,wpot,womega,wzeta
  use definitions
  use filename_module
  use seren_sim_module
  use particle_module
  use type_module
  implicit none

  character(len=256) :: store_file   ! parameters output
#if defined(DEBUG_KERNEL)
  integer :: i                       ! Aux. counter
  real(kind=PR) :: s                 ! Aux. real variable
#endif

! Read in and process command line arguments
  incomplete_params = .false.
  call read_arguments

! Set default parameters in case we read-in an incomplete parameters file
  call default_parameters

! Reading parameter file
  call read_parameters(param_file)

! Initialise some variables using parameters before IC file is read in
  call initialize_seren_variables_1
  call set_default_particle_types

! Checking compiler flags and parameter values before proceeding further
  call sanitycheck

! Setting up scaling units for simulation
  call units

#if defined(USE_MPI)
! Root task will decompose domain into new files.  Then data is read.
  call mpi_setup_decomposition
#else
! Reading input snapshot data file
  call read_data(in_file,in_file_form)
#endif

MPI_ROOT
! Writing compiler flags and parameters to file for record
  store_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".params"
  call paramstore(store_file)
MPI_END

! Setting up variables for different particle types
  call types

! Converting all important physical variables to dimensionless code units
  call convert_to_code_units_1(.FALSE.)
  call convert_to_code_units_2

! Initialising SPH kernel tables
#if defined(KERNEL_TABLES)
  call tabulate_kernel_functions
#endif

! Output kernel tables to file "kernel.dat"
#if defined(DEBUG_KERNEL)
  write(6,*) "Writing kernel tables to file"
  open(1, file="kernel.dat", status="unknown", form="formatted")
  do i=0,2*KERNTOT
     s = KERNRANGE*real(i,PR)/real(KERNTOT,PR)
     write(1,'(8E15.7)') s,w0(s),w1(s),w2(s),&
          &wgrav(s),wpot(s),womega(s),wzeta(s)
  end do
  close(1)
#endif

! Read in stellar model table
#if defined(STELLAR_FEEDBACK)
  call read_stellar_model_table
#endif

! Read in opacity tables for radiation transport
#if defined(HYDRO) && defined(SELF_GRAVITY) && defined(RAD_WS)
  call read_cooling_table_ws
#endif

! Create Ewald correction table
#if defined(GRAVITY) && defined(EWALD)
  call ewald_init
#endif

! Calculate COM of system, and if required, change to COM frame
  call COM

! Initialise all other particles arrays
  call initialize_seren_variables_2

! Initialise turbulent forcing
#if defined(TURBULENT_FORCING)
  call turb_force_init
#endif
  
  return
END SUBROUTINE seren_setup
