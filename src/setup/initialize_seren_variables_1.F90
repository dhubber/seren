! INITIALIZE_SEREN_VARIABLES_1.F90
! D. A. Hubber - 1/10/2007
! Initializes various variables in Seren before the main simulation begins.  
! Should be called immediatly after the parameters have been set/read-in.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_seren_variables_1
  use neighbour_module
  use filename_module
  use type_module
  use particle_module
  use sink_module
  use hydro_module
  use timing_module
  use periodic_module
  use diagnostics_module
  use time_module
  use tree_module
  use turbulence_module
#ifdef USE_MPI
  use mpi
  use mpi_communication_module, only : MPI_ext, numtasks
#endif
  implicit none

  logical :: flag                       ! flag if restart file exists
  character(len=256) :: restart_file    ! snapshot from restart log
#ifdef USE_MPI
  integer :: ierr
#endif

  debug1("First initialisation of seren variables [initialize_seren_variables_1.F90]")


! Calculate average neighbour number for 'grad-h' SPH scheme
! ----------------------------------------------------------------------------
#if defined(GRAD_H_SPH)
#if NDIM==1
  pp_gather = int(2.0_PR*KERNRANGE*h_fac)
#elif NDIM==2
  pp_gather = int(PI*(KERNRANGE*h_fac)**2)
#elif NDIM==3
  pp_gather = int((4.0_PR*PI/3.0_PR)*(KERNRANGE*h_fac)**3)
#endif
#endif
  pp_limit = 2*pp_gather


! Set filename variables
! ----------------------------------------------------------------------------
  if (out_file_form=="dragon_form" .or. out_file_form=="df") then
     fileform_ext = ".df"
  else if (out_file_form=="dragon_unform" .or. out_file_form=="du") then
     fileform_ext = ".du"
  else if (out_file_form=="seren_form" .or. out_file_form=="sf") then
     fileform_ext = ".sf"
  else if (out_file_form=="seren_unform" .or. out_file_form=="su") then
     fileform_ext = ".su"
  else if (out_file_form=="ascii" .or. out_file_form=="column") then
     fileform_ext = ".as"
  end if

  run_dir = trim(adjustl(run_dir))//"/"
#if defined(USE_MPI)
  out_init  = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
       &trim(adjustl(fileform_ext))//".ini.MPI."//trim(adjustl(MPI_ext))
  out_final = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
       &trim(adjustl(fileform_ext))//".fin.MPI."//trim(adjustl(MPI_ext))
  out_temp  = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
       &trim(adjustl(fileform_ext))//".tmp.MPI."//trim(adjustl(MPI_ext))
#else
  out_init  = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
       &trim(adjustl(fileform_ext))//".ini"
  out_final = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
       &trim(adjustl(fileform_ext))//".fin"
#endif
  out_temp1 = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
       &trim(adjustl(fileform_ext))//".tmp1"
  out_temp2 = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
       &trim(adjustl(fileform_ext))//".tmp2"
  restart_log = trim(adjustl(run_dir))//trim(adjustl(run_id))//".restart"
  error_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".error"


! Check if a temp restart file has been generated.  If yes, use this 
! instead of the file named in the params file.  
! ----------------------------------------------------------------------------
MPI_ROOT
  inquire(file=restart_log,exist=flag)
  if (flag) then
     open(1,file=restart_log,status='unknown',form='formatted')
     read(1,'(a)') restart_file
     read(1,'(a)') in_file_form
     in_file = restart_file
     write(6,*) "Restart file   : ",trim(restart_file)
     write(6,*) "Restart format : ",trim(in_file_form)
     if (restart_file == out_temp1) then
        ntemp = 2
     else if (restart_file == out_temp2) then
        ntemp = 1
     else
        ntemp = 1
     end if
     close(1)
     restart = .true.
     inifile = .true.
  else
     ntemp = 1
     inifile = .false.
  end if
MPI_END
#if defined(USE_MPI)
  call MPI_BCAST(ntemp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(restart, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(inifile, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

#if defined(TURBULENT_FORCING)
  nturbtemp = 1 ! If this is a restart, we correct this later
#endif

! Initialise particle counters
! ----------------------------------------------------------------------------
  pboundary = 0
  picm      = 0
  pgas      = 0
  pcdm      = 0
  pdust     = 0
  pion      = 0
  stot      = 0


! Initialize all timing variables
! ----------------------------------------------------------------------------
  nsteps    = 0
  time      = 0.0_DP
  nextsnap  = firstsnap
  lastsnap  = 0.0_DP
  snapshot  = 0
  ntempnext = 0_ILP
  ndiagnext = 0_ILP
  nsnapnext = 0_ILP
  nsinknext = 0_ILP
#if defined(TIMING)
  itime    = 0
  rtime    = 0.0_DP
  last_id  = 0
  mark_tot = 0
  marker_id(1:NBLOCKS) = " "
  ngravcomp = 0_ILP
  nhydrocomp = 0_ILP
  nsphcomp = 0_ILP
  iblock(1:NBLOCKS) = 0
  rblock(1:NBLOCKS) = 0.0_DP
#endif


! Initialize periodic variables
! ----------------------------------------------------------------------------
  periodic_size = periodic_max - periodic_min
  periodic_half = 0.5_PR*periodic_size
  ghost_search_min = periodic_min
  ghost_search_max = periodic_max
#if defined(USE_MPI) && defined(PERIODIC)
  periodic_mpi_min = periodic_min
  periodic_mpi_max = periodic_max
  if (numtasks > 0) then
#if defined(PERIODIC_X)
     periodic_mpi_min(1) = periodic_min(1) - MPI_PERIODIC_OVERLAP*periodic_size(1)
     periodic_mpi_max(1) = periodic_max(1) + MPI_PERIODIC_OVERLAP*periodic_size(1)
#endif
#if defined(PERIODIC_Y)
     periodic_mpi_min(2) = periodic_min(2) - MPI_PERIODIC_OVERLAP*periodic_size(2)
     periodic_mpi_max(2) = periodic_max(2) + MPI_PERIODIC_OVERLAP*periodic_size(2)
#endif
#if defined(PERIODIC_Z)
     periodic_mpi_min(3) = periodic_min(3) - MPI_PERIODIC_OVERLAP*periodic_size(3)
     periodic_mpi_max(3) = periodic_max(3) + MPI_PERIODIC_OVERLAP*periodic_size(3)
#endif
  end if
#endif

#if defined(PERIODIC)
  periodic_half_minval = minval(periodic_half)
  leftwall = .FALSE.
  rightwall = .FALSE.
! Force 'root' of domain tree to take the size of the periodic box
#if defined(PERIODIC_X)
  leftwall(1) = .TRUE.
  rightwall(1) = .TRUE.
#endif
#if defined(WALL_X_LHS)
  leftwall(1) = .TRUE.
#endif
#if defined(WALL_X_RHS)
  rightwall(1) = .TRUE.
#endif
#if NDIM==2 || NDIM==3
#if defined(PERIODIC_Y)
  leftwall(2) = .TRUE.
  rightwall(2) = .TRUE.
#endif
#if defined(WALL_Y_LHS)
  leftwall(2) = .TRUE.
#endif
#if defined(WALL_Y_RHS)
  rightwall(2) = .TRUE.
#endif
#endif
#if NDIM==3
#if defined(PERIODIC_Z)
  leftwall(3) = .TRUE.
  rightwall(3) = .TRUE.
#endif
#if defined(WALL_Z_LHS)
  leftwall(3) = .TRUE.
#endif
#if defined(WALL_Z_RHS)
  rightwall(3) = .TRUE.
#endif
#endif
#endif
  
! Initialize misc. variables
! ----------------------------------------------------------------------------
  call random_seed(rseed)
  rzero(1:NDIM) = 0.0_PR
  gammaone = gamma - 1.0_PR
  etot0 = 2.0_DP*BIG_NUMBER_DP
  build_tree = .TRUE.
#if defined(ARTIFICIAL_CONDUCTIVITY)
  alpha_cond = 1.0_PR
#endif

  return
END SUBROUTINE initialize_seren_variables_1
