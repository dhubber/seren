! INITIALIZE_SEREN_VARIABLES_2.F90
! D. A. Hubber - 1/10/2007
! Initializes various variables in Seren before the main simulation begins.  
! Should be called before SPH/hybrid/N-body simulation has started.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_seren_variables_2
  use particle_module
  use time_module
  use filename_module
  use type_module
  use neighbour_module
  use sink_module
#if defined(USE_MPI)
  use mpi_communication_module
  use mpi
#endif
  implicit none

  integer :: p                      ! Particle counter
#if defined(USE_MPI)
  integer :: ptot_tasks(0:lastrank) ! Number of particles in each task
  integer :: porig_slot             ! Slot to start allocating porig
  integer :: ierr                   ! MPI error value
#if defined(FIXED_HMULT_SINKRAD)
  integer :: totalpgas              ! Total pgas in simulation
#endif
#endif

  debug2("Initializing variables [initialize_seren_variables_2.F90]")


! Record original particle identifiers if not a restart
! ----------------------------------------------------------------------------
  if (.not. restart) then
#if defined(USE_MPI)
     call MPI_ALLGATHER(ptot, 1, MPI_INTEGER, ptot_tasks, 1, MPI_INTEGER, &
                       &MPI_COMM_WORLD, ierr)
     if (rank == 0) then
        porig_slot = 0
     else
        porig_slot = sum(ptot_tasks(0:rank-1))
     end if
     do p=1,ptot
        sph(p)%porig = p + porig_slot
     end do
#else
     do p=1,ptot
        sph(p)%porig = p
     end do
#endif
  end if

! Initialise all time counters and variables
! ----------------------------------------------------------------------------
  if (.not. restart) then
     nsteps    = 0
     time      = 0.0_DP
     nextsnap  = firstsnap
     lastsnap  = 0.0_DP
     snapshot  = 0
     ntempnext = ntempstep
     ndiagnext = ndiagstep
     nsnapnext = nsnapstep
     nsinknext = nsinkstep
     noutput   = noutputstep
  else
     if (time < firstsnap) then
        nextsnap = firstsnap
     else
        nextsnap = lastsnap + snaptime
        if (nextsnap < time) nextsnap = time + snaptime
     end if
     if (ntempnext <= nsteps) ntempnext = nsteps + ntempstep
     if (ndiagnext <= nsteps) ndiagnext = nsteps + ndiagstep
     if (nsnapnext <= nsteps) nsnapnext = nsteps + nsnapstep
     if (nsinknext <= nsteps) nsinknext = nsteps + nsinkstep
     if (noutput   <= nsteps) noutput   = nsteps + noutputstep
  end if
  n         = 0
  nresync   = n
  timestep  = 0.0_DP
  nmaxsteps = -1_ILP
  synchronise_all = .false.


! Misc. variables
! ----------------------------------------------------------------------------
  mgas = 0.0_DP
  do p=pgravitystart,ptot
     mgas = mgas + real(sph(p)%m,DP)
  end do
#if defined(USE_MPI)
  call MPI_ALLREDUCE(MPI_IN_PLACE, mgas, 1, &
     & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  if (.not. restart) mgas_orig = mgas
  if (.not. restart) pgas_orig = pgas
!#if defined(SINKS) && defined(SMOOTH_ACCRETION)
  if (pgas_orig > 0) mmean = real(mgas_orig,DP) / real(pgas_orig,DP)
!#endif


! Set sink radius for fixed (constant) multiple of h
#if defined(FIXED_HMULT_SINKRAD)
#if defined(USE_MPI)
  call MPI_ALLREDUCE(pgas, totalpgas, 1, &
     & MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  sinkrad = sinkrad*INVKERNRANGE*&
     &((3.0_PR*pp_gather*mgas)/(4.0_PR*PI*totalpgas*rhosink))**(ONETHIRD)
#else
  sinkrad = sinkrad*INVKERNRANGE*&
       &((3.0_PR*pp_gather*mgas)/(4.0_PR*PI*pgas*rhosink))**(ONETHIRD)
#endif
#endif

! Override sink radii for non-HMULT_SINKRAD options (important when
! using DRAGON files)
#if defined(FIXED_HMULT_SINKRAD) || defined(FIXED_ABSOLUTE_SINKRAD)
  sink(1:stot)%radius = sinkrad
#endif

! Set minimum smoothing length to sink smoothing length for smooth accretion
#if defined(SINKS) && defined(SMOOTH_ACCRETION) && defined(MINIMUM_H)
#if defined(HMULT_SINKRAD)
  if (stot > 0) then
     hmin = INVKERNRANGE*sink(1)%radius
     !hmin = INVKERNRANGE*&
     !     &((3.0_PR*pp_gather*mgas)/(4.0_PR*PI*pgas*rhosink))**(ONETHIRD)
     !hmin = 0.0_PR
  else
     hmin = 0.0_PR
  end if
#else
  hmin = INVKERNRANGE*sinkrad
#endif
#endif


  return
END SUBROUTINE initialize_seren_variables_2
