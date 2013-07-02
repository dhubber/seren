! SEREN.F90
! C. P. Batty, D. A. Hubber, A. McLeod & A. P. Whitworth - 8/12/2006 
! Main Seren program.
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM seren
  use definitions
  implicit none

#if defined(DEBUG1) && !defined(USE_MPI)
  write(6,*) "**********************************************************************"
  write(6,*) "*            ****     ******    *****     ******   *     *           *"
  write(6,*) "*           *    *    *         *    *    *        **    *           *"
  write(6,*) "*           *         *         *    *    *        * *   *           *"
  write(6,*) "*            ****     *****     *****     ******   *  *  *           *"
  write(6,*) "*                *    *         *    *    *        *   * *           *"
  write(6,*) "*           *    *    *         *    *    *        *    **           *"
  write(6,*) "*            ****     ******    *    *    ******   *     *           *"
  write(6,*) "*                                                                    *"
  write(6,*) "*                            Version 1.5.0                           *"
  write(6,*) "*                              15/01/2013                            *"
  write(6,*) "*                                                                    *"
  write(6,*) "*   Coders           : David Hubber, Chris Batty & Andrew McLeod     *"
  write(6,*) "*   Contributions by : Thomas Bisbas, Simon Goodwin,                 *"
  write(6,*) "*                      Krisada Rawiraswattana, Dimitri Stamatellos,  *"
  write(6,*) "*                      Stefanie Walch, Anthony Whitworth             *"
  write(6,*) "*                                                                    *"
  write(6,*) "*               http://www.astro.group.shef.ac.uk/seren              *"
  write(6,*) "**********************************************************************"
#endif

! Call MPI initialization routines
#if defined(USE_MPI)
  call mpi_start_seren
#endif

! Call all important routines to set-up up Seren prior to individual sims.
  call seren_setup

! SPH simulation
! ============================================================================
#if defined(SPH_SIMULATION)
  call sph_simulation
#endif

! Hybrid SPH-Nbody simulation
! ============================================================================
#if defined(NBODY_SPH_SIMULATION)
  call nbody_sph_simulation
#endif

! N-body simulation
! ============================================================================
#if defined(NBODY_SIMULATION)
  call nbody_simulation
#endif

! Cleaning up memory
  call clean_up

! End MPI processes (print walltimes)
#if defined(USE_MPI)
  call mpi_finish(.false.)
#endif

  stop
END PROGRAM seren
