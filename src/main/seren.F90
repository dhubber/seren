! ============================================================================
! SEREN.F90
! 
! Copyright (C) 2006 - 2013  C. P. Batty, D. A. Hubber, A. McLeod,
!                            A. P. Whitworth, T. Bisbas, S. P. Goodwin,
!                            K. Rawiraswattana, D. Stamatellos, S. Walch
!
! This file is part of SEREN.
!
! SEREN is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
!
! Foobar is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with SEREN.  If not, see <http://www.gnu.org/licenses/>.
!
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
  write(6,*) "*                            Version 1.6.0                           *"
  write(6,*) "*                              17/07/2013                            *"
  write(6,*) "*                                                                    *"
  write(6,*) "*   Coders           : David Hubber, Chris Batty & Andrew McLeod     *"
  write(6,*) "*   Contributions by : Thomas Bisbas, Simon Goodwin,                 *"
  write(6,*) "*                      Krisada Rawiraswattana, Dimitri Stamatellos,  *"
  write(6,*) "*                      Stefanie Walch, Anthony Whitworth             *"
  write(6,*) "*                                                                    *"
  write(6,*) "*                   https://github.com/dhubber/seren                 *"
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
