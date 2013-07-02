! DEFINITIONS.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Definition of precision kind variables
! ============================================================================

#include "macros.h"

! ============================================================================
MODULE definitions
#ifdef USE_MPI
  use mpi
#endif

  integer, parameter :: DP = selected_real_kind(p=15) ! double precision
  integer, parameter :: QP = selected_real_kind(p=33) ! quadruple precision
  integer, parameter :: SP = selected_real_kind(p=6)  ! single precision
  integer, parameter :: ILP = selected_int_kind(r=15) ! integer long precision

  integer, parameter :: LARGEST_INT = huge(0)         ! Largest default integer
  
#if defined(QUADRUPLE_PRECISION)
  integer, parameter :: PR = QP
#elif defined(DOUBLE_PRECISION)
  integer, parameter :: PR = DP
#else
  integer, parameter :: PR = SP                       ! default = single
#endif

#ifdef USE_MPI
  integer :: rank                    ! Rank (id) of mpi task
#ifdef DOUBLE_PRECISION
  integer, parameter :: MPI_REAL_PR = MPI_DOUBLE_PRECISION
                                               ! MPI_REAL - double precision
  integer, parameter :: MPI_COMPLEX_PR = MPI_DOUBLE_COMPLEX
                                               ! MPI_COMPLEX - double precision
#else
  integer, parameter :: MPI_REAL_PR = MPI_REAL ! MPI_REAL - single precision
  integer, parameter :: MPI_COMPLEX_PR = MPI_COMPLEX 
                                               ! MPI_COMPLEX - single precision
#endif
  integer            :: MPI_INTEGER_ILP        ! MPI_INTEGER - long precision
#endif

END MODULE definitions


! ============================================================================
