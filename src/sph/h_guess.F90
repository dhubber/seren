! H_GUESS.F90
! C. P. Batty & D. A. Hubber - 11/5/2007
! Guesses average global value of h and assigns to all particles (used when 
! the tree is not employed).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE h_guess
  use interface_module, only : bounding_box
  use neighbour_module, only : pp_gather
  use particle_module
  use hydro_module
#if defined(USE_MPI)
  use mpi
#endif
  implicit none

  integer :: p                           ! ..
  real(kind=PR) :: hguess                ! Global guess of h
  real(kind=PR) :: vol                   ! Volume of bounding box
  real(kind=PR), allocatable :: r(:,:)   ! ..
#if defined(USE_MPI)
  integer       :: ierr                  ! MPI error value
#endif

  debug2("Calculating global guess for smoothing length [h_guess.F90]")

! Get max and min of x,y,z coordinates
  allocate(r(1:NDIM,1:ptot))
  do p=1,ptot
     r(1:NDIM,p) = sph(p)%r(1:NDIM)
  end do
  call bounding_box(ptot,r,rmax(1:NDIM),rmin(1:NDIM))
  deallocate(r)
#if defined(USE_MPI)
   call MPI_ALLREDUCE(MPI_IN_PLACE,rmin,NDIM,MPI_REAL_PR,&
      MPI_MIN,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE,rmax,NDIM,MPI_REAL_PR,&
      MPI_MAX,MPI_COMM_WORLD,ierr)
#endif
#if defined(H_RHO)
  rextent = maxval(rmax(1:NDIM) - rmin(1:NDIM))
#endif

! Make guess of h assuming constant density depending on dimensionality
! (For 'grad-h' method, pp_gather is approximated in parameters.F90)
! ----------------------------------------------------------------------------
#if NDIM==1
  vol = rmax(1) - rmin(1)
  hguess = (vol*real(pp_gather,PR))/(4.0_PR*real(ptot,PR))
! ----------------------------------------------------------------------------
#elif NDIM==2
  vol = (rmax(1) - rmin(1))*(rmax(2) - rmin(2))
  hguess = sqrt((real(pp_gather,PR)*vol)/(4.0_PR*PI*real(ptot,PR)))
! ----------------------------------------------------------------------------
#elif NDIM==3
  vol = (rmax(1) - rmin(1))*(rmax(2) - rmin(2))*(rmax(3) - rmin(3))
  hguess = ((3.0_PR*real(pp_gather,PR)*vol) / &
       & (32.0_PR*PI*real(ptot,PR)))**(ONETHIRD)
#endif
! ----------------------------------------------------------------------------
  
! Store guess of h in main arrays
  sph(1:ptot)%h = hguess

#if defined(DEBUG_H_GUESS)
  write(6,*) "Estimated smoothing length =", hguess,vol,rmax,rmin
  if (hguess < 0.0_PR) stop "Fatal error: hguess is negative"
#endif


  return
END SUBROUTINE h_guess
