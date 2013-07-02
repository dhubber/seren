! MPI_H_GATHER
! A. McLeod - 18/03/2011
! Wrapper for first-pass h_gather routines in MPI
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE mpi_h_gather(p,acctot,h_not_done,h_list,h_old,typemask)
   use interface_module, only : h_gather,h_rho_iteration
   use particle_module
   use mpi_communication_module
   use type_module
   use domain_comparisons
   implicit none

   integer, intent(in) :: p              ! particle id
   integer, intent(in) :: acctot         ! Length of arrays
   integer, intent(inout)       :: h_not_done    ! Thread h_not_done
   integer, intent(inout)       :: h_list(1:acctot) ! Thread h_list
   real(kind=PR), intent(inout) :: h_old(1:acctot) !  Thread h_old list
   logical, optional, intent(in) :: typemask(1:ntypes) ! part. types to include?
   
   real(kind=PR)                :: h_old_p    ! Old value of h
   real(kind=PR)                :: hrangehguess  ! KERNRANGE * hguess
   real(kind=PR)                :: minaux(1:NDIM), maxaux(1:NDIM) ! temporary storage
   
   h_old_p = sph(p)%h

#if defined(H_RHO)
   call h_rho_iteration(p,sph(p)%h,typemask)
#elif defined(HGATHER) || defined(HMASS)
   call h_gather(p,sph(p)%h,sph(p)%r(1:NDIM),typemask)
#endif

   hrangehguess = KERNRANGE * sph(p)%h
   minaux = sph(p)%r - hrangehguess
   maxaux = sph(p)%r + hrangehguess

   ! Test if particle's smoothing length overlaps boundary of non-overlapped box
   ! Note that this also tests for particles that overlap periodic boundary
   if (.NOT. inside_box(minaux, maxaux, local_min, local_max)) then
      h_not_done = h_not_done + 1
      h_list(h_not_done) = p
      ! We are not done, save old smoothing length for later
      h_old(h_not_done) = h_old_p
   end if

   return
END SUBROUTINE mpi_h_gather
