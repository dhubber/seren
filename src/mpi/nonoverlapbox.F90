! NONOVERLAPBOX.F90
! A. McLeod - 03/07/08
! Find a box which is not overlapped by other domains, starting with the local
! domain box
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine nonoverlapbox()
   use mpi
   use mpi_communication_module
   use particle_module
   use domain_comparisons
   implicit none

   integer :: d                            ! Loop counter over domains
   integer :: ind(1:2)                     ! Index of lowest value
   integer :: i                            ! Loop counter over ghost domains
   integer :: iloop                        ! Range of loop over ghost domains
   integer :: k                            ! Loop counter over dimensions
   integer :: p                            ! Loop counter over particles
   integer :: xyz(1:NDIM)                  ! Loop counters
   integer :: x, y, z                      ! Loop counters
   integer :: n_ghosts                     ! Number of domain ghosts
   real(kind=PR) :: A(1:NDIM)              ! Area of domain
   real(kind=PR) :: aux                    ! Temporary variable
   real(kind=PR) :: ds(1:NDIM)             ! Width of domain
   real(kind=PR) :: d_bbmin(1:NDIM)        ! Domain box minima
   real(kind=PR) :: d_bbmax(1:NDIM)        ! Domain box maxima
   real(kind=PR) :: d_ghost_bbmin(1:NDIM,27) ! All possible ghosts of domain box
   real(kind=PR) :: d_ghost_bbmax(1:NDIM,27) ! All possible ghosts of domain box
   real(kind=PR) :: localcut(2,1:NDIM)     ! Amount by which a cut is needed
   real(kind=PR), allocatable :: r(:,:)    ! Particle positions
   real(kind=PR) :: r_max(1:NDIM)          ! maximum positions
   real(kind=PR) :: r_min(1:NDIM)          ! minimum positions
   real(kind=PR) :: V                      ! Volume of domain

   debug2("Finding non-overlapped boxes [nonoverlapbox.F90]")
   debug_timing("NONOVERLAPBOX")

#ifdef DEBUG_BOXES
   write (6,*) "********* NON_OVERLAP_BOX *********"
#endif

   localcut = BIG_NUMBER

! Find local domain box which is unique (does not overlap other domain boxes)
   local_min = domain_bbmin(1:NDIM,rank)
   local_max = domain_bbmax(1:NDIM,rank)
#if defined(PERIODIC)
#if defined(PERIODIC_X)
   local_min(1) = max(local_min(1), periodic_min(1))
   local_max(1) = min(local_max(1), periodic_max(1))
#endif
#if defined(PERIODIC_Y) && (NDIM == 2 || NDIM == 3)
   local_min(2) = max(local_min(2), periodic_min(2))
   local_max(2) = min(local_max(2), periodic_max(2))
#endif
#if defined(PERIODIC_Z) && (NDIM == 3)
   local_min(3) = max(local_min(3), periodic_min(3))
   local_max(3) = min(local_max(3), periodic_max(3))
#endif
#endif
   allocate(r(1:NDIM,1:ptot))
   do p=1,ptot
      r(1:NDIM,p) = sph(p)%r(1:NDIM)
   end do
   call bounding_box(ptot,r,r_max(1:NDIM),r_min(1:NDIM))
   deallocate(r)
#ifdef DEBUG_BOXES
   write (6,*) "bounding_box r_min = ", r_min
   write (6,*) "bounding_box r_max = ", r_max
#endif

  ! ==========================================================================
   do d=0,lastrank
      if (d==rank) cycle
#ifdef DEBUG_NONOVERLAP
      write (6,*) "testing to rank ", d
#endif

      n_ghosts = 1 ! One ghost for original task
      d_ghost_bbmin(1:NDIM,1) = domain_bbmin(1:NDIM,d)
      d_ghost_bbmax(1:NDIM,1) = domain_bbmax(1:NDIM,d)
      
#if defined(PERIODIC)
#if defined(PERIODIC_X)
      i = 1
      if (d_ghost_bbmin(1,i) < periodic_min(1)) then
         n_ghosts = n_ghosts + 1
         d_ghost_bbmin(1:NDIM,n_ghosts) = d_ghost_bbmin(1:NDIM,i)
         d_ghost_bbmax(1:NDIM,n_ghosts) = d_ghost_bbmax(1:NDIM,i)
         d_ghost_bbmin(1,n_ghosts) = d_ghost_bbmin(1,n_ghosts) + periodic_size(1)
         d_ghost_bbmax(1,n_ghosts) = d_ghost_bbmax(1,n_ghosts) + periodic_size(1)
      end if
      if (d_ghost_bbmax(1,i) > periodic_max(1)) then
         n_ghosts = n_ghosts + 1
         d_ghost_bbmin(1:NDIM,n_ghosts) = d_ghost_bbmin(1:NDIM,i)
         d_ghost_bbmax(1:NDIM,n_ghosts) = d_ghost_bbmax(1:NDIM,i)
         d_ghost_bbmin(1,n_ghosts) = d_ghost_bbmin(1,n_ghosts) - periodic_size(1)
         d_ghost_bbmax(1,n_ghosts) = d_ghost_bbmax(1,n_ghosts) - periodic_size(1)
      end if
#endif
#if defined(PERIODIC_Y) && (NDIM == 2 || NDIM == 3)
      iloop = n_ghosts
      do i=1,iloop
         if (d_ghost_bbmin(2,i) < periodic_min(2)) then
            n_ghosts = n_ghosts + 1
            d_ghost_bbmin(1:NDIM,n_ghosts) = d_ghost_bbmin(1:NDIM,i)
            d_ghost_bbmax(1:NDIM,n_ghosts) = d_ghost_bbmax(1:NDIM,i)
            d_ghost_bbmin(2,n_ghosts) = d_ghost_bbmin(2,n_ghosts) + periodic_size(2)
            d_ghost_bbmax(2,n_ghosts) = d_ghost_bbmax(2,n_ghosts) + periodic_size(2)
         end if
         if (d_ghost_bbmax(2,i) > periodic_max(2)) then
            n_ghosts = n_ghosts + 1
            d_ghost_bbmin(1:NDIM,n_ghosts) = d_ghost_bbmin(1:NDIM,i)
            d_ghost_bbmax(1:NDIM,n_ghosts) = d_ghost_bbmax(1:NDIM,i)
            d_ghost_bbmin(2,n_ghosts) = d_ghost_bbmin(2,n_ghosts) - periodic_size(2)
            d_ghost_bbmax(2,n_ghosts) = d_ghost_bbmax(2,n_ghosts) - periodic_size(2)
         end if
      end do
#endif
#if defined(PERIODIC_Z) && (NDIM == 3)
      iloop = n_ghosts
      do i=1,iloop
         if (d_ghost_bbmin(3,i) < periodic_min(3)) then
            n_ghosts = n_ghosts + 1
            d_ghost_bbmin(1:NDIM,n_ghosts) = d_ghost_bbmin(1:NDIM,i)
            d_ghost_bbmax(1:NDIM,n_ghosts) = d_ghost_bbmax(1:NDIM,i)
            d_ghost_bbmin(3,n_ghosts) = d_ghost_bbmin(3,n_ghosts) + periodic_size(3)
            d_ghost_bbmax(3,n_ghosts) = d_ghost_bbmax(3,n_ghosts) + periodic_size(3)
         end if
         if (d_ghost_bbmax(3,i) > periodic_max(3)) then
            n_ghosts = n_ghosts + 1
            d_ghost_bbmin(1:NDIM,n_ghosts) = d_ghost_bbmin(1:NDIM,i)
            d_ghost_bbmax(1:NDIM,n_ghosts) = d_ghost_bbmax(1:NDIM,i)
            d_ghost_bbmin(3,n_ghosts) = d_ghost_bbmin(3,n_ghosts) - periodic_size(3)
            d_ghost_bbmax(3,n_ghosts) = d_ghost_bbmax(3,n_ghosts) - periodic_size(3)
         end if
      end do
#endif
#endif
      
#ifdef DEBUG_NONOVERLAP
      write (6,*) "Created ", n_ghosts, " ghost domains"
#endif

      do i=1,n_ghosts
      
         d_bbmin = d_ghost_bbmin(1:NDIM,i)
         d_bbmax = d_ghost_bbmax(1:NDIM,i)
#ifdef DEBUG_NONOVERLAP
         write (6,*) "Ghost number ", i
         write (6,*) "d_bbmin = ", d_bbmin
         write (6,*) "d_bbmax = ", d_bbmax
#endif
      
         ! ----------------------------------------------------------------------
         if (overlap_box(local_min, local_max, d_bbmin, d_bbmax)) then
#ifdef DEBUG_NONOVERLAP
            write (6,*) "rank ", d, " overlaps!"
            write (6,'(A,3F10.7)') "domain_bbmin(d) = ", domain_bbmin(1:NDIM,d)
            write (6,'(A,3F10.7)') "domain_bbmin(d) = ", domain_bbmax(1:NDIM,d)
#endif

            ! Compare with min/max particle positions
            r_max = min(r_max,local_max)
            r_min = max(r_min,local_min)
#ifdef DEBUG_NONOVERLAP
            write (6,*) "update r_min/r_max"
            write (6,*) "r_min = ", r_min
            write (6,*) "r_max = ", r_max
#endif

            ds = r_max - r_min
#ifdef DEBUG_NONOVERLAP
            write (6,*) "ds = ", ds
#endif
            V = product(ds)
#ifdef DEBUG_NONOVERLAP
            write (6,*) "V = ", V
#endif
            if (V == 0.0_PR) exit ! Box is already of zero size
            do k=1,NDIM
               ! Calculate the current projection area of the box for each axis
               if (V > 0._PR) then
                  A(k) = V / ds(k)
               else
                  A(k) = 0.0_PR
               end if
            end do
#ifdef DEBUG_NONOVERLAP
            write (6,*) "A = ", A
#endif

            ! -------------------------------------------------------------------
            do k=1,NDIM
#ifdef DEBUG_NONOVERLAP
               write (6,*) "k = ", k
#endif
               if (d_bbmin(k) > local_min(k) .AND. d_bbmax(k) < local_max(k)) then
#ifdef DEBUG_NONOVERLAP
                  write (6,*) "Cut from either side"
#endif
                  localcut(1,k) = (r_max(k) - d_bbmin(k)) * A(k)
                  localcut(2,k) = (d_bbmax(k) - r_min(k)) * A(k)
               else if (d_bbmin(k) > local_min(k)) then
#ifdef DEBUG_NONOVERLAP
                  write (6,*) "Cut from right side"
#endif
                  localcut(1,k) = (r_max(k) - d_bbmin(k)) * A(k)
                  localcut(2,k) = V
               else if (d_bbmax(k) < local_max(k)) then
#ifdef DEBUG_NONOVERLAP
                  write (6,*) "Cut from left side"
#endif
                  localcut(1,k) = V
                  localcut(2,k) = (d_bbmax(k) - r_min(k)) * A(k)
               else
#ifdef DEBUG_NONOVERLAP
                  write (6,*) "Destroyed!"
#endif
                  ! Totally ensconced
                  localcut(1,k) = V
                  localcut(2,k) = V
               end if
            end do
            ! -------------------------------------------------------------------

            ind = minloc(localcut)
            k = ind(2) ! Dimension to cut along
#ifdef DEBUG_NONOVERLAP
            write (6,*) "Cutting in dimension ", k
            write (6,*) "cutting type: ", ind(1)
#endif
            if (ind(1)==1) then
               local_max(k) = d_bbmin(k)
            else if (ind(1)==2) then
               local_min(k) = d_bbmax(k)
            end if
#ifdef DEBUG_NONOVERLAP
            write (6,*) "New local min/max"
            write (6,*) "local_min = ", local_min
            write (6,*) "local_max = ", local_max
#endif

            local_min(k) = min(local_min(k),local_max(k))

         end if
      ! ----------------------------------------------------------------------

      end do
   end do
   ! ==========================================================================

! Write local boxes to screen
#ifdef DEBUG_BOXES
   write(6,'(A,I0,A,3(F0.12,X))') "Rank : ",rank," local min = ", local_min(1:NDIM)
   write(6,'(A,I0,A,3(F0.12,X))') "Rank : ",rank," local max = ", local_max(1:NDIM)
#endif

#ifdef DEBUG2
   do d=0,lastrank
      if (d==rank) cycle
      if (overlap_domain_box(local_min, local_max, d)) then
         if (minval(abs(local_max-local_min))>SMALL_NUMBER) then
            write (6,*) "Still an overlap!!!"
            write (6,*) "Overlap in task: ", rank
            write (6,*) "with domain box of task: ", d
            write (6,*) "local_min=", local_min
            write (6,*) "local_max=", local_max
            write (6,*) "domain_bbmin=", domain_bbmin(:,d)
            write (6,*) "domain_bbmax=", domain_bbmax(:,d)
            stop
         end if
      end if
   end do
#endif

   return
end subroutine nonoverlapbox
