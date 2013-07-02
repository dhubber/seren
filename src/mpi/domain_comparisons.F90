! DOMAIN_COMPARISONS.F90
! A. McLeod - 27/10/2008
! Contains routines comparing points and/or particles with domain boundaries
! and domain bounding boxes
! ============================================================================

#include "macros.h"

! ============================================================================
module domain_comparisons
   use definitions
   use mpi_communication_module
   use periodic_module
#if defined(OPENMP)
   use omp_lib
#endif
   implicit none

   contains


   ! =========================================================================
   ! FIND_BOX_CENTRE
   ! Find the centre of box
   ! =========================================================================
   function find_box_centre(minbox, maxbox)
      real(kind=PR) :: find_box_centre(1:NDIM)     ! Central point of box
      real(kind=PR), intent(in) :: minbox(1:NDIM)   ! ..
      real(kind=PR), intent(in) :: maxbox(1:NDIM)   ! ..

      find_box_centre = 0.5_PR*(minbox + maxbox)

      return
   end function find_box_centre


   ! =========================================================================
   ! INSIDE_BOUNDING_SMOOTHING_BOX
   ! TRUE if particle in smoothing kernel bounding box
   ! =========================================================================
   function inside_bounding_smoothing_box(rp,d)
      logical :: inside_bounding_smoothing_box   ! ..
      integer, intent(in) :: d                   ! Domain of bounding box
      real(kind=PR), intent(in) :: rp(1:NDIM)     ! Particle position

      inside_bounding_smoothing_box = .NOT. &
         & any(rp(1:NDIM) < bbmin(1:NDIM,d))&
         & .OR. any(rp(1:NDIM) > bbmax(1:NDIM,d))
      
      return
   end function inside_bounding_smoothing_box


   ! =========================================================================
   ! OVERLAP_SMOOTHING_BOX
   ! TRUE if box overlaps smoothing kernel bounding box
   ! =========================================================================
   function overlap_smoothing_box(minbox, maxbox, d)
      logical                    :: overlap_smoothing_box
      real (kind=PR), intent(in) :: minbox(1:NDIM)   ! Minimum box position
      real (kind=PR), intent(in) :: maxbox(1:NDIM)   ! Maximum box position
      integer, intent(in)        :: d                ! Domain of bounding box

      overlap_smoothing_box = .NOT. &
         & (any(maxbox(1:NDIM) < bbmin(1:NDIM,d)) &
         & .OR. any(minbox(1:NDIM) > bbmax(1:NDIM,d)))
      
      return
   end function overlap_smoothing_box


   ! =========================================================================
   ! OVERLAP_GATHER_SCATTER
   ! TRUE if either
   !    box overlaps smoothing kernel bounding box
   ! OR
   !    box plus maxh overlaps active domain box
   ! =========================================================================
   function overlap_gather_scatter(minbox, maxbox, maxh, d)
      logical :: overlap_gather_scatter            ! ..
      real(kind=PR), intent(in) :: minbox(1:NDIM)   ! Minimum box position
      real(kind=PR), intent(in) :: maxbox(1:NDIM)   ! Maximum box position
      real(kind=PR), intent(in) :: maxh             ! Maximum smoothing length extent of box
      integer, intent(in)  :: d               ! Domain of bounding box
      logical :: overlap_gather              ! Does the box overlap the smoothing kernel box?
      logical :: overlap_scatter             ! Does the box plus max overlap the active box?
      real(kind=PR) :: minhbox(1:NDIM)       ! Minimum box position minus maxh
      real(kind=PR) :: maxhbox(1:NDIM)       ! Maximum box position plus maxh

      minhbox = minbox - maxh
      maxhbox = maxbox + maxh

      overlap_gather = .NOT. &
         & (any(maxbox(1:NDIM) < bbmin(1:NDIM,d)) &
         & .OR. any(minbox(1:NDIM) > bbmax(1:NDIM,d)))
      
      overlap_scatter = .NOT. &
         & (any(maxhbox(1:NDIM) < activemin(1:NDIM,d)) &
         & .OR. any(minhbox(1:NDIM) > activemax(1:NDIM,d)))
      
      overlap_gather_scatter = overlap_gather .OR. overlap_scatter
      
      return
   end function overlap_gather_scatter


   ! =========================================================================
   ! OVERLAP_DOMAIN_BOX
   ! TRUE if box overlaps domain box
   ! =========================================================================
   function overlap_domain_box(minbox, maxbox, d)
      logical                   :: overlap_domain_box
      real(kind=PR), intent(in) :: minbox(1:NDIM)   ! Minimum box position
      real(kind=PR), intent(in) :: maxbox(1:NDIM)   ! Maximum box position
      integer, intent(in)       :: d                ! Domain of bounding box

      overlap_domain_box = .NOT. &
         & (any(minbox(1:NDIM) >= domain_bbmax(1:NDIM,d)) &
         & .OR. any(maxbox(1:NDIM) <= domain_bbmin(1:NDIM,d)))
      
      return
   end function overlap_domain_box


   ! =========================================================================
   ! OVERLAP_BOX
   ! TRUE if box overlaps box
   ! =========================================================================
   function overlap_box(minbox1, maxbox1, minbox2, maxbox2)
      logical                   :: overlap_box
      real(kind=PR), intent(in) :: minbox1(1:NDIM)   ! Minimum box position
      real(kind=PR), intent(in) :: maxbox1(1:NDIM)   ! Maximum box position
      real(kind=PR), intent(in) :: minbox2(1:NDIM)   ! Minimum box position
      real(kind=PR), intent(in) :: maxbox2(1:NDIM)   ! Maximum box position

      overlap_box = .NOT. &
         & (any(minbox1(1:NDIM) >= maxbox2(1:NDIM)) &
         & .OR. any(maxbox1(1:NDIM) <= minbox2(1:NDIM)))
      
      return
   end function overlap_box


   ! =========================================================================
   ! INSIDE_BOX
   ! TRUE if box 1 entirely inside box 2
   ! =========================================================================
   function inside_box(minbox1, maxbox1, minbox2, maxbox2)
      logical                   :: inside_box
      real(kind=PR), intent(in) :: minbox1(1:NDIM)   ! Minimum box 1 position
      real(kind=PR), intent(in) :: maxbox1(1:NDIM)   ! Maximum box 1 position
      real(kind=PR), intent(in) :: minbox2(1:NDIM)   ! Minimum box 2 position
      real(kind=PR), intent(in) :: maxbox2(1:NDIM)   ! Maximum box 2 position

      inside_box = .NOT. (any(minbox1 <= minbox2(1:NDIM)) .OR. &
                        & any(maxbox1 >= maxbox2(1:NDIM)))

      return
   end function inside_box


   ! =========================================================================
   ! INSIDE_DOMAIM_BOX
   ! TRUE if point is inside domain box
   ! =========================================================================
   function inside_domain_box(rp, d)
      logical                   :: inside_domain_box
      real(kind=PR), intent(in) :: rp(1:NDIM)       ! Point position
      integer, intent(in)       :: d                ! Domain of bounding box

      inside_domain_box = .NOT. &
         & (any(rp(1:NDIM) >= domain_bbmax(1:NDIM,d)) &
         & .OR. any(rp(1:NDIM) < domain_bbmin(1:NDIM,d)))

      return
   end function inside_domain_box


   ! =========================================================================
   ! INSIDE_GATHER_SCATTER
   ! TRUE if either
   !   the particle's smoothing kernel overlaps the remote domain's active box,
   ! OR
   !   the particle lies within the other domain's smoothing kernel bounding box
   ! =========================================================================
   function inside_gather_scatter(rp,hsearch,d)
      logical :: inside_gather_scatter           ! ..
      real(kind=PR), intent(in) :: rp(1:NDIM)    ! Particle position
      real(kind=PR), intent(in) :: hsearch       ! Size of smoothing kernel
      integer, intent(in)       :: d             ! Domain of bounding box

      inside_gather_scatter = ( any(rp - activemax(1:NDIM,d) > hsearch) &
      & .OR. any(activemin(1:NDIM,d) - rp > hsearch) ) &
      & .AND. ( any(rp(1:NDIM) < bbmin(1:NDIM,d)) &
      & .OR. any(rp(1:NDIM) > bbmax(1:NDIM,d)) )

      inside_gather_scatter = .NOT. inside_gather_scatter

      return
   end function inside_gather_scatter


   ! =========================================================================
   ! MAXMIN_HYDRO_FILTER
   ! Gets maximum and minimum of r values for hydro with filtering options
   ! =========================================================================
   subroutine maxmin_hydro_filter(rmaxaux,rminaux,filter,nfilter,use_smoo)
      use particle_module
      use periodic_module
      use mpi_communication_module, only : rank
      implicit none

      real(kind=PR), intent(out) :: rmaxaux(1:NDIM)      ! maximum r value
      real(kind=PR), intent(out) :: rminaux(1:NDIM)      ! minimum r value
      integer, intent(in)        :: filter(1:nfilter) ! Filter list
      integer, intent(in)        :: nfilter           ! Length of filter list
      logical, intent(in), optional :: use_smoo       ! Use smoothing lengths?

      logical :: use_smoo_local    ! ..
      integer :: p, i              ! counters to loop over particles
#if defined(OPENMP)
      integer :: chunksize         ! Loop chunksize for OMP
#endif

      debug2("Finding maximum extent of particles [maxmin_hydro_filter.F90]")
      debug_timing("MAXMIN_HYDRO_FILTER")
      
#if defined(OPENMP)
      chunksize = int(CHUNKFRAC*real(nfilter,PR)/&
         &real(omp_get_max_threads(),PR)) + 1
#endif

      rmaxaux(1:NDIM) = -BIG_NUMBER
      rminaux(1:NDIM) =  BIG_NUMBER

      use_smoo_local = .FALSE.
      if (present(use_smoo)) then
         if (use_smoo) use_smoo_local = .TRUE.
      end if

      if (use_smoo_local) then
         !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) &
         !$OMP DEFAULT(SHARED) PRIVATE(p) &
         !$OMP REDUCTION(MIN:rminaux) REDUCTION(MAX:rmaxaux)
         do i=1,nfilter
            p = filter(i)
            rminaux = min(rminaux, sph(p)%r - KERNRANGE*sph(p)%h)
            rmaxaux = max(rmaxaux, sph(p)%r + KERNRANGE*sph(p)%h)
         end do
         !$OMP END PARALLEL DO
      else
         !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) &
         !$OMP DEFAULT(SHARED) PRIVATE(p) &
         !$OMP REDUCTION(MIN:rminaux) REDUCTION(MAX:rmaxaux)
         do i=1,nfilter
            p = filter(i)
            rmaxaux = max(rmaxaux, sph(p)%r)
            rminaux = min(rminaux, sph(p)%r)
         end do
         !$OMP END PARALLEL DO
      end if
   
      return
   end subroutine maxmin_hydro_filter

end module domain_comparisons
