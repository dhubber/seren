! SPH_UPDATE.F90
! C. P. Batty & D. A. Hubber - 18/1/2007
! Control subroutine for calculating all SPH properties, 
! i.e. smoothing lengths, neighbour lists and densities, of all particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE mpi_sph_update
  use interface_module, only : active_particle_list,all_sph,&
       &BHhydro_stock,bounding_box,get_neib,h_gather,h_rho_iteration
  use particle_module
  use neighbour_module
  use type_module
  use timing_module
  use type_module
#if defined(OPENMP)
  use omp_lib
#endif
  use mpi_communication_module
  use domain_comparisons
  use tree_module
  use mpi
  implicit none

  integer :: acctot                   ! No. of particles on acc. step
  integer :: acchydrotot              ! No. of active hydro particles
  integer :: i                        ! Aux. particle counter
  integer :: p                        ! particle counter
  integer, allocatable :: acclist(:)  ! List of particles on acc. step
  real(kind=PR), allocatable :: r(:,:)     ! ..
#if defined(GRAD_H_SPH)
  real(kind=PR) :: rmax(1:NDIM)       ! ..
  real(kind=PR) :: rmin(1:NDIM)       ! ..
#endif
#if defined(OPENMP)
  integer :: chunksize                ! Loop chunksize for OMP
#endif
  integer :: h_not_done                  ! Number of particles who need second pass
  integer, allocatable :: h_list(:)  ! List of particles not done
  real(kind=PR), allocatable :: h_old(:)  ! List of particles not done old smoothing length
  integer, allocatable :: t_h_list(:)! Thread list of particles not done
  real(kind=PR), allocatable :: t_h_old(:)  ! List of particles not done old smoothing length
  integer :: t_h_not_done            ! Thread number of particles who need second pass
  integer :: h_chunksize             ! Chunksize for second pass of h_gather

  debug2("Calculating SPH quantities [sph_update.F90]")


! Prepare list of particles which are on an acceleration step
! ----------------------------------------------------------------------------
  allocate(acclist(1:ptot))
  call active_particle_list(acctot,acclist,sphmask)
#if defined(OPENMP)
  chunksize = int(CHUNKFRAC*real(acctot,PR)/real(omp_get_max_threads(),PR)) + 1
#endif
#if defined(TIMING)
  nsphcomp = nsphcomp + int(acctot,ILP)
#endif
  sum_acctot = sum_acctot + acctot

! Set up some variables
  h_not_done = 0           ! This is incremented if any particle requires ghosts
  pghost = -1

! Set initial bounding box so that it is 'negative' and no particle will fall within it
  bbmin(1:NDIM,rank) = BIG_NUMBER
  bbmax(1:NDIM,rank) = -BIG_NUMBER

! Call all subroutines only when there are active particles
! ============================================================================
  if (acctot > 0) then

     ! Find the size of the non-overlapped box
     call nonoverlapbox()

#ifdef GRAD_H_SPH
     allocate(r(1:NDIM,1:ptot))
     do p=1,ptot
        r(1:NDIM,p) = sph(p)%r(1:NDIM)
     end do
     call bounding_box(ptot,r,rmax(1:NDIM),rmin(1:NDIM))
#endif
     allocate(h_list(1:acctot))
     allocate(h_old(1:acctot))

     ! h_gather pass 1
     ! -----------------------------------------------------------------------

!$OMP PARALLEL PRIVATE(i, p, t_h_not_done, t_h_list, t_h_old)
     allocate(t_h_list(1:acctot))
     allocate(t_h_old(1:acctot))
     t_h_not_done = 0
!$OMP DO SCHEDULE(DYNAMIC,chunksize)
     do i=1,acctot
        p = acclist(i)
        call mpi_h_gather(p,acctot,t_h_not_done,t_h_list,t_h_old,typeinfo(sph(p)%ptype)%h)
     end do
!$OMP END DO NOWAIT
!$OMP CRITICAL
     h_list(h_not_done+1:h_not_done+t_h_not_done) = t_h_list(1:t_h_not_done)
     h_old(h_not_done+1:h_not_done+t_h_not_done) = t_h_old(1:t_h_not_done)
     h_not_done = h_not_done + t_h_not_done
!$OMP END CRITICAL
     deallocate(t_h_list)
     deallocate(t_h_old)
!$OMP END PARALLEL
     if (h_not_done > 0) then
        call maxmin_hydro_periodic(bbmax(1:NDIM,rank),bbmin(1:NDIM,rank),&
           & domain_centre,ptot,h_list(1:h_not_done),h_not_done,.TRUE.)
     end if
#if defined(OPENMP)
     h_chunksize = int(CHUNKFRAC*real(h_not_done,PR)) + 1
#endif
     sum_second_h = sum_second_h + h_not_done

#ifdef PERIODIC
     call make_minmax_periodic(bbmin(1:NDIM,rank),bbmax(1:NDIM,rank))

     ! For wrapped components, if 99.999% of domain is already covered in a
     ! dimension, just unwrap and use periodic boundaries
     where (bbmin(1:NDIM,rank)>bbmax(1:NDIM,rank) .AND. bbmin(1:NDIM,rank)-bbmax(1:NDIM,rank) <= tiny_gap)
        bbmin(1:NDIM,rank) = periodic_min(1:NDIM)
        bbmax(1:NDIM,rank) = periodic_max(1:NDIM)
     end where
#endif

     debug_timing("BROADCAST_BOUNDING_BOXES 1")
     call broadcastboundingboxes()

     ! obtain ghosts
     ! -----------------------------------------------------------------------

! Walk tree (or direct-sum search) on current processor using bounding boxes
     call get_pot_list_bb

#if defined(BH_TREE)
     if (pghost > 0) then
        ! Build neighbour tree of ghost particles
        call BHghost_build

        ! Stock ghost tree
        call BHhydro_stock(cmax_ghost,ctot_ghost,&
           &first_cell_ghost,last_cell_ghost,BHghost)
     end if
#endif

     ! h_gather pass 2
     ! -----------------------------------------------------------------------
     
! Reset initial bounding box so that it is 'negative' and
! no particle will fall within it
     bbmin(1:NDIM,rank) = BIG_NUMBER
     bbmax(1:NDIM,rank) = -BIG_NUMBER

! Calculate smoothing lengths for all overlapping particles
      if (h_not_done > 0) then
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,h_chunksize) PRIVATE(p)
         do i=1,h_not_done
            p = h_list(i)
            sph(p)%h = h_old(i)
#ifdef GRAD_H_SPH
            call h_rho_iteration(p,typeinfo(sph(p)%ptype)%h)
#else
            call h_gather(p,sph(p)%h,sph(p)%r(1:NDIM),typeinfo(sph(p)%ptype)%h)
#endif
         end do
!$OMP END PARALLEL DO
      end if
      
      deallocate(h_list)
      deallocate(h_old)

! Update hmax in hydro cells now smoothing lengths have been calculated
#ifdef BH_TREE
      call BHhydro_update_hmax
#endif

     ! expand for hydro
     ! -----------------------------------------------------------------------

! We now need to expand this box to include all particles,
! not just those on an acceleration step
#ifdef HYDRO
#ifdef BH_TREE
      call BHhydrowalkreverse_bb(1,ptot)
#else
      debug_timing("EXPAND_BBMIN_BBMAX")
      bbmin(1:NDIM,rank) = BIG_NUMBER
      bbmax(1:NDIM,rank) = -BIG_NUMBER
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p)&
!$OMP& REDUCTION(MIN:bbmin) REDUCTION(MAX:bbmax)
      do p=1,ptot
        bbmin(1:NDIM,rank) = min(bbmin(1:NDIM,rank),parray(1:NDIM,p)-KERNRANGE*parray(SMOO,p))
        bbmax(1:NDIM,rank) = max(bbmax(1:NDIM,rank),parray(1:NDIM,p)+KERNRANGE*parray(SMOO,p))
      end do
!$OMP END PARALLEL DO
#endif
#endif

#ifdef PERIODIC
      call make_minmax_periodic(bbmin(1:NDIM,rank),bbmax(1:NDIM,rank))

      ! For wrapped components, if 99.999% of domain is already covered in a
      ! dimension, just unwrap and use periodic boundaries
      where (bbmin(1:NDIM,rank)>bbmax(1:NDIM,rank) .AND. bbmin(1:NDIM,rank)-bbmax(1:NDIM,rank) <= tiny_gap)
         bbmin(1:NDIM,rank) = periodic_min(1:NDIM)
         bbmax(1:NDIM,rank) = periodic_max(1:NDIM)
      end where
#endif

! Broadcast this new bounding box which is used to find potential hydro
! particle neighbours
      debug_timing("BROADCAST_BOUNDING_BOXES 2")
      call broadcastboundingboxes()

     ! Getting neighbour lists for all particles 
#if defined(NEIGHBOUR_LISTS) && defined(HYDRO)
     debug2("Calculating neighbour lists for all particles")
     debug_timing("GET_NEIB")

     ! Reconstruct active particle list for hydro particles only
     call active_particle_list(acchydrotot,acclist,hydromask)
#if defined(OPENMP)
     chunksize = int(CHUNKFRAC*real(acchydrotot,PR)/&
          &real(omp_get_max_threads(),PR)) + 1
#endif

     if (acchydrotot > 0) then
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) &
        !$OMP DEFAULT(SHARED) PRIVATE(p)
        do i=1,acchydrotot
           p = acclist(i)
           call get_neib(p)
        end do
        !$OMP END PARALLEL DO
     end if
#endif

  end if
! ============================================================================

  if (allocated(acclist)) deallocate(acclist)

#ifdef BH_TREE
!   if (allocated(BHghost)) then
!      deallocate(BHghost)
! !      deallocate(BHghost_cell)
!   end if
#endif
  pghost = -1

  return
END SUBROUTINE mpi_sph_update
