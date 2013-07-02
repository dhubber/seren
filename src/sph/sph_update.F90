! SPH_UPDATE.F90
! C. P. Batty & D. A. Hubber - 18/1/2007
! Control subroutine for calculating all SPH properties, 
! i.e. smoothing lengths, neighbour lists and densities, of all particles.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_update
  use interface_module, only : active_particle_list,all_sph,&
       &bounding_box,get_neib,h_gather,h_rho_iteration
  use particle_module
  use neighbour_module
  use type_module
  use timing_module
  use tree_module
#if defined(OPENMP)
  use omp_lib
#endif
  implicit none

  integer :: acctot                     ! No. of particles on acc. step
  integer :: acchydrotot                ! No. of active hydro particles
  integer :: i                          ! Aux. particle counter
  integer :: p                          ! particle counter
  integer, allocatable :: acclist(:)    ! List of particles on acc. step
#if defined(H_RHO)
  real(kind=PR), allocatable :: r(:,:)  ! List of particle positions
#endif
#if defined(OPENMP)
  integer :: chunksize                  ! Loop chunksize for OMP
#endif

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


! Call all subroutines only when there are active particles
! ============================================================================
  if (acctot > 0) then
     
     ! h-rho iteration method for calculating h (cf. grad-h SPH)
     ! -----------------------------------------------------------------------
#if defined(H_RHO)
     debug2("Obtaining bounding box containing all particles")
     ! Get max and min of x,y,z coordinates
     allocate(r(1:NDIM,1:ptot))
     do p=1,ptot
        r(1:NDIM,p) = sph(p)%r(1:NDIM)
     end do
     call bounding_box(ptot,r,rmax(1:NDIM),rmin(1:NDIM))
     rextent = maxval(rmax(1:NDIM) - rmin(1:NDIM))
     deallocate(r)
    
     debug2("Calculating densities and smoothing lengths for all particles")
     debug_timing("H_RHO_ITERATION")
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
        call h_rho_iteration(p,sph(p)%h,typeinfo(sph(p)%ptype)%h)
     end do
     !$OMP END PARALLEL DO


     ! Standard SPH
     ! -----------------------------------------------------------------------
#elif defined(HGATHER) || defined(HMASS)
     debug2("Calculating smoothing lengths for all particles")
     debug_timing("H_GATHER")

     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
        call h_gather(p,sph(p)%h,sph(p)%r(1:NDIM),typeinfo(sph(p)%ptype)%h)
     end do
     !$OMP END PARALLEL DO

     ! Calculate second smoothing length determined by second fluid
#if defined(TWO_FLUIDS)
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
        call h_gather(p,sph(p)%h2,sph(p)%r(1:NDIM),typeinfo(sph(p)%ptype)%h2)
     end do
     !$OMP END PARALLEL DO
#endif

#endif
     ! -----------------------------------------------------------------------


     ! Calculate all SPH gather quantities
     ! -----------------------------------------------------------------------
     debug2("Calculating SPH quantities for all particles")
     debug_timing("ALL_SPH")

     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
        call all_sph(p,typeinfo(sph(p)%ptype)%h)
     end do
     !$OMP END PARALLEL DO

     
     ! Update hmax in hydro cells now smoothing lengths have been calculated
     ! -----------------------------------------------------------------------
#if defined(BH_TREE)
     debug_timing("BH_TREE")
     call BHhydro_update_hmax(cmax_hydro,ctot_hydro,&
          &ltot_hydro,first_cell_hydro,last_cell_hydro,BHhydro)
#endif


     ! Update ghost particle properties
     ! -----------------------------------------------------------------------
#if defined(GHOST_PARTICLES) && defined(PERIODIC)
     if (pghost > 0) then
        call copy_data_to_ghosts
#if defined(BH_TREE)
        call BHhydro_update_hmax(cmax_ghost,ctot_ghost,&
             &ltot_ghost,first_cell_ghost,last_cell_ghost,BHghost)
#endif
     end if
#endif


     ! Getting neighbour lists for all particles
     ! -----------------------------------------------------------------------
#if defined(NEIGHBOUR_LISTS)
     debug2("Calculating neighbour lists for all particles")
     debug_timing("GET_NEIB")

     ! Reconstruct active particle list for hydro particles only
     call active_particle_list(acchydrotot,acclist,sphmask) !hydromask)
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
     ! -----------------------------------------------------------------------


     ! Mean-h zeta for self-gravity
     ! -----------------------------------------------------------------------
#if defined(GRAD_H_SPH) && defined(SELF_GRAVITY) && defined(MEANH_GRAVITY)
     debug2("Calculating meanh zeta for all particles")
     debug_timing("MEANH_ZETA")

     !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p)
     do i=1,acctot
        p = acclist(i)
        call meanh_zeta(p)
     end do
     !$OMP END PARALLEL DO
#endif
     ! -----------------------------------------------------------------------

  end if
! ============================================================================

  if (allocated(acclist)) deallocate(acclist)

  return
END SUBROUTINE sph_update
