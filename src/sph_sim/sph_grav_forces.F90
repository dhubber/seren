! SPH_GRAV_FORCES.F90
! C. P. Batty & D. A. Hubber - 11/1/2007
! Calculates gravitational accelerations for all particles 
! (including sink particles) - controls calls to relevant subroutines.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sph_grav_forces
  use interface_module, only : active_particle_list,&
       &add_external_gravitational_force,BHgrav_accel,&
       &direct_sph_gravity,direct_sink_gravity
  use type_module
  use particle_module
  use time_module
  use timing_module
  use type_module
#if defined(OPENMP)
  use omp_lib
#endif
#if defined(CELL_WALK)
  use tree_module
#endif
#if defined(USE_MPI)
  use mpi_communication_module
  use mpi
#endif
  implicit none

#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  integer :: acctot                   ! No. of particles on acc. step
  integer :: i                        ! Aux. particle counter
  integer :: p                        ! particle counter
  integer, allocatable :: acclist(:)  ! List of particles on acc. step
  real(kind=DP) :: agravp(1:NDIM)     ! Grav. acceleration of particle p
  real(kind=DP) :: atemp(1:NDIM)      ! Aux. grav acceleration vector
  real(kind=DP) :: potp               ! Grav. potential of particle p
#if defined(OPENMP) || defined(USE_MPI)
  integer :: chunksize                ! Data packet size for dynamic OpenMP
#endif
#if defined(CELL_WALK)
  logical :: activecell               ! Does cell contain active particles?
  integer :: c                        ! Cell counter
  integer :: pcounter                 ! Aux. particle counter
  integer :: pstart                   ! First particle id
#endif
#if defined(EXTERNAL_FORCE)
  real(kind=DP) :: adottemp(1:NDIM)   ! 'Jerk' aux. variable
#endif
#if defined(USE_MPI)
  integer :: ierr                     ! Return code
  integer :: grav_send(0:endranklist) ! Array of send requests (one for each OTHER task)
  integer :: grav_recv(0:endranklist) ! Array of recv requests (one for each OTHER task)
#endif

  debug2("Calculating grav. forces for all particles [sph_grav_forces.F90]")

  allocate(acclist(1:ptot))

! Calculate gravitational forces of groups using tree
! ============================================================================
#if defined(GRAVITY) && defined(SELF_GRAVITY) && defined(BH_TREE) && defined(CELL_WALK)
  debug_timing("GROUPED_GRAVITY")

  acctot = 0
  pstart = pgravitystart
  pcounter = pstart - 1
  c = 0

  ! Walk gravity tree and find leaf cells containing active particles
  do
     if (BHgrav(c)%leaf > 0) then
        activecell = .false.
        do i=1,BHgrav(c)%leaf
           p = BHgrav(c)%plist(i)
           if (sph(p)%accdo) activecell = .true.
        end do
        if (activecell) then
           acctot = acctot + 1
           acclist(acctot) = c
        end if
        c = BHgrav(c)%nextcell
     else if (BHgrav(c)%leaf == 0) then
        c = BHgrav(c)%ifopen
     else
        c = BHgrav(c)%nextcell
     end if
     if (c > ctot_grav) exit
  end do

  if (acctot > 0) then
     !$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC,1) PRIVATE(c)
     do i=1,acctot
        c = acclist(i)
        call BHgrav_grouped_walk(c)
     end do
     !$OMP END PARALLEL DO
  end if


! Loop over all self-gravitating particles and calculate forces
! ============================================================================
#elif defined(GRAVITY)

! For multiple particle timesteps, first make a list of all self-gravitating
! SPH particles on an acceleration step, and then parallelize over that list.
  call active_particle_list(acctot,acclist,gravmask)
#if defined(OPENMP)
  chunksize = int(CHUNKFRAC*real(acctot,PR)/real(omp_get_max_threads(),PR)) + 1
#endif
#if defined(TIMING)
  ngravcomp = ngravcomp + int(acctot,ILP)
#endif

#if defined(USE_MPI) && defined(SELF_GRAVITY)
  ! Export particles to other domains
  call gravity_export_to(grav_send,grav_recv,acclist,acctot,chunksize)

! Do gravity calculations for exported particles
  call do_export_gravity(grav_recv)

  ! Return exported particles, and receive gravity forces for particles
  ! exported from this domain; merge them with the local calculations
  call gravity_export_return(grav_send)
#endif

  ! --------------------------------------------------------------------------
  if (acctot > 0) then
     debug_timing("GRAVITY_FORCES")

     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) &
     !$OMP PRIVATE(agravp,atemp,p,potp)
     do i=1,acctot
        p = acclist(i)
#if defined(USE_MPI)
#if !defined(LOW_MEM)
        sph(p)%gravity_calcs = 0._PR
#endif
        agravp(1:NDIM) = sph(p)%a      ! This should be the 'remote' component
#else
        agravp(1:NDIM) = 0.0_DP
#endif
        atemp(1:NDIM)  = 0.0_DP

        ! SPH self-gravity
        ! --------------------------------------------------------------------
#if defined(SELF_GRAVITY) && defined(BH_TREE) && defined(GRAD_H_SPH)
        call BHgrav_accel(p,1.0_PR/sph(p)%h,sph(p)%zo,&
             &sph(p)%agravmag,sph(p)%r(1:NDIM),atemp(1:NDIM),potp)
#elif defined(SELF_GRAVITY) && defined(BH_TREE)
        call BHgrav_accel(p,1.0_PR/sph(p)%h,0.0_PR,&
             &sph(p)%agravmag,sph(p)%r(1:NDIM),atemp(1:NDIM),potp)
#elif defined(SELF_GRAVITY) && defined(BINARY_TREE)
        call binary_gravacc(p,1.0_PR/sph(p)%h,&
             &sph(p)%r(1:NDIM),atemp(1:NDIM),potp)
#elif defined(SELF_GRAVITY) && defined(GRAD_H_SPH) 
        call direct_sph_gravity(p,1.0_PR/sph(p)%h,sph(p)%zo,&
             &sph(p)%r(1:NDIM),atemp(1:NDIM),potp)
#elif defined(SELF_GRAVITY)   
        call direct_sph_gravity(p,1.0_PR/sph(p)%h,0.0_PR,&
             &sph(p)%r(1:NDIM),atemp(1:NDIM),potp)
#endif

#if defined(USE_MPI) && defined(SELF_GRAVITY)
        agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),PR)
        sph(p)%gpot = sph(p)%gpot + real(potp,PR) ! Adding to remote component
#else
        agravp(1:NDIM) = atemp(1:NDIM)
        sph(p)%gpot = real(potp,PR)
#endif
        ! Store gravitational potential at this point ONLY due to SPH particles
#if defined(RAD_WS)
        sph(p)%sphgpot = sph(p)%gpot
#endif

        ! Sink gravity contribution
        ! --------------------------------------------------------------------
#if defined(SINKS)
        call direct_sink_gravity(p,sph(p)%h,&
             &sph(p)%r(1:NDIM),atemp(1:NDIM),potp)
        agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),PR)
        sph(p)%gpot = sph(p)%gpot + real(potp,PR)
#endif

        ! Compute values for tree MAC
        ! --------------------------------------------------------------------
#if defined(BH_TREE) && defined(SELF_GRAVITY) && defined(EIGEN_MAC)
        sph(p)%agravmag = sph(p)%gpot
#elif defined(BH_TREE) && defined(SELF_GRAVITY) && !defined(GEOMETRIC_MAC)
        sph(p)%agravmag = real(sqrt(dot_product(agravp(1:NDIM),&
             &agravp(1:NDIM))),PR)
#endif

        ! External gravitational force contribution
        ! --------------------------------------------------------------------
#if defined(EXTERNAL_FORCE)
        call add_external_gravitational_force(real(sph(p)%r(1:NDIM),DP),&
             &real(sph(p)%v(1:NDIM),DP),atemp(1:NDIM),adottemp(1:NDIM),potp)
        agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),PR)
        sph(p)%gpot = sph(p)%gpot + real(potp,PR)
#endif

#if defined(USE_MPI)
        sph(p)%a(1:NDIM) = real(agravp(1:NDIM),PR) ! Don't double-count remote
                                                   ! gravity contribution
#else
        sph(p)%a(1:NDIM) = sph(p)%a(1:NDIM) + real(agravp(1:NDIM),PR)
#endif
#if defined(DEBUG_FORCES)
        sph(p)%a_grav(1:NDIM) = real(agravp(1:NDIM),PR)
#endif
     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------

  end if
  ! --------------------------------------------------------------------------

#endif
! ============================================================================


  deallocate(acclist)

! Special routine for ionization simulations
#if defined(HII_NO_INFALL)
  call hii_no_infall
#endif
#endif

! ----------------------------------------------------------------------------
  ! Add costs for MPI loadbalancing
#if defined(USE_MPI) && defined(GRAVITY) && !defined(LOW_MEM)
  sum_costs = sum_costs + COST_GRAVITY * &
     &sum(sph(pgravitystart:ptot)%gravity_calcs,&
     &mask=sph(pgravitystart:ptot)%accdo)
#endif

  return
END SUBROUTINE sph_grav_forces
