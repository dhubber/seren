! BHGRAV_ACCEL.F90
! D. A. Hubber - 23/01/2008
! Computes gravitational force exerted on particle p due to all other 
! particles by walking the BH tree.  Starting on the first level, each 
! cell is checked with one the opening criterion
! i) Opening angle
! ii) Octupole moment error term (GADGET MAC)
! iii) Maximum absolute multipole moment error (Salmon & Warren MAC)
! iv) Maximum quadrupole error using eigenvalues
! If the criterion is satisfied, the cell is not opened and the 
! gravitational force is given by the centre of mass of the cell plus 
! quadrupole moment correction terms.  If the above criterion is not 
! satisfied, then the cell is opened and the sub-cells are checked.  
! If a leaf cell is opened, then the contributions due to the particles 
! within the cell are added directly. 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHgrav_accel(p,invhp,zo_p,agravmag_p,rp,agravp,potp)
  use interface_module, only : distance3,ewald_force,&
       &gravity_gradh,gravity_nbody,gravity_sph,gravity_meanh,wpot
  use definitions
  use tree_module
  use particle_module
#if defined(SINKS)
  use sink_module
#endif
#if defined(USE_MPI)
  use mpi_communication_module
#endif
  implicit none

  integer, intent(in) :: p                      ! Id of current particle
  real(kind=PR), intent(in) :: invhp            ! Smoothing length of p
  real(kind=PR), intent(in) :: agravmag_p       ! Mag. of grav accel.
  real(kind=PR), intent(in) :: rp(1:NDIM)       ! Position of particle p
  real(kind=DP), intent(out) :: agravp(1:NDIM)  ! Gravitational accelertation
  real(kind=DP), intent(out) :: potp            ! Gravitational potential
  real(kind=PR), optional, intent(in) :: zo_p   ! zeta/omega for grad-h grav.

  integer :: c                      ! Cell counter
  integer :: i                      ! Auxilary counter
  integer :: nlist                  ! Number of particles in plist
  integer :: pp                     ! Second particle identifier
  integer :: plist(1:GLISTSIZE)     ! List of nearby particle ids
  real(kind=PR) :: atemp(1:NDIM)    ! Grav acceleration contribution vector
  real(kind=PR) :: dr(1:NDIM)       ! Relative displacement vector
  real(kind=PR) :: dpot             ! Grav. potential contribution
  real(kind=PR) :: drsqd            ! Distance squared
  real(kind=PR) :: invdrmag         ! ( 1 / drmag )
  real(kind=PR) :: invdrsqd         ! ( 1 / drsqd )
#if defined(EWALD)
  real(kind=PR) :: eaccel(1:NDIM)   ! Ewald grav. acceleration
#endif
#if defined(MEANH_GRAVITY)
  real(kind=PR) :: hp               ! Smoothing length
#endif
#ifndef GEOMETRIC_MAC
  real(kind=PR) :: afactor          ! 'old' acceleration factor for MACs
#endif
#if defined(USE_MPI) && !defined(LOW_MEM)
  real(kind=PR) :: gravity_calcs    ! Temporary value for sph(p)%gravity_calcs
#endif

  debug3("Calculating gravity forces [BHgrav_accel.F90] for particle ", p)

! Zero arrays and variables
  agravp(1:NDIM) = 0.0_DP
  potp = 0.0_DP
  if (p > 0) potp = sph(p)%m*invhp*wpot(0.0_PR)
  nlist = 0
#if defined(MEANH_GRAVITY)
  hp = 1.0_PR / invhp
#endif
#if defined(USE_MPI) && !defined(LOW_MEM)
  gravity_calcs = 0.0_PR
#endif

! Prepare various MAC variables for tree walk
#ifndef GEOMETRIC_MAC
  afactor = agravmag_p
  if (afactor > SMALL_NUMBER) then
#if defined(GADGET_MAC)
     afactor = afactor**(-ONETHIRD)
#elif defined(GADGET2_MAC)
     afactor = 1.0_PR/sqrt(afactor)
#elif defined(EIGEN_MAC)
     afactor = (1.0_PR/afactor)**(2.0_PR*ONETHIRD)
#endif
  else
     afactor = 0.0_PR
  end if
#endif

! Start on cell 0 and work our way down
  c = 0

! Walk gravity tree until we reach end cell pointer
! ============================================================================
  do

     ! If cell is leaf cell with just one particle, more efficient to 
     ! simply record particle id in list straight away
     ! -----------------------------------------------------------------------
     if (BHgrav(c)%leaf == 1) then
        pp = BHgrav(c)%plist(1)
        if (pp /= p) then
           nlist = nlist + 1
           plist(nlist) = pp
        end if

        ! Point to next cell in list
        c = BHgrav(c)%nextcell
        
     ! Else check the opening angle of the cell
     ! -----------------------------------------------------------------------
     else

        ! Calculate distance depending on whether we use Ewald gravity
#if defined(EWALD)
        call distance3(BHgrav(c)%r(1:NDIM),rp(1:NDIM),dr(1:NDIM),drsqd)
#else
        dr(1:NDIM) = rp(1:NDIM) - BHgrav(c)%r(1:NDIM)
        drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif

        ! If distance between p and c is greater than min distance, and the 
        ! cell is not a leaf cell containing one particle, then calculate 
        ! gravitational acceleration due to cell
        ! --------------------------------------------------------------------
#if defined(GEOMETRIC_MAC)
        if (drsqd > BHgrav(c)%dminsqd) then
#else 
        if (drsqd > BHgrav(c)%mac*afactor .and. drsqd > BHgrav(c)%dminsqd) then
#endif        

           ! Calculate gravitational contribution from cell c
           call BHgrav_node_accel(BHgrav(c),rp(1:NDIM),atemp(1:NDIM),dpot)

           ! Add contribution due to cell c to summation vector
           agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),DP)
           potp = potp + real(dpot,DP)

           ! Move to next cell
           c = BHgrav(c)%nextcell

        ! If distance is smaller, open cell
        ! --------------------------------------------------------------------
        else
	
           ! Check if cell is actually a leaf cell
           if (BHgrav(c)%leaf > 0) then

              ! If so, loop over all particles in leaf cell
              do i=1,BHgrav(c)%leaf
                 pp = BHgrav(c)%plist(i)
                 if (pp == p) cycle
                 nlist = nlist + 1
                 plist(nlist) = pp
              end do

              ! Point to next cell in list
              c = BHgrav(c)%nextcell

           ! If it's not a leaf cell, open cell.  If it's a dead cell 
           ! (e.g. due to particle accretion), point to next cell in list.
           else if (BHgrav(c)%leaf == 0) then
              c = BHgrav(c)%ifopen
           else
              c = BHgrav(c)%nextcell
           end if
        end if
        ! --------------------------------------------------------------------

     end if
     ! -----------------------------------------------------------------------

     ! Loop over particle list and add contributions to grav. accel.
     ! This is done when i) the particle id buffer has been filled up, 
     ! and/or ii) we have finished traversing the tree. 
     ! -----------------------------------------------------------------------
     if ((nlist > GLISTSIZE - LEAFMAX .or. &
          &(nlist > 0 .and. c > ctot_grav))) then
        do i=1,nlist
           pp = plist(i)
#if defined(N_BODY)
           call gravity_nbody(sph(pp)%m,rp(1:NDIM),&
                &sph(pp)%r(1:NDIM),atemp(1:NDIM),dpot)
#elif defined(MEANH_GRAVITY) && defined(GRAD_H_SPH)
           if (p < 0) then
              call gravity_meanh(0.5_PR*(hp + sph(pp)%h),sph(pp)%m,&
                   &rp(1:NDIM),sph(pp)%r(1:NDIM),atemp(1:NDIM),dpot)
           else
              call gravity_gradh_meanh(hp,sph(pp)%h,sph(pp)%m,&
                   &rp(1:NDIM),sph(pp)%r(1:NDIM),zo_p,sph(pp)%zo,&
                   &atemp(1:NDIM),dpot)
           end if
#elif defined(GRAD_H_SPH)
           if (p < 0) then
              call gravity_sph(invhp,sph(pp)%h,sph(pp)%m,&
                   &rp(1:NDIM),sph(pp)%r(1:NDIM),atemp(1:NDIM),dpot)
           else
              call gravity_gradh(invhp,sph(pp)%h,sph(pp)%m,rp(1:NDIM)&
                   &,sph(pp)%r(1:NDIM),zo_p,sph(pp)%zo,atemp(1:NDIM),dpot)
           end if
#elif defined(MEANH_GRAVITY)
           call gravity_meanh(0.5_PR*(hp + sph(pp)%h),sph(pp)%m,&
                &rp(1:NDIM),sph(pp)%r(1:NDIM),atemp(1:NDIM),dpot)
#else
           call gravity_sph(invhp,sph(pp)%h,sph(pp)%m,&
                &rp(1:NDIM),sph(pp)%r(1:NDIM),atemp(1:NDIM),dpot)
#endif
           agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),DP)
           potp = potp + real(dpot,DP)
           
#if defined(USE_MPI) && !defined(LOW_MEM)
           gravity_calcs = gravity_calcs + COST_MONOPOLE
#endif

        end do
        nlist = 0
     end if
     ! -----------------------------------------------------------------------

     ! Exit loop if we have finished traversing the tree
     if (c > ctot_grav) exit

  end do
! ============================================================================

! Add to particle gravity_calcs for MPI loadbalancing instrumentation
#if defined(USE_MPI) && !defined(LOW_MEM)
  if (p > 0) sph(p)%gravity_calcs = sph(p)%gravity_calcs + gravity_calcs
#endif

  return
END SUBROUTINE BHgrav_accel
