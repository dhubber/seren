! BH_TREE_REMOTE_GRAV
! D. A. Hubber - 23/01/2008
! Computes gravitational force exerted on particle p due to all other
! particles by walking the remote BH tree.  Starting on the first level,
! each cell is checked with the opening criterion
! theta = size / distance < thetamax
! If the above criterion is satisfied, the cell is not opened and the
! gravitational force is given by the centre of mass of the cell plus
! quadrupole moment correction terms.
! If the above criterion is not satisfied, then the cell is opened and
! the sub-cells are checked.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHtree_remote_grav(p,d,rp,agravp,potp,export)
  use interface_module, only : distance3,ewald_force
  use definitions
  use particle_module
!  use tree_module, only : ctot_grav, BHgrav
  use tree_module
  use mpi_communication_module
#ifdef SINKS
  use sink_module
#endif
  implicit none

  integer, intent(in) :: p                      ! Id of current particle
  integer, intent(in) :: d                      ! Domain to use
  real(kind=PR), intent(in) :: rp(1:NDIM)       ! Position of particle p
  real(kind=DP), intent(out) :: agravp(1:NDIM)  ! Gravitational acceleration
  real(kind=DP), intent(out) :: potp            ! Gravitational potential
  logical, intent(out)       :: export          ! If we need to export the particle

  integer :: c                      ! Cell counter
  integer :: i, cc                  ! Loop counters
  integer :: nlist                  ! Number of nodes in nodelist
  integer :: nodelist(1:GLISTSIZE)  ! List of nearby nodes

  real(kind=PR) :: atemp(1:NDIM)    ! Grav acceleration contribution vector
  real(kind=PR) :: dpot             ! Grav. potential contribution
  real(kind=PR) :: dr(1:NDIM)       ! Relative displacement vector
  real(kind=PR) :: drsqd            ! Distance squared

#ifdef EWALD
  real(kind=PR) :: eaccel(1:NDIM)   ! Ewald grav. acceleration
#endif
#ifndef GEOMETRIC_MAC
  real(kind=PR) :: afactor          ! 'old' acceleration factor for MACs
#endif

  debug3("Calculating remote gravity tree forces [BHtree_remote_grav.F90] for particle ", p)
  
! Zero arrays and variables
  agravp(1:NDIM) = 0._DP
  potp = 0.d0
  export = .FALSE.
  nlist = 0

! Prepare various MAC variables for tree walk
#ifndef GEOMETRIC_MAC
  afactor = 0.0_PR
#ifdef SINKS
  if (p < 0) afactor = sink(-p)%agravmag
#endif
  if (p > 0) afactor = sph(p)%agravmag
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

     ! If this is a dead cell, skip
     if (BHremote_grav(c,d)%leaf == -1) then
        c = BHremote_grav(c,d)%nextcell

     ! Else check the opening angle of the cell
     ! -----------------------------------------------------------------------
     else
        ! Calculate distance depending on whether we use Ewald gravity
#ifdef EWALD
        call distance3(BHremote_grav(c,d)%r(1:NDIM),rp(1:NDIM),dr(1:NDIM),drsqd)
#else
        dr(1:NDIM) = rp(1:NDIM) - BHgrav(c)%r(1:NDIM)
        drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif

        ! If distance between p and c is greater than min distance, and the 
        ! cell is not a leaf cell containing one particle, then calculate 
        ! gravitational acceleration due to cell
        ! --------------------------------------------------------------------
#if defined(GEOMETRIC_MAC)
        if (drsqd > BHremote_grav(c,d)%dminsqd) then
#else 
        if (drsqd > BHremote_grav(c,d)%mac*afactor .and. drsqd > BHremote_grav(c,d)%dminsqd) then
#endif

           nlist = nlist + 1
           nodelist(nlist) = c

! Move to next cell
           c = BHremote_grav(c,d)%nextcell

! If distance is smaller, open cell
! ----------------------------------------------------------------------------
        else
           if (BHremote_grav(c,d)%ifopen==-1) then
              ! If c==0, we need to export particle
              export = .TRUE.
              !if (rank==0) write (6,*) "Particle ", p, " must be exported to domain ", d
              return
           else
              c = BHremote_grav(c,d)%ifopen
           end if
        end if
     end if
! ----------------------------------------------------------------------------
     ! Loop over node list and add contributions to grav. accel.
     ! This is done when i) the node id buffer has been filled up,
     ! and/or ii) we have finished traversing the tree.
     ! -----------------------------------------------------------------------
     if (nlist >= GLISTSIZE .OR. c > cmax_remotegrav) then
        do i=1,nlist
           cc = nodelist(i)

           call BHgrav_node_accel(BHremote_grav(cc,d),rp,atemp,dpot)
           
           ! Add contribution due to cell c to summation vector
           agravp(1:NDIM) = agravp(1:NDIM) + real(atemp(1:NDIM),DP)
           potp = potp + real(dpot,DP)

        end do
#if !defined(LOW_MEM)
#if defined (QUADRUPOLE)
        if (p>0) sph(p)%gravity_calcs = sph(p)%gravity_calcs + (nlist * COST_QUADRUPOLE)
#elif defined (OCTUPOLE)
        if (p>0) sph(p)%gravity_calcs = sph(p)%gravity_calcs + (nlist * COST_OCTUPOLE)
#else
        if (p>0) sph(p)%gravity_calcs = sph(p)%gravity_calcs + (nlist * COST_MONOPOLE)
#endif
#endif
        nlist = 0
     end if

! Exit loop if we have finished traversing the tree
     !if (c > cmax_remotegrav.AND.rank==0) write (6,*)&
     !& "Particle ", p, " does not need exporting to domain ", d
     if (c > cmax_remotegrav) exit

  end do
! ----------------------------------------------------------------------------

  return
END SUBROUTINE BHtree_remote_grav
