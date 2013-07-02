! BHGRAV_ACCEL_JERK.F90
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
SUBROUTINE BHgrav_accel_jerk(s,hs,rs,vs,agravs,adots,pots)
  use interface_module, only : distance3_dp,ewald_force,gravity_hermite4_meanh
  use particle_module, only : sph,ptot
  use definitions
  use tree_module
  use Nbody_module
  implicit none

  integer, intent(in) :: s                      ! Id of star
  real(kind=DP), intent(in) :: hs               ! Smoothing length of p
  real(kind=DP), intent(in) :: rs(1:NDIM)       ! Position of particle p
  real(kind=DP), intent(in) :: vs(1:NDIM)       ! Velocity of particle p
  real(kind=DP), intent(out) :: agravs(1:NDIM)  ! Gravitational acceleration
  real(kind=DP), intent(out) :: adots(1:NDIM)   ! Gravitational 'jerk'
  real(kind=DP), intent(out) :: pots            ! Gravitational potential

  integer :: c                       ! Cell counter
  integer :: i                       ! Auxilary counter
  integer :: nlist                   ! Number of particles in treelist
  integer :: pp                      ! Second particle identifier
  integer :: ss                      ! Secondary star counter
  integer :: treelist(1:GLISTSIZE)   ! List of nearby particle ids
  real(kind=DP) :: adottemp(1:NDIM)  ! ..
  real(kind=DP) :: atemp(1:NDIM)     ! Grav acceleration contribution vector
  real(kind=DP) :: dr(1:NDIM)        ! Relative displacement vector
  real(kind=DP) :: dpot              ! Grav. potential contribution
  real(kind=DP) :: drsqd             ! Distance squared
  real(kind=DP) :: drdt              ! Rate of change of radial distance
  real(kind=DP) :: dv(1:NDIM)        ! Relative velocity
  real(kind=DP) :: invdrmag          ! ( 1 / drmag )
  real(kind=DP) :: invdrsqd          ! ( 1 / drsqd )
#if defined(QUADRUPOLE)
  real(kind=DP) :: invdr5            ! ( 1 / drmag^5 )
  real(kind=DP) :: qscalar           ! Inner product of quad tensor
#endif
#if defined(OCTUPOLE)
  real(kind=DP) :: sscalar           ! Scalar component of octupole terms
#endif
#if defined(EWALD)
  real(kind=DP) :: eaccel(1:NDIM)    ! Ewald grav. acceleration
#endif
#ifndef GEOMETRIC_MAC
  real(kind=DP) :: afactor           ! 'old' acceleration factor for MACs
#endif

  debug3("Calculating gravity forces [BHtreegravity.F90] for particle ", p)

! Zero arrays and variables
  agravs(1:NDIM) = 0.0_DP
  adots(1:NDIM) = 0.0_DP
  pots = 0.0_DP
  nlist = 0

! Prepare various MAC variables for tree walk
#ifndef GEOMETRIC_MAC
  afactor = star(s)%agravmag
  if (afactor <= SMALL_NUMBER_DP) afactor = 0.0_DP
  if (afactor > SMALL_NUMBER) then
#if defined(GADGET_MAC)
     afactor = afactor**(-ONETHIRD_DP)
#elif defined(GADGET2_MAC)
     afactor = 1.0_DP/sqrt(afactor)
#elif defined(EIGEN_MAC)
     afactor = (1.0_DP/afactor)**(2.0_DP*ONETHIRD_DP)
#endif
  else
     afactor = 0.0_DP
  end if
#endif

! Start on cell 0 and work our way down
  c = 0

! Walk gravity tree until we reach end cell pointer
! ============================================================================
  do

     ! If cell is leaf cell with just one particle, more efficient to 
     ! simply record particle id in list straight away
     if (BHgrav(c)%leaf == 1) then
        pp = BHgrav(c)%plist(1)
!        if (pp /= p) then
           nlist = nlist + 1
           treelist(nlist) = pp
!        end if

        ! Point to next cell in list
        c = BHgrav(c)%nextcell
        
     ! Else check the opening angle of the cell
     ! -----------------------------------------------------------------------
     else

        ! Calculate distance depending on whether we use Ewald gravity
#if defined(EWALD)
        call distance3_dp(real(BHgrav(c)%r(1:NDIM),DP),&
             &rs(1:NDIM),dr(1:NDIM),drsqd)
#else
        dr(1:NDIM) = rs(1:NDIM) - BHgrav(c)%r(1:NDIM)
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
           invdrsqd = 1.0_DP / (drsqd + SMALL_NUMBER_DP)
           invdrmag = sqrt(invdrsqd)
#if defined(CELL_VELOCITIES)
           dv(1:NDIM) = vs(1:NDIM) - BHgrav(c)%v(1:NDIM)
           drdt = dot_product(dv(1:NDIM),dr(1:NDIM))*invdrmag
           !adottemp(1:NDIM) = -BHgrav(c)%m*invdrsqd*invdrmag*dv(1:NDIM) + &
           !     & 3.0_DP*BHgrav(c)%m*drdt*invdrsqd*invdrsqd*dr(1:NDIM)
#endif
           atemp(1:NDIM) = -BHgrav(c)%m*invdrsqd*invdrmag*dr(1:NDIM)
           adottemp(1:NDIM) = 0.0_DP
           dpot = BHgrav(c)%m*invdrmag

           ! Add quadrupole moment correction terms
#if defined(QUADRUPOLE)
           invdr5 = invdrsqd*invdrsqd*invdrmag
#if NDIM==2
           qscalar = BHgrav(c)%q(1)*dr(1)*dr(1) + BHgrav(c)%q(3)*dr(2)*dr(2) &
                    & + 2.0_DP*BHgrav(c)%q(2)*dr(1)*dr(2)
           atemp(1) = atemp(1) + (BHgrav(c)%q(1)*dr(1) + BHgrav(c)%q(2)*dr(2))&
                     &*invdr5 - 2.5_DP*qscalar*dr(1)*invdr5*invdrsqd
           atemp(2) = atemp(2) + (BHgrav(c)%q(2)*dr(1) + BHgrav(c)%q(3)*dr(2))&
                     &*invdr5 - 2.5_DP*qscalar*dr(2)*invdr5*invdrsqd
#elif NDIM==3
           qscalar = BHgrav(c)%q(1)*dr(1)*dr(1) + BHgrav(c)%q(3)*dr(2)*dr(2) &
                & - (BHgrav(c)%q(1) + BHgrav(c)%q(3))*dr(3)*dr(3) &
                & + 2.0_DP*(BHgrav(c)%q(2)*dr(1)*dr(2) &
                & + BHgrav(c)%q(4)*dr(1)*dr(3) + BHgrav(c)%q(5)*dr(2)*dr(3))
           atemp(1) = atemp(1) + (BHgrav(c)%q(1)*dr(1) + BHgrav(c)%q(2)*dr(2)&
                & + BHgrav(c)%q(4)*dr(3))*invdr5 &
                & - 2.5_DP*qscalar*dr(1)*invdr5*invdrsqd
           atemp(2) = atemp(2) + (BHgrav(c)%q(2)*dr(1)+BHgrav(c)%q(3)*dr(2)&
                & + BHgrav(c)%q(5)*dr(3))*invdr5 &
                & - 2.5_DP*qscalar*dr(2)*invdr5*invdrsqd
           atemp(3) = atemp(3) + (BHgrav(c)%q(4)*dr(1) + BHgrav(c)%q(5)*dr(2)&
                & - (BHgrav(c)%q(1)+BHgrav(c)%q(3))*dr(3))*invdr5 &
                & - 2.5_DP*qscalar*dr(3)*invdr5*invdrsqd
#endif
           dpot = dpot + 0.5_DP*qscalar*invdr5 
#endif
           ! Add octupole moment correction terms
#if defined(OCTUPOLE)
           sscalar = BHgrav(c)%s(1)*dr(1)*dr(1)*dr(1) + &
                   & BHgrav(c)%s(2)*dr(1)*dr(1)*dr(2) + &
                   & BHgrav(c)%s(3)*dr(2)*dr(2)*dr(1) + &
                   & BHgrav(c)%s(4)*dr(2)*dr(2)*dr(2) + &
                   & BHgrav(c)%s(5)*dr(3)*dr(3)*dr(1) + &
                   & BHgrav(c)%s(6)*dr(3)*dr(3)*dr(2) + &
                   & BHgrav(c)%s(7)*dr(3)*dr(3)*dr(3) + &
                   & BHgrav(c)%s(8)*dr(1)*dr(1)*dr(3) + &
                   & BHgrav(c)%s(9)*dr(2)*dr(2)*dr(3) + &
                   & BHgrav(c)%s(10)*dr(1)*dr(2)*dr(3)
           atemp(1) = atemp(1) + 0.5_DP*(3.0_DP*BHgrav(c)%s(1)*dr(1)*dr(1) + &
                & 2.0_DP*BHgrav(c)%s(2)*dr(1)*dr(2)+2.0_DP*BHgrav(c)%s(8)*dr(1)*dr(3) &
                & + BHgrav(c)%s(3)*dr(2)*dr(2) + BHgrav(c)%s(5)*dr(3)*dr(3) &
                & + BHgrav(c)%s(10)*dr(2)*dr(3))*invdr5*invdrsqd - &
                & 3.5_DP*sscalar*dr(1)*invdr5*invdrsqd*invdrsqd
           atemp(2) = atemp(2) + 0.5_DP*(3.0_DP*BHgrav(c)%s(4)*dr(2)*dr(2) + &
                & 2.0_DP*BHgrav(c)%s(3)*dr(1)*dr(2)+2.0_DP*BHgrav(c)%s(9)*dr(2)*dr(3) &
                & + BHgrav(c)%s(2)*dr(1)*dr(1) + BHgrav(c)%s(6)*dr(3)*dr(3) &
                & + BHgrav(c)%s(10)*dr(1)*dr(3))*invdr5*invdrsqd - &
                & 3.5*sscalar*dr(2)*invdr5*invdrsqd*invdrsqd
           atemp(3) = atemp(3) + 0.5*(3.*BHgrav(c)%s(7)*dr(3)*dr(3) + &
                & 2.*BHgrav(c)%s(5)*dr(1)*dr(3)+2.*BHgrav(c)%s(6)*dr(2)*dr(3) &
                & + BHgrav(c)%s(8)*dr(1)*dr(1) + BHgrav(c)%s(9)*dr(2)*dr(2) &
                & + BHgrav(c)%s(10)*dr(1)*dr(2))*invdr5*invdrsqd - &
                & 3.5*sscalar*dr(3)*invdr5*invdrsqd*invdrsqd
           dpot = dpot + 0.5*sscalar*invdr5*invdrsqd
#endif

           ! Add contribution due to cell c to summation vector
           agravs(1:NDIM) = agravs(1:NDIM) + atemp(1:NDIM)
           adots(1:NDIM) = adots(1:NDIM) + adottemp(1:NDIM)
           pots = pots + dpot

           ! Add Ewald correction force if required
#if defined(EWALD)
           call ewald_force(dr(1:NDIM),BHgrav(c)%m,eaccel(1:NDIM))
           agravs(1:NDIM) = agravs(1:NDIM) + real(eaccel(1:NDIM),DP)
#endif

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
!                 if (p == pp) cycle
                 nlist = nlist + 1
                 treelist(nlist) = pp
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
     if ((nlist > GLISTSIZE - LEAFMAX .or. (nlist > 0 .and. c > ctot_grav))) then
        do i=1,nlist
           pp = treelist(i)
#if defined(NBODY_HERMITE4) && defined(MEANH_GRAVITY)
           call gravity_hermite4_meanh(0.5_DP*(hs + real(sph(pp)%h,DP)),&
                &real(sph(pp)%m,DP),rs(1:NDIM),&
                &real(sph(pp)%r(1:NDIM),DP),vs(1:NDIM),&
                &real(sph(pp)%v(1:NDIM),DP),atemp(1:NDIM),&
                &adottemp(1:NDIM),dpot)
#elif defined(NBODY_HERMITE4)
           call gravity_hermite4(1.0_DP/hs,real(sph(pp)%h,DP),&
                &real(sph(pp)%m,DP),rs(1:NDIM),real(sph(pp)%r(1:NDIM),DP),&
                &vs(1:NDIM),real(sph(pp)%v(1:NDIM),DP),&
                &atemp(1:NDIM),adottemp(1:NDIM),dpotp)
#endif
           agravs(1:NDIM) = agravs(1:NDIM) + atemp(1:NDIM)
           adots(1:NDIM) = adots(1:NDIM) + adottemp(1:NDIM)
           pots = pots + dpot
        end do
        nlist = 0
     end if
     ! -----------------------------------------------------------------------

     ! Exit loop if we have finished traversing the tree
     if (c > ctot_grav) exit

  end do
! ============================================================================


  return
END SUBROUTINE BHgrav_accel_jerk
