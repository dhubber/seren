! BHGRAV_STOCK.F90
! D. A. Hubber - 22/01/2008
! Stocks cells in BH tree with various important properties, i.e. 
! Total mass, position of centre of mass, quadrupole moment terms, ...
! Calculates the opening criteria for various different MACs.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHgrav_stock
  use interface_module, only : eigenvalue_mac
  use particle_module
  use tree_module
  implicit none

  integer :: c                   ! Cell counter
  integer :: cc                  ! Child cell counter
  integer :: cfirst              ! First child
  integer :: cend                ! Final child marker
  integer :: i                   ! Auxilary counter
  integer :: k                   ! Dimension counter
  integer :: l                   ! Level counter
  integer :: nalive              ! No. of 'alive' child cells
  integer :: npart               ! Number of particles in leaf cell c
  integer :: p                   ! Particle id
  real(kind=PR) :: daux          ! Aux. variable for calculating min dist
  real(kind=PR) :: diagonalsqd   ! Diagonal distance squared
  real(kind=PR) :: dr(1:NDIM)    ! Relative displacement vector
  real(kind=PR) :: drsqd         ! Distance squared
  real(kind=PR) :: mc            ! Mass in cell c
  real(kind=PR) :: mcc           ! Mass in child cell cc
  real(kind=PR) :: mp            ! Mass of particle p
  real(kind=PR) :: rc(1:NDIM)    ! Position of centre of mass of cell c
  real(kind=PR) :: rsize         ! 'Size' of cell
  real(kind=PR) :: rp(1:NDIM)    ! Position of particle p
  real(kind=PR) :: qc(1:NQUAD)   ! Quadrupole moment terms of cell c
  real(kind=PR) :: sc(1:NOCT)    ! Octupole moment terms of cell c
#if defined(CELL_VELOCITIES)
  real(kind=PR) :: vc(1:NDIM)    ! Cell COM velocity
#endif

  debug2("Stocking cells with important quantities [BHgrav_stock.F90]")
  debug_timing("BH_TREE")

! Initialise all tree-stock variables
  BHgrav(0:ctot_grav)%dminsqd = 0.0_PR  
  do c=0,ctot_grav
     BHstock(c)%hmax = 0.0_PR
     BHstock(c)%bbmin(1:NDIM) = BIG_NUMBER
     BHstock(c)%bbmax(1:NDIM) = -BIG_NUMBER
  end do

  
! Loop over all levels in gravity tree
! ============================================================================
  do l=ltot_grav,0,-1

#if defined(DEBUG_BHTREESTOCK)
     write(6,*) "Stocking cells on level ",l
     write(6,*) "First cell :",first_cell_grav(l),&
          &"  Last cell :",last_cell_grav(l)
#endif

     ! Loop over all cells on current level
     ! -----------------------------------------------------------------------
     !$OMP PARALLEL DO IF (l > 1) DEFAULT(SHARED) &
     !$OMP PRIVATE(cc,cend,cfirst,daux,diagonalsqd,dr,drsqd,i,k,mc) &
     !$OMP PRIVATE(mcc,mp,nalive,npart,p,qc,rc,rp,rsize,sc)
     do c=first_cell_grav(l),last_cell_grav(l)

#if defined(DEBUG_BHTREESTOCK)
        write(6,*) "Stocking cell ",c," with leaf :",BHgrav(c)%leaf
#endif

        ! Zero variables for current cell
        mc = 0.0_PR
        nalive = 0
        rc(1:NDIM) = 0.0_PR
        rsize = 0.0_PR
#if defined(CELL_VELOCITIES)
        vc(1:NDIM) = 0.0_PR
#endif

        ! If a dead cell, set all quantities to skip cells
        ! --------------------------------------------------------------------
        if (BHgrav(c)%leaf == -1) then
           BHgrav(c)%r(1:NDIM) = 0.0_PR
           BHgrav(c)%m = 0.0_PR
           BHgrav(c)%dminsqd = 0.0_PR
#if defined(QUADRUPOLE)
           BHgrav(c)%q(1:NQUAD) = 0.0_PR
#endif
#if defined(OCTUPOLE)
           BHgrav(c)%s(1:NOCT) = 0.0_PR
#endif
#if !defined(GEOMETRIC_MAC)
           BHgrav(c)%mac = 0.0_PR
#endif
#if defined(CELL_VELOCITIES)
           BHgrav(c)%v(1:NDIM) = 0.0_PR
#endif
           BHstock(c)%hmax = 0.0_PR
           BHstock(c)%bbmin(1:NDIM) = -0.0_PR
           BHstock(c)%bbmax(1:NDIM) = 0.0_PR
           cycle

        ! If a leaf cell with only 1 particle, calculate non-zero quantities
        ! --------------------------------------------------------------------
        else if (BHgrav(c)%leaf == 1) then
           p = BHgrav(c)%plist(1)
           rc(1:NDIM) = sph(p)%r(1:NDIM)
           BHgrav(c)%r(1:NDIM) = rc(1:NDIM)
           BHgrav(c)%m = sph(p)%m
#if defined(QUADRUPOLE)
           BHgrav(c)%q(1:NQUAD) = 0.0_PR
#endif
#if defined(OCTUPOLE)
           BHgrav(c)%s(1:NOCT) = 0.0_PR
#endif
#if defined(CELL_VELOCITIES)
           BHgrav(c)%v(1:NDIM) = sph(p)%v(1:NDIM)
#endif
           BHstock(c)%hmax = sph(p)%h
           BHstock(c)%bbmin(1:NDIM) = rc(1:NDIM)
           BHstock(c)%bbmax(1:NDIM) = rc(1:NDIM)
           

        ! If a leaf cell, add contributions due to all particles
        ! --------------------------------------------------------------------
        else if (BHgrav(c)%leaf > 0) then
           npart = BHgrav(c)%leaf

           ! Loop over all particles in leaf cell and calculate the total 
           ! mass, centre of mass and the bounding box.
           do i=1,npart
              p = BHgrav(c)%plist(i)
              rp(1:NDIM) = sph(p)%r(1:NDIM)
              mp = sph(p)%m
              mc = mc + mp
              BHstock(c)%hmax = max(BHstock(c)%hmax,sph(p)%h)
              rc(1:NDIM) = rc(1:NDIM) + mp*rp(1:NDIM)
#if defined(CELL_VELOCITIES)
              vc(1:NDIM) = vc(1:NDIM) + mp*sph(p)%v(1:NDIM)
#endif
              do k=1,NDIM
                 BHstock(c)%bbmax(k) = max(BHstock(c)%bbmax(k),rp(k))
                 BHstock(c)%bbmin(k) = min(BHstock(c)%bbmin(k),rp(k))
              end do
           end do
           rc(1:NDIM) = rc(1:NDIM) / mc
#if defined(CELL_VELOCITIES)
           vc(1:NDIM) = vc(1:NDIM) / mc
#endif

           ! Compute quadrupole moment tensor for cell
#if defined(OCTUPOLE)
           sc(1:NOCT) = 0.0_PR
#endif
#if defined(QUADRUPOLE)
           qc(1:NQUAD) = 0.0_PR      
           do i=1,npart
              p = BHgrav(c)%plist(i)
              dr(1:NDIM) = sph(p)%r(1:NDIM) - rc(1:NDIM)
              mp = sph(p)%m              
#if NDIM==2
              drsqd = dr(1)*dr(1) + dr(2)*dr(2)
              qc(1) = qc(1) + mp*(3.0_PR*dr(1)*dr(1) - drsqd)
              qc(2) = qc(2) + mp*(3.0_PR*dr(1)*dr(2))
              qc(3) = qc(3) + mp*(3.0_PR*dr(2)*dr(2) - drsqd)
#elif NDIM==3  
              drsqd = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
              qc(1) = qc(1) + mp*(3.0_PR*dr(1)*dr(1) - drsqd)
              qc(2) = qc(2) + mp*(3.0_PR*dr(1)*dr(2))
              qc(3) = qc(3) + mp*(3.0_PR*dr(2)*dr(2) - drsqd)
              qc(4) = qc(4) + mp*(3.0_PR*dr(3)*dr(1))
              qc(5) = qc(5) + mp*(3.0_PR*dr(3)*dr(2))
#endif
#if defined(OCTUPOLE)
              sc(1) = sc(1) + mp*dr(1)*(5.0_PR*dr(1)*dr(1) - 3.0_PR*drsqd)
              sc(2) = sc(2) + mp*dr(2)*(15.0_PR*dr(1)*dr(1) - 3.0_PR*drsqd)
              sc(3) = sc(3) + mp*dr(1)*(15.0_PR*dr(2)*dr(2) - 3.0_PR*drsqd)
              sc(4) = sc(4) + mp*dr(2)*(5.0_PR*dr(2)*dr(2) - 3.0_PR*drsqd)
              sc(5) = sc(5) + mp*dr(1)*(15.0_PR*dr(3)*dr(3) - 3.0_PR*drsqd)
              sc(6) = sc(6) + mp*dr(2)*(15.0_PR*dr(3)*dr(3) - 3.0_PR*drsqd)
              sc(7) = sc(7) + mp*dr(3)*(5.0_PR*dr(3)*dr(3) - 3.0_PR*drsqd)
              sc(8) = sc(8) + mp*dr(3)*(15.0_PR*dr(1)*dr(1) - 3.0_PR*drsqd)
              sc(9) = sc(9) + mp*dr(3)*(15.0_PR*dr(2)*dr(2) - 3.0_PR*drsqd)
              sc(10) = sc(10) + 15.0_PR*mp*dr(1)*dr(2)*dr(3)
#endif
           end do
#endif

           ! Record mass, COM and quadrupole moments in main tree arrays
           BHgrav(c)%m = mc
           BHgrav(c)%r(1:NDIM) = rc(1:NDIM)
#if defined(QUADRUPOLE)
           BHgrav(c)%q(1:NQUAD) = qc(1:NQUAD)
#endif
#if defined(OCTUPOLE)
           BHgrav(c)%s(1:NOCT) = sc(1:NOCT)
#endif
#if defined(CELL_VELOCITIES)
           BHgrav(c)%v(1:NDIM) = vc(1:NDIM)
#endif


        ! If it's a branch cell, add contributions due to child cells
        ! -------------------------------------------------------------------- 
        else if (BHgrav(c)%leaf == 0) then
	   
           ! Set limits for walking through child cells
           cfirst = BHgrav(c)%ifopen
           cend = BHgrav(c)%nextcell

           cc = cfirst
           do
              if (BHgrav(cc)%leaf >= 0) then
                 nalive = nalive + 1
                 mc = mc + BHgrav(cc)%m
                 rc(1:NDIM) = rc(1:NDIM) + BHgrav(cc)%m*BHgrav(cc)%r(1:NDIM)
#if defined(CELL_VELOCITIES)
                 vc(1:NDIM) = vc(1:NDIM) + BHgrav(cc)%m*BHgrav(cc)%v(1:NDIM)
#endif
                 BHstock(c)%hmax = max(BHstock(cc)%hmax,BHstock(c)%hmax)
                 do k=1,NDIM
                    BHstock(c)%bbmax(k) = &
                         &max(BHstock(cc)%bbmax(k),BHstock(c)%bbmax(k))
                    BHstock(c)%bbmin(k) = &
                         &min(BHstock(cc)%bbmin(k),BHstock(c)%bbmin(k))
                 end do
              end if
              cc = BHgrav(cc)%nextcell
              if (cc == cend) exit
           end do

           ! If c happens to contain only dead child cells, mark also as dead
           if (nalive == 0) then
              BHgrav(c)%leaf = -1
              BHgrav(c)%r(1:NDIM) = BIG_NUMBER
              BHgrav(c)%m = 0.0_PR
#if defined(CELL_VELOCITIES)
              BHgrav(c)%v = 0.0_PR
#endif
#if defined(QUADRUPOLE)
              BHgrav(c)%q(1:NQUAD) = 0.0_PR
#endif
#if defined(OCTUPOLE)
              BHgrav(c)%s(1:NOCT) = 0.0_PR
#endif
              BHstock(c)%hmax = BIG_NUMBER
              BHstock(c)%bbmin(1:NDIM) = -BIG_NUMBER
              BHstock(c)%bbmax(1:NDIM) = BIG_NUMBER  
              cycle
           end if

           ! Normalise centre of mass of cell c
           if (mc > 0.0_PR) rc(1:NDIM) = rc(1:NDIM) / mc
#if defined(CELL_VELOCITIES)
           if (mc > 0.0_PR) vc(1:NDIM) = vc(1:NDIM) / mc
#endif

           ! Compute quadrupole moment tensor for cell
#if defined(OCTUPOLE)
           sc(1:NOCT) = 0.0_PR
#endif
#if defined(QUADRUPOLE)
           qc(1:NQUAD) = 0.0_PR
           cc = cfirst
           do
              if (BHgrav(cc)%leaf == -1) then
                 cc = BHgrav(c)%nextcell
                 if (cc == cend) exit
                 cycle
              end if

              mcc = BHgrav(cc)%m 
              qc(1:NQUAD) = qc(1:NQUAD) + BHgrav(cc)%q(1:NQUAD)
              dr(1:NDIM) = BHgrav(cc)%r(1:NDIM) - rc(1:NDIM)
#if NDIM==2
              drsqd = dr(1)*dr(1) + dr(2)*dr(2)
              qc(1) = qc(1) + mcc*(3.0_PR*dr(1)*dr(1) - drsqd)
              qc(2) = qc(2) + mcc*(3.0_PR*dr(1)*dr(2))
              qc(3) = qc(3) + mcc*(3.0_PR*dr(2)*dr(2) - drsqd)
#elif NDIM==3  
              drsqd = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
              qc(1) = qc(1) + mcc*(3.0_PR*dr(1)*dr(1) - drsqd)
              qc(2) = qc(2) + mcc*(3.0_PR*dr(1)*dr(2))
              qc(3) = qc(3) + mcc*(3.0_PR*dr(2)*dr(2) - drsqd)
              qc(4) = qc(4) + mcc*(3.0_PR*dr(3)*dr(1))
              qc(5) = qc(5) + mcc*(3.0_PR*dr(3)*dr(2))
#endif
#if defined(OCTUPOLE)
              sc(1:NOCT) = sc(1:NOCT) + BHgrav(cc)%s(1:NOCT)
              sc(1) = sc(1) + mcc*dr(1)*(5.0_PR*dr(1)*dr(1) - 3.0_PR*drsqd) + &
                   & 1.5_PR*dr(1)*BHgrav(cc)%q(1) - dr(2)*BHgrav(cc)%q(2) - &
                   & dr(3)*BHgrav(cc)%q(4)
              sc(2) = sc(2) + mcc*dr(2)*(15.0_PR*dr(1)*dr(1) - 3.0_PR*drsqd) +&
                   & 4.0_PR*dr(1)*BHgrav(cc)%q(2) + 2.5*dr(2)*BHgrav(cc)%q(1) &
                   & - dr(2)*BHgrav(cc)%q(3) - dr(3)*BHgrav(cc)%q(5)
              sc(3) = sc(3) + mcc*dr(1)*(15.0_PR*dr(2)*dr(2) - 3.0_PR*drsqd) +&
                   & 4.0_PR*dr(2)*BHgrav(cc)%q(2) + 2.5*dr(1)*BHgrav(cc)%q(3) &
                   & - dr(1)*BHgrav(cc)%q(1) - dr(3)*BHgrav(cc)%q(4)
              sc(4) = sc(4) + mcc*dr(2)*(5.0_PR*dr(2)*dr(2) - 3.0_PR*drsqd) + &
                   & 1.5_PR*dr(2)*BHgrav(cc)%q(3) - dr(1)*BHgrav(cc)%q(2) - &
                   & dr(3)*BHgrav(cc)%q(5)
              sc(5) = sc(5) + mcc*dr(1)*(15.0_PR*dr(3)*dr(3) - 3.0_PR*drsqd) +&
                   & 4.0_PR*dr(3)*BHgrav(cc)%q(4) - &
                   & 2.5*dr(1)*(BHgrav(cc)%q(1) + BHgrav(cc)%q(3)) - &
                   & dr(1)*BHgrav(cc)%q(1) - dr(2)*BHgrav(cc)%q(2)
              sc(6) = sc(6) + mcc*dr(2)*(15.0_PR*dr(3)*dr(3) - 3.0_PR*drsqd) +&
                   & 4.0_PR*dr(3)*BHgrav(cc)%q(5) - &
                   & 2.5*dr(2)*(BHgrav(cc)%q(1) + BHgrav(cc)%q(3)) - &
                   & dr(1)*BHgrav(cc)%q(2) - dr(2)*BHgrav(cc)%q(3)
              sc(7) = sc(7) + mcc*dr(3)*(5.0_PR*dr(3)*dr(3) - 3.0_PR*drsqd) - &
                   & 1.5_PR*dr(3)*(BHgrav(cc)%q(1) + BHgrav(cc)%q(3)) - &
                   & dr(1)*BHgrav(cc)%q(4) - dr(2)*BHgrav(cc)%q(5)
              sc(8) = sc(8) + mcc*dr(3)*(15.0_PR*dr(1)*dr(1) - 3.0_PR*drsqd) +&
                   & 4.0_PR*dr(1)*BHgrav(cc)%q(4) + 2.5*dr(3)*BHgrav(cc)%q(1) &
                   & - dr(2)*BHgrav(cc)%q(5) + &
                   & dr(3)*(BHgrav(cc)%q(1) + BHgrav(cc)%q(3))
              sc(9) = sc(9) + mcc*dr(3)*(15.0_PR*dr(2)*dr(2) - 3.0_PR*drsqd) +&
                   & 4.0_PR*dr(2)*BHgrav(cc)%q(5) + 2.5*dr(3)*BHgrav(cc)%q(3) &
                   & - dr(1)*BHgrav(cc)%q(4) + &
                   & dr(3)*(BHgrav(cc)%q(1) + BHgrav(cc)%q(3))
              sc(10) = sc(10) + 15.0_PR*mcc*dr(1)*dr(2)*dr(3) + &
                   & 25.0_PR*(dr(1)*BHgrav(cc)%q(5) + dr(2)*BHgrav(cc)%q(4) + &
                   &      dr(3)*BHgrav(cc)%q(2))
#endif
              cc = BHgrav(cc)%nextcell
              if (cc == cend) exit
           end do
#endif

           ! Record mass, COM and quadrupole moments in main tree arrays
           BHgrav(c)%m = mc
           BHgrav(c)%r(1:NDIM) = rc(1:NDIM)
#if defined(CELL_VELOCITIES)
           BHgrav(c)%v(1:NDIM) = vc(1:NDIM)
#endif
#if defined(QUADRUPOLE)
           BHgrav(c)%q(1:NQUAD) = qc(1:NQUAD)
#endif
#if defined(OCTUPOLE)
           BHgrav(c)%s(1:NOCT) = sc(1:NOCT)
#endif
        end if

        ! Determine MAC quantities for cell
        ! --------------------------------------------------------------------
        do k=1,NDIM
           dr(k) = max(BHstock(c)%bbmax(k) - rc(k),rc(k) - BHstock(c)%bbmin(k))
        end do
        diagonalsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
        rsize = sqrt(diagonalsqd) + KERNRANGE*BHstock(c)%hmax

#if defined(GEOMETRIC_MAC)
        if (thetamaxsqd > 0.0_PR) then
           BHgrav(c)%dminsqd = diagonalsqd/thetamaxsqd
        else 
           BHgrav(c)%dminsqd = BIG_NUMBER 
        end if
#elif defined(GADGET_MAC)
        BHgrav(c)%dminsqd = diagonalsqd/thetamaxsqd
        BHgrav(c)%mac = (mc*diagonalsqd*diagonalsqd/abserror)**(ONETHIRD)
#elif defined(GADGET2_MAC)
        BHgrav(c)%dminsqd = diagonalsqd/thetamaxsqd
        BHgrav(c)%mac = sqrt(mc*diagonalsqd/abserror)
#elif defined(EIGEN_MAC)
        BHgrav(c)%dminsqd = diagonalsqd/thetamaxsqd
        call eigenvalue_mac(real(BHgrav(c)%q(1:NQUAD),DP),BHgrav(c)%mac)
#endif

        ! Ensure all cells including (scatter) neighbours are opened
#ifndef N_BODY
        BHgrav(c)%dminsqd = max(BHgrav(c)%dminsqd,rsize*rsize)
#endif

        ! Ensure cell is opened if too close due to precision problem 
        ! when calculating quadrupole moment terms (i.e. invdr5)
        daux = 10.0_PR*(huge(diagonalsqd))**(-0.2_PR)
        if (BHgrav(c)%dminsqd < daux) BHgrav(c)%dminsqd = daux

     end do
     !$OMP END PARALLEL DO
     ! -----------------------------------------------------------------------

  end do
! ============================================================================

  return 
END SUBROUTINE BHgrav_stock
