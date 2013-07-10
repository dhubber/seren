! EWALD_INIT.F90
! A. McLeod & D. A. Hubber - 21/01/2008
! Calculates look-up tables for Ewald periodic forces
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE ewald_init
  use definitions
  use ewald_module
  use periodic_module
  implicit none

  integer :: grid(1:NDIM)                ! ..
  integer :: h2(1:NDIM)                  ! ..
  integer :: hx                          ! ..
  integer :: hy                          ! ..
  integer :: hz                          ! ..
  integer :: i                           ! ..
  integer :: j                           ! ..
  integer :: k                           ! ..
  integer :: n(1:NDIM)                   ! ..
  integer :: nx                          ! ..
  integer :: ny                          ! ..
  integer :: nz                          ! ..
  integer, parameter :: nrange=4         ! ..
  integer, parameter :: hrange=4         ! ..

  real(kind=PR) :: alphae                ! ..
  real(kind=PR) :: dist                  ! ..
  real(kind=PR) :: ds(1:NDIM)            ! ..
  real(kind=PR) :: dr(1:NDIM)            ! ..
  real(kind=PR) :: dr2(1:NDIM)           ! ..
  real(kind=PR) :: f(1:NDIM)             ! ..
  real(kind=PR) :: k2(1:NDIM)            ! ..
  real(kind=PR) :: kr                    ! ..
  real(kind=PR) :: ksqrd                 ! ..
  real(kind=PR) :: x                     ! ..

  debug2("Initializing Ewald tables [ewald_init.F90]")

#if NDIM==1
  write (6,*) "One dimensional Ewald gravity not supported in this release"
  stop
#endif

! ----------------------------------------------------------------------------
#if NDIM==2
  debug1("Creating Ewald correction force table")

  allocate(fcorr(1:2,1:ewsize(1),1:ewsize(2),1))
  L = (/periodic_size(1),periodic_size(2)/)

! For a grid -L/2 < x,y,z < L/2 - we have quasi-periodic gravity
  ewsizeil = 2.0_PR*real(ewsize - 1,PR) / L    

! crazy scaling factor between real (close) and fourier (far) components
  alphae = 2.0 / minval(L) 

! Grid cell size/spacing
  ds = 1.0 / ewsizeil 

! ============================================================================
   do i=1,ewsize(1)

     ! =======================================================================
     do j=1,ewsize(2)
        
        ! quasi-periodic
        if (i == 1 .AND. j == 1) then 
           ! We are at dr=0, force is zero, skip
           cycle
        end if
        
        grid = (/i,j/)
        
        ! Zero correction term before recalculating for this grid point
        f = 0.0_PR 
        
        ! Grid POINT (top left corner of a grid cell, 
        ! except the last points)
        dr = real(grid - 1,PR)*ds 
        
        ! Compute first sum in ewald summation
        do nx=-nrange,nrange
           do ny=-nrange,nrange
              n(1:2)=(/nx,ny/)
              
              ! Implicit in this is that dr is actually -dr
              dr2 = dr + L*n 
              
              ! Distance to periodic replica
              dist = sqrt(sum(dr2**2))           
              x = erfc(alphae*dist) + (2*alphae*dist * SQRT(INVPI) &
                   & * EXP((-1.0_PR)*(alphae**2)*(dist**2)))
              
              ! First part of correction term
              f = f - dr2*x/(dist**3)       
           end do
        end do
        
        ! Now compute second sum in k-space
        do hx=-hrange,hrange
           do hy=-hrange,hrange
              if (hx == 0 .AND. hy == 0) cycle
              h2(1:2) = (/hx,hy/)
              
              ! Convert to k-space (k=2*PI*h/L)
              k2 = REAL(h2,PR)*2.0*PI/L               
              ksqrd = sum(k2**2)
              
              ! Because dr is actually -dr
              kr = dot_product(k2,-dr)           
              x = 4.0_PR*PI*exp((-1.0_PR)*ksqrd/(4.0_PR*(alphae**2)))&
                   &*sin(kr)/ksqrd
              
              ! Second part of correction term
              f = f - k2*x/product(2*L)      
              
           end do
        end do
        
        ! Subtract force from original particle
        dist = sqrt(sum(dr**2))
        f = f + dr/(dist**3)
           
        ! Record force in table
        fcorr(1:2,i,j,1)=real(f(1:NDIM),PR)
        
     end do
  end do
  
  ! At dr=0, fcorr = 0
  fcorr(1:NDIM,1,1,1) = 0.0_PR 
  
  open(1, file="ewald.dat", status="replace", form="formatted")
  
  do i=1,ewsize(1)
     do j=1,ewsize(2)
        grid = (/i,j/)
        dr = real(grid - 1,PR)*ds
        write (1,*) dr(1:2), fcorr(1:2,i,j,1)
     end do
  end do

  
  close (1)


! ----------------------------------------------------------------------------
#elif NDIM==3

  debug1("Creating Ewald correction force table")

  allocate(fcorr(1:3,1:ewsize(1),1:ewsize(2),1:ewsize(3)))
  L = (/periodic_size(1),periodic_size(2),periodic_size(3)/)

! For a grid -L/2 < x,y,z < L/2 - we have quasi-periodic gravity
  ewsizeil(1:NDIM) = 2.0_PR*real(ewsize - 1,PR) / L(1:NDIM)

! crazy scaling factor between real (close) and fourier (far) components
  alphae = 2.0 / minval(L) 

! Grid cell size/spacing
  ds = 1.0 / ewsizeil 

! ============================================================================
  do i=1,ewsize(1)

     ! =======================================================================
     do j=1,ewsize(2)

        ! Iterate over all grid POINTS (which define the corners of grid CELLS
        ! =====================================================================
        !$OMP PARALLEL DO PRIVATE(grid,x,f,n,nx,ny,nz,hx,hy,hz,h2,dist) &
        !$OMP PRIVATE (dr,dr2,k2,ksqrd,kr)
        do k=1,ewsize(3) 
           
           ! We are at dr=0, force is zero, skip
           if (i == 1 .and. j == 1 .and. k == 1) cycle

           grid = (/i,j,k/)
           
           ! Zero correction term before recalculating for this grid point
           f = 0.0_PR 
           
           ! Grid POINT (top left corner of a grid cell, 
           ! except the last points)
           dr(1:NDIM) = real(grid - 1,PR)*ds(1:NDIM)
           
           ! Compute first sum in ewald summation
           do nx=-nrange,nrange
              do ny=-nrange,nrange
                 do nz=-nrange,nrange
                    n(1:3)=(/nx,ny,nz/)

                    ! Implicit in this is that dr is actually -dr
                    dr2 = dr + L*n 
                    
                    ! Distance to periodic replica
                    dist = sqrt(sum(dr2**2))           
                    x = erfc(alphae*dist) + (2*alphae*dist * SQRT(INVPI) &
                         & * EXP((-1.0_PR)*(alphae**2)*(dist**2)))
                    
                    ! First part of correction term
                    f = f - dr2*x/(dist**3)       
                 end do
              end do
           end do
           
           ! Now compute second sum in k-space
           do hx=-hrange,hrange
              do hy=-hrange,hrange
                 do hz=-hrange,hrange
                    if (hx == 0 .AND. hy == 0 .AND. hz == 0) cycle
                    h2(1:3) = (/hx,hy,hz/)

                    ! Convert to k-space (k=2*PI*h/L)
                    k2 = REAL(h2,PR)*2.0*PI/L               
                    ksqrd = sum(k2**2)
                    
                    ! Because dr is actually -dr
                    kr = dot_product(k2,-dr)           
                    x = 4.0_PR*PI*exp((-1.0_PR)*ksqrd/(4.0_PR*(alphae**2)))&
                         &*sin(kr)/ksqrd
                    
                    ! Second part of correction term
                    f = f - k2*x/product(2*L)      

                 end do
              end do
           end do
           
           ! Subtract force from original particle
           dist = sqrt(sum(dr**2))
           f = f + dr/(dist**3)
           
           ! Record force in table
           fcorr(1:3,i,j,k)=real(f(1:NDIM),PR)

        end do
        !$OMP END PARALLEL DO

     end do
  end do
  
! At dr=0, fcorr = 0
  fcorr(1:NDIM,1,1,1) = 0.0_PR 
  
  open(1, file="ewald.dat", status="replace", form="formatted")
  
  do i=1,ewsize(1)
     do j=1,ewsize(2)
        do k=1,ewsize(3)
           grid = (/i,j,k/)
           dr = real(grid-1,PR)*ds
           write (1,*) dr(1:3), fcorr(1:3,i,j,k)
        end do
     end do
  end do
  
  close (1)

#endif
! ----------------------------------------------------------------------------

     
  return
END SUBROUTINE ewald_init
