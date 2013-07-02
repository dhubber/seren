! LANE_EMDEN.F90
! Solves Lane-Emden equation (or isothermal equation) to obtain the
! Lane-Emden functions for various Polytropes.  The results are tabulated 
! and outputted to file
! D. A. Hubber - 30/09/2008
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM lane_emden
   use definitions
   implicit none

   character(len=256) :: out_file    ! Name of output file

   logical :: isocloud               ! Flag if isothermal cloud or not
   integer :: i                      ! integration counter
   integer :: ibound                 ! 
   integer :: nmax                   ! no. of integration points
   real(kind=PR) :: a0sqd            ! Isothermal sound speed squared
   real(kind=PR) :: delta_xi         ! Step-size of xi
   real(kind=PR) :: etapoly
   real(kind=PR) :: Kpoly            ! Polytropic constant
   real(kind=PR) :: mcloud           ! Mass of cloud
   real(kind=PR) :: mu               ! Dimensionless mass variable
   real(kind=PR) :: npoly            ! Polytropic index
   real(kind=PR) :: phi              ! d(psi)/d(xi) or d(theta)/d(xi)
   real(kind=PR) :: psi              ! (log_e) dimensionless density
   real(kind=PR) :: rcloud           ! Physical radius of cloud
   real(kind=PR) :: remainder        ! aux. interpolation varibale
   real(kind=PR) :: R0               ! 
   real(kind=PR) :: rho0             ! central density
   real(kind=PR) :: theta            ! LE density variable
   real(kind=PR) :: xi               ! Dimensionless distance
   real(kind=PR) :: xi_bound         ! Dimensionless radius/boundary
   real(kind=PR) :: k1_phi, k1_psi, k1_theta    ! 4th order Runge-Kutta values
   real(kind=PR) :: k2_phi, k2_psi, k2_theta    ! for phi, psi and theta
   real(kind=PR) :: k3_phi, k3_psi, k3_theta    ! 
   real(kind=PR) :: k4_phi, k4_psi, k4_theta    !
   real(kind=PR), allocatable :: xi_array(:)    ! 
   real(kind=PR), allocatable :: psi_array(:)   ! 
   real(kind=PR), allocatable :: phi_array(:)   !
   real(kind=PR), allocatable :: mu_array(:)    ! 
   real(kind=PR), allocatable :: pressure(:)    ! 
   real(kind=PR), allocatable :: density(:)     ! 
   real(kind=PR), allocatable :: theta_array(:) ! 
   real(kind=PR), allocatable :: mass_array(:)  ! 


! READ initial conditions set-up file  
! ----------------------------------------------------------------------------
   open(unit=1,file='polytrope.dat',status='old')
   read(1,*) isocloud         ! 
   read(1,*) npoly            ! 
   read(1,*) xi_bound         ! dimensionless cloud boundary
   read(1,*) delta_xi         ! integration step
   read(1,*) rho0             !
   read(1,*) Kpoly            ! 
   read(1,'(a)') out_file     ! output filename
   close(1)
   out_file = trim(adjustl(out_file))

! Add 2 for round-off error and end points
   nmax = int(xi_bound/delta_xi) + 2
   write (6,*) "nmax - ",nmax

! Allocate arrays now number of array points is known
   allocate(xi_array(1:nmax))
   allocate(psi_array(1:nmax))
   allocate(phi_array(1:nmax))
   allocate(mu_array(1:nmax))
   allocate(pressure(1:nmax))
   allocate(density(1:nmax))
   allocate(theta_array(1:nmax))
   allocate(mass_array(1:nmax))
   ibound = nmax

! Isothermal sphere solution
! ============================================================================
   if (isocloud) then

      write(6,*) "Solving for isothermal sphere"

! Boundary Conditions (at centre) of Isothermal Sphere -----------------------
      xi  = 0.0
      phi = 0.0
      psi = 0.0

! Tabulate central values
      xi_array(1)  = xi
      psi_array(1) = psi
      phi_array(1) = phi
      mu_array(1)  = phi*(xi**2)
      density(1)   = exp(-psi_array(1))

! Use first few terms of series solution for first step (due to singularity 
! in differential equation at xi = 0)
      xi = delta_xi
      psi = (1/6.0)*(xi**2) - (1/120.0)*(xi**4) + (1/1890.0)*(xi**6)
      phi = (1/3.0)*(xi) -  (1/30.0)*(xi**3) + (1/315.0)*(xi**5)

      xi_array(2)  = xi
      psi_array(2) = psi
      phi_array(2) = phi
      mu_array(2)  = phi*(xi**2)
      density(2)   = exp(-psi_array(2))


! Now loop over all over integration points
! ----------------------------------------------------------------------------
      do i = 3,nmax

       ! Solve using 4th order Runge-Kutta method
         k1_phi = delta_xi*(exp(-psi) - 2*phi/xi)
         k1_psi = delta_xi*phi
         
         k2_phi = delta_xi*(exp(-psi - 0.5*k1_psi) - &
              &2.*(phi + 0.5*k1_phi)/(xi + delta_xi/2))
         k2_psi = delta_xi*(phi + k1_phi/2)
         
         k3_phi = delta_xi*(exp(-psi - 0.5*k2_psi) - &
              &2.*(phi + 0.5*k2_phi)/(xi + delta_xi/2))
         k3_psi = delta_xi*(phi + k2_phi/2)
         
         k4_phi = delta_xi*(exp(-psi - k3_psi) - &
              &2*(phi + k3_phi)/(xi + delta_xi))
         k4_psi = delta_xi*(phi + k3_phi)
         
         phi = phi + ONESIXTH*(k1_phi + k4_phi) + ONETHIRD*(k2_phi + k3_phi)
         psi = psi + ONESIXTH*(k1_psi + k4_psi) + ONETHIRD*(k2_psi + k3_psi)

         xi = real(i - 1,PR)*delta_xi

         ! Tabulate values
         xi_array(i)  = xi
         psi_array(i) = psi
         phi_array(i) = phi
         mu_array(i)  = phi*(xi**2)
         density(i)   = exp(-psi_array(i))
         if (density(i) < 0. .AND. i > ibound) ibound = i
         
      end do
! ----------------------------------------------------------------------------

! Find boundary of sphere in tables
      i = int(xi_bound/delta_xi) + 1
      remainder = (xi_bound - xi_array(i))/delta_xi
 
      xi = xi_array(i) + (xi_array(i+1) - xi_array(i))*remainder
      mu = mu_array(i) + (mu_array(i+1) - mu_array(i))*remainder
      
!      rcloud = (xi_bound*mcloud)/(Kpoly*mu)
!      rho0   = (Kpoly**3)*(mu**2)/(4.*PI*mcloud**2)

      density(1:nmax) = rho0*(theta_array(1:nmax)**(npoly))
      mass_array(1:nmax) = mcloud*mu_array(1:nmax)/mu
      pressure(1:nmax) = a0sqd*density(1:nmax)

      write(6,*) "Mass of cloud   :",mcloud
      write(6,*) "Radius of cloud :",rcloud
      write(6,*) "Central density :",rho0
   

! General polytrope solution
! ============================================================================
   else

      write(6,*) "Solving Lane-Emden equation for polytrope"
      etapoly = 1 + 1/npoly

! Boundary Conditions (at centre) of Isothermal Sphere -----------------------
      xi    = 0.0
      phi   = 0.0
      theta = 1.0

! Tabulate central values
      xi_array(1)    = xi
      theta_array(1) = theta
      phi_array(1)   = phi
      mu_array(1)    = -phi*(xi**2)
      density(1)     = theta**(npoly)

! Use first few terms of series solution for first step (due to singularity 
! in differential equation at xi = 0)
      xi = delta_xi
      theta = 1. - (1./6.)*(xi**2) + (npoly/120.)*(xi**4)
      phi = (-1./3.)*(xi) + (npoly/30.)*(xi**3) 

      xi_array(2)    = xi
      theta_array(2) = theta
      phi_array(2)   = phi
      mu_array(2)    = -phi*(xi**2)
      density(2)     = theta**(npoly)


! Now loop over all over integration points
! ----------------------------------------------------------------------------
      do i = 3,nmax
         ibound = i - 1

       ! Solve using 4th order Runge-Kutta method
         k1_phi = delta_xi*(-theta**(npoly) - 2*phi/xi)
         k1_theta = delta_xi*phi
         
         k2_phi = delta_xi*(-(theta + 0.5*k1_theta)**(npoly) - &
              &2.*(phi + 0.5*k1_phi)/(xi + delta_xi/2))
         k2_theta = delta_xi*(phi + k1_phi/2)
         
         k3_phi = delta_xi*(-(theta + 0.5*k2_theta)**(npoly) - &
              &2.*(phi + 0.5*k2_phi)/(xi + delta_xi/2))
         k3_theta = delta_xi*(phi + k2_phi/2)
         
         k4_phi = delta_xi*(-(theta + k3_theta)**(npoly) - &
              &2*(phi + k3_phi)/(xi + delta_xi))
         k4_theta = delta_xi*(phi + k3_phi)
         
         phi = phi + ONESIXTH*(k1_phi + k4_phi) + ONETHIRD*(k2_phi + k3_phi)
         theta = theta + ONESIXTH*(k1_theta + k4_theta) + &
              &ONETHIRD*(k2_theta + k3_theta)

         xi = real(i - 1,PR)*delta_xi

         ! Tabulate values
         xi_array(i)    = xi
         theta_array(i) = theta
         phi_array(i)   = phi
         mu_array(i)    = -phi*(xi**2)
         if (.NOT.(theta > 0. .AND. theta < BIG_NUMBER)) exit
         density(i)     = theta**(npoly)
         
      end do
! ----------------------------------------------------------------------------

      write (6,*) "Finished integrating loop "

! Find where sphere ends i.e. xi_bound as confirmation
!      i = int(xi_bound/delta_xi) + 1
!      remainder = (xi_bound - xi_array(i))/delta_xi
!      remainder = theta_array(ibound-1) / &
!           &(theta_array(ibound - 1) - theta_array(ibound))
!      write(6,*) "remainder :",remainder,theta_array(ibound-1),theta_array(ibound)
 
!      xi = xi_array(i) + (xi_array(ibound) - xi_array(ibound-1))*remainder
!      mu = mu_array(i) + (mu_array(ibound) - mu_array(ibound-1))*remainder
      xi = xi_array(ibound)
      mu = mu_array(ibound)
      
!      rcloud = (xi_bound*mcloud)/(Kpoly*mu)
!      rho0   = (Kpoly**3)*(mu**2)/(4.*PI*mcloud**2)
      R0 = sqrt((Kpoly*(npoly + 1)*rho0**(1./npoly - 1.))/(4.*PI))
      rcloud = xi*R0
      mcloud = 4.*PI*mu*R0**3

      density(1:nmax) = rho0*(theta_array(1:nmax)**(npoly))
      mass_array(1:nmax) = mcloud*mu_array(1:nmax)/mu
      pressure(1:nmax) = Kpoly*(density(1:nmax))**etapoly

      write(6,*) "R0              :",R0
      write(6,*) "Radius of cloud :",rcloud
      write(6,*) "Mass of cloud   :",mcloud
      write(6,*) "Central density :",rho0   

   end if
! ============================================================================


! Write tables to file
   out_file = trim(adjustl(out_file))
   open(unit=1,file=out_file,status='unknown')
   do i=1,ibound
      write (1,'(9E15.7)') xi_array(i),psi_array(i),phi_array(i),&
           &mu_array(i),theta_array(i),xi_array(i)*R0,&
           &density(i),pressure(i)
   end do
   close(1)

#define RSPH_OUTPUT
#ifdef RSPH_OUTPUT
   out_file = 'startupRSPH3D.dat'
   open(unit=1,file=out_file,status='unknown')
   write(1,'(I6)') ibound
   do i=1,ibound
      write (1,'(4E15.7)') xi_array(i)*R0, density(i), &
           & pressure(i), 0.
   end do
   if (ibound < nmax) then
      do i=ibound+1,nmax
         write (1,'(4E15.7)') real(i - 1,PR)*delta_xi*R0,0.,0.,0.
      end do
   end if
   close(1)
#endif

! Deallocate arrays
   deallocate(mass_array)
   deallocate(theta_array)
   deallocate(density)
   deallocate(pressure)
   deallocate(mu_array)
   deallocate(phi_array)
   deallocate(psi_array)
   deallocate(xi_array)

   stop
 END PROGRAM lane_emden


