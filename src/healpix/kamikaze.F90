! KAMIKAZE.F90
! J. E. Dale, D. A. Hubber, J. Ngoumou - 12/10/2011
! Code for simulating pressure-driven O-star winds by injecting hot
! gas particles. Before this code is called, the momentum winds
! code must first be called to establish the location of the wind
! 'working face', so that the ISM particles which the wind will
! strike can be identified. The momentum winds code is also used
! to simulate the action of the free wind at the inner wind shock.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE kamikaze(isource)
  use particle_module
  use HP_module
  use type_module
  use hydro_module
  use scaling_module
  use time_module
  implicit none

  integer, intent(in) :: isource          ! HP source id
  
  integer :: i                            ! ..
  integer :: p                            ! ..
  integer :: iwind                        ! loop index
  integer :: iwinner                      ! ISM particle which will be struck
  integer :: ninject                      ! no. wind particles actually put in
  integer :: nwind                        ! no. injected particles 

  real(kind=DP) :: cuprob                 ! cumulative prob of hitting surface
  real(kind=PR) :: dput                   ! distance fr source to put new part.
  real(kind=PR) :: dr(1:NDIM)             ! vect. between source & target part
  real(kind=PR) :: drandv                 ! size of random vector
  real(kind=PR) :: dsource                ! ...
  real(kind=PR) :: frac_four_pi           ! frac of 4pi taken up by ISM
  real(kind=PR) :: mdotwind               ! mass loss rate
  real(kind=PR) :: mwind                  ! mass of each wind particle
  real(kind=DP) :: prob                   ! prob. of striking particle
  real(kind=PR) :: random                 ! random number
  real(kind=PR) :: randtheta              ! random angle from random vector
  real(kind=PR) :: randv(1:NDIM)          ! random unit vector
  real(kind=PR) :: rput(1:NDIM)           ! location to put new part.
  real(kind=PR) :: rwsource(1:NDIM)       ! location of wind source
  real(kind=PR) :: vput(1:NDIM)           ! velocity of new part.
  real(kind=PR) :: vwind                  ! v_infinity
  real(kind=PR) :: thetatarget            ! angle subtended by target particle
  real(kind=PR) :: Twind                  ! shocked wind temp.    
  real(kind=PR) :: uwind                  ! shocked wind internal energy
  real(kind=PR) :: wtstep                 ! wind timestep
  real(kind=PR) :: xrand                  ! ..
  real(kind=PR) :: yrand                  ! ..
  real(kind=PR) :: zrand                  ! ..
  real(kind=PR) :: mwratio                ! Ratio of wind to sph particle mass 
  real(kind=PR) :: mpart 

  type(sph_particle) :: windpart          ! struct for new wind particle
  
  real(kind=PR) :: time0                  ! time for the 1st injection of ptcls
  real(kind=PR) :: mwind0                 ! mass of wind ptcls injected at time0 
  integer :: nwind0                       ! number of particles injected at time0
  real(kind=PR) :: Rshell0                ! radius of shell at time0 
  integer :: numbinj                      ! flag for first injection
  real(kind=PR) :: tol_param              ! Tolerance parameter for injection radius
  real(kind=PR) :: s                      ! r/h
  real(kind=PR) :: w0                     ! kernel

  debug2("[kamikaze.F90]")

! Obtain location of wind source
  rwsource(1:NDIM) = HPsource(isource)%r(1:NDIM)
  nnewwind = 0
  ninject = 0
  
! Get information on target particle
  do p=1,ptot
     if (windsurface(p) .gt. -1)then
        iwinner = p
        exit
     end if
  end do
   

! Obtain mass flux and vinfinity - NB these must agree with mass
! and vinfinity being used by the momentum winds code
  mdotwind = HPsource(isource)%M_loss
  vwind = HPsource(isource)%v_wind
  numbinj = HPsource(isource)%numb_inj
 
  mpart = 0.01_PR
  mwratio = 0.01_PR
  mwind = mwratio*mpart
  wtstep = 0.0_PR
  tol_param = 1.0_PR
  
  nwind0 = 100
  mwind0 = 0.001_PR*mmean

  
! Calculate the time for first particle injection and first Radius of injection
  time0 = nwind0*mwind0/mdotwind    
  Rshell0 = sph(iwinner)%h 
  !Rshell0 = ((125.0_PR*mdotwind*vwind**2.0_PR*time0**(3.0_PR))/&
  ! (308.0_PR*PI*sph(iwinner)%rho))**(0.2_PR)
  
! Set number of particles injected depending on the timestep
  if (HPsource(isource)%tlastwind .lt. 0.0_DP)then
     HPsource(isource)%tlastwind = time
     print*,'NO INJECTION'
     return
  else if (HPsource(isource)%tlastwind .lt. time0)then
     HPsource(isource)%tlastwind = time
     print*,'time less than time0', time0, time
     return
  else 
     if (numbinj .eq. 0)then    
        wtstep = time
        HPsource(isource)%tlastwind = time
        nwind = nwind0
        mwind = mwind0
        dput = 1.0_PR*Rshell0
        print *, "OOOOHH!!!! Prepare for first injection! ",wtstep 
        numbinj = numbinj+1
        HPsource(isource)%numb_inj = numbinj
        print*, 'nwind =', nwind,'injection number: ', numbinj
        synchronise_all = .true.
     else
        wtstep = time - HPsource(isource)%tlastwind 
        mwind = mwind0
        nwind = int(mdotwind*wtstep/(mwind)) 
        nwind=0
        if (nwind .lt. 1)then
           print*,'nwind =',nwind,' returning...'
           return
        end if
        HPsource(isource)%tlastwind = time
        nwind = int(mdotwind*wtstep/(mwind))
        mwind = mwind0 
        dput = tol_param*sph(iwinner)%h*nwind**(1.0_PR/3.0_PR)
        print*, 'NOT THE FIRST INJECTION'
        numbinj = numbinj+1
        HPsource(isource)%numb_inj = numbinj
        print*, 'nwind =', nwind,'injection number: ', numbinj
     end if
  end if
  print*, 'RSHELL0', Rshell0*rscale, 'dput',dput*rscale,&
       & 'target', sqrt(dot_product(sph(iwinner)%r(1:NDIM),&
       &sph(iwinner)%r(1:NDIM)))*rscale       
  

! compute temperature of shocked wind (see Lamers and Cassinelli, 
! eqn 12.3)
  Twind = 1.4e5_PR*(vwind/(100.0_PR*1000.0_PR/vscale/v_SI))**2
  write(*,*)'mdotwind, vwind, mwind, Twind:'
  write(*,*) mdotwind,vwind,mwind,Twind
! convert temperature to internal energy if required
  uwind = Pconst*Twind/(gamma - 1.0_PR)
  write(*,*)'wind int. energy: ',uwind

! mwratio = mwind/mpart

! windsurface(1:ptot) is an array valid for each wind source where -1 indicates
! the particle is not on the wind working face, otherwise indicates the 
! healpix level, l, on which that particle is, which determines the solid 
! angle subtended at the source


! debugging - print out locations of surface particles
!  open(30,file='surface.out')
!  do p=1,ptot
!     if(windsurface(p).gt.-1)then
!        dr(1:NDIM)=sph(p)%r(1:NDIM)-rwsource(1:NDIM)
!        dsource=sqrt(dot_product(dr,dr))
!        write(30,*)sph(p)%r(1:NDIM),sph(p)%h,dsource,windsurface(p)
!     endif
!  end do
!  close(30)
!  open(30,file='wind.out')
  
!Get information on target particle
!  do p=1,ptot
!        if (windsurface(p) .gt. -1)then
!        	iwinner = p
!        	exit
!        end if
!   end do


!----------------------------------- INJECTION --------------------------------
  do iwind=1,nwind
     
     positions: do
        ! generate random coordinates
        ! ---------------------------------------------------------------------
        call random_number(xrand)
        call random_number(yrand)
        call random_number(zrand)
        
        randv(1) = 2.0_PR*xrand - 1.0_PR
        randv(2) = 2.0_PR*yrand - 1.0_PR
        randv(3) = 2.0_PR*zrand - 1.0_PR        
        drandv = dot_product(randv(1:3),randv(1:3))
        
        ! throw away random vectors outside the unit sphere   
        if (drandv .gt. 1.0_PR) cycle positions
        
        drandv = sqrt(drandv)
        !randv(1:NDIM) = randv(1:NDIM)/drandv
        exit positions
        
     end do positions
     
     ! scale coordinates
     ! ------------------------------------------------------------------------
     
     rput(1:NDIM) = rwsource(1:NDIM) + dput*(randv(1:NDIM)) 
     vput(1:NDIM) = 0.0_PR    
     write(30,*) rput(1:NDIM)
     
     ! give wind particle position and other attributes
     ! give new wind particle the same h as the target
     windpart%h = 0.5_PR*sph(iwinner)%h
     windpart%r(1:NDIM) = rput(1:NDIM)
     windpart%v(1:NDIM) = vput(1:NDIM)
     windpart%m = mwind
     windpart%temp = Twind
     
     s= 2.0_PR*sqrt(dot_product(randv(1:NDIM),randv(1:NDIM)))
     if (s < 1.0_PR) then
        w0 = (1.0_PR - 1.5_PR*s*s + 0.75_PR*s**3)
     else if (s < 2.0_PR) then
        w0 = 0.25_PR*(2.0_PR - s)**3
     else
        w0 = 0.0_PR
     end if
     windpart%u = uwind*w0
     windpart%rho = (1.2_PR/windpart%h)**3*mwind*(nwind/60)**(1.0_PR/3.0_PR)
     if (windpart%u .lt. sph(iwinner)%u*sph(iwinner)%rho/windpart%rho)then
        windpart%u = sph(iwinner)%u*sph(iwinner)%rho/windpart%rho  
     end if
     
     ninject = ninject + 1
     
     !write(6,*) "NEW WIND PARTICLE : ",iwind,windpart%r(1:NDIM)*rscale,&
     !     &windpart%v(1:NDIM)*vscale,windpart%m*mscale,windpart%h*rscale,&
     !     &windpart%temp,windpart%u*uscale
     !write(6,*) "MASSES : ",windpart%m,sph(1)%m,windpart%m/sph(1)%m
     
     ! insert new wind particle into main arrays
     call create_new_sph_particle(gasid,ptot + iwind,windpart)
     
  end do

  !call write_injection_wind_data(isource)
  
  write(*,*) 'injected ',ninject,' of ',nwind,&
       &' wind particles.',' injection nr: ',numbinj
  nnewwind = nnewwind + nwind

  return
END SUBROUTINE kamikaze
