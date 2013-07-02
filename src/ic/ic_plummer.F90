! IC_PLUMMER.F90
! D. A. Hubber & R. Smith - 11/06/2010
! Creates a plummer sphere.
! Reads in a parameters file 'plummer.dat' with the following data : 
! out_file      : Output filename
! out_file_form : Format of output file
! runit         : Length unit
! munit         : Mass unit
! rhounit       : Density unit
! vunit         : Velocity unit
! rplummer      : 'Plummer' radius
! rmax          : Maximum radius of particles from centre of core
! mplummer      : Total mass in plummer sphere
! ptot          : Total number of particles inside plummer sphere
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_plummer
  use interface_module, only : insertion_sort_real,paramerror,write_data
  use particle_module
  use filename_module
  use scaling_module
  use time_module
  use hydro_module
  use type_module
  use neighbour_module
  use constant_module
  use sink_module
  implicit none
  
  character(len=256) :: out_file            ! Name of output file
  logical :: istar                          ! Are particles 'stars'?
  integer :: i                              ! Integration counter
  integer :: ntot                           ! ..
  integer :: p                              ! Particle counter
  integer :: s                              ! ..
  real(kind=PR) :: dr(1:NDIM)               ! Relative position
  real(kind=PR) :: drmag                    ! Distance
  real(kind=PR) :: drsqd                    ! Distance squared
  real(kind=PR) :: mplummer                 ! Mass of plummer sphere
  real(kind=PR) :: mp                       ! Mass of particle p
  real(kind=PR) :: radmax                   ! Radius ...
  real(kind=PR) :: rcentre(1:NDIM)          ! ..
  real(kind=PR) :: rplummer                 ! ..
  real(kind=PR) :: rstar                    ! ..
  real(kind=PR) :: cdmfrac,gasfrac,starfrac ! ..
  real(kind=PR) :: raux                     ! ..

  integer, allocatable :: porder(:)         ! ..
  real(kind=PR), allocatable :: radsqd(:)   ! ..

  integer :: idum
  real(kind=PR) :: rpl,rlim,tcr,g,gp,rpold
  real(kind=PR) :: mrpl,mrlim
  real(kind=PR) :: x1,x2,x3,x4,x5,x6,x7
  real(kind=PR) :: r,vm,ve,t1,t2,w,z

  integer :: mcount
  mcount = 0
  
  debug1("Creating plummer sphere [ic_plummer.F90]")

! Read parameters file
  call default_parameters

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types
  
! Read polytrope set-up file  
  open(unit=1,file='plummer.dat',status='old')
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) out_file     
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) runit
  read(1,*) munit
  read(1,*) rhounit
  read(1,*) vunit
  read(1,*) rplummer
  read(1,*) radmax
  read(1,*) mplummer
  read(1,*) rstar
  read(1,*) pcdm
  read(1,*) pgas
  read(1,*) stot
  read(1,*) cdmfrac
  read(1,*) gasfrac
  read(1,*) starfrac
  read(1,*) gamma
  read(1,*) idum
  close(1)

  if (pcdm + pgas + stot <= 0) call paramerror("ntot <= 0")
  if (rplummer <= 0.0_PR) call paramerror("rplummer <= 0")
  if (rplummer > radmax) call paramerror("rplummer > radmax")
#if !defined(SINKS)
  if (stot > 0) call paramerror("Sink particles not activated")
#endif
#if !defined(INTERNAL_ENERGY)
  if (pgas > 0) call paramerror("Energy equation not activiated")
#endif

  in_file       = trim(adjustl(in_file))
  in_file_form  = trim(adjustl(in_file_form))
  out_file      = trim(adjustl(out_file))
  out_file_form = trim(adjustl(out_file_form))
  restart       = .false.
  run_id        = ''
  rcentre(1:3)  = 0.0_PR
  raux = cdmfrac + gasfrac + starfrac
  cdmfrac = cdmfrac / raux
  gasfrac = gasfrac / raux
  starfrac = starfrac / raux

! Calculate scaling variables
  call units

! Scale input parameters
  mplummer = mplummer / mscale
  rplummer = rplummer / rscale
  radmax   = radmax / rscale
  rstar    = rstar / rscale
  
  pboundary = 0
  picm      = 0
  pdust     = 0
  pion      = 0
  ptot      = pcdm + pgas
  ntot      = pcdm + pgas + stot
  call allocate_memory(.FALSE.)
  call types

  write(6,*) "pcdm : ",pcdm,"    pgas : ",pgas
  write(6,*) "ptot : ",ptot,"    stot : ",stot


! Initialize random number generator :
  x1    = ran3(idum) 
  
  write (6,*) 'Generating Plummer sphere'

! ============================================================================
  do i=1,ntot

!     write(6,*) "Creating particle :",p

10   continue      

     ! Position of particle
     ! -----------------------------------------------------------------------
     x1 = ran3(idum) 
     x2 = ran3(idum)
     x3 = ran3(idum) 

     if (x1 == 0.0_PR .and. x2 == 0.0_PR .and. x3 == 0.0_PR) goto 10

     r = 1.0_PR / sqrt(x1**(-2.0_PR/3.0_PR) - 1.0_PR)
!     write(6,*) x1,x2,x3,r
     
     if (r .gt. radmax/rplummer) goto 10
     
     z = (1.0_PR - 2.0_PR*x2)*r

     if (i > stot .and. i <= pgas + stot) then
        p = i - stot
        sph(p)%r(3) = z
        sph(p)%r(1) = sqrt(r*r - z*z) * cos(TWOPI*x3)
        sph(p)%r(2) = sqrt(r*r - z*z) * sin(TWOPI*x3)
        sph(p)%m = gasfrac / real(pgas,PR)
     else if (i > stot) then
        p = i - stot
        sph(p)%r(3) = z
        sph(p)%r(1) = sqrt(r*r - z*z) * cos(TWOPI*x3)
        sph(p)%r(2) = sqrt(r*r - z*z) * sin(TWOPI*x3)
        sph(p)%m = cdmfrac / real(pcdm,PR)
#if defined(SINKS)
     else
        s = i
        sink(s)%r(3) = z
        sink(s)%r(1) = sqrt(r*r - z*z) * cos(TWOPI*x3)
        sink(s)%r(2) = sqrt(r*r - z*z) * sin(TWOPI*x3)
        sink(s)%m = starfrac / real(stot,PR)
#endif
     end if

     ! Maximum velocity for this distance 
     ve = sqrt(2.0_PR / sqrt(1.0_PR + r*r))
     
20   continue
     

     ! Velocity of particle :
     ! -----------------------------------------------------------------------
     x4 = ran3(idum) 
     x5 = ran3(idum)       
     
     t1 = 0.1_PR*x5
     t2 = x4*x4*(1.0_PR - x4*x4)**(7.0_PR/2.0_PR)
     
     if (t1 .gt. t2) goto 20
     
     vm = ve*x4
     x6 = ran3(idum) 
     x7 = ran3(idum)
     
     w = (1.0_PR - 2.0_PR*x6)*vm

     if (i > stot .and. i <= pgas + stot) then
        p = i - stot
        sph(p)%v(1:3) = 0.0_PR
        sph(p)%sound = sqrt(0.1666666_PR / sqrt(1.0_PR + r*r))
        sph(p)%rho = 1.0_PR
#if defined(INTERNAL_ENERGY)
        sph(p)%u = sph(p)%sound*sph(p)%sound/(gamma - 1.0_PR)
        call thermal(p)
#endif
     else if (i > stot) then
        p = i - stot
        sph(p)%v(3) = w
        sph(p)%v(1) = sqrt(vm*vm - w*w)*cos(TWOPI*x7)
        sph(p)%v(2) = sqrt(vm*vm - w*w)*sin(TWOPI*x7)
#if defined(SINKS)
     else
        s = i
        sink(s)%v(3) = w
        sink(s)%v(1) = sqrt(vm*vm - w*w)*cos(TWOPI*x7)
        sink(s)%v(2) = sqrt(vm*vm - w*w)*sin(TWOPI*x7)
#endif
     end if


  end do
! ============================================================================


! Instanly move to COM
  call diagnostics

  com_frame = .true.
  call COM

  write(6,*) "Generated plummer sphere!!"


! Calculate mass-distribution as a function of radius and write to file
! ----------------------------------------------------------------------------
  allocate(porder(1:ntot))
  allocate(radsqd(1:ntot))

  do p=1,ptot
     porder(p) = p
     radsqd(p) = sqrt(dot_product(sph(p)%r(1:NDIM),sph(p)%r(1:NDIM)))
  end do
#if defined(SINKS)
  do i=1,stot
     p = i + ptot
     porder(p) = p
     radsqd(p) = sqrt(dot_product(sink(i)%r(1:NDIM),sink(i)%r(1:NDIM)))
  end do
#endif

  call insertion_sort_real(ptot,porder(1:ptot),radsqd(1:ptot))

  open(1,file="plummer-radmass.dat",status="unknown",form="formatted")
  do i=1,ptot
     p = porder(i)
     write(1,'(2E18.10)') radsqd(i),real(i,PR)/real(ntot,PR)
  end do
  close(1)

  deallocate(radsqd)
  deallocate(porder)


! Now scale variables to required physical size
  do p=1,ptot
     sph(p)%r(1:NDIM) = sph(p)%r(1:NDIM)*rplummer
     sph(p)%m = sph(p)%m*mplummer
     sph(p)%v(1:VDIM) = sph(p)%v(1:VDIM)*sqrt(mplummer/rplummer)
#if defined(INTERNAL_ENERGY)
     if (p <= phydroend) sph(p)%u = sph(p)%u*(mplummer/rplummer)
     if (p <= phydroend) sph(p)%temp = (gamma - 1.0_PR)*sph(p)%u/Pconst
#endif
  end do
#if defined(SINKS)
  do i=1,stot
     sink(i)%r(1:NDIM) = sink(i)%r(1:NDIM)*rplummer
     sink(i)%m         = sink(i)%m*mplummer
     sink(i)%v(1:VDIM) = sink(i)%v(1:NDIM)*sqrt(mplummer/rplummer)
     sink(i)%radius    = rstar
     sink(i)%h         = INVKERNRANGE*sink(i)%radius
  end do
#endif

! Setting up for different particle types
  call types

! Write particle data to file
  call diagnostics

  restart = .false.
  call initialize_seren_variables_2

  out_file = trim(adjustl(out_file))
  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  call write_data_debug("ICPLUMMER.debug.dat",rcentre(1:NDIM))
#endif

! Clean-up all memory
  call clean_up
  
  stop

contains




! ******************************************************************
      function  ran1(idum)
!     Random number generator (Numerical Recepies 2nd edition)
!     DEFAULT IST INTEGER*8 (=INTEGER)
!     siehe (max-funktion,variable HILF)
! ******************************************************************


      integer*4 idum,ia,im,iq,ir,ntab,ndiv,HILF
      real*4    ran1,am,eps,rnmax

      parameter (ia = 16807 , im = 2147483647 , am = 1./im ,&
     &     iq = 127773 , ir = 2836 , ntab = 32 , &
     &     ndiv = 1+(im-1)/ntab , eps = 1.2e-7 , rnmax = 1.-eps)

      integer*4 j,k,iv(ntab),iy

      SAVE iv,iy
      DATA iv /ntab*0/, iy /0/

      HILF = 1

      if (idum.le.0.or.iy.eq.0) then

!       ...........
!         idum = max(-idum,1) : idum hat type integer*4,und1 integer*8
         idum = max(-idum,HILF)
!       .............

         do 11 j = ntab+8,1,-1
            k    = idum/iq
            idum = ia*(idum-k*iq)-ir*k
            if (idum.lt.0) idum  = idum+im
            if (j.le.ntab) iv(j) = idum
 11      continue

         iy = iv(1)

      endif

      k    = idum/iq
      idum = ia*(idum-k*iq)-ir*k

      if (idum.lt.0) idum = idum+im

      j     = 1+iy/ndiv
      iy    = iv(j)
      iv(j) = idum
      ran1  = min(am*iy,rnmax)

      return


 end function ran1

! ******************************************************************
!     end of file
! ******************************************************************



! ==================================================================
! ==================================================================
!
!     THE FUNCTION SUBROUTINE RAN3.        
        DOUBLE PRECISION FUNCTION RAN3(IDUM)                                   
        IMPLICIT NONE
        INTEGER idum,iff,k,i,ii,inext,inextp
        DOUBLE PRECISION FAC,mbig,mseed,mz,ma,mj,mk
        PARAMETER (MBIG=4000000., MSEED=1618033., MZ=0., FAC=1./MBIG)
        DIMENSION MA(55)
              SAVE
            DATA IFF/0/
          IF(IDUM<0 .OR. IFF==0)THEN
             IFF=1
             MJ=MSEED-IABS(IDUM)
             MJ=MOD(MJ,MBIG)
             MA(55)=MJ
             MK=1
            DO 100 I=1,54
              II=MOD(21*I,55)
              MA(II)=MK
              MK=MJ-MK
              IF(MK.LT.MZ)MK=MK+MBIG
              MJ=MA(II)
100         CONTINUE
            DO 200 K=1,4
              DO 300 I=1,55
                MA(I)=MA(I)-MA(1+MOD(I+30,55))
                IF(MA(I)<MZ)MA(I)=MA(I)+MBIG
300           CONTINUE
200         CONTINUE
             INEXT=0
             INEXTP=31
             IDUM=1
         ENDIF
            INEXT=INEXT+1
            IF(INEXT==56)INEXT=1
            INEXTP=INEXTP+1
            IF(INEXTP==56)INEXTP=1
            MJ=MA(INEXT)-MA(INEXTP)
            IF(MJ<MZ)MJ=MJ+MBIG
            MA(INEXT)=MJ
            RAN3=MJ*FAC
             RETURN
             END FUNCTION ran3




END PROGRAM ic_plummer
