! IC_KEPDISC.F90
! D. H. FORGAN - 14/7/2011
! Creates initial conditions for a Keplerian, Self-Gravitating Disc
! Assumes a constant Toomre Q Parameter
! Uses the Numerical Recipes ran2 random number generator
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_kepdisc
  use interface_module, only:comperror,paramerror,&
	& read_data, write_data
  use particle_module
  use hydro_module
  use constant_module
  use filename_module
  use scaling_module
  use sink_module
  use type_module
  use time_module
  use neighbour_module, only : pp_gather
  use Tprof_module, only : temp_inf
  implicit none

  character(len=256) :: out_file    ! Name of output file
  integer :: k                      ! Dimension counter
  integer :: p                      ! Particle counter
  integer :: iseed		    ! Random Number Seed
  integer, parameter :: ienc = 5000 ! Enclosed Mass Table Counter
  integer :: itable
  real(kind=PR) :: rin		    ! Inner Disc Radius
  real(kind=PR) :: rout		    ! Outer Disc Radius
  real(kind=PR) :: Mstar	    ! Sink Mass
  real(kind=PR) :: massratio	    ! Disc-to-Sink Mass Ratio
  real(kind=PR) :: mdisc	    ! Disc Mass
  real(kind=PR) :: phi		    ! cylindrical azimuth
  real(kind=PR) :: Qmin		    ! Minimum Toomre Q Parameter
  real(kind=PR) :: p_index	    ! Surface Density Power Law Index
  real(kind=PR) :: q_index  	    ! Sound Speed Power Law Index
  real(kind=PR) :: sigma_0	    ! Surface Density at r=1
  real(kind=PR) :: cs_0		    ! Sound Speed at r=1
  real(kind=PR) :: omega_rot	    ! Angular Velocity
  real(kind=PR) :: cs		    ! Sound Speed
  real(kind=PR) :: H		    ! Scale Height
  real(kind=PR) :: rtest	    ! Cylindrical Radius
  real(kind=PR) :: ztest	    ! Cartesian Z Coordinate
  real(kind=PR) :: rsphere	    ! Spherical Radius
  real(kind=PR) :: randtest	    ! Random number for testing
  real(kind=PR) :: f		    !	"	"
  real(kind=PR) :: f_max 	    ! 	"	"
  real(kind=PR) :: dr		    !	Annulus width
  real(kind=PR) :: vk 		    ! Keplerian tangential velocity
  real(kind=PR) :: cmx, cmy,cmz     ! Position of system centre of mass
  real(kind=PR) :: vcmx, vcmy, vcmz ! Velocity of system centre of mass
  real(kind=PR), dimension(ienc) :: enc_m	! Enclosed Mass
  real(kind=PR), dimension(ienc) :: R_an	! Annular Radius
  real(kind=PR), dimension(ienc) :: r_c		! Centre of Annulus


  debug1("Create Keplerian self-gravitating disc [ic_kepdisc.F90]")

! Set default parameters
  call default_parameters
	
! Set up other initial variables
  call initialize_seren_variables_1

!  Read in kepdiscparams.dat
  write(6,*) "------------------------------"
  write(6,*) "         ic_kepdisc           "
  write(6,*) "------------------------------"
  open(unit=1,file="kepdiscparams.dat",status="old")
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)    
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) pgas
  read(1,*) iseed
  read(1,*) rin
  read(1,*) rout
  read(1,*) Mstar
  read(1,*) massratio
  read(1,*) Qmin
  read(1,*) p_index
  read(1,*) q_index
  read(1,*) gamma  
  read(1,*) munit
  read(1,*) runit
  read(1,*) vunit
  read(1,*) uunit
read(1,*) dudtunit
  close(1)

#if !defined(SINKS)
  call comperror("SINKS not switched on")
#endif


  out_file = trim(adjustl(out_file))
  out_file_form = trim(adjustl(out_file_form))

  restart = .false.
  run_id = ''

!  Scale input variables

  rin = rin/rscale
  rout = rout/rscale
  Mstar = Mstar/mscale

  mdisc = Mstar*massratio

  ptot = pgas
  picm = 0
  pboundary = 0
  pcdm = 0

  write(6,*) "Inner Radius: ", rin
  write(6,*) "Outer Radius: ", rout
  write(6,*) "Central Object Mass: ", Mstar
  write(6,*) "Disc/Object Mass Ratio: ", massratio
  write(6,*) "Minimum Toomre Q: ", Qmin
  write(6,*) "Surface Density Profile: ", p_index
  write(6,*) "Sound Speed Profile: ", q_index
  write(6,*) "Disc Mass: ", mdisc
  write(6,*) "SPH Particle No.: ",pgas

! Allocate arrays for storage
  call allocate_memory

! Setting up scaling units for simulation 
  call units

! Set up particle types
  call types

! Calculate surface density at r=1, sigma_0
    
  if (p_index==2.0_PR) then
     sigma_0 = mdisc / (2.0_PR*pi*log(rout/rin))
  else
     rtest = (rout**(2.0_PR-p_index)-rin**(2.0_PR-p_index))
     sigma_0 = mdisc*((2.0_PR-p_index)/(2.0_PR*pi*rtest))
  end if
  
! Calculate sound speed at r=1, cs_0

      if (p_index-q_index-1.50_PR<0.0_PR) then         
            cs_0 = Qmin*pi*sigma_0*rout**(1.5_PR + q_index - p_index)         
      else         
            cs_0 = Qmin*pi*sigma_0*rin**(1.5_PR + q_index - p_index)     
      endif  

  sigma_0 = sigma_0*(mscale*m_SI)/(rscale*r_SI)**2.0_PR     
  cs_0 = cs_0*vscale*v_SI

  write(6,*) 'sigma_0: ', sigma_0
  write (6,*) 'cs_0:  ',cs_0

! Set up central sink at origin

   stot = 1
  
   sink(stot)%m = Mstar  
   sink(stot)%r(1:NDIM) = 0.0_PR
   sink(stot)%v(1:NDIM) = 0.0_PR
   sink(stot)%h = 1.0_PR ! Give sink small h to avoid timestepping issues  
   sink(stot)%tcreate = 0.0_PR
   sink(stot)%macc(1:DMDT_RANGE) = 0.0_PR
   sink(stot)%tacc(1:DMDT_RANGE) = 0.0_PR
   sink(stot)%angmom(1:3) = 0.0_PR
   sink(stot)%angmomnet(1:3) = 0.0_PR

   write(6,*) 'Sink placed at system origin'
   write(6,*) 'Mass: ', sink(1)%m

!	Calculate enclosed mass as a function of r
!	We will use this for calculating rotation speeds

	enc_m(:) = 0.0_PR

   do p = 1, ienc
         dr     = (rout - rin)/real(ienc-1)
         R_an(p)= rin + real(p-1)*dr
         r_c(p) = rin + real(p-0.5_PR)*dr
   if (p_index==2.0_PR) then 
     enc_m(p) = Mstar*(1.0_PR + massratio*(log(R_an(p)/rin)/ &
                  log(rout/rin)))
   else
     enc_m(p) = Mstar*(1.0_PR + massratio*((R_an(p)**(2.0_PR-p_index) - rin**(2.0_PR-p_index)) / &
                  (rout**(2.0_PR-p_index) - rin**(2.0_PR-p_index))))
   endif
         
   enddo

!  Set up maximum sigma values for accept-reject

   if (p_index<=1.0_PR) then
       f_max = rout**(1.0_PR-p_index)
   else
       f_max = rin**(1.0_PR-p_index) 
   endif

  print*, 'Placing gas particles'

  do p=1,pgas

! Set up particle masses, other properties

  sph(p)%m = mdisc/REAL(pgas)
  sph(p)%h = 0.0_PR

! Set up randomly positioned particles in (r,phi,z) coordinates
! Assume phi is uniformly distributed
! r placed according to surface density profile (accept-reject)

  phi = 2.0_PR*pi*ran2(iseed)

  randtest = 1.0_PR
  f = 0.0_PR

  DO WHILE(randtest>f)

  rtest = rin + (rout-rin)*ran2(iseed)
  f = rtest**(1.0_PR-p_index)

  randtest = f_max*ran2(iseed)

  ENDDO

! Convert to x,y

  sph(p)%r(1) = rtest*COS(phi)
  sph(p)%r(2) = rtest*SIN(phi)

!  print*, sph(p)%r(1), sph(p)%r(2), rtest,phi
! Give each particle a sound speed based on r and profile

  cs = cs_0*(rtest**(-q_index))

! Convert this to an internal energy using gamma
  sph(p)%temp = cs*cs*real((rho_SI*rhoscale)/(Pscale*P_SI),PR)/(gamma*Pconst)

#if defined(INTERNAL_ENERGY)
  sph(p)%u = cs*cs/(gamma*(gamma-1.0_PR))
#endif

! Calculate the local angular velocity
! Assume orbits are initially Keplerian
  itable = nint((rtest-r_c(1))/(r_c(2)-r_c(1))) + 1

! Calculation in code units
  omega_rot = sqrt(enc_m(itable)/rtest**3)

! Calculate the local scale height (in code units)
  H = cs/(vscale*v_SI*omega_rot)

! z placed assuming hydrostatic equilibrium: rho(z) = rho_0 e^{-z*z/(2*H*H)}
! (accept-reject)

  randtest = 1.0_PR
  f = 0.0_PR

  DO WHILE(randtest>f)
  
  ztest = -3.0_PR*H + 6.0_PR*H*ran2(iseed)
  f = exp(-ztest*ztest/(2.0_PR*H*H))
  randtest = ran2(iseed)

  ENDDO

  sph(p)%r(3) = ztest  
  sph(p)%rho = 0.0_PR
  
!  Now calculate velocities - first spherical separation

  rsphere = SQRT(rtest*rtest + ztest*ztest)

  vk = enc_m(itable)/ &
	(rsphere*rsphere*rsphere)
   vk = sqrt(vk)*rtest

! Set up particle velocities in x,y

  sph(p)%v(1) = -vk*SIN(phi)
  sph(p)%v(2) = vk*COS(phi)
  sph(p)%v(3) = 0.0_PR

  ENDDO	   ! End loop over gas particles

! Set system into the centre of mass frame

  cmx = 0.0
  cmy = 0.0	
  cmz = 0.0

  vcmx = 0.0
  vcmy = 0.0	
  vcmz = 0.0
  
  DO p =1,pgas
	
  cmx = cmx + sph(p)%r(1)*sph(p)%m
  cmy = cmy + sph(p)%r(2)*sph(p)%m
  cmz = cmz + sph(p)%r(3)*sph(p)%m

  vcmx = vcmx + sph(p)%v(1)*sph(p)%m
  vcmy = vcmy + sph(p)%v(2)*sph(p)%m
  vcmz = vcmz + sph(p)%v(3)*sph(p)%m

  ENDDO

  cmx = cmx/mdisc
  cmy = cmy/mdisc
  cmz = cmz/mdisc

  vcmx = vcmx/mdisc
  vcmy = vcmy/mdisc
  vcmz = vcmz/mdisc

  sph(:)%r(1) = sph(:)%r(1) - cmx
  sph(:)%r(2) = sph(:)%r(2) - cmy
  sph(:)%r(3) = sph(:)%r(3) - cmz

  sph(:)%v(1) = sph(:)%v(1) - vcmx
  sph(:)%v(2) = sph(:)%v(2) - vcmy
  sph(:)%v(3) = sph(:)%v(3) - vcmz

! Zero time variables

  nsteps = 0
  n = 0

! Write data to file
  call write_data(out_file,out_file_form)

! Clean up memory
  call clean_up

  STOP

CONTAINS
  REAL FUNCTION ran2(idum)
  INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1

          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END FUNCTION


END PROGRAM ic_kepdisc

