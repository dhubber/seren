! Bondi.F90
! A. P. Whitworth
! ~~~~~~~~~
! Programme to generate parameters of the Bondi accretion flow. 
! In the output file "BONDI.DAT", 
! the first column gives the ID, n; 
! the second, the dimensionless radius, x(n); 
! the third, the dimensionless mass, z(n);
! the fourth, the dimensionless density, y(n); 
! the fifth, the dimensionless velocity, w(n). 
! ============================================================================
SUBROUTINE bondi_isothermal_solution(ntable,lw,lx,ly,lz)
  use definitions
  IMPLICIT NONE

  integer, intent(in) :: ntable
  REAL(KIND=DP), intent(inout) :: lw(0:1000)     ! tabulated log[w]-values
  REAL(KIND=DP), intent(inout) :: lx(0:1000)     ! tabulated log[x]-values
  REAL(KIND=DP), intent(inout) :: ly(0:1000)     ! tabulated log[y]-values
  REAL(KIND=DP), intent(inout) :: lz(0:1000)     ! tabulated log[z]-values

  REAL(KIND=DP)                  :: asympt         ! asymptotic value of yx^2
  REAL(KIND=DP)                  :: disc           ! convergence on asymptotic solution
  REAL(KIND=DP)                  :: dlnx(0:1000)   ! fractional displacement
  REAL(KIND=DP)                  :: dx             ! integration stepsize
  REAL(KIND=DP)                  :: epsilon        ! step size
  REAL(KIND=DP)                  :: f1,f2,f3,f4,f5 ! dwdx-values for single R-K step
  REAL(KIND=DP)                  :: g1,g2,g3,g4,g5 ! dzdx-values for single R-K step
  REAL(KIND=DP)                  :: hdx            ! half stepsize
  INTEGER                        :: n              ! ID of tabulated value
  REAL(KIND=DP)                  :: w(0:1000)      ! tabulated w-values
  REAL(KIND=DP)                  :: w0             ! asymptotic w
  REAL(KIND=DP)                  :: w1,w2,w3,w4    ! w=values for single R-K step
  REAL(KIND=DP)                  :: x(0:1000)      ! tabulated x-values
  REAL(KIND=DP)                  :: x0             ! asymptotic x
  REAL(KIND=DP)                  :: x1,x2,x3,x4    ! x-values for single R-K step
  REAL(KIND=DP)                  :: y(0:1000)      ! tabulated y-values
  REAL(KIND=DP)                  :: y0             ! asymptotic y
  REAL(KIND=DP)                  :: z(0:1000)      ! tabulated z-values
  REAL(KIND=DP)                  :: z0             ! asymptotic z
  REAL(KIND=DP)                  :: z1             ! z-value

  write(6,*) "ntable : ",ntable

                                                ! FORMAT STATEMENTS
                                                ! ~~~~~~~~~~~~~~~~~
602 FORMAT (/,X,'INDEX:',3X,'RADIUS(x):',24X, &
 &          'VELOCITY(w):',17X,'DENSITY(y):', &
 &          18X,'MASS(z):')                     ! format for header
603 FORMAT (2X,I5,3X,E10.3,F7.3,2X,E10.3,4X,     &
 &          E11.4,F7.3,F7.3,4X,E10.3,F7.3,    &
 &          F7.3,3X,E14.7)                      ! format for printout

                                                ! INITIALISATION
                                                ! ~~~~~~~~~~~~~~
epsilon=0.001                                   ! set accuracy of integration step
WRITE (6,602)                                   ! write header
asympt=EXP(1.5)                                 ! evaluate asympt=e^3/2 
x(700)=1.                                       ! set x(700) to 1 (sonic radius)
lx(700)=0.                                      ! set lx(700) to 0
w(700)=1.                                       ! set w(700) to 1 (sonic inflow velocity)
lw(700)=0.                                      ! set lw(700) to 0
y(700)=asympt                                   ! set y(700) to e^3/2
ly(700)=LOG10(y(700))                           ! compute ly(700)
z(700)=2.4102434440                             ! set z(700) to 0
lz(700)=LOG10(z(700))                           ! compute lz(700)

                                                ! OUTWARD INTEGRATION
                                                ! ~~~~~~~~~~~~~~~~~~~
x1=1.+epsilon                                   ! advance outward from sonic point
w1=1.-epsilon                                   ! adjust w
z1=2.4102434440+asympt*epsilon                  ! adjust z
n=701                                           ! set n=701
DO WHILE (n<1001)                               ! loop out to n=1000
  disc=100.*LOG10(x1)+700.-DBLE(n)              !   discriminator for printout
  IF (disc>0.) THEN                             !   [IF] printout 
    lx(n)=0.01*DBLE(n-700)                      !     compute log[x(n)]
    x(n)=10.**lx(n)                             !     compute x(n)
    dlnx(n)=(x(n)-x1)/x1                        !     compute fractional displacement
    w(n)=w1+f5*(x(n)-x1)                        !     compute w(n)
    lw(n)=LOG10(w(n))                           !     compute log[w(n)]
    y(n)=asympt/(x(n)*x(n)*w(n))                !     compute y(n)
    ly(n)=LOG10(y(n))                           !     compute log[y(n)]
    z(n)=z1+g5*(x(n)-x1)                        !     compute z(n)
    lz(n)=LOG10(z(n))                           !     compute log[w(n)]
    n=n+1                                       !     advance n
  ENDIF                                         !   [ENDIF]
  dx=x1*epsilon                                 !   compute step
  hdx=0.5*dx                                    !   compute half-step
  f1=2.*((1./x1)-(1./(x1*x1)))/(w1-(1./w1))     !   compute f1
  g1=asympt/w1                                  !   compute g1
  x2=x1+hdx                                     !   advance x
  w2=w1+f1*hdx                                  !   advance w
  f2=2.*((1./x2)-(1./(x2*x2)))/(w2-(1./w2))     !   compute f2
  g2=asympt/w2                                  !   compute g2
  w3=w1+f2*hdx                                  !   advance w
  f3=2.*((1./x2)-(1./(x2*x2)))/(w3-(1./w3))     !   compute f3
  g3=asympt/w3                                  !   compute g3
  x4=x1+dx                                      !   advance x
  w4=w1+f3*dx                                   !   advance w
  f4=2.*((1./x4)-(1./(x4*x4)))/(w4-(1./w4))     !   compute f4
  g4=asympt/w4                                  !   compute g4
  x1=x4                                         !   advance x fully
  f5=(f1+2.*f2+2.*f3+f4)/6.                     !   compute f5
  w1=w1+f5*dx                                   !   advance w fully
  g5=(g1+2.*g2+2.*g3+g4)/6.                     !   compute g5
  z1=z1+g5*dx                                   !   advance z fully
ENDDO                                           ! [ENDIF]

                                                ! INWARD INTEGRATION
                                                ! ~~~~~~~~~~~~~~~~~~
epsilon=-epsilon                                ! reverse epsilon
x1=1.+epsilon                                   ! advance inward
w1=1.-epsilon                                   ! adjust w
z1=2.4102434440+asympt*epsilon                  ! adjust z
n=699                                           ! set n=699
DO WHILE (n>-1)                                 ! loop inwards to n=0
  disc=DBLE(n)-100.*LOG10(x1)-700.              !   discriminator for printout
  IF (disc>0.) THEN                             !   [IF] printout
    lx(n)=0.01*DBLE(n-700)                      !     compute log[x(n)]
    x(n)=10.**lx(n)                             !     compute x(n)
    dlnx(n)=(x(n)-x1)/x1                        !     compute fractional displacement
    w(n)=w1+f5*(x(n)-x1)                        !     compute w(n)
    lw(n)=LOG10(w(n))                           !     compute log[w(n)]
    y(n)=asympt/(x(n)*x(n)*w(n))                !     compute y(n)
    ly(n)=LOG10(y(n))                           !     compute log[y(n)]
    z(n)=z1+g5*(x(n)-x1)                        !     compute z(n)
    lz(n)=-8.                                   !     set lz(n) to -8. (precautionary)
    IF (z(n)>0.1E-07) lz(n)=LOG10(z(n))         !     compute log[z(n)] if >-8.
    n=n-1                                       !     advance n
  ENDIF                                         !   [ENDIF]
  dx=x1*epsilon                                 !   compute step
  hdx=0.5*dx                                    !   compute half-step
  f1=2.*((1./x1)-(1./(x1*x1)))/(w1-(1./w1))     !   compute f1
  g1=asympt/w1                                  !   compute g1
  x2=x1+hdx                                     !   advance x
  w2=w1+f1*hdx                                  !   advance w
  f2=2.*((1./x2)-(1./(x2*x2)))/(w2-(1./w2))     !   compute f2
  g2=asympt/w2                                  !   compute g2
  w3=w1+f2*hdx                                  !   advance w
  f3=2.*((1./x2)-(1./(x2*x2)))/(w3-(1./w3))     !   compute f3
  g3=asympt/w3                                  !   compute g3
  x4=x1+dx                                      !   advance x
  w4=w1+f3*dx                                   !   advance w
  f4=2.*((1./x4)-(1./(x4*x4)))/(w4-(1./w4))     !   compute f4
  g4=asympt/w4                                  !   compute g4
  x1=x4                                         !   advance x fully
  f5=(f1+2.*f2+2.*f3+f4)/6.                     !   compute f5
  w1=w1+f5*dx                                   !   advance w fully
  g5=(g1+2.*g2+2.*g3+g4)/6.                     !   compute g5
  z1=z1+g5*dx                                   !   advance z fully
ENDDO                                           ! [ENDDO]

                                                ! PRINT TO SCREEN
                                                ! ~~~~~~~~~~~~~~~
DO n=0,1000,10                                  ! loop over array
  IF (n<700) THEN                               !   [IF] supersonic, compute freefall asymptotes
    w0=0.30103-0.5*lx(n)                        !     freefall velocity
    y0=1.5*(0.43429-lx(n))-0.30103              !     freefall density
    z0=1.5*(0.43429+lx(n))-0.47712              !     freefall mass
  ELSE                                          !   [ELSE] compute subsonic merging asymptotes
    w0=0.65144-2.*lx(n)                         !     faraway velocity
    y0=1.                                       !     faraway density
    z0=3.*lx(n)                                 !     faraway mass
  ENDIF                                         !   [ENDIF]
  WRITE (6,603) n,x(n),lx(n),dlnx(n),w(n),    &
 &              lw(n),w0,y(n),ly(n),y0,z(n)     !   printout
ENDDO                                           ! end loop over values near the sonic point

                                                ! WRITE TO FILE
                                                ! ~~~~~~~~~~~~~
111 FORMAT (I8,4(F12.5))                        ! file format
OPEN (UNIT=1,FILE="BONDI.DAT",STATUS='REPLACE') ! open file
DO n=0,1000                                     ! loop over IDs
  WRITE (1,111) n,lx(n),lz(n),ly(n),lw(n)       !   record ID, position and velocity
ENDDO                                           ! end loop over IDs
CLOSE(1)                                        ! close fle

                                                ! END PROGRAMME
                                                ! ~~~~~~~~~~~~~
RETURN
END SUBROUTINE bondi_isothermal_solution        ! end
