      SUBROUTINE VELFIELD(nindex,iseed,iseed2,ngrid,vel)
!
!==============================================================
!  Generates a random turbulent velocity field which is
!  divergence-free.  Follows Dubinski, Narayan & Phillips 1985.
!  Based on original code for Zeldovich shift from Volker Bromm.
!  Velocity field code written by M. R. Bate (21/11/2000).
!
!==============================================================
!
      IMPLICIT NONE
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: &
&                                phix,phiy,phiz,pow,ampx,ampy,ampz
      INTEGER :: iseed,iseed2,ngrid,kx,ky,kz,ii,jj,kk
      INTEGER :: i,j,k
      DOUBLE PRECISION :: powsum,RTOT,velx,vely,velz
      DOUBLE PRECISION :: zdmjk,ydmjk,xdmjk,contrib,AA
      DOUBLE PRECISION :: ran3,ran4,rayldev,sigma2,pi
      DOUBLE PRECISION :: nindex,kmod,kdotq,LF
      REAL :: vel(4,2*NGRID,2*NGRID,2*NGRID)
      EXTERNAL ran3,ran4,rayldev
!
      WRITE(6,*) 'getting velocity field'
!
      PI=3.1415927D0
      RTOT=1.0d0
      LF=2.0d0*RTOT
      AA=0.5d-1/32.52d0
!
      WRITE(6,*) 'nindex', nindex
      WRITE(6,*) 'nrands', iseed, iseed2
      WRITE(6,*) 'ngrid', ngrid
!
      ALLOCATE(phix(1:2*NGRID,1:2*NGRID,1:2*NGRID))
      ALLOCATE(phiy(1:2*NGRID,1:2*NGRID,1:2*NGRID))
      ALLOCATE(phiz(1:2*NGRID,1:2*NGRID,1:2*NGRID))
      ALLOCATE(pow(1:2*NGRID,1:2*NGRID,1:2*NGRID))
      ALLOCATE(ampx(1:2*NGRID,1:2*NGRID,1:2*NGRID))
      ALLOCATE(ampy(1:2*NGRID,1:2*NGRID,1:2*NGRID))
      ALLOCATE(ampz(1:2*NGRID,1:2*NGRID,1:2*NGRID))
!
!--Gives P(k)=k^{nindex) power spectrum
!     Note: rayldev returns the square root of the -Log of a random number
      powsum=0.d0
      DO k=1,NGRID
        DO j=1,2*NGRID
          DO i=1,2*NGRID
            phix(i,j,k)=(-PI + 2.d0*PI*RAN3(iseed))
            phiy(i,j,k)=(-PI + 2.d0*PI*RAN3(iseed))
            phiz(i,j,k)=(-PI + 2.d0*PI*RAN3(iseed))
            kmod=SQRT(REAL((i-NGRID)**2+(j-NGRID)**2+k*k))
            pow(i,j,k)=AA*(kmod)**(nindex)
            powsum=powsum+pow(i,j,k)
            ampx(i,j,k)=rayldev(iseed2)*SQRT(pow(i,j,k))
            ampy(i,j,k)=rayldev(iseed2)*SQRT(pow(i,j,k))
            ampz(i,j,k)=rayldev(iseed2)*SQRT(pow(i,j,k))
          ENDDO
        ENDDO
      ENDDO
!
      WRITE(6,*) 'done first loop'
!
      sigma2=0.d0
      DO kk = 1, 2*NGRID
         DO jj = 1, 2*NGRID
            DO ii = 1, 2*NGRID
               xdmjk = (REAL(ii-NGRID)-0.5)/REAL(NGRID)
               ydmjk = (REAL(jj-NGRID)-0.5)/REAL(NGRID)
               zdmjk = (REAL(kk-NGRID)-0.5)/REAL(NGRID)
               velx = 0.
               vely = 0.
               velz = 0.

!$OMP PARALLEL DO default(NONE) shared(ngrid,pi) &
!$OMP shared(LF,xdmjk,ydmjk,zdmjk,ampx,ampy,ampz) &
!$OMP shared(phix,phiy,phiz) &
!$OMP PRIVATE(i,j,k,kx,ky,kz,kdotq,contrib) &
!$OMP reduction(+:velx,vely,velz)
               DO k=1,NGRID
                  DO j=1,2*NGRID
                     DO i=1,2*NGRID
                        kx = i-NGRID
                        ky = j-NGRID
                        kz = k
                        kdotq=(2.d0*PI/LF)*(REAL(kx)*XDMJK+          &
     &                        REAL(ky)*YDMJK+REAL(kz)*ZDMJK)
                        contrib = ampz(i,j,k)*(2.d0*PI/LF)*REAL(ky)* &
     &                       2.d0*SIN(kdotq+phiz(i,j,k))             &
     &                       - ampy(i,j,k)*(2.d0*PI/LF)*REAL(kz)*    &
     &                       2.d0*SIN(kdotq+phiy(i,j,k))
                        velx = velx + contrib 
                        contrib = ampx(i,j,k)*(2.d0*PI/LF)*REAL(kz)* &
     &                       2.d0*SIN(kdotq+phix(i,j,k))             &
     &                       - ampz(i,j,k)*(2.d0*PI/LF)*REAL(kx)*    &
     &                       2.d0*SIN(kdotq+phiz(i,j,k))
                        vely = vely + contrib
                        contrib = ampy(i,j,k)*(2.d0*PI/LF)*REAL(kx)* &
     &                       2.d0*SIN(kdotq+phiy(i,j,k))             &
     &                       - ampx(i,j,k)*(2.d0*PI/LF)*REAL(ky)*    &
     &                       2.d0*SIN(kdotq+phix(i,j,k))
                        velz = velz + contrib
                     ENDDO
                  ENDDO
               ENDDO
!$OMP END PARALLEL DO
!
               vel(1,ii,jj,kk) = velx
               vel(2,ii,jj,kk) = vely
               vel(3,ii,jj,kk) = velz
            ENDDO
         ENDDO
      ENDDO
!
      WRITE(6,*) 'done second loop'
!
      DO kk = 1, 2*NGRID
         DO jj = 1, 2*NGRID
            DO ii = 1, 2*NGRID
               vel(4,ii,jj,kk) = SQRT(vel(1,ii,jj,kk)**2 +     &
     &              vel(2,ii,jj,kk)**2 + vel(3,ii,jj,kk)**2)
            ENDDO
         ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE velfield
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
! ===================================================================
      DOUBLE PRECISION FUNCTION rayldev(idum)
      INTEGER :: idum
      DOUBLE PRECISION :: ran4
      EXTERNAL ran4
      rayldev=dsqrt(-LOG(ran4(idum)))
      RETURN
      END
! ===================================================================
        DOUBLE PRECISION FUNCTION RAN4(IDUM)                                   
!
        IMPLICIT NONE
        DOUBLE PRECISION :: FAC,mbig,mseed,mz,mj,ma,mk
        INTEGER :: idum,i,ii,k,iff,inext,inextp
        PARAMETER (MBIG=4000000., MSEED=1618033., MZ=0., FAC=1./MBIG)
          DIMENSION MA(55)
              SAVE
            DATA IFF/0/
          IF(IDUM<0.OR.IFF==0)THEN
             IFF=1
             MJ=MSEED-IABS(IDUM)
             MJ=MOD(MJ,MBIG)
             MA(55)=MJ
             MK=1
            DO 100 I=1,54
              II=MOD(21*I,55)
              MA(II)=MK
              MK=MJ-MK
              IF(MK<MZ)MK=MK+MBIG
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
            IF(MJ.LT.MZ)MJ=MJ+MBIG
            MA(INEXT)=MJ
            RAN4=MJ*FAC
             RETURN
             END



