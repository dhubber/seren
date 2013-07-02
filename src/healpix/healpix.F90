! HEALPIX.F90
! ????
! Contains all HEALPix routines for the ionization algorithm.  
! ============================================================================



  !=======================================================================
  subroutine vec2ang(vector, theta, phi)
    !=======================================================================
    !     renders the angles theta, phi corresponding to vector (x,y,z)
    !     theta (co-latitude measured from North pole, in [0,Pi] radians)
    !     and phi (longitude measured eastward, in [0,2Pi[ radians)
    !     North pole is (x,y,z)=(0,0,1)
    !     added by EH, Feb 2000
    !=======================================================================
    use healpix_types

    REAL(KIND=DP), INTENT(IN), dimension(1:3) :: vector
    REAL(KIND=DP), INTENT(OUT) :: theta, phi

    REAL(KIND=DP) :: dnorm, z
    !=======================================================================

    dnorm = SQRT(vector(1)**2+vector(2)**2+vector(3)**2)

    z = vector(3) / dnorm
    theta = ACOS(z)

    phi = 0.0_dp
    if (vector(1) /= 0.0_dp .or. vector(2) /= 0.0_dp) &
         &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]
    if (phi < 0.0)     phi = phi + twopi ! phi in [0,2pi[

    return
  end subroutine vec2ang



  !=======================================================================
  subroutine mk_pix2xy(pix2x,pix2y)
    !=======================================================================
    !     constructs the array giving x and y in the face from pixel number
    !     for the nested (quad-cube like) ordering of pixels
    !
    !     the bits corresponding to x and y are interleaved in the pixel number
    !     one breaks up the pixel number by even and odd bits
    !=======================================================================
    use healpix_types

    INTEGER(KIND=I4B), INTENT(INOUT), DIMENSION(0:1023) :: pix2x, pix2y

    INTEGER(KIND=I4B) ::  kpix, jpix, ix, iy, ip, id

    !cc cf block data      data      pix2x(1023) /0/
    !-----------------------------------------------------------------------
    !      print *, 'initiate pix2xy'
    do kpix=0,1023          ! pixel number
       jpix = kpix
       IX = 0
       IY = 0
       IP = 1               ! bit position (in x and y)
!        do while (jpix/=0) ! go through all the bits
       do
          if (jpix == 0) exit ! go through all the bits
          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in ix
          jpix = jpix/2
          IX = ID*IP+IX

          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in iy
          jpix = jpix/2
          IY = ID*IP+IY

          IP = 2*IP         ! next bit (in x and y)
       enddo
       pix2x(kpix) = IX     ! in 0,31
       pix2y(kpix) = IY     ! in 0,31
    enddo

    return
  end subroutine mk_pix2xy



  !=======================================================================
  subroutine mk_xy2pix(x2pix,y2pix)
    !=======================================================================
    !     sets the array giving the number of the pixel lying in (x,y)
    !     x and y are in {1,128}
    !     the pixel number is in {0,128**2-1}
    !
    !     if  i-1 = sum_p=0  b_p * 2^p
    !     then ix = sum_p=0  b_p * 4^p
    !          iy = 2*ix
    !     ix + iy in {0, 128**2 -1}
    !=======================================================================
    use healpix_types

    INTEGER(KIND=I4B), INTENT(INOUT), DIMENSION(128) :: x2pix, y2pix

    INTEGER(KIND=I4B):: k,ip,i,j,id
    !=======================================================================

    do i = 1,128           !for converting x,y into
       j  = i-1            !pixel numbers
       k  = 0
       ip = 1

       do
          if (j==0) then
             x2pix(i) = k
             y2pix(i) = 2*k
             exit
          else
             id = MODULO(J,2)
             j  = j/2
             k  = ip*id+k
             ip = ip*4
          endif
       enddo

    enddo

    RETURN
  END subroutine mk_xy2pix
  !=======================================================================



  !=======================================================================
  subroutine ang2pix_nest(nside, theta, phi, ipix, x2pix, y2pix)
    !=======================================================================
    ! renders the pixel number ipix (NESTED scheme) for a pixel which contains
    ! a point on a sphere at coordinates theta and phi, given the map
    ! resolution parametr nside
    !
    ! the computation is made to the highest resolution available (nside=8192)
    ! and then degraded to that required (by integer division)
    ! this doesn't cost more, and it makes sure
    ! that the treatement of round-off will be consistent
    ! for every resolution
    !=======================================================================
    use healpix_types

    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=I4B), INTENT(OUT) :: ipix
    INTEGER(KIND=I4B), INTENT(IN), DIMENSION(128) :: x2pix, y2pix
    REAL(KIND=DP), INTENT(IN) ::  theta, phi

    REAL(KIND=DP) ::  z, za, tt, tp, tmp
    INTEGER(KIND=I4B) :: jp, jm, ifp, ifm, face_num, &
         &     ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    if (theta<0.0_dp .or. theta>pi)  then
       print*,"ANG2PIX_NEST: theta : ",theta," is out of range [0,Pi]"
       call fatal_error("Oops!!")
    endif
    !!if (x2pix(128) <= 0) call mk_xy2pix()

    z  = COS(theta)
    za = ABS(z)
    tt = MODULO(phi, twopi) / halfpi  ! in [0,4[

    if (za <= twothird) then ! equatorial region

       !        (the index of edge lines increase when the longitude=phi goes up)
       jp = INT(ns_max*(0.5_dp + tt - z*0.75_dp)) !  ascending edge line index
       jm = INT(ns_max*(0.5_dp + tt + z*0.75_dp)) ! descending edge line index

       !        finds the face
       ifp = jp / ns_max  ! in {0,4}
       ifm = jm / ns_max
       if (ifp == ifm) then          ! faces 4 to 7
          face_num = MODULO(ifp,4) + 4
       else if (ifp < ifm) then     ! (half-)faces 0 to 3
          face_num = MODULO(ifp,4)
       else                            ! (half-)faces 8 to 11
          face_num = MODULO(ifm,4) + 8
       endif

       ix = MODULO(jm, ns_max)
       iy = ns_max - MODULO(jp, ns_max) - 1

    else ! polar region, za > 2/3

       ntt = INT(tt)
       if (ntt >= 4) ntt = 3
       tp = tt - ntt
       tmp = SQRT( 3.0_dp*(1.0_dp - za) )  ! in ]0,1]

       !        (the index of edge lines increase when distance from the closest pole goes up)
       jp = INT( ns_max * tp          * tmp ) ! line going toward the pole as phi increases
       jm = INT( ns_max * (1.0_dp - tp) * tmp ) ! that one goes away of the closest pole
       jp = MIN(ns_max-1, jp) ! for points too close to the boundary
       jm = MIN(ns_max-1, jm)

       !        finds the face and pixel's (x,y)
       if (z >= 0) then
          face_num = ntt  ! in {0,3}
          ix = ns_max - jm - 1
          iy = ns_max - jp - 1
       else
          face_num = ntt + 8 ! in {8,11}
          ix =  jp
          iy =  jm
       endif

       !         print*,z,face_num,ix,iy
    endif

    ix_low = MODULO(ix,128)
    ix_hi  =     ix/128
    iy_low = MODULO(iy,128)
    iy_hi  =     iy/128

    ipf =  (x2pix(ix_hi +1)+y2pix(iy_hi +1)) * (128 * 128) &
         &     + (x2pix(ix_low+1)+y2pix(iy_low+1))

    ipf = ipf / ( ns_max/nside ) **2  ! in {0, nside**2 - 1}

    ipix = ipf + face_num* nside **2    ! in {0, 12*nside**2 - 1}

    return
  end subroutine ang2pix_nest



  !=======================================================================
  subroutine pix2vec_nest(nside, ipix, pix2x, pix2y, vector)
    !=======================================================================
    !     renders vector (x,y,z) coordinates of the nominal pixel center
    !     for the pixel number ipix (NESTED scheme)
    !     given the map resolution parameter nside
    !     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
    !     in the order N,W,S,E
    !=======================================================================
    use healpix_types

    INTEGER(KIND=I4B), INTENT(IN) :: nside, ipix
    REAL(KIND=DP), INTENT(OUT), dimension(1:3) :: vector
!    REAL(KIND=DP),     INTENT(OUT),dimension(1:3,1:4), optional :: vertex
    INTEGER(KIND=I4B), INTENT(INOUT), DIMENSION(0:1023) :: pix2x, pix2y

    INTEGER(KIND=I4B) :: npix, npface, &
         &     ipf, ip_low, ip_trunc, ip_med, ip_hi, &
         &     jrt, jr, nr, jpt, jp, kshift, nl4
    REAL(KIND=DP) :: z, fn, fact1, fact2, sth, phi

    INTEGER(KIND=I4B) ::  ix, iy, face_num
!     common /xy_nest/ ix, iy, face_num ! can be useful to calling routine

    ! coordinate of the lowest corner of each face
    INTEGER(KIND=I4B), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(KIND=I4B), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2

!    real(kind=DP) :: phi_nv, phi_wv, phi_sv, phi_up, phi_ev, phi_dn
!    real(kind=DP) :: z_nv, z_sv, sth_nv, sth_sv
!    real(kind=DP) :: hdelta_phi
!    integer(kind=I4B) :: iphi_mod, iphi_rat
!    logical(kind=LGT) :: do_vertex
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    npix = 12 * nside**2
    if (ipix <0 .or. ipix>npix-1) call fatal_error("ipix out of range")

    !     initiates the array for the pixel number -> (x,y) mapping
    !!if (pix2x(1023) <= 0) call mk_pix2xy()

    fn = real(nside,kind=dp)
    fact1 = 1.0_dp/(3.0_dp*fn*fn)
    fact2 = 2.0_dp/(3.0_dp*fn)
    nl4   = 4*nside

!    do_vertex = .false.

    !     finds the face, and the number in the face
    npface = nside**2

    face_num = ipix/npface  ! face number in {0,11}
    ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    ip_low = MODULO(ipf,1024)       ! content of the last 10 bits
    ip_trunc =   ipf/1024        ! truncation of the last 10 bits
    ip_med = MODULO(ip_trunc,1024)  ! content of the next 10 bits
    ip_hi  =     ip_trunc/1024   ! content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}

    nr = nside                  ! equatorial region (the most frequent)
    z  = (2*nside-jr)*fact2
    kshift = MODULO(jr - nside, 2)

    if (jr < nside) then     ! north pole region
       nr = jr
       z = 1.0_dp - nr*nr*fact1
       kshift = 0
    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       z = - 1.0_dp + nr*nr*fact1
       kshift = 0
    endif

    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4

    phi = (jp - (kshift+1)*0.5_dp) * (halfpi / nr)

    sth = SQRT((1.0_dp-z)*(1.0_dp+z))
    vector(1) = sth * COS(phi)
    vector(2) = sth * SIN(phi)
    vector(3) = z


    return
  end subroutine pix2vec_nest
  !=======================================================================



  subroutine fatal_error (msg)
    character(len=*), intent(in) :: msg
       print *,'Fatal error: ', trim(msg)
       stop
  end subroutine fatal_error
