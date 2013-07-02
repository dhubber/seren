! SINK_SEARCH.F90
! D. A. Hubber - 12/6/2007; 9/6/2010
! Searches through all gas particles for possible sink candidates.  
! Creates a sink particle once the following (optional) conditions are met:
! (i)   the density of the particle exceeds the sink density (rhosink).
! (ii)  the particle lies at the bottom of its local gravitational potential 
!       well (i.e. no SPH neighbour has a more negative potential).
! (iii) the candidate particle is not overlapping with any other sinks.
! (iv)  the Hill sphere of the particle does not overlap the Hill sphere of 
!       any existing sinks.
! (v)   the local acceleration divergence is negative, i.e. particles
!       are not accelerating away from local maximum.  This would otherwise
!       indicate a transient object that is 'bouncing' back out possibly
!       due to tidal forces of a nearby dense object/sink particle.
! (vi)  the local velocity divergence is negative, i.e. particles are on 
!       on average converging.
! (vii) the particle is unlikely to be accreted by an existing sink before 
!       it has time to condense out to a new object.
! Currently creates only one sink in any one timestep for simplicity.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE sink_search(psink)
  use interface_module, only : create_sink,distance2,distance3,&
       &gather_neib_on_fly,insertion_sort_real,w0,w2
  use type_module
  use particle_module
  use hydro_module
  use sink_module
  use scaling_module
  use constant_module
  use kernel_module
  implicit none

  integer, intent(out) :: psink          ! i.d. of candidate sink particle
  
  logical, allocatable :: testflag(:)    ! Has ptcl passed all previous tests?
  integer :: i                           ! Auxilary loop counter
  integer :: iaux                        ! Aux. variable
  integer :: j                           ! Aux. loop variable
  integer :: jj                          ! Aux. loop variable
  integer :: p                           ! Particle counter
  integer :: pcandidate                  ! Counts number of particles in list
  integer :: pp                          ! Neighbour particle id
  integer :: ppp                         ! Aux. particle id
  integer :: pp_max                      ! Max. no. of potential neighbours
  integer :: pp_pot                      ! No. of potential neighbours
  integer :: s                           ! Sink counter
  integer, allocatable :: pp_potlist(:)  ! List of potential neighbours
  integer, allocatable :: rholist(:)     ! List of ptcls with rho > rhosink
  real(kind=PR) :: agravp(1:NDIM)        ! Grav. accel of particle p
  real(kind=PR) :: ap(1:VDIM)            ! Total accel of particle p
  real(kind=PR) :: arad                  ! Radial acceleration
  real(kind=PR) :: da(1:NDIM)            ! Relative acceleration
  real(kind=PR) :: div_a_p               ! Divergence of acceleration for p
  real(kind=PR) :: div_v_p               ! Divergence of velocity for p
  real(kind=PR) :: drad                  ! Radial distance
  real(kind=PR) :: dr(1:NDIM)            ! Relative position vector 
  real(kind=PR) :: dr_unit(1:NDIM)       ! Relative unit vector
  real(kind=PR) :: drmag                 ! Distance
  real(kind=PR) :: drsqd                 ! Distance squared
  real(kind=PR) :: dv(1:NDIM)            ! Relative velocity
  real(kind=PR) :: dvdr                  ! Scalar product of dv and dr
  real(kind=PR) :: gpetot                ! Total grav. pot. energy of clump
  real(kind=PR) :: hfactor               ! invhp^(NDIM)
  real(kind=PR) :: hrange                ! Tree search radius
  real(kind=PR) :: hrangesqd             ! hrange*hrange
  real(kind=PR) :: invdrmag              ! 1 / drmag
  real(kind=PR) :: invhp                 ! Inverse of smoothing length
  real(kind=PR) :: ketot                 ! Total k.e. of clump
  real(kind=PR) :: mpp                   ! Mass of neighbour pp
  real(kind=PR) :: dpotp                 ! Grav. potential between 2 ptcls
  real(kind=PR) :: raux                  ! Aux. real variable
  real(kind=PR) :: rho_hill              ! Density of gas inside Hill sphere
  real(kind=PR) :: rhotemp               ! Aux. density summation variable
  real(kind=PR) :: rhop                  ! Density of clump
  real(kind=PR) :: RKH                   ! Kelvin-Helmholz contraction rate
  real(kind=PR) :: rotketot              ! Total rot. k.e. of clump
  real(kind=PR) :: rp(1:NDIM)            ! Position of particle p
  real(kind=PR) :: rs(1:NDIM)            ! Position of sink particle s
  real(kind=PR) :: taccrete              ! Accretion timescale
  real(kind=PR) :: tcollapse             ! Collapse timescale
  real(kind=PR) :: utot                  ! Total thermal energy of clump
  real(kind=PR) :: vp(1:VDIM)            ! Velocity of particle p
  real(kind=PR) :: vrad                  ! Radial velocity
  real(kind=PR) :: vtansqd               ! Relative tangential velocity
  real(kind=PR), allocatable :: rhovalues(:)  ! rho for candidate particles

  real(kind=PR) :: Mpoly
  real(kind=PR) :: Rpoly
  real(kind=PR) :: thetafac
  real(kind=PR) :: angmom(1:3)
  real(kind=PR) :: Eorbit,aorbit,ecc,rperi,Lsqd,reduced_mass

#if defined(USE_MPI)
  integer :: densest_task                ! MPI task with densest particle
#endif

  debug2("Searching for candidate sink particles [sink_search.F90]")
  debug_timing("SINK_SEARCH")

#if defined(USE_MPI)
     if (div_a_search .OR. div_v_search .OR. &
         &thermal_search .OR. timescale_search) then
        write (6,*) "div_a / div_v / thermal_search / timescale_search "//&
                   &"not implemented for MPI!"
     end if
#endif

! Initialize some values
  pcandidate = 0
  psink = -1
  iaux = 0

! First, find no. of particles with rho > rhosink which lie at the bottom of 
! their local potential well.
  do p=pgasstart,pgasend
     if (rho_search .and. sph(p)%rho < rhosink) cycle
     if (potmin_search .and. (.not. sph(p)%ispotmin)) cycle
     pcandidate = pcandidate + 1
  end do

! If there are no sink candidates, return from subroutine immediatly.
  if (pcandidate == 0) then
#if defined(USE_MPI)
     call mpi_find_densest_particle(0._PR,densest_task)
#endif
     return
  end if

! First, allocate required memory for sinks, and then create a list of all 
! particles with density greater than the sink density that also lie at the 
! bottom of their local potential minimum.
  allocate(rholist(1:pcandidate))
  allocate(rhovalues(1:pcandidate))
  allocate(testflag(1:pcandidate))
  do p=pgasstart,pgasend
     if (rho_search .and. sph(p)%rho < rhosink) cycle
     if (potmin_search .and. (.not. sph(p)%ispotmin)) cycle
     iaux = iaux + 1
     rholist(iaux) = p
     rhovalues(iaux) = sph(p)%rho
     testflag(iaux) = .true.
  end do

! Order rholist and rhovalues in ascending order of rho
  call insertion_sort_real(pcandidate,&
       &rholist(1:pcandidate),rhovalues(1:pcandidate))

#if defined(DEBUG_SINK_SEARCH)
  write(6,*) "Found ",pcandidate," sink candidate particles"
  if (pcandidate > 0) then
     do i=1,pcandidate
        p = rholist(i)
        write(6,*) p,sph(p)%porig,&
             &sph(p)%rho*rhoscale*rhocgs,sph(p)%ispotmin,sph(p)%gpot
     end do
  end if
#endif


! Determine if any of the candidate sinks overlap existing sinks, and if not, 
! see if the Hills sphere of any particles overlap those of other nearby sinks.
! ============================================================================
  if (stot > 0) then

     ! Loop over list of candidate particles
     ! -----------------------------------------------------------------------
     do i=1,pcandidate
        p = rholist(i)

        ! Set variables for candidate sink
        rhop           = sph(p)%rho
        ap(1:NDIM)     = sph(p)%a(1:NDIM)
        rp(1:NDIM)     = sph(p)%r(1:NDIM)
        agravp(1:NDIM) = sph(p)%a_grav(1:NDIM)

        ! Now loop over all existing sinks
        ! --------------------------------------------------------------------
        do s=1,stot
           call distance3(rp(1:NDIM),sink(s)%r(1:NDIM),dr(1:NDIM),drsqd)

           ! Sink-overlap criterion
           if (drsqd < (NEW_SINK_RMAX*2.0_PR*sink(s)%radius)**2) &
                & testflag(i) = .false.
#if defined(DEBUG_SINK_SEARCH)
           if (drsqd < (NEW_SINK_RMAX*2.0_PR*sink(s)%radius)**2) write(6,*) &
                &"Overlap between sink ",s," and particle ",p,sph(p)%porig
#endif
           ! Hill-sphere search
           if (hill_sphere_search) then
              da(1:NDIM) = sink(s)%agrav(1:NDIM) - agravp(1:NDIM)
              rho_hill = -(3.0_PR*dot_product(da(1:NDIM),dr(1:NDIM))) &
                   & / (4.0_PR*PI*drsqd)
              if (rhop < 3.0_PR*rho_hill) testflag(i) = .false.
#if defined(DEBUG_SINK_SEARCH)
              if (rhop < 3.0_PR*rho_hill) write(6,*) &
                   &"Hill criteria violated : ",&
                   &p,sph(p)%porig,s,rhop,3.0_PR*rho_hill
#endif
           end if

        end do
        ! --------------------------------------------------------------------

     end do
     ! -----------------------------------------------------------------------

  end if
! ============================================================================


! Record new shortened list of candidate sink particles
  iaux = 0
  do i=1,pcandidate
     if (testflag(i)) then
        iaux = iaux + 1
        rholist(iaux) = rholist(i)
        testflag(iaux) = testflag(i)
     end if
  end do
  pcandidate = iaux
  
! If there are no sink candidates, return from subroutine immediatly
  if (pcandidate == 0) then
     if (allocated(testflag)) deallocate(testflag)
     if (allocated(rhovalues)) deallocate(rhovalues)
     if (allocated(rholist)) deallocate(rholist)
#if defined(USE_MPI)
     call mpi_find_densest_particle(0._PR,densest_task)
#endif
     return
  end if

#if !defined(USE_MPI)
! Now loop over all candidate particles and calculate relevant SPH quantities 
! with enlarged smoothing kernel (if selected in params file).
! ============================================================================
  do i=1,pcandidate
     p = rholist(i)

     rp(1:NDIM) = sph(p)%r(1:NDIM)
     vp(1:NDIM) = sph(p)%v(1:NDIM)
     ap(1:NDIM) = sph(p)%a(1:NDIM)

     ! Find all potential gather neighbours within required radius
#if defined(HMULT_SINKRAD)
     hrange = max(sinkrad*sph(p)%h,KERNRANGE*sph(p)%h)
#else
     hrange = max(sinkrad,KERNRANGE*sph(p)%h)
#endif
     hrangesqd = hrange*hrange
     invhp = KERNRANGE / hrange
     hfactor = invhp**(NDIM)
     pp_max = max(LISTSIZE,ptot)
     if (allocated(pp_potlist)) deallocate(pp_potlist)
     call gather_neib_on_fly(p,pp_max,pp_pot,pp_potlist,rp(1:NDIM),hrange)

     ! First, add self-contribution to various SPH quantities
     rhotemp  = 0.0_PR
     div_v_p  = 0.0_PR
     div_a_p  = 0.0_PR
     ketot    = 0.0_PR
     gpetot   = 0.0_PR
     utot     = 0.0_PR
     rotketot = 0.0_PR

     ! Now loop over all gather neighbours and add contributions to various
     ! SPH quantities, such as div_v and div_a.
     ! -----------------------------------------------------------------------
     do j=1,pp_pot
        pp = pp_potlist(j)
        call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
        if (drsqd > hrangesqd) cycle
        drmag = sqrt(drsqd) + SMALL_NUMBER
        invdrmag = 1.0_PR / drmag
        mpp = sph(pp)%m
        dv(1:NDIM) = sph(pp)%v(1:NDIM) - vp(1:NDIM)
        dvdr = dot_product(dv(1:NDIM),dr(1:NDIM))
        rhotemp = rhotemp + mpp*hfactor*w0(drmag*invhp)
        div_v_p = div_v_p - mpp*dvdr*w2(drmag*invhp)*hfactor*invhp*invdrmag
        div_a_p = div_a_p - mpp*w2(drmag*invhp)*hfactor*invhp/&
             &dot_product(sph(pp)%a(1:NDIM) - ap(1:NDIM),dr(1:NDIM))*invdrmag
        
        ! Now compute all energy contributions
        ketot = ketot + 0.5_PR*sph(pp)%m*dot_product(dv(1:VDIM),dv(1:VDIM))
        rotketot = rotketot + &
             &0.5_PR*mpp*dot_product(dv(1:VDIM),dv(1:VDIM)) - &
             &0.5_PR*mpp*dvdr*dvdr*invdrmag*invdrmag
        utot = utot + mpp*sph(pp)%u
#if defined(SELF_GRAVITY)
        do jj=j,pp_pot
           ppp = pp_potlist(jj)
           call gravity_sph(invhp,sph(ppp)%h,sph(ppp)%m,rp(1:NDIM),&
                &sph(ppp)%r(1:NDIM),agravp(1:NDIM),dpotp)
           gpetot = gpetot + sph(ppp)%m*dpotp
        end do
#endif
     end do
     ! -----------------------------------------------------------------------

     div_v_p = div_v_p / rhotemp
     div_a_p = div_a_p / rhotemp

     ! Perform div_v, div_a and energy checks here if selected
     if (div_a_search .and. div_a_p > 0.0_PR) testflag(i) = .false.
     if (div_v_search .and. div_v_p > 0.0_PR) testflag(i) = .false.
#if defined(SELF_GRAVITY)
     if (thermal_search .and. utot > 0.5_PR*gpetot) testflag(i) = .false.
     if (thermal_search .and. utot + rotketot > gpetot) testflag(i) = .false.
#endif

     ! Now loop over all existing sinks
     ! -----------------------------------------------------------------------
     do s=1,stot
        rs(1:NDIM) = sink(s)%r(1:NDIM)
        call distance3(rp(1:NDIM),rs(1:NDIM),dr(1:NDIM),drsqd)
        dv(1:NDIM) = sink(s)%v(1:NDIM) - vp(1:NDIM)
        RKH = 1.0_PR/(sph(p)%rho*rhoscale*rhocgs/1.0E-13)/&
             &(1000.0_PR*yr/tscale/t_SI)
        taccrete = -0.25_PR*drsqd/dot_product(dv(1:NDIM),dr(1:NDIM))
        if (abs(RKH - div_v_p) > SMALL_NUMBER) then
           tcollapse = -0.5_PR*(RKH - div_v_p)/(RKH*div_v_p)
        else
           tcollapse = BIG_NUMBER
        end if
        
#if defined(DEBUG_SINK_SEARCH)
        write(6,*) "Timescale check : ",p,s,"   rho : ",sph(p)%rho*rhoscale*rhocgs,"g cm^-3"
        write(6,*) "RKH : ",RKH/tscale,"    div_v : ",div_v_p/tscale
        write(6,*) "tcollapse : ",tcollapse*tscale,trim(tunit)
        write(6,*) "taccrete  : ",taccrete*tscale,trim(tunit)
#endif

        ! 
        if (timescale_search .and. tcollapse > taccrete) testflag(i) = .false.
        
     end do
     ! -----------------------------------------------------------------------


#if defined(DEBUG_SINK_SEARCH)
!      write(6,*) "Smoothed quantities for : ",i,p,sph(p)%porig
!      write(6,*) "rho                     : ",sinktest(i)%rho*rhoscale*rhocgs,&
!           &sph(p)%rho*rhoscale*rhocgs,rhosink*rhoscale*rhocgs
!      write(6,*) "div_v                   : ",sinktest(i)%div_v
!      write(6,*) "div_a                   : ",sinktest(i)%div_a
!      write(6,*) "alpha        : ",utot/gpetot
!      write(6,*) "alpha + beta : ",(utot + ketot)/gpetot
#endif

  end do
! ============================================================================


! Record new shortened list of candidate sink particles
  iaux = 0
  do i=1,pcandidate
     if (testflag(i)) then
        iaux = iaux + 1
        rholist(iaux) = rholist(i)
        testflag(iaux) = testflag(i)
     end if
  end do
  pcandidate = iaux
#endif

! Find the densest particle that has passed all the above tests, if any 
! exist at all.
! ----------------------------------------------------------------------------
  do i=pcandidate,1,-1
     if (.not. testflag(i)) cycle

     ! If all tests passed, record id of candidate and exit loop
     p = rholist(i)
     psink = p
     exit
  end do

#if defined(USE_MPI)
! Find out which task has the highest density particle
  if (psink > 0) then
     call mpi_find_densest_particle(sph(psink)%rho,densest_task)
     if (densest_task /= rank) psink = -1
  else
     call mpi_find_densest_particle(0._PR,densest_task)
  end if
#endif

! Deallocate arrays
! ----------------------------------------------------------------------------
  if (allocated(pp_potlist)) deallocate(pp_potlist)
  if (allocated(testflag)) deallocate(testflag)
  if (allocated(rhovalues)) deallocate(rhovalues)
  if (allocated(rholist)) deallocate(rholist)


! Look at polytropic quantities to check sink doesn't overlap
! ----------------------------------------------------------------------------
  if (gamma > 1.33333333_PR .and. stot > 999999999) then

     thetafac = 100.0_PR
     Rpoly = sqrt(thetafac)*(abs(sph(psink)%gpot)*(gamma - 1.0_PR) &
          & - sqrt(gamma/(gamma - 1.0_PR))*sph(psink)%press/sph(psink)%rho)&
          & / sqrt(4.0_PR*PI*sph(psink)%press)
     Mpoly = (abs(sph(psink)%gpot) - gamma*sph(psink)%press/&
          &((gamma - 1.0_PR)*sph(psink)%rho))*Rpoly

     write(6,*) "Rpoly : ",Rpoly*rscale,trim(runit)
     write(6,*) "Mpoly : ",Mpoly*mscale,trim(munit)
     write(6,*) "Rpoly/rsink : ",Rpoly/2.0_PR/sph(psink)%h
     write(6,*) "Mpoly/mp    : ",Mpoly/58/sph(psink)%m

     ! Now compute orbital properties and Hill radius to all other sinks
     do s=1,stot
        call distance2(sink(s)%r(1:NDIM),psink,dr(1:NDIM),drsqd)
        dv(1:NDIM) = sph(psink)%v(1:NDIM) - sink(s)%v(1:NDIM)
        reduced_mass = Mpoly*sink(s)%m/(Mpoly + sink(s)%m)

        angmom(1) = dr(2)*dv(3) - dr(3)*dv(2)
        angmom(2) = dr(3)*dv(1) - dr(1)*dv(3)
        angmom(3) = dr(1)*dv(2) - dr(2)*dv(1)
        Lsqd = dot_product(angmom(1:3),angmom(1:3))

        Eorbit = 0.5_PR*dot_product(dv(1:NDIM),dv(1:NDIM)) - &
             & reduced_mass/sqrt(drsqd)

        if (Eorbit < 0.0_PR) then
           aorbit = -0.5_PR*reduced_mass/Eorbit
           ecc = sqrt(1.0_PR - Lsqd/reduced_mass/aorbit)
           rperi = (1.0_PR - ecc)*aorbit

           write(6,*) "Checking sink : ",s,Eorbit,reduced_mass
           write(6,*) "Dist   : ",sqrt(drsqd),sqrt(drsqd)/sink(s)%radius,Lsqd
           write(6,*) "ecc    : ",ecc
           write(6,*) "aorbit : ",aorbit*rscale,trim(runit)
           write(6,*) "rperi  : ",rperi*rscale,trim(runit)
           write(6,*) "rnew?  : ",NEW_SINK_RMAX*(Rpoly + sink(s)%radius),&
                &trim(runit)

           if (rperi < NEW_SINK_RMAX*(Rpoly + sink(s)%radius)) then
              write(6,*) "Potential object overlap detected : ",&
                   &sqrt(drsqd),Rpoly,sink(s)%radius
              psink = -1
              exit
           end if
        end if
     end do
     
     !do s=1,stot
     !   call distance2(sink(s)%r(1:NDIM),psink,dr(1:NDIM),drsqd)
     !   if (drsqd <= (NEW_SINK_RMAX*(Rpoly + sink(s)%radius))**2) then
     !      write(6,*) "Potential object overlap detected : ",&
     !           &sqrt(drsqd),Rpoly,sink(s)%radius
     !      psink = -1
     !      exit
     !   end if
     !end do

  end if


! If all conditions are met, then create sink particle 
! ----------------------------------------------------------------------------

  if (psink > 0) then
     if (stot >= SMAX) then
        stop 'Maximum number of sinks exceeded. Need to increase SMAX.'
     end if
  end if


  return
END SUBROUTINE sink_search
