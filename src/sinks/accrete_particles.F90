! ACCRETE_PARTICLES.F90
! D. A. Hubber - 17/6/2007
! Removes particles from computational domain that have been accreted by 
! sink particles.  First, we identify all particles that are within any sink 
! boundaries.  Next, we determine if any of these particles are bound to any 
! of the sinks.  If accreted, the mass, momentum and angular momentum of the 
! particle is added to the sink.  The particle is then removed from the 
! simulation.  Finally, the particles are reordered so all live particles 
! are contiguous in the main data arrays (i.e. no holes in memory with dead 
! particles).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE accrete_particles
  use particle_module
  use sink_module
  use type_module
  use kernel_module
  use scaling_module
  use time_module, only : nsteps,nstock
#if defined(USE_MPI)
  use mpi_communication_module
  use mpi
#endif
  implicit none

  integer :: i                           ! Auxilary counter
  integer :: ndead                       ! Total number of dead particles
  integer :: p                           ! Particle counter
  integer :: pp_tot                      ! Number of particles accreted by s
  integer :: s                           ! Sink counter
  integer :: ss                          ! Aux sink i
  integer, allocatable :: deadlist(:)    ! List of dead (accreted) particles
  integer, allocatable :: pp_templist(:) ! Temp. list for particle ids
  integer, allocatable :: sinkflag(:)    ! Flag sink-particle overlap
  real(kind=DP) :: acom_temp(1:VDIM)     ! Accel. of particle-sink system
  real(kind=DP) :: angmomsink(1:3)       ! Angular momentum of sink s
  real(kind=DP) :: ap(1:VDIM)            ! Acceleration of particle p
  real(kind=DP) :: as(1:VDIM)            ! Acceleration of sink s
  real(kind=PR) :: atemp(1:NDIM)         ! Aux. acceleration vector
  real(kind=PR) :: dpotp                 ! Aux. grav potential variable
  real(kind=DP) :: dr(1:NDIM)            ! Relative position vector
  real(kind=DP) :: dv(1:VDIM)            ! Relative velocity vector
  real(kind=DP) :: drsqd                 ! Distance squared
  real(kind=DP) :: drsqd2                ! Second distance squared variable
  real(kind=PR) :: invhp                 ! 1.0 / hp
  real(kind=PR) :: invhs                 ! 1.0 / hsink
  real(kind=DP) :: mp                    ! Mass of particle p
  real(kind=DP) :: ms                    ! Mass of sink s
  real(kind=DP) :: mtot_temp             ! Mass of particle-sink system
  real(kind=DP) :: rads                  ! Radius of sink particle s
  real(kind=DP) :: rcom_temp(1:NDIM)     ! COM of particle-sink system
  real(kind=PR) :: reduced_mass          ! Reduced mass of 2-body system
  real(kind=PR) :: rp(1:NDIM)            ! Position of particle p
  real(kind=DP) :: rp_dp(1:NDIM)         ! Position of particle p (DP precision)
  real(kind=PR) :: rs(1:NDIM)            ! Position of sink s
  real(kind=DP) :: rs_dp(1:NDIM)         ! Position of sink s (DP precision)
  real(kind=DP) :: total_energy          ! Total energy of sink-particle system
  real(kind=DP) :: vcom_temp(1:VDIM)     ! COM velocity of particle-sink system
  real(kind=DP) :: vp(1:VDIM)            ! Velocity of particle p
  real(kind=DP) :: vs(1:VDIM)            ! Velocity of sink s
#if defined(ACCRETION_RATE)
  real(kind=DP) :: maccreted             ! Total mass of accreted particles
#endif
#if defined(SINK_REMOVE_ANGMOM)
  real(kind=DP) :: cmean                 ! Mean sound speed of accreted ptcls
#endif
#if defined(DEBUG_ACCRETE)
  logical :: accreteflag                 ! Flag if any particles are accreted
#endif
#if defined(USE_MPI)
  integer :: ierr                        ! MPI error value
#endif

  debug2("Accreting SPH particles to sinks [accrete_particles.F90]") 
  debug_timing("ACCRETE_PARTICLES")

! Initialize sink accretion arrays
  allocate(deadlist(1:ptot))
  allocate(sinkflag(1:ptot))
  allocate(pp_templist(1:ptot))
  sinkflag(1:ptot) = -1
  ndead = 0

! Find id of sink to which particle is closest and within the sink radius
! ----------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dr,drsqd,drsqd2,rads,rs_dp,s,ss)
  do p=pgasstart,pgasend
     do s=1,stot
        if (.NOT. sink(s)%accrete) cycle
        
        rs_dp(1:NDIM) = real(sink(s)%r(1:NDIM),DP)
        rads = real(sink(s)%radius,DP)
        call distance2_dp(rs_dp(1:NDIM),p,dr(1:NDIM),drsqd)

        if (drsqd < rads*rads) then
           if (sinkflag(p) == -1) then
              sinkflag(p) = s
           else
              ss = sinkflag(p)
              rs_dp(1:NDIM) = sink(ss)%r(1:NDIM)
              call distance2_dp(rs_dp(1:NDIM),p,dr(1:NDIM),drsqd2)
              if (drsqd < drsqd2) sinkflag(p) = s
           end if
#if defined(DEBUG_ACCRETE)
           write(6,*) "Check accretion : ",p,sinkflag(p),&
                &sqrt(drsqd)*rscale,rads*rscale
#endif
        end if

     end do
  end do
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------------

! Write some information to screen for debugging purposes
#if defined(DEBUG_ACCRETE)
  do p=1,ptot
     if (sinkflag(p) /= -1) ndead = ndead + 1
  end do
  if (ndead > 0) then
     write(6,*) "ptot (before accretion) : ",ptot
     write(6,*) "sink candidates : ",ndead
     accreteflag = .true.
#if !defined(USE_MPI)
     call diagnostics ! Can't call diagnostics for single task with MPI
#endif
  else
     accreteflag = .false.
  end if
  ndead = 0
#endif


! If a particle is within the sink radius of and bound to a sink,
! accrete particle to the sink.
! ============================================================================
  do s=1,stot
#if defined(USE_MPI)
     if (sink(s)%domain /= rank) cycle
#endif

#if defined(DEBUG_ACCRETE)
     write(6,*) "Accreting sink? : ",s,sink(s)%accrete
#endif

     ! Skip sink if it's not flagged as accreting
     if (.NOT. sink(s)%accrete) cycle

     ! Create local copies of old sink properties
     rs(1:NDIM) = sink(s)%r(1:NDIM)
     rs_dp(1:NDIM) = real(rs(1:NDIM),DP)
     vs(1:VDIM) = real(sink(s)%v(1:VDIM),DP)
     as(1:VDIM) = real(sink(s)%a(1:VDIM),DP)
     ms    = real(sink(s)%m,DP)
     rads  = real(sink(s)%radius,DP)
     invhs = 1.0_PR / sink(s)%h

     ! Temporary array to store new sink properties
     rcom_temp(1:NDIM) = ms*rs_dp(1:NDIM)
     vcom_temp(1:VDIM) = ms*vs(1:VDIM)
     acom_temp(1:VDIM) = ms*as(1:VDIM)
     mtot_temp = ms
#if defined(DEBUG_FORCES)
     sink(s)%agrav(1:NDIM) = real(ms,PR)*sink(s)%agrav(1:NDIM)
     sink(s)%ahydro(1:NDIM) = real(ms,PR)*sink(s)%ahydro(1:NDIM)
#endif
#if defined(ACCRETION_RATE)
     maccreted = 0.0_DP
#endif
     pp_tot = 0


     ! Search through all gravitating particles and find those bound to sink s
     ! -----------------------------------------------------------------------
     do p=pgasstart,pgasend

        ! If the particle is within radius of sink s (and only sink s), then
        ! check if particle p is bound to sink s
        if (sinkflag(p) /= s) cycle

        ! Make local copy of particle quantities
        invhp      = 1.0_PR / sph(p)%h
        mp         = real(sph(p)%m,DP)
        rp(1:NDIM) = sph(p)%r(1:NDIM)
        ap(1:VDIM) = real(sph(p)%a(1:VDIM),DP)
        vp(1:VDIM) = real(sph(p)%v(1:VDIM),DP)
        dv(1:VDIM) = vs(1:VDIM) - vp(1:VDIM)
        reduced_mass = ms*mp/(ms + mp)

        ! Only perform energy check if option chosen in params file
        if (energy_accrete) then
#if defined(N_BODY)
           call gravity_nbody(real(ms,PR),rp(1:NDIM),&
                &rs(1:NDIM),atemp(1:NDIM),dpotp)
#elif defined(GRAD_H_SPH)
           call gravity_gradh(invhp,invhs,real(ms,PR),&
                &rp(1:NDIM),rs(1:NDIM),&
                &0.0_PR,0.0_PR,atemp(1:NDIM),dpotp)
#else
           call gravity_sph(invhp,invhp,real(ms,PR),&
                &rp(1:NDIM),rs(1:NDIM),atemp(1:NDIM),dpotp)
#endif
           total_energy = 0.5_DP*reduced_mass*&
                &dot_product(dv(1:NDIM),dv(1:NDIM)) - real(dpotp,DP)*mp
        else
           total_energy = -BIG_NUMBER_DP
        end if

#if defined(DEBUG_ACCRETE)
        write(6,*) "candidate : ",p,total_energy
#endif

        ! If bound, add particle properties to the sink and then add 
        ! particle id to list of 'dead' (i.e. accreted) particles.
        if (total_energy < 0.0_DP) then
           mtot_temp = mtot_temp + mp
           rcom_temp(1:NDIM) = rcom_temp(1:NDIM) + mp*real(rp(1:NDIM),DP)
           vcom_temp(1:VDIM) = vcom_temp(1:VDIM) + mp*vp(1:VDIM)
           acom_temp(1:VDIM) = acom_temp(1:VDIM) + mp*ap(1:VDIM)
#if defined(DEBUG_FORCES)
           sink(s)%agrav(1:NDIM) = sink(s)%agrav(1:NDIM) + &
                &real(mp,PR)*sph(p)%a_grav(1:NDIM)
           sink(s)%ahydro(1:NDIM) = sink(s)%ahydro(1:NDIM) + &
                &real(mp,PR)*sph(p)%a_hydro(1:NDIM)
#endif
           ndead  = ndead + 1
           pp_tot = pp_tot + 1
           deadlist(ndead) = p
#if defined(ACCRETION_RATE)
           maccreted = maccreted + mp
#endif
        end if
     end do

     ! Store new sink properties if any particles have been accreted
     ! -----------------------------------------------------------------------
     if (pp_tot > 0) then

        ! First calculate centres of mass, velocity and acceleration
        rcom_temp(1:NDIM) = rcom_temp(1:NDIM) / mtot_temp
        vcom_temp(1:VDIM) = vcom_temp(1:VDIM) / mtot_temp
        acom_temp(1:VDIM) = acom_temp(1:VDIM) / mtot_temp

        ! Angular momentum of old COM around new COM
        angmomsink(1:3) = 0.0_DP
        call distance3_dp(rcom_temp(1:NDIM),rs_dp(1:NDIM),dr(1:NDIM),drsqd)
        dv(1:VDIM) = vs(1:VDIM) - vcom_temp(1:VDIM)
#if NDIM==3
        angmomsink(1) = angmomsink(1) + ms*dr(2)*dv(3) - ms*dr(3)*dv(2)
        angmomsink(2) = angmomsink(2) + ms*dr(3)*dv(1) - ms*dr(1)*dv(3)
#endif
        angmomsink(3) = angmomsink(3) + ms*dr(1)*dv(2) - ms*dr(2)*dv(1)

        ! Now add angular momentum of accreted particles to new COM
        if (ndead > 0) then
           do i=1,ndead
              p = deadlist(i)
              
              ! Ensure we only include particles accreted by this sink
              if (sinkflag(p) /= s) cycle
              mp            = real(sph(p)%m,DP)
              rp_dp(1:NDIM) = real(sph(p)%r(1:NDIM),DP)
              dv(1:VDIM)    = vcom_temp(1:VDIM) - real(sph(p)%v(1:VDIM),DP)
              call distance3_dp(rp_dp(1:NDIM),rcom_temp(1:NDIM),&
                               &dr(1:NDIM),drsqd)
#if NDIM==3
              angmomsink(1) = angmomsink(1) + mp*dr(2)*dv(3) - mp*dr(3)*dv(2)
              angmomsink(2) = angmomsink(2) + mp*dr(3)*dv(1) - mp*dr(1)*dv(3)
#endif
              angmomsink(3) = angmomsink(3) + mp*dr(1)*dv(2) - mp*dr(2)*dv(1)
           end do
        end if


        ! Now store all sink properties in main data structures
        sink(s)%m              = real(mtot_temp,PR)
        sink(s)%r(1:NDIM)      = real(rcom_temp(1:NDIM),PR)
        sink(s)%v(1:VDIM)      = real(vcom_temp(1:VDIM),PR)
        sink(s)%a(1:VDIM)      = real(acom_temp(1:VDIM),PR)
        sink(s)%rold(1:NDIM)   = real(rcom_temp(1:NDIM),PR)
        sink(s)%vold(1:VDIM)   = real(vcom_temp(1:VDIM),PR)
        sink(s)%angmom(1:3)    = sink(s)%angmom(1:3) + &
             &real(angmomsink(1:3),PR)
        sink(s)%angmomnet(1:3) = sink(s)%angmomnet(1:3) + &
             &real(angmomsink(1:3),PR)
#if defined(LEAPFROG_KDK)
        sink(s)%aold(1:VDIM)   = real(acom_temp(1:VDIM),PR)
#endif
#if defined(DEBUG_FORCES)
        sink(s)%agrav(1:NDIM)  = sink(s)%agrav(1:NDIM) / real(mtot_temp,PR)
        sink(s)%ahydro(1:NDIM) = sink(s)%ahydro(1:NDIM) / real(mtot_temp,PR)
#endif

     end if
     ! -----------------------------------------------------------------------

     ! Calculate accretion rate and protostar properties if required
#if defined(ACCRETION_RATE)
     call sink_accretion_properties(s,maccreted)
#endif

  end do
! ============================================================================


! Record accretion history for each sink
! ----------------------------------------------------------------------------
  do s=1,stot
#if defined(USE_MPI)
     if (sink(s)%domain /= rank) cycle
#endif
     pp_tot = 0
     if (ndead > 0) then
        do i=1,ndead
           p = deadlist(i)
           if (sinkflag(p) /= s) cycle
           pp_tot = pp_tot + 1
           pp_templist(pp_tot) = p
        end do
        if (pp_tot > 0) then
           call write_accreted_particles(s,pp_tot,pp_templist(1:pp_tot))
        end if
     end if
  end do

! If particles have been accreted, need to reorder arrays 
  if (ndead > 0) then
     call insertion_sort_int(ndead,deadlist(1:ndead))
     call remove_from_list(ndead,deadlist(1:ndead))
  end if

! Redistribute angular momentum contained in sinks to neighbouring particles
#if defined(SINK_REMOVE_ANGMOM)
  call redistribute_sink_angmom
#endif

! Free memory
  deallocate(pp_templist)
  deallocate(sinkflag)
  deallocate(deadlist)

! Re-stock tree immediately in next step
  nstock = nsteps

#if defined(DEBUG_ACCRETE)
  if (accreteflag) then
     write(6,*) "ptot (after accretion) : ",ptot
#if !defined(USE_MPI)
     call diagnostics ! cannot call diagnostics with single MPI task
#endif
  end if
#endif

! Calculate fraction of gas accreted by sinks in case of switching 
! to the N-body integrator
#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  mgas = 0.0_DP
  do p=pgasstart,pgasend
     mgas = mgas + real(sph(p)%m,DP)
  end do
#if defined(USE_MPI)
  WAIT_TIME_MACRO
  call MPI_ALLREDUCE(MPI_IN_PLACE, mgas, 1, &
     & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  CALC_TIME_MACRO
#endif
  sink_frac = (mgas_orig - mgas) / mgas_orig
#endif

  return
END SUBROUTINE accrete_particles
