! BINARY_SEARCH.F90
! D. A. Hubber - 28/01/2008
! Identifies bound binary systems from the ensemble of sink particles 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE binary_search
  use interface_module, only : binary_energy,binary_properties,&
       &distance3_dp,gravity_hermite4,gravity_hermite4_meanh,insertion_sort_dp
  use definitions
  use sink_module
  use Nbody_module
  use scaling_module
  use constant_module
  use time_module, only : n,nsteps,time
  use filename_module, only : run_dir,run_id
  implicit none

  character(len=256) :: filename1             ! Output filename
  character(len=256) :: filename2             ! Output filename
  character(len=256) :: filename3             ! Output filename
  logical :: ex                               ! Does file exist?
  integer :: i                                ! Auxilary loop counter
  integer :: nsingles                         ! Number of single stars
  integer :: s                                ! Sink counter
  integer :: ss                               ! Secondary sink counter
  integer :: binflag(1:2*SMAX)                ! Status flag of stars
  integer :: most_bound(1:2*SMAX,1:2*SMAX)    ! List of most bound members
  integer :: sinkid(1:2*SMAX)                 ! Sink id list
  real(kind=DP) :: binen                      ! Two-body energy
  real(kind=DP) :: energy(1:2*SMAX,1:2*SMAX)  ! Energy of all system pairs
  real(kind=DP) :: energy_s(1:2*SMAX)         ! Energy of pairs for system s

  real(kind=DP) :: atemp(1:NDIM)              ! Grav. accel of star s
  real(kind=DP) :: adot(1:NDIM)               ! Jerk of star s
  real(kind=DP) :: dpotp                      ! Potential of star s
  real(kind=DP) :: da(1:NDIM)                 ! ..
  real(kind=DP) :: apert(1:NDIM)              ! ..
  real(kind=DP) :: gammapert                  ! ..


! Immediately return if there are not enough sinks/stars
  if (stot <= 1) return

  debug2("Search for bound multiple systems [binary_search.F90]")

! Initialize variables and arrays
  nbin = 0
  nsingles = stot
  do s=1,2*SMAX
     binflag(s) = 0
     do ss=1,2*SMAX
        energy(s,ss) = BIG_NUMBER_DP
        most_bound(s,ss) = -1
     end do
  end do

#if defined(DEBUG_BINARY_SEARCH)
  write(6,*) "Finished initialising variables [binary_search.F90]"
#endif


! Calculate the energy of all sink pairs
! ----------------------------------------------------------------------------
  do s=1,stot
     do ss=1,stot
        if (s == ss) cycle
        call binary_energy(s,ss,star(s)%h,star(ss)%h,star(s)%m,star(ss)%m,&
             &star(s)%r(1:NDIM),star(ss)%r(1:NDIM),star(s)%v(1:NDIM),&
             &star(ss)%v(1:NDIM),binen)
        energy(s,ss) = binen
     end do
  end do

#if defined(DEBUG_BINARY_SEARCH)
  write(6,*) "Sort energies into ascending order [binary_search.F90]"
#endif


! Sort energies into ascending order and store for later analysis
! ----------------------------------------------------------------------------
  do s=1,stot
     i = 0
     do ss=1,stot
        if (s == ss) cycle        
        i = i + 1
        energy_s(i) = energy(s,ss)
        sinkid(i) = ss
     end do

     ! Arrange in ascending order and store in main array
     call insertion_sort_dp(i,sinkid(1:i),energy_s(1:i))
     most_bound(s,1:i) = sinkid(1:i)
  end do


! Identify all mutually bound pairs of stars as binaries and calculate 
! and store binary properties
! ----------------------------------------------------------------------------
  do s=1,stot-1
     do ss=s+1,stot
        if (most_bound(s,1) == ss .and.  most_bound(ss,1) == s .and. &
             &energy(s,ss) < 0.0_DP .and. energy(ss,s) < 0.0_DP) then

           ! Now calculate perturbing force on binary due to other 
           ! nearby stars
#if defined(NBODY_HERMITE4) && defined(MEANH_GRAVITY)
           call gravity_hermite4_meanh(0.5_DP*(star(s)%h + star(ss)%h),&
                &star(ss)%m,star(s)%r(1:NDIM),star(ss)%r(1:NDIM),&
                &star(s)%v(1:NDIM),star(ss)%v(1:NDIM),&
                &atemp(1:NDIM),adot(1:NDIM),dpotp)
#elif defined(NBODY_HERMITE4)
           call gravity_hermite4(1.0_DP/star(s)%h,star(ss)%h,star(ss)%m,&
                &star(s)%r(1:NDIM),star(ss)%r(1:NDIM),star(s)%v(1:NDIM),&
                &star(ss)%v(1:NDIM),atemp(1:NDIM),adot(1:NDIM),dpotp)
#endif
           atemp(1:NDIM) = (star(s)%m + star(ss)%m)*atemp(1:NDIM)/star(ss)%m
           da(1:NDIM) = star(s)%a(1:NDIM) - &
                & star(ss)%a(1:NDIM)
           apert(1:NDIM) = da(1:NDIM) - atemp(1:NDIM)

           gammapert = sqrt(dot_product(apert(1:NDIM),apert(1:NDIM)) / &
                dot_product(atemp(1:NDIM),atemp(1:NDIM)))


           if (gammapert < 1.E-03) then
              write(6,*) "Checking perturbation : ",nbin,s,ss,gammapert,&
                   &dot_product(da(1:NDIM),da(1:NDIM)),&
                   &dot_product(apert(1:NDIM),apert(1:NDIM)),&
                   &dot_product(atemp(1:NDIM),atemp(1:NDIM))
           end if


           ! If perturbation is too large, ignore this pair as a 
           ! potential binary star
           if (gammapert > gammapertmax) cycle


           ! Adjust binary/single counters, and flag stars as part of binary
           nbin        = nbin + 1
           nsingles    = nsingles - 2
           binflag(s)  = stot + nbin
           binflag(ss) = stot + nbin
           
           ! Calculate all properties of binary and store in memory
           call binary_properties(s,ss,energy(s,ss),sink(s)%m,&
                &sink(ss)%m,sink(s)%r(1:NDIM),sink(ss)%r(1:NDIM),&
                &sink(s)%v(1:NDIM),sink(ss)%v(1:NDIM))

        end if
     end do
  end do

#if defined(DEBUG_BINARY_SEARCH)
  write(6,*) "Search for multiple systerms [binary_search.F90]"
#endif



! Write properties to file
! ============================================================================
  if (nbin > 0) then
     filename1 = trim(adjustl(run_dir))//trim(adjustl(run_id))//".binstats"
     filename2 = trim(adjustl(run_dir))//trim(adjustl(run_id))//".finbinstats"
     inquire(file=filename1,exist=ex)
     if (.not. ex) then
        filename3 = trim(adjustl(run_dir))//trim(adjustl(run_id))&
             &//".inibinstats"
        open(unit=3,file=filename3,status="unknown")
     end if
     open(unit=1,file=filename1,status="unknown",position="append")
     open(unit=2,file=filename2,status="unknown")
     
#if NDIM==2
     ! -----------------------------------------------------------------------
10   format(2i10,2X,E12.4,2X,3i10,11E12.4)
     
     do s=1,nbin
        write(1,10) n,nsteps,time*tscale,&
             &binary(s)%id, &
             &binary(s)%s1, &
             &binary(s)%s2, &
             &binary(s)%r(1)*rscale, &
             &binary(s)%r(2)*rscale, &
             &binary(s)%v(1)*vscale, &
             &binary(s)%v(2)*vscale, &
             &binary(s)%m*mscale, &
             &binary(s)%angmom(3)*angmomscale, &
             &binary(s)%ecc, &
             &binary(s)%period*tscale, &
             &binary(s)%q, &
             &binary(s)%sma*rscale, &
             &binary(s)%drmag*rscale
        write(2,10) n,nsteps,time*tscale,&
             &binary(s)%id, &
             &binary(s)%s1, &
             &binary(s)%s2, &
             &binary(s)%r(1)*rscale, &
             &binary(s)%r(2)*rscale, &
             &binary(s)%v(1)*vscale, &
             &binary(s)%v(2)*vscale, &
             &binary(s)%m*mscale, &
             &binary(s)%angmom(3)*angmomscale, &
             &binary(s)%ecc, &
             &binary(s)%period*tscale, &
             &binary(s)%q, &
             &binary(s)%sma*rscale, &
             &binary(s)%drmag*rscale
        if (.not. ex) then
           write(2,10) n,nsteps,time*tscale,&
                &binary(s)%id, &
                &binary(s)%s1, &
                &binary(s)%s2, &
                &binary(s)%r(1)*rscale, &
                &binary(s)%r(2)*rscale, &
                &binary(s)%v(1)*vscale, &
                &binary(s)%v(2)*vscale, &
                &binary(s)%m*mscale, &
                &binary(s)%angmom(3)*angmomscale, &
                &binary(s)%ecc, &
                &binary(s)%period*tscale, &
                &binary(s)%q, &
                &binary(s)%sma*rscale, &
                &binary(s)%drmag*rscale
        end if
     end do

#elif NDIM==3
     ! -----------------------------------------------------------------------
10   format(2i10,2X,E12.4,2X,3i10,15E12.4)
     
     do s=1,nbin
        write(1,10) n,nsteps,time*tscale,&
             &binary(s)%id, &
             &binary(s)%s1, &
             &binary(s)%s2, &
             &binary(s)%r(1)*rscale, &
             &binary(s)%r(2)*rscale, &
             &binary(s)%r(3)*rscale, &
             &binary(s)%v(1)*vscale, &
             &binary(s)%v(2)*vscale, &
             &binary(s)%v(3)*vscale, &
             &binary(s)%m*mscale, &
             &binary(s)%angmom(1)*angmomscale, &
             &binary(s)%angmom(2)*angmomscale, &
             &binary(s)%angmom(3)*angmomscale, &
             &binary(s)%ecc, &
             &binary(s)%period*tscale, &
             &binary(s)%q, &
             &binary(s)%sma*rscale, &
             &binary(s)%drmag*rscale
        write(2,10) n,nsteps,time*tscale,&
             &binary(s)%id, &
             &binary(s)%s1, &
             &binary(s)%s2, &
             &binary(s)%r(1)*rscale, &
             &binary(s)%r(2)*rscale, &
             &binary(s)%r(3)*rscale, &
             &binary(s)%v(1)*vscale, &
             &binary(s)%v(2)*vscale, &
             &binary(s)%v(3)*vscale, &
             &binary(s)%m*mscale, &
             &binary(s)%angmom(1)*angmomscale, &
             &binary(s)%angmom(2)*angmomscale, &
             &binary(s)%angmom(3)*angmomscale, &
             &binary(s)%ecc, &
             &binary(s)%period*tscale, &
             &binary(s)%q, &
             &binary(s)%sma*rscale, &
             &binary(s)%drmag*rscale
        if (.not. ex) then
           write(2,10) n,nsteps,time*tscale,&
                &binary(s)%id, &
                &binary(s)%s1, &
                &binary(s)%s2, &
                &binary(s)%r(1)*rscale, &
                &binary(s)%r(2)*rscale, &
                &binary(s)%r(3)*rscale, &
                &binary(s)%v(1)*vscale, &
                &binary(s)%v(2)*vscale, &
                &binary(s)%v(3)*vscale, &
                &binary(s)%m*mscale, &
                &binary(s)%angmom(1)*angmomscale, &
                &binary(s)%angmom(2)*angmomscale, &
                &binary(s)%angmom(3)*angmomscale, &
                &binary(s)%ecc, &
                &binary(s)%period*tscale, &
                &binary(s)%q, &
                &binary(s)%sma*rscale, &
                &binary(s)%drmag*rscale
        end if
     end do
#endif
     ! -----------------------------------------------------------------------

     ! Close all files for binary properties
     if (.not. ex) close(3)
     close(2)
     close(1)

  end if
! ============================================================================

  return
END SUBROUTINE binary_search
