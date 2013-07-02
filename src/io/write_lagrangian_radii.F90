! WRITE_LAGRANGIAN_RADII.F90
! C. P. Batty & D. A. Hubber - 12/12/2006
! Writes data to ASCII file for debugging purposes. 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_lagrangian_radii(rorigin)
  use interface_module, only : distance2,distance3
  use particle_module
  use hydro_module
  use scaling_module
  use filename_module
  use time_module
  use type_module
  use sink_module
  implicit none

  real(kind=PR), intent(inout) :: rorigin(1:NDIM)  ! Origin for output

  character(len=256) :: debug_file                 ! Debug filename
  character(len=12) :: fs = '(11E15.7)'                ! ..
  integer, parameter :: nradii=10                 ! ..
  integer :: i                                     ! .. 
  integer :: p                                     ! ..
  integer :: s                                     ! ..
  integer :: iaux                                  ! ..
  integer, allocatable :: rad_ids(:)               ! ..
  real(kind=PR) :: dr(1:NDIM)                      ! ..
  real(kind=PR) :: drsqd                           ! ..
  real(kind=PR) :: mtot_aux                        ! ..
  real(kind=DP) :: mtot_temp                       ! ..
  real(kind=DP) :: rmass(1:nradii)                 ! ..
  real(kind=PR), allocatable:: radius(:)           ! ..

  debug2("[write_lagrangian_radii.F90]")

  allocate(rad_ids(1:ptot+stot))
  allocate(radius(1:ptot+stot))

! Calculate Lagrangian radii for SPH particles only
! -----------------------------------------------------------------------------
#if defined(SPH_SIMULATION) || defined(NBODY_SPH_SIMULATION)
  if (ptot > 0) then
     debug2("Computing Lagrangian radii for gas particles [write_lagrangian_radii.F90]")
     mtot_temp = 0.0_PR
     do p=1,ptot
        mtot_temp = mtot_temp + sph(p)%m
        call distance2(rorigin(1:NDIM),p,dr(1:NDIM),drsqd)
        rad_ids(p) = p
        radius(p)  = sqrt(drsqd) + SMALL_NUMBER
     end do
     
     ! Order particles in order of increasing radii
     call insertion_sort_real(ptot,rad_ids(1:ptot),radius(1:ptot))
     
     mtot_aux = 0.0_DP
     rmass = 0.0_PR
     iaux = 1
     do i=1,ptot
        p = rad_ids(i)
        mtot_aux = mtot_aux + real(sph(p)%m,DP)
        if (mtot_aux >= 0.9999999999_DP*mtot_temp*&
             &real(iaux,DP)/real(nradii,DP)) then
           rmass(iaux) = radius(i)
           iaux = iaux + 1
        end if
     end do
     
     debug_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".gasradius"
     open(1,file=debug_file,status="unknown",&
          &form="formatted",position="append")
     write(1,fs) time, rmass(1:nradii)
     !do i=1,10
     !   write(1,*) 0.1_PR*real(i,PR)*mtot_temp,rmass(i)
     !end do
     close(1)
  end if
#endif
! -----------------------------------------------------------------------------


! Calculate Lagrangian radii for sinks/stars only
! -----------------------------------------------------------------------------
#if defined(SINKS)
  if (stot > 0) then
     debug2("Computing Lagrangian radii for stars [write_lagrangian_radii.F90]")
     mtot_temp = 0.0_PR
     do s=1,stot
        mtot_temp = mtot_temp + sink(s)%m
        call distance3(rorigin(1:NDIM),sink(s)%r(1:NDIM),dr(1:NDIM),drsqd)
        rad_ids(s) = s
        radius(s)  = sqrt(drsqd) + SMALL_NUMBER
     end do
     
     ! Order particles in order of increasing radii
     call insertion_sort_real(stot,rad_ids(1:stot),radius(1:stot))
     
     mtot_aux = 0.0_DP
     rmass = 0.0_PR
     iaux = 1
     do i=1,stot
        s = rad_ids(i)
        mtot_aux = mtot_aux + real(sink(s)%m,DP)
        if (mtot_aux >= 0.99999999999_DP*mtot_temp*real(iaux,DP)/real(nradii,DP)) then
           rmass(iaux) = radius(i)
           iaux = iaux + 1
        end if
     end do
     
     debug_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".starradius"
     open(1,file=debug_file,status="unknown",&
          &form="formatted",position="append")
     write(1,'(101E15.7)') time, rmass(1:nradii)
     !do i=1,10
     !   write(1,*) 0.1_PR*real(i,PR)*mtot_temp,rmass(i)
     !end do
     close(1)
  end if
#endif
! -----------------------------------------------------------------------------


! Calculate Lagrangian radii for complete distribution
! -----------------------------------------------------------------------------
#if defined(SPH_SIMULATION) || defined(NBODY_SPH_SIMULATION)
#if defined(SINKS)
  if (ptot > 0 .and. stot > 0) then
     debug2("Computing Lagrangian radii for all particles [write_lagrangian_radii.F90]")
     mtot_temp = 0.0_PR
     do p=1,ptot
        mtot_temp = mtot_temp + sph(p)%m
        call distance2(rorigin(1:NDIM),p,dr(1:NDIM),drsqd)
        rad_ids(p) = p
        radius(p)  = sqrt(drsqd) + SMALL_NUMBER
     end do
     do s=1,stot
        mtot_temp = mtot_temp + sink(s)%m
        call distance3(rorigin(1:NDIM),sink(s)%r(1:NDIM),dr(1:NDIM),drsqd)
        rad_ids(s+ptot) = ptot + s
        radius(s+ptot)  = sqrt(drsqd) + SMALL_NUMBER
     end do
     
     ! Order particles in order of increasing radii
     call insertion_sort_real(ptot+stot,&
          &rad_ids(1:ptot+stot),radius(1:ptot+stot))
     
     mtot_aux = 0.0_DP
     rmass = 0.0_PR
     iaux = 1
     do i=1,ptot+stot
        p = rad_ids(i)
        if (p <= ptot) mtot_aux = mtot_aux + real(sph(p)%m,DP)
        if (p > ptot) mtot_aux = mtot_aux + real(sink(p-ptot)%m,DP)
        if (mtot_aux >= 0.999999999_DP*mtot_temp*real(iaux,DP)/real(nradii,DP)) then
           rmass(iaux) = radius(i)
           iaux = iaux + 1
        end if
     end do
     
     debug_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".allradius"
     open(1,file=debug_file,status="unknown",&
          &form="formatted",position="append")
     write(1,'(101E15.7)') time, rmass(1:nradii)
     !do i=1,10
     !   write(1,*) 0.1_PR*real(i,PR)*mtot_temp,rmass(i)
     !end do
     close(1)
  end if
#endif
#endif
! -----------------------------------------------------------------------------


  deallocate(radius)
  deallocate(rad_ids)


#if defined(DEBUG_FREEFALL) || defined(FREEFALL_TEST)
  debug_file = "freefall.dat"
  open(1,file=debug_file,status="unknown",form="formatted")
  do i=1,500
     drsqd = min(0.01_PR*real(i,PR),0.99999999999_PR)
     write(1,'(2E15.7)') 1.0_PR - 2.0_PR*INVPI*(asin(sqrt(drsqd)) - &
          &sqrt(drsqd)*sqrt(1 - drsqd)),drsqd
  end do
  close(1)
#endif


  return
END SUBROUTINE write_lagrangian_radii
