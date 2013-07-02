! WRITE_DATA_DEBUG.F90
! C. P. Batty & D. A. Hubber - 12/12/2006
! Writes data to ASCII file for debugging purposes. 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_data_debug(out_file,rorigin)
  use particle_module
  use hydro_module
  use neighbour_module
  use scaling_module
  use filename_module
  use time_module
  use type_module
  use sink_module
#if defined(RAD_WS)
  use Eos_module
#endif
  implicit none

  character(len=*), intent(in) :: out_file         ! Name of debug file
  real(kind=PR), intent(inout) :: rorigin(1:NDIM)  ! Origin for output

  integer :: k                            ! Dimension counter
  integer :: nelements                    ! No. of elements in alldata
  integer :: p                            ! Counter to loop over particles
  real(kind=PR) :: alldata(1:100)         ! Array with all output data
#if defined(SINKS)
  integer :: s                            ! Sink counter
#endif
#if defined(SINKS) || defined(GHOST_PARTICLES) || defined(DEBUG_FREEFALL) || defined(FREEFALL_TEST) || defined(PLUMMER_TEST)
  character(len=256) :: debug_file        ! Debug filename
#endif
#if defined(RAD_WS) && defined(DEBUG_RAD)
  integer :: prhomax                      ! id of particle with max. density
  real(kind=PR) :: rhomax                 ! max. density of all particles
#endif
#if defined(DEBUG_FREEFALL) || defined(FREEFALL_TEST) || defined(PLUMMER_TEST)
  integer :: i                            ! ..
  integer :: iaux                         ! ..
  integer, allocatable :: rad_ids(:)      ! ..
  real(kind=PR) :: dr(1:NDIM)             ! ..
  real(kind=PR) :: drsqd                  ! ..
  real(kind=DP) :: mtot_temp              ! ..
  real(kind=DP) :: rmass(1:100)           ! ..
  real(kind=PR), allocatable:: radius(:)  ! ..
#endif

  debug2("Writing data to debug file [write_data_debug.F90]")

  open(1, file=out_file, status="unknown", form="formatted")
#if defined(DEBUG1)
  write(6,*) "Debugging output file :",trim(out_file),"   (debug file)"
#endif

! Find id and position of densest particle
#if defined(RAD_WS) && defined(DEBUG_RAD)
  prhomax = 1
  rhomax = SMALL_NUMBER
  do p=pgravitystart,ptot
     if (sph(p)%rho > rhomax) then
        prhomax = p
        rhomax = sph(p)%rho
     end if
  end do
  rorigin(1:NDIM) = sph(prhomax)%r(1:NDIM)
#endif

! Loop over all particles, and scan through all possible arrays that 
! can go in alldata.
! ----------------------------------------------------------------------------
  do p=1,ptot

     ! Record all data for particle p in array 'alldata'
     call record_particle_data(p,nelements,alldata,rorigin)

     ! Now write all values to the file
#if defined(USE_MPI)
     write(1,'(I10,X,100E15.7)') sph(p)%porig, (alldata(k),k=1,nelements)
#else
     write(1,'(100E15.7)') (alldata(k),k=1,nelements)
#endif

  end do
! ----------------------------------------------------------------------------

  close(1)


! Write ghost particle information if required
! ----------------------------------------------------------------------------
#if defined(GHOST_PARTICLES) && defined(PERIODIC)
  if (pperiodic > 0) then
     debug_file = trim(adjustl(out_file))//".ghosts"
     open(1, file=debug_file, status="unknown", form="formatted")
     write(6,*) "Debugging output file: ",trim(debug_file)

     do p=ptot+1,ptot+pperiodic
        
        ! Record all data for particle p in array 'alldata'
        call record_particle_data(p,nelements,alldata,rorigin)
        
        ! Now write all values to the file
#if defined(USE_MPI)
        write(1,'(100E15.7)') (alldata(k),k=1,nelements)
#else
        write(1,'(I10,X,100E15.7)') sph(p)%porig, (alldata(k),k=1,nelements)
#endif
        
     end do

     close(1)

  end if  
#endif
! ----------------------------------------------------------------------------



! Write sink information if required
! ----------------------------------------------------------------------------
#if defined(SINKS)
  if (stot > 0) then
     debug_file = trim(adjustl(out_file))//".sinks"
     open(1, file=debug_file, status="unknown", form="formatted")
     write(6,*) "Debugging output file: ",trim(debug_file)

     do s=1,stot
        nelements = 0
        alldata = 0.0_PR
        
        ! Positions 
        do k=1,NDIM
           nelements = nelements + 1
           alldata(nelements) = sink(s)%r(k)*real(rscale,PR)
        end do

        ! Velocities 
        do k=1,NDIM
           nelements = nelements + 1
           alldata(nelements) = sink(s)%v(k)*real(vscale,PR)
        end do

        ! Accelerations 
        do k=1,NDIM
           nelements = nelements + 1
           alldata(nelements) = sink(s)%a(k)*real(ascale,PR)
        end do

        nelements = nelements + 1
        alldata(nelements) = sink(s)%h*real(rscale,PR)

        ! Now write all values to the file
        write(1,'(10E15.7)') (alldata(k),k=1,nelements)

     end do

     close(1)

  end if  
#endif
! ----------------------------------------------------------------------------


! ----------------------------------------------------------------------------
!#if defined(DEBUG_FREEFALL) || defined(FREEFALL_TEST) || defined(PLUMMER_TEST)
!  allocate(rad_ids(1:ptot+stot))
!  allocate(radius(1:ptot+stot))
!
!  do p=1,ptot
!     call distance2(rorigin(1:NDIM),p,dr(1:NDIM),drsqd)
!     rad_ids(p) = p
!     radius(p)  = sqrt(drsqd) + SMALL_NUMBER
!  end do
!#if defined(SINKS)
!  do s=1,stot
!     call distance3(rorigin(1:NDIM),sink(s)%r(1:NDIM),dr(1:NDIM),drsqd)
!     rad_ids(s+ptot) = s+ptot
!     radius(s+ptot)  = sqrt(drsqd) + SMALL_NUMBER
!  end do
!#endif
!
!! Order particles in order of increasing radii
!  call insertion_sort_real(ptot+stot,rad_ids(1:ptot+stot),radius(1:ptot+stot))
!
!! Write to first file
!#if defined(DEBUG_FREEFALL) || defined(FREEFALL_TEST)
!  debug_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".freefall"
!#elif defined(PLUMMER_TEST)
!  debug_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".plummer"
!#endif
!  open(1,file=debug_file,status="unknown",form="formatted",position="append")
!
!  mtot_temp = 0.0_DP
!  iaux = 1
!  do i=1,ptot
!     p = rad_ids(i)
!     if (p <= ptot) mtot_temp = mtot_temp + real(sph(p)%m,DP)
!#if defined(SINKS)
!     if (p > ptot) mtot_temp = mtot_temp + real(sink(p - ptot)%m,DP)
!#endif
!     if (mtot_temp >= 0.0099999999_DP*real(iaux,DP)) then
!        rmass(iaux) = radius(i)
!        iaux = iaux + 1
!     end if
!  end do
!#if defined(DEBUG_FREEFALL) || defined(FREEFALL_TEST)
!  write(1,'(101E15.7)') time/1.1107207, rmass(1:100)
!#elif defined(PLUMMER_TEST)
!  write(1,'(101E15.7)') time, rmass(1:100)
!#endif
!  close(1)
!
!#if defined(DEBUG_FREEFALL) || defined(FREEFALL_TEST)
!  debug_file = "freefall.dat"
!  open(1,file=debug_file,status="unknown",form="formatted")
!  do i=1,500
!     drsqd = min(0.01_PR*real(i,PR),0.9999999999_PR)
!     write(1,'(2E15.7)') 1.0_PR - 2.0_PR*INVPI*(asin(sqrt(drsqd)) - &
!          &sqrt(drsqd)*sqrt(1 - drsqd)),drsqd
!  end do
!  close(1)
!#endif
!
!  deallocate(radius)
!  deallocate(rad_ids)
!
!#endif
!! ----------------------------------------------------------------------------

  return
END SUBROUTINE write_data_debug
