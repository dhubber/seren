! CREATE_HP_SOURCE.F90
! D. A. Hubber - 18/10/2009
! Create a source of ionizing UV radiation at location rsource.  If the 
! source is a dynamically moving sink particles, associate source with sink 
! id so the position and any other properties can be updated later on.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE create_HP_source(ssink,rsource)
  use interface_module, only : distance2,heapsort_real
  use particle_module
  use HP_module
  use sink_module
  use filename_module, only : restart
  implicit none

  integer, intent(in) :: ssink                 ! Associated sink id
  real(kind=PR), intent(in) :: rsource(1:NDIM) ! Position of source

  integer :: p                                 ! Particle counter
  integer, allocatable :: distids(:)           ! List of particles (unsorted)
  real(kind=PR) :: dr(1:NDIM)                  ! Relative displacement
  real(kind=PR) :: drsqd                       ! Distance squared
  real(kind=PR), allocatable :: distsqd(:)     ! Distance from source squared

  debug2("Creating a HP ionizing radiation source [create_HP_source.F90]")


! Set main ionizing source properties
! ----------------------------------------------------------------------------
  HPtot = HPtot + 1
  HPsource(HPtot)%sinkid = ssink
  HPsource(HPtot)%r(1:NDIM) = rsource(1:NDIM)
  if (ssink /= -1) sink(ssink)%HPid = HPtot
#if defined(SINGLE_SINK_SOURCE)
  if (ssink == 1) HPsource(1)%r(1:NDIM) = sink(1)%r(1:NDIM)
#endif
#if defined(IONIZING_UV_RADIATION)
  HPsource(HPtot)%N_LyC  = N_LyC
  HPsource(HPtot)%intmax = real(intmax,PR)
#endif
#if defined(STELLAR_WIND)
  HPsource(HPtot)%M_loss = M_loss
  HPsource(HPtot)%v_wind = v_wind
  HPsource(HPtot)%tlastwind = -1.0_PR
#endif
#if defined(PARTICLE_INJECTION_WINDS)
  if (.not. restart) HPsource(HPtot)%numb_inj = 0
  !print*,'SETTING TO ZERO AGAIN'
#endif
#if defined(STELLAR_LUMINOSITY)
  !HPsource(HPtot)%L_star = L_star
#endif

#if defined(DEBUG_CREATE_HP_SOURCE)
  write(6,*) "Created new HP source;  HPtot = ",HPtot
  write(6,*) "ssink : ",ssink,"   rsource : ",rsource(1:NDIM)
#endif


! Now allocate memory for list and order first list using (parallel) heapsort
! ----------------------------------------------------------------------------
  allocate(HPsource(HPtot)%distorder(1:pmax))
  allocate(distsqd(1:pmax))
  allocate(distids(1:pmax))

  do p=1,pmax
     HPsource(HPtot)%distorder(p) = p
  end do
  
  do p=1,ptot
     call distance2(HPsource(HPtot)%r(1:NDIM),p,dr(1:NDIM),drsqd)
     distsqd(p) = drsqd
     distids(p) = p
     HPsource(HPtot)%distorder(p) = p
  end do

#if defined(DEBUG_CREATE_HP_SOURCE)
  write(6,*) "Sorting initial particle list for source ",HPtot
#endif
  debug_timing("HP_ORDERING")

! Sort lists in order of increasing distance from source
  call heapsort_real(ptot,distsqd(1:ptot),distids(1:ptot))
  HPsource(HPtot)%distorder(1:ptot) = distids(1:ptot)

! Syncronise all timesteps now we have new HP source
!  call syncronise_timesteps

  deallocate(distids)
  deallocate(distsqd)

  return
END SUBROUTINE create_HP_source
