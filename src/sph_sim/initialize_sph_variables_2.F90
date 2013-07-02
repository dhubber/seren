! INITIALIZE_SPH_VARIABLES_2.F90
! D. A. Hubber - 1/10/2007
! Sets values for particular variables that need to be initialized AFTER 
! the first force calculations in setup.  Other variables (before the first 
! force calculation) are initialized in initialize_variables_1.F90.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_sph_variables_2
  use interface_module, only : find_equilibrium_temp_ws
  use filename_module, only : restart
  use particle_module
  use hydro_module
  use scaling_module
  use type_module
  use time_module
  use neighbour_module
  use timing_module
#if defined(SINKS)
  use sink_module
#endif
#if defined(RAD_WS)
  use Eos_module
#endif
#if defined(USE_MPI)
  use mpi_communication_module
  use mpi
#endif
  implicit none

  integer :: p               ! Particle counter
#if defined(SINKS)
  integer :: s               ! Sink counter
#endif
#if defined(USE_MPI)
  integer :: ierr            ! MPI error value
#endif

  debug2("Initializing variables for sph simulation [initialize_sph_variables_2.F90]")


! Calculate flux-limited diffusion terms
! ----------------------------------------------------------------------------
#if defined(RAD_WS)
  do p=pgasstart,pgasend
     call find_idens(sph(p)%rho,sph(p)%idens)
     call find_itemp(sph(p)%temp,sph(p)%itemp)
#if defined(RAD_WS_SINK_POT)
     sph(p)%column2 = (fcolumn**2)*sph(p)%gpot*sph(p)%rho
#else
     sph(p)%column2 = (fcolumn**2)*sph(p)%sphgpot*sph(p)%rho
#endif
#if defined(DIFFUSION)
     sph(p)%du_dt_diff = 0.0_PR
#endif
     call find_equilibrium_temp_ws(p)
 end do
#endif
! ----------------------------------------------------------------------------


! Record 'old' particle properties for integration scheme
! ----------------------------------------------------------------------------
  do p=1,ptot
     sph(p)%r_old(1:NDIM)  = sph(p)%r(1:NDIM)
     sph(p)%v_old(1:VDIM)  = sph(p)%v(1:VDIM)
#if defined(ENERGY_EQN)
     sph(p)%u_old          = sph(p)%u
#elif defined(ENTROPIC_FUNCTION) && defined(ENTROPY_EQN)
     sph(p)%Aold           = sph(p)%Aent
#endif
#if defined(RUNGE_KUTTA)
     sph(p)%v_half(1:VDIM) = sph(p)%v(1:VDIM)
#endif
#if defined(LEAPFROG_KDK)
     sph(p)%a_old(1:VDIM)  = sph(p)%a(1:VDIM)
#endif
#if defined(ENTROPIC_FUNCTION) && defined(ENTROPY_EQN) && defined(LEAPFROG_KDK)
     sph(p)%Ahalf          = sph(p)%Aent
#elif defined(ENERGY_EQN) && defined(LEAPFROG_KDK)
     sph(p)%dudt_old       = sph(p)%dudt
#endif
  end do


! Initialise sink particles info not contained in IC file
! ----------------------------------------------------------------------------
#if defined(SINKS)
  do s=1,stot
     sink(s)%rold(1:NDIM) = sink(s)%r(1:NDIM)
     sink(s)%vold(1:VDIM) = sink(s)%v(1:VDIM)
     sink(s)%menc = real(sink(s)%m,DP)
#if defined(RUNGE_KUTTA)
     sink(s)%vhalf(1:VDIM) = sink(s)%v(1:VDIM)
#endif
#if defined(LEAPFROG_KDK)
     sink(s)%aold(1:VDIM) = sink(s)%a(1:VDIM)
#endif
#if defined(SMOOTH_ACCRETION)
     !sink(s)%mmax = 4.0_PR*PI/3.0_PR*rhosink*sink(s)%radius**3
     sink(s)%trot = TWOPI_DP*sqrt(real(sink(s)%radius,DP)**3/sink(s)%menc)
#endif
  end do

  if (.not. restart) then
     do s=1,stot
        sink(s)%macc(1:DMDT_RANGE) = 0.0_DP
        sink(s)%tacc(1:DMDT_RANGE) = 0.0_DP
        sink(s)%dmdt               = 0.0_DP
        sink(s)%angmom(1:3)        = 0.0_DP
        sink(s)%angmomnet(1:3)     = 0.0_DP
     end do
  end if
#endif

! Set accdo to FALSE since we have already calculated the initial 
! accelerations for integration schemes which require them.
  sph(1:ptot)%accdo = .false.

! Set accdo_sinks depending on integration scheme used
#if defined(SINKS)
  nsearchnext = nsteps + nsearchstep
#if defined(RUNGE_KUTTA2) || defined(EULER) || defined(LEAPFROG_KDK)
  accdo_sinks = .true.
#elif defined(LEAPFROG_DKD)
  accdo_sinks = .false.
#endif
#endif

! Set time for next tree build and stock, and next ionizing HEALPix walk
#if defined(LEAPFROG_KDK)
  nbuild = nsteps + 2
  nstock = nsteps + 2
  nionize = nsteps + 2
#else
  nbuild = nsteps + 1
  nstock = nsteps + 1
  nionize = nsteps + 1
#endif
  nionall = nsteps  + nionallstep

! Set some MPI timing variables
#if defined(USE_MPI)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  calcstart = MPI_WTIME()
  calctime = 0._DP
  exporttime = 0._DP
  waittime = 0._DP
  sum_acctot = 0
  predict_acctot = real(ptot,PR) * real(loadbalance_nsteps,PR)
  sum_second_h = 0
  sum_costs = 0._PR
  non_start_time = MPI_WTIME()
#endif

  return
END SUBROUTINE initialize_sph_variables_2
