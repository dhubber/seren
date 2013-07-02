! INITIALIZE_NBODY_SPH_VARIABLES_2.F90
! D. A. Hubber - 1/10/2007
! Sets values for particular variables that need to be initialized AFTER 
! the first force calculations in setup.  Other variables (before the first 
! force calculation) are initialized in initialize_variables_1.F90
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_nbody_sph_variables_2
  use particle_module
  use hydro_module
  use scaling_module
  use type_module
  use time_module
  use neighbour_module
  use filename_module, only : restart
  use sink_module
  use Nbody_module
#if defined(RAD_WS)
  use Eos_module
#endif
  implicit none

  integer :: p               ! Particle counter
#if defined(SINKS)
  integer :: s               ! Sink counter
#endif

  debug2("Initializing variables [nbody_sph_initialize_variables_2.F90]")


! Set initial column densities
! ----------------------------------------------------------------------------
#if defined(RAD_WS)
  do p=pgasstart,pgasend
#if defined(RAD_WS_SINK_POT)
     sph(p)%column2 = (fcolumn**2)*sph(p)%gpot*sph(p)%rho
#else
     sph(p)%column2 = (fcolumn**2)*sph(p)%sphgpot*sph(p)%rho
#endif
     call find_equilibrium_temp_ws(p)
  end do
#endif


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


! Initialise star particles info not contained in IC file
! ----------------------------------------------------------------------------
  do s=1,stot
     star(s)%rold(1:NDIM) = star(s)%r(1:NDIM)
     star(s)%vold(1:VDIM) = star(s)%v(1:VDIM)
  end do

  if (.not. restart) then
     do s=1,stot
        star(s)%macc(1:DMDT_RANGE) = 0.0_DP
        star(s)%tacc(1:DMDT_RANGE) = 0.0_DP
        star(s)%dmdt               = 0.0_DP
        star(s)%angmom(1:3)        = 0.0_DP
     end do
  end if

! Set accdo to FALSE since we have already calculated the initial 
! accelerations for integration schemes which require it.
  sph(1:ptot)%accdo = .false.
  star(1:stot)%accdo = .false.

! Set time for next tree build and stock
#if defined(LEAPFROG_KDK)
  nbuild = nsteps + 2
  nstock = nsteps + 2
  nionize = nsteps + 2
#else
  nbuild = nsteps + 1
  nstock = nsteps + 1
  nionize = nsteps + 1
#endif
  nionall = nsteps + nionallstep

! Set sink radius for fixed (constant) multiple of h
#if defined(FIXED_HMULT_SINKRAD)
  sinkrad = sinkrad*((3.0_PR*pp_gather*mgas) / &
       & (32.0_PR*PI*pgas*rhosink))**(ONETHIRD)
#endif

#if defined(SMOOTH_ACCRETION) && defined(MINIMUM_H)
  if (stot > 0) then
     hmin = ((3.0_PR*pp_gather*mgas)/(32.0_PR*PI*pgas*rhosink))**(ONETHIRD)
  else
     hmin = 0.0_PR
  end if
#endif


  return
END SUBROUTINE initialize_nbody_sph_variables_2
