! HP_UPDATE.F90
! D. A. Hubber - 24/06/2011
! Update any HEALPix quantities.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE HP_update
  use interface_module, only : HP_walk_all_rays
  use particle_module
  use HP_module
  use time_module
  use scaling_module
  use sink_module
  use type_module
  implicit none

  integer :: i                    ! HEALPix source counter
  integer :: p                    ! Particle counter

  debug2("HP_update.F90")

! Calculate thermal properties of all (active) particles if on an ionization 
! step due to ionizing radiation from sources.
! ----------------------------------------------------------------------------
  if (nionall == nsteps) then

     nnewwind = 0

     ! Initialise all particle variables before HEALPix walk
     do p=1,ptot
#if defined(IONIZING_UV_RADIATION)
        sph(p)%newtemp = 0
        sph(p)%tempmin = 0.0_PR
        sph(p)%tempaux = 0.0_PR
#endif
#if defined(STELLAR_WIND)
        sph(p)%a_wind(1:NDIM) = 0.0_PR
#endif
#if defined(STELLAR_LUMINOSITY)
        sph(p)%dudt_rad = 0.0_PR
#endif
#if defined(DEBUG_HP_WALK_ALL_RAYS)
        whichHPlevel(p) = HP_LEVELS + 1
#endif
     end do

#if defined(SINGLE_SINK_SOURCE)
     if (stot > 0) then
        HPsource(1)%r(1:NDIM) = sink(1)%r(1:NDIM)
     else
        stop 'Error : SINGLE_SINK_SOURCE flag acivated for HEALPIX, &
             &but not sinks present'
     end if
#endif

     ! If there are any sources, call ionization/wind routines
     if (HPtot > 0) then
        do i=1,HPtot
           call HP_walk_all_rays(i)
        end do
#if defined(IONIZING_UV_RADIATION)
        if (nsteps /= 0) call write_ionization_data
#endif
     end if
     
     ! Update time counters
     if (nionall == nsteps) nionall = nionall + nionallstep

#if defined(PARTICLE_INJECTION_WINDS)
     if (nnewwind > 0) then
        nbuild = nsteps
        nstock = nsteps
        pgas = pgas + nnewwind
        ptot = ptot + nnewwind
        call types
        call tree_update(nbuild,nstock)
        !write(6,*) "HP_update"
        !write(6,*) "ptot : ",ptot,"    nnewwind : ",nnewwind
        !if (nnewwind > 0) call diagnostics
     end if
#endif

  end if
! ----------------------------------------------------------------------------

  return
END SUBROUTINE HP_update
