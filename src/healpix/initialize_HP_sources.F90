! INITIALIZE_HP_SOURCES.F90
! D. A. Hubber - 18/10/2009
! Create and initialize all HEALPix sources depending on Makefile options.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE initialize_HP_sources
  use interface_module, only : create_HP_source
  use filename_module
  use particle_module
  use hydro_module
  use type_module
  use sink_module
  use HP_module
  implicit none

  debug2("[initialize_HP_sources.F90]")

  HPtot = 0

! Construct arrays for HEALPix functions
  call mk_xy2pix(x2pix,y2pix)
  call mk_pix2xy(pix2x,pix2y)


! ----------------------------------------------------------------------------
#if defined(SINGLE_STATIC_SOURCE) || defined(SINGLE_SINK_SOURCE)

  call create_HP_source(-1,rstatic(1:NDIM))

! ----------------------------------------------------------------------------
#elif defined(MULTIPLE_SINK_SOURCES)

  sink(1:SMAX)%HPid = -1
  call calculate_sink_properties

#endif
! ----------------------------------------------------------------------------

  return
END SUBROUTINE initialize_HP_sources
