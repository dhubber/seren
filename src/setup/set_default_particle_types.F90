! SET_DEFAULT_PARTICLE_TYPES
! D. A. Hubber & A.McLeod  - 14/03/2011 
! Set up all properties for default particle types, including masks for 
! various SPH routines.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE set_default_particle_types
  use type_module
  implicit none
  
  integer :: i                  ! Type counter
  
  debug2("Initialise all variables for particle types [set_default_particle_types.F90]")
  
  ! Set order of particle types in main SPH arrays here.
  typeorder(1) = boundaryid
  typeorder(2) = icmid
  typeorder(3) = gasid
  typeorder(4) = cdmid
  typeorder(5) = dustid
  typeorder(6) = ionid
  ntypes = nmaxtypes
  
  ! Set false values to all flags for now, and set true flags below.
  do i=1,nmaxtypes
     typeinfo(i)%static = .false.
     typeinfo(i)%h      = .false.
     typeinfo(i)%h2     = .false.
     typeinfo(i)%hydro  = .false.
     typeinfo(i)%grav   = .false.
     typeinfo(i)%drag   = .false.
  end do
  sphmask(1:nmaxtypes) = .false.
  hydromask(1:nmaxtypes) = .false.
  gravmask(1:nmaxtypes) = .false.
  dragmask(1:nmaxtypes) = .false.
  
  ! Boundary particles
  ! --------------------------------------------------------------------------
  typeinfo(boundaryid)%name = "boundary"
  typeinfo(boundaryid)%dragonid = 6
  typeinfo(boundaryid)%static = .true.
  typeinfo(boundaryid)%h(boundaryid) = .true.
  typeinfo(boundaryid)%h(icmid) = .true.
  typeinfo(boundaryid)%h(gasid) = .true.
  typeinfo(boundaryid)%h(cdmid) = .true.
#if defined(TWO_FLUIDS)
  typeinfo(boundaryid)%h2(dustid) = .true.
#endif
  
  ! icm particles
  ! --------------------------------------------------------------------------
  typeinfo(icmid)%name = "icm"
  typeinfo(icmid)%dragonid = 9
  typeinfo(icmid)%h(boundaryid) = .true.
  typeinfo(icmid)%h(icmid) = .true.
  typeinfo(icmid)%h(gasid) = .true.
  typeinfo(icmid)%h(cdmid) = .true.
  typeinfo(icmid)%hydro(boundaryid) = .true.
  typeinfo(icmid)%hydro(icmid) = .true.
  typeinfo(icmid)%hydro(gasid) = .true.
#if defined(TWO_FLUIDS)
  typeinfo(icmid)%h2(dustid) = .true.
#endif
  
  ! gas particles
  ! --------------------------------------------------------------------------
  typeinfo(gasid)%name = "gas"
  typeinfo(gasid)%dragonid = 1
  typeinfo(gasid)%h(boundaryid) = .true.
  typeinfo(gasid)%h(icmid) = .true.
  typeinfo(gasid)%h(gasid) = .true.
  typeinfo(gasid)%h(cdmid) = .true.
  typeinfo(gasid)%hydro(boundaryid) = .true.
  typeinfo(gasid)%hydro(icmid) = .true.
  typeinfo(gasid)%hydro(gasid) = .true.
  typeinfo(gasid)%grav(gasid) = .true.
  typeinfo(gasid)%grav(cdmid) = .true.
#if defined(TWO_FLUIDS)
  typeinfo(gasid)%h2(dustid) = .true.
#endif
  
  ! cdm particles
  ! --------------------------------------------------------------------------
  typeinfo(cdmid)%name = "cdm"
  typeinfo(cdmid)%eos = "none"
  typeinfo(cdmid)%dragonid = 10
  typeinfo(cdmid)%h(boundaryid) = .true.
  typeinfo(cdmid)%h(icmid) = .true.
  typeinfo(cdmid)%h(gasid) = .true.
  typeinfo(cdmid)%h(cdmid) = .true.
  typeinfo(cdmid)%grav(gasid) = .true.
  typeinfo(cdmid)%grav(cdmid) = .true.
  
  ! dust particles
  ! --------------------------------------------------------------------------
  typeinfo(dustid)%name = "dust"
  typeinfo(dustid)%dragonid = 11
#if defined(TWO_FLUIDS)
  typeinfo(dustid)%h(boundaryid) = .true.
  typeinfo(dustid)%h(icmid) = .true.
  typeinfo(dustid)%h(gasid) = .true.
  typeinfo(dustid)%h2(dustid) = .true.
#endif
  
  ! ion particles
  ! --------------------------------------------------------------------------
  typeinfo(ionid)%name = "ion"
  typeinfo(ionid)%dragonid = 12
  typeinfo(ionid)%h(ionid) = .true.
  
  
  ! Now set-up mask arrays for calculating SPH quantities and forces.
  ! ==========================================================================
  sphmask(boundaryid) = .true.
  sphmask(icmid) = .true.
  sphmask(gasid) = .true.
  sphmask(cdmid) = .true.
  
  hydromask(icmid) = .true.
  hydromask(gasid) = .true.
  
  gravmask(gasid) = .true.
  gravmask(cdmid) = .true.
  
  dragmask(dustid) = .true.
  
  
  return
END SUBROUTINE set_default_particle_types

