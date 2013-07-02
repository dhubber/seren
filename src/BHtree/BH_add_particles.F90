! BH_ADD_PARTICLES.F90
! A. McLeod - 20/09/2012
! Adds particles to the tree. The tree must then be restocked.
! Includes routines for hydro and grav trees, this wouldn't be necessary
! with F2003...
! ============================================================================

#include "macros.h"

! ============================================================================

!#define DEBUG_BH_ADD_PARTICLES 1

#define BH_ADD_PART_SUB_NAME BH_add_particles_hydro
#define BH_ADD_PART_TYPE BHhydro_node
#define BH_ADD_PART_DIST BHtree(c)%rmax**2
#define BH_ADD_PART_GRAV 0

#include "BH_add_particles.h"

#if defined(SELF_GRAVITY)

#undef BH_ADD_PART_SUB_NAME
#undef BH_ADD_PART_TYPE
#undef BH_ADD_PART_DIST
#undef BH_ADD_PART_GRAV

#define BH_ADD_PART_SUB_NAME BH_add_particles_grav
#define BH_ADD_PART_TYPE BHgrav_node
#define BH_ADD_PART_DIST 0.0
#define BH_ADD_PART_GRAV 1

#include "BH_add_particles.h"

#endif