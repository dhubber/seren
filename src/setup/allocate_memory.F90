! ALLOCATE_MEMORY.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Allocates main data arrays.  These arrays are deallocated in reverse order 
! when the simulation terminates (in clean_up.F90).  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE allocate_memory(minimal)
  use particle_module
  use neighbour_module
  use kernel_module
  use type_module
  use sink_module
#if defined(BH_TREE)
  use tree_module  
#endif
#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  use Nbody_module
#endif
#if defined(RAD_WS)
  use Eos_module
#endif
#if defined(HEALPIX)
  use HP_module
#endif
#if defined(USE_MPI)
  use mpi_communication_module
#endif

  logical, intent(in) :: minimal     ! Use minimal memory for MPI decomposition

#if defined(USE_MPI)
#if defined(SELF_GRAVITY)
  integer :: i               ! Loop counter
  integer :: j               ! Summing variable
#endif
#endif

  debug1("Allocating memory for main simulation arrays [allocate_memory.F90]")


! Kernel arrays
! ----------------------------------------------------------------------------
#if defined(KERNEL_TABLES)
  debug2("Allocating kernel arrays...")
  allocate(w0table(0:KERNTOT))
  allocate(w1table(0:KERNTOT))
  allocate(w2table(0:KERNTOT))
  allocate(wgravtable(0:KERNTOT))
  allocate(wpottable(0:KERNTOT))
  allocate(womegatable(0:KERNTOT))
  allocate(wzetatable(0:KERNTOT))
#endif


! If this is a 'minimal' run, only allocate minimal SPH particles and
! sinks/stars
! ============================================================================
  if (minimal) then
     debug2("Allocating minimal particle data arrays...")
     pmax = 0
     allocate(minimal_sph(1:ptot))
     
! Only allocate SPH particles when there are some SPH particles
! (in case simulation is N-body only)
! ============================================================================
  else if (ptot > 0) then

     ! Set arbitrary maximum no. of particles here
     pmax = max(int(real(PMAXMULT,DP)*real(ptot,DP)),ptot)
     write(6,*) "PMAX : ",PMAXMULT,ptot,pmax

     ! Main particle data
     ! -----------------------------------------------------------------------
     debug2("Allocating main particle data arrays...")
     allocate(sph(1:pmax))
     !allocate(parray(1:DATATOT,1:pmax))

     ! Neighbour arrays
     ! -----------------------------------------------------------------------
#if defined(NEIGHBOUR_LISTS)
     debug2("Allocating neighbour arrays...")
     listtot = pmax
     pp_limit = 2*pp_gather
     allocate(pptot(1:listtot))
     allocate(pplist(1:pp_limit,1:listtot))
#endif

     ! BH tree arrays
     ! -----------------------------------------------------------------------
#if defined(BH_TREE)
     debug2("Allocating tree arrays...")
     if (LEAFMAX > 8) then
        cmax_grav  = (2*(pgas + pcdm)*PMAXMULT) / (LEAFMAX - 4)
        cmax_hydro = (2*pmax) / (LEAFMAX - 4)
     else
        cmax_grav  = 2*(pgas + pcdm)*PMAXMULT
        cmax_hydro = 2*pmax
     end if
     cmax_skeleton = cmax_hydro

#if defined(SELF_GRAVITY)
     allocate(BHgrav(0:cmax_grav))
#endif
#if defined(GHOST_PARTICLES)
     cmax_ghost = cmax_hydro
     allocate(BHghost(0:cmax_ghost))
#endif
     allocate(BHhydro(0:cmax_hydro))
     allocate(BHskeleton(0:cmax_skeleton))

     allocate(cellof(1:pmax))
     allocate(BHnextptcl(1:pmax))
     allocate(whichchild(1:pmax))
     allocate(BHstock(0:cmax_hydro))
#if defined(REORDER_TREE)
     allocate(hydroorder(0:cmax_grav))
#if defined(SELF_GRAVITY)
     allocate(gravorder(0:cmax_grav))
#endif
#endif

#if defined(USE_MPI) && defined(SELF_GRAVITY)
     cmax_remotegrav = 1
     j = 1
     do i=1,remotetreedepth
        j = j * NCHILD
        cmax_remotegrav = cmax_remotegrav + j
     end do
     allocate(BHlocal_grav(0:cmax_remotegrav))
     allocate(BHremote_grav(0:cmax_remotegrav,0:lastrank))
     allocate(treesend(0:endranklist))
     allocate(treerecv(0:endranklist))
#endif
#endif

     ! Radiative transport approximation arrays
     ! -----------------------------------------------------------------------
#if defined(RAD_WS) && defined(DEBUG_RAD)
     allocate(rad_info(1:10,ptot))
#endif

     ! HEALPix arrays
     ! -----------------------------------------------------------------------
#if defined(HEALPIX)
     imax = 100*ptot
     write(6,*) "imax : ",imax
     allocate(HPray(1:imax))
     !allocate(newtemp(1:pmax))
     !allocate(ionizedo(1:pmax))
#if defined(IONIZING_UV_RADIATION)
     !allocate(temp_min(1:pmax))
     !allocate(temp_aux(1:pmax))
#endif
#if defined(STELLAR_WIND)
     allocate(windsurface(1:pmax))
#endif
#if defined(DEBUG_HP_WALK_ALL_RAYS)
     allocate(whichHPlevel(1:pmax))
#endif
#endif

  end if
! ============================================================================


! Sink particle arrays
! ----------------------------------------------------------------------------
#if defined(SINKS) || defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  allocate(sink(1:SMAX))
#endif

! Star particle arrays
! ----------------------------------------------------------------------------
#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  allocate(star(1:SMAX))
#if defined(BINARY_STATS)
  allocate(binary(1:SMAX))
#endif
#endif


  return
END SUBROUTINE allocate_memory
