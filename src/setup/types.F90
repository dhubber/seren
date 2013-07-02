! TYPES.F90
! C. P. Batty & D. A. Hubber - 11/9/2007
! Sets loop counters and limits for different particle types.  
! Requires that particle types are ordered as follows : 
! 1 - Boundary particles
! 2 - Inter-cloud medium particles
! 3 - (Self-gravitating) gas particles
! 4 - Cold dark matter (CDM) particles
! 5 - Dust particles
! 6 - Ion particles
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE types
  implicit none
  
  debug2("Setting particle type counters [types.F90]")
  
  call types_1        ! Allocate type variables
  call types_2        ! Allocate ptype to sph variables
  
  return
END SUBROUTINE types

SUBROUTINE types_1
  use type_module
  use particle_module
  implicit none

  integer :: i                         ! ..
  integer :: iprev                     ! ..
  integer :: j                         ! ..

! Set numbers of particles for each type
  typeinfo(boundaryid)%N = pboundary
  typeinfo(icmid)%N = picm
  typeinfo(gasid)%N = pgas
  typeinfo(cdmid)%N = pcdm
  typeinfo(dustid)%N = pdust
  typeinfo(ionid)%N = pion

! Calculate particle id ranges of all types
  i = typeorder(1)
  typeinfo(i)%pfirst = min(1,typeinfo(i)%N)
  typeinfo(i)%plast = typeinfo(i)%N
  do j=2,ntypes
     i = typeorder(j)
     iprev = typeorder(j-1)
     typeinfo(i)%pfirst = typeinfo(iprev)%plast + min(1,typeinfo(i)%N)
     typeinfo(i)%plast = typeinfo(iprev)%plast + typeinfo(i)%N
  end do

! Start and end points for particle types in main arrays
  pboundarystart = min(1,pboundary)
  pboundaryend   = pboundary
  picmstart      = pboundary + min(1,picm)
  picmend        = pboundary + picm
  pgasstart      = pboundary + picm + min(1,pgas)
  pgasend        = pboundary + picm + pgas
  pcdmstart      = pboundary + picm + pgas + min(1,pcdm)
  pcdmend        = pboundary + picm + pgas + pcdm
  pduststart     = pboundary + picm + pgas + pcdm + min(1,pdust)
  pdustend       = pboundary + picm + pgas + pcdm + pdust
  pionstart      = pboundary + picm + pgas + pcdm + pdust + min(1,pion)
  pionend        = pboundary + picm + pgas + pcdm + pdust + pion


! Start and end points for calculating hydro and gravitational forces in 
! main arrays
  phydrostart    = pboundary + 1
  phydroend      = pboundary + picm + pgas
  pgravitystart  = pboundary + picm + 1
  pgravityend    = pboundary + picm + pgas + pcdm

! Start and end points in main arrays for building various trees
  phtreestart    = 1
  phtreeend      = pboundary + picm + pgas + pcdm
  pgtreestart    = pboundary + picm + 1
  pgtreeend      = pboundary + picm + pgas + pcdm
!  psftreestart   = 0
!  psftreeend     = 0

! For the moment, allow all particles to be ionized by UV radiation
  pionizestart   = 1

#if defined(DEBUG_TYPES)
  write(6,*) "pboundary : ",pboundary,pboundarystart,pboundaryend
  write(6,*) "picm      : ",picm,picmstart,picmend
  write(6,*) "pgas      : ",pgas,pgasstart,pgasend
  write(6,*) "pcdm      : ",pcdm,pcdmstart,pcdmend
  write(6,*) "pdust     : ",pdust,pduststart,pdustend
  write(6,*) "pion      : ",pion,pionstart,pionend
  write(6,*) "phydro    : ",phydrostart,phydroend
  write(6,*) "pgravity  : ",pgravitystart,pgravityend
  write(6,*) "phtree    : ",phtreestart,phtreeend
  write(6,*) "pgtree    : ",pgtreestart,pgtreeend
#endif

  return
END SUBROUTINE types_1

SUBROUTINE types_2
  use type_module
  use particle_module
  implicit none

  integer :: i                         ! ..
  integer :: p                         ! ..

! Set ptypes for all particles
  do i=1,ntypes
     if (typeinfo(i)%N == 0) cycle
     do p=typeinfo(i)%pfirst,typeinfo(i)%plast
        sph(p)%ptype = i
     end do
  end do

  return
END SUBROUTINE types_2