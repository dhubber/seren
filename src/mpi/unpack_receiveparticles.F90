! UNPACK_RECEIVEPARTICLES.F90
! A. McLeod - 1/08/08
! Unpacks particles from the receiveparticles array
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine unpack_receiveparticles(pcount)
   use particle_module
   use hydro_module
   use type_module
   use time_module
   use mpi_communication_module
   use eos_module
   use type_module
#ifdef DEBUG_TRANSFER
   use filename_module
#endif
   implicit none

   integer, intent(in) :: pcount(0:endranklist) ! No. of ptcls in diff sections

   integer :: ptype                      ! Particle type
   integer :: ptypes(1:ntypes)           ! Particle type identifier
   integer :: i, j                       ! Loop counters
   integer :: p                          ! Particle counter
   integer :: t                          ! Type counter
   integer :: t_ord                      ! Order of type in memory
   integer :: typeslot                   ! Type slot
   integer :: summove                    ! Sum variable
   integer :: numcopy(1:ntypes)          ! Number to copy
   integer :: shift(1:ntypes)            ! Amount to shift by
   integer :: copyfrom(1:ntypes)         ! Copy from here
   integer :: loadhere(1:ntypes)         ! Where to load particles
#ifdef DEBUG_TRANSFER
   character(len=256) :: out_file
#endif

   debug_timing("UNPACK_RECEIVE_PARTICLES")

#ifdef DEBUG_TRANSFER
   out_file = trim(adjustl(run_id))//".receiving."//trim(adjustl(MPI_ext))
   open(2,file=out_file,position='APPEND',form='formatted')
   write (2,*) " ==================== ", nsteps, " ===================="
#endif

   ! Determine number of each particle types
   ! -------------------------------------------------------------------------
   ptypes = 0
   do i=0,lastrank-1
      do p=1,pcount(i)
         ptype = receiveparticles(i)%data(p)%ptype
         ptypes(ptype) = ptypes(ptype) + 1
      end do
   end do
!    if (sum(ptypes) > 0) write (6,*) "ptypes = ", ptypes
!    if (sum(ptypes) > 0) write (6,*) "ptot = ", ptot

   ! Make some room for that many particle of each type
   ! -------------------------------------------------------------------------
   summove = 0
   t_ord = typeorder(1)
   typeslot = typeinfo(t_ord)%N + 1
   loadhere(1) = typeslot
!    if (sum(ptypes) > 0) write (6,*) "summove, typeslot, loadhere(1): ", summove, typeslot, loadhere
   do t=2,ntypes
      t_ord = typeorder(t)
      summove = summove + ptypes(t-1)
      numcopy(t) = min(summove,typeinfo(t_ord)%N)
      shift(t) = max(summove,typeinfo(t_ord)%N)
      copyfrom(t) = typeslot
      typeslot = typeslot + typeinfo(t_ord)%N
      loadhere(t) = typeslot + summove
   end do
!    if (sum(ptypes) > 0) write (6,*) "Generated the following:"
!    if (sum(ptypes) > 0) write (6,*) "numcopy = ", numcopy
!    if (sum(ptypes) > 0) write (6,*) "copyfrom = ", copyfrom
!    if (sum(ptypes) > 0) write (6,*) "shift = ", shift
   do t=ntypes,2,-1
!       if (sum(ptypes) > 0) write (6,*) "Doing for type ", t
      p = copyfrom(t)
      do i=1,numcopy(t)
!          if (sum(ptypes) > 0) write (6,*) "i = ", i
         !call copy_particle_data(p+shift(t),p)
         sph(p+shift(t)) = sph(p)
         p = p + 1
      end do
   end do

   ! Update numtypes
   ! -------------------------------------------------------------------------
   do t=1,ntypes
      t_ord = typeorder(t)
      typeinfo(t_ord)%N = typeinfo(t_ord)%N + ptypes(t)
   end do
   ! Types uses these to set typeinfo, so we must set them here
   pboundary = typeinfo(boundaryid)%N
   picm = typeinfo(icmid)%N
   pgas = typeinfo(gasid)%N
   pcdm = typeinfo(cdmid)%N
   pdust = typeinfo(dustid)%N
   pion = typeinfo(ionid)%N

   ! Update ptot
   ptot = ptot + sum(ptypes)

   ! Load data from receiveparticles into arrays
   ! -------------------------------------------------------------------------
   do i=0,lastrank-1
      do j=1,pcount(i)
         ptype = receiveparticles(i)%data(j)%ptype
         p = loadhere(ptype)
         loadhere(ptype) = p + 1
         sph(p) = receiveparticles(i)%data(j)
#ifdef DEBUG_TRANSFER
         write (2,*) "p ", p, "; porig ", sph(p)%porig
         write (2,*) receiveparticles(i)%data(j)
#endif
      end do
   end do
   ! -------------------------------------------------------------------------
#ifdef DEBUG_TRANSFER
   write (2,*) "------------------------------"
   close(2)
#endif

   ! Call types to update non-pointer type integers (phydrostart,pgravitystart)
   call types

   return
end subroutine unpack_receiveparticles
