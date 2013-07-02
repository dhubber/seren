! REMOVE_PARTICLES.F90
! A. McLeod - 31/07/08
! Removes particles from arrays, and shuffles down particles as required
! =========================================================================

#include "macros.h"

! =========================================================================
subroutine remove_particles(pp_tot,plist,remove_from_tree)
   use time_module, only : nbuild,nstock,nsteps
   use type_module
   use particle_module
   use mpi_communication_module
   use interface_module, only : BH_remove_particles
   implicit none

   integer, intent(in) :: pp_tot            ! Number of ptcls to delete
   integer, intent(in) :: plist(1:pp_tot)   ! List of ptcls to delete
   logical, intent(in) :: remove_from_tree  ! Remove particles from tree?
   
   integer :: curslot                    ! Current slot
   integer :: first_in_type(1:ntypes)    ! plist id of first deleted ptcl
   integer :: last_in_type(1:ntypes)     ! plist id of last deleted ptcl
   integer :: move                       ! Number of particles to move
   integer :: new_numtypes(1:ntypes)     ! New no. ptcls. after deleting
   integer :: num_del_type(1:ntypes)     ! Number of each type deleted
   integer :: numdel                     ! Number deleted summing variable
   integer :: newslots(1:ntypes)         ! New slots after deleting ptcls
   integer :: p                          ! Particle counter
   integer :: pp                         ! Particle counter
   integer :: plive                      ! Live particle counter
   integer :: plistdead                  ! Dead particle counter
   integer :: t                          ! Type counter
   integer :: t_ord                      ! Order of type in memory
   integer :: typeslots(1:ntypes)        ! Slots at which types start
   integer :: typeslotsend(1:ntypes)     ! Last particle in type
#if defined(BH_TREE)
   integer :: origid(1:ptot)             ! Old position of particle
   integer :: newid(1:ptot)              ! For removing from tree
#endif

   ! ids relating to p and the arrays of deleted AND non-deleted particles

   ! ids relating to plist and the deleted particles list

   ! Throughout we are using types numbered 1 to ntypes in the order they
   ! appear in memory, so typeinfo must be referenced by typeinfo(typeorder( ))

#ifdef DEBUG_TRANSFER
   write (1,*) "ptot = ", ptot
   write (1,*) "Removing ", pp_tot, " particles: ", plist(1:pp_tot)
#endif

#ifdef BH_TREE
   ! Prepare the list of new particle ids for use in BH_remove_particles
   if (remove_from_tree) then
      do p=1,ptot
         origid(p) = p
      end do
   end if
#endif

   ! Find slot at which types currently start and end
   typeslots(1) = 1
   typeslotsend(1:ntypes) = ptot
   do t=2,ntypes
      t_ord = typeorder(t-1)
      typeslots(t) = typeslots(t-1) + typeinfo(t_ord)%N
   end do
   do t=1,ntypes-1
      typeslotsend(t) = typeslots(t+1) - 1
   end do

   ! Find out which is the first particle of each type
   num_del_type = 0
   curslot = 1
   first_in_type(1) = 1

   ! -------------------------------------------------------------------------
   do p=1,pp_tot
      ! Assume particles are correctly ordered in memory
      if (plist(p) >= typeslots(curslot+1)) then
         ! We are in a new type, find which one in case we skip some
         do t=curslot+1,ntypes
            if (plist(p) >= typeslots(t)) then
               ! Particle is above typeslots(t)
               curslot = t
               first_in_type(curslot) = p
            else
               ! Particle is below typeslots(t)
               exit
            end if
         end do
      end if

      if (curslot == ntypes) then
         num_del_type(curslot) = pp_tot - p + 1
         exit ! We are done, the last type is the rest of the particles
      else
         ! Increment number of deleted particles in this type
         num_del_type(curslot) = num_del_type(curslot) + 1
      end if
   end do
   ! -------------------------------------------------------------------------

   ! Find last deleted particle in each type
   do t=1,ntypes-1
      last_in_type(t) = first_in_type(t) + num_del_type(t) - 1
   end do
   last_in_type(ntypes) = pp_tot

   ! Now remove particles from each type by looping over types
   ! and shuffling particles down from the end
   ! -------------------------------------------------------------------------
   bigdo: do t=1,ntypes
      if (num_del_type(t) == 0) cycle ! If this type has no particles to delete

      ! Algorithm for removing particles
      p = first_in_type(t)
      plive = typeslotsend(t)
      plistdead = last_in_type(t)
      do
         ! Find last live particle
         if (plive == plist(plistdead)) then
            plive = plive - 1
            plistdead = plistdead - 1
            if (plistdead < first_in_type(t)) then
               ! We have finished
               exit
            end if
            cycle
         else
            ! We have a live particle, we can delete by copying
            if (plist(p) > plive) then
               ! We have finished
               exit
            end if
            sph(plist(p)) = sph(plive)
#if defined(BH_TREE)
            ! plive goes to plist(p) which is now dead
            if (remove_from_tree) then
               origid(plist(p)) = origid(plive)
            end if
#endif
            plive = plive - 1
            p = p + 1
            if (p > pp_tot) then
               ! We have finished with everything
               exit bigdo
            end if
         end if
      end do

   end do bigdo
   ! -------------------------------------------------------------------------

! Now shuffle types down, by taking particles from the end of each type
! and moving it to where the type should start

   ! Translate numtypes(-1:9) to new_numtypes(1:ntypes)
   ! Update this and also find the CURRENT end of live particles in types
   do t=1,ntypes
      t_ord = typeorder(t)
      new_numtypes(t) = typeinfo(t_ord)%N - num_del_type(t)
      typeslotsend(t) = typeslotsend(t) - num_del_type(t)
   end do

   ! Calculate where the new types will start
   newslots(1) = 1
   do t=2,ntypes
      newslots(t) = newslots(t-1) + new_numtypes(t-1)
   end do

   numdel = num_del_type(1)

   ! Do the work of copying particles down into their new slots
   ! -------------------------------------------------------------------------
   do t=2,ntypes

      move = min(new_numtypes(t),numdel)
      plive = typeslotsend(t) - move + 1 ! Last live particle of type - move+1
      pp = newslots(t)                   ! Slot where first particle should go

      do p=1,move
         sph(pp) = sph(plive)
#if defined(BH_TREE)
         ! plive goes to pp which is now dead
         if (remove_from_tree) then
            origid(pp) = origid(pp)
         end if
#endif
         plive = plive + 1
         pp = pp + 1
      end do

      numdel = numdel + num_del_type(t)

   end do
   ! -------------------------------------------------------------------------
   
#if defined(BH_TREE)
   ! Remove dead particles from tree
   if (remove_from_tree) then
      newid = -2
      do p=1,pp_tot
         newid(plist(p)) = -1
      end do
      do p=1,ptot-numdel
         pp = origid(p)
         newid(pp) = p
      end do
      call BH_remove_particles(newid)
   end if
#endif

   ! Update ptot and numtypes
   ptot = ptot - numdel
   do t=1,ntypes
      t_ord = typeorder(t)
      typeinfo(t_ord)%N = new_numtypes(t)
   end do

   ! Types uses these to set typeinfo, so we must set them here
   pboundary = typeinfo(boundaryid)%N
   picm = typeinfo(icmid)%N
   pgas = typeinfo(gasid)%N
   pcdm = typeinfo(cdmid)%N
   pdust = typeinfo(dustid)%N
   pion = typeinfo(ionid)%N
   ! Call types to update non-pointer type integers (phydrostart, pgravitystart)
   call types

   ! Need to restock trees before any tree search is done

#ifdef DEBUG_TRANSFER
   write (1,*) "ptot = ", ptot
#endif

   return
end subroutine remove_particles
