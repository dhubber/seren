! SINK_TRANSFER_PARTICLES.F90
! A. McLeod - 31/07/08
! ...
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine sink_transfer_particles(psink)
   use neighbour_module, only : pp_gather
   use particle_module
   use sink_module
   use mpi_communication_module
   use time_module
   use mpi
!   use domain_comparisons
   use tree_module
   use type_module
! #ifdef DEBUG_TRANSFER
!   use filename_module
! #endif
   use interface_module, only : BH_add_particles
   implicit none

   integer, intent(in) :: psink             ! Candidate sink particle (or 0)
   integer :: i                             ! Loop counters
   integer :: d                             ! Domain counter
   integer :: p                             ! Particle counter
   integer :: s, ss                         ! Sink counter
   integer :: sstart                        ! Start for sink loop
   integer :: ierr                          ! MPI error value
   integer :: pleft                         ! Number of particles left
   integer :: old_ptot                      ! ptot before new particle arrival
   integer :: numsendto(0:lastrank)         ! Number sent to domain
   integer :: numdel                        ! Number to delete
   integer :: pcount(0:endranklist)         ! Length of message
   integer :: recv(0:endranklist)           ! Array of receive requests
   integer :: slots(0:lastrank)             ! Slots
   integer :: send(0:endranklist)           ! Array of send requests
   integer :: stat(MPI_STATUS_SIZE)         ! MPI Status variable
   integer, allocatable :: deletelist(:)    ! List of ptcls to delete
   integer, allocatable :: newlist(:)       ! List of ptcls to add
   integer, allocatable :: sinkflag(:)      ! Flag sink-particle overlap
   real(kind=DP) :: rs_dp(1:NDIM)           ! Position of sink s
   real(kind=DP) :: rads                    ! New sink radius
   real(kind=DP) :: new_sink_info(1:NDIM+2) ! New sink information for sharing
   real(kind=DP) :: dr(1:NDIM)              ! Relative position vector
   real(kind=DP) :: drsqd                   ! Distance squared
   real(kind=DP) :: drsqd2                  ! Second distance squared variable

   ! If there's only one MPI task, then return immediately
   if (numtasks==1) return
     
   ! Share the position and radius of any new sink particle
   if (psink > 0) then
      ! We, and in theory we alone, have a new sink to share
#if defined(HMULT_SINKRAD)
      rads = real(sinkrad*sph(psink)%h,DP)
#else
      rads = real(sinkrad,DP)
#endif
      new_sink_info(1:NDIM) = real(sph(psink)%r(1:NDIM),DP)
      new_sink_info(NDIM + 1) = rads
      new_sink_info(NDIM + 2) = real(rank,DP) + 1.0001_DP
                 ! Has to be rank+1 so we can tell if there is a sink on rank 0
   else
      new_sink_info = 0.0_DP
   end if
   
   if (psink == -2) then
      ! If psink == -2, indicates there cannot be a new sink
      sstart = 1
   else
      ! Check for new sink (by receiving information)
      WAIT_TIME_MACRO
      call MPI_ALLREDUCE(MPI_IN_PLACE, new_sink_info, NDIM+2, &
                        & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALC_TIME_MACRO
   
      if (new_sink_info(NDIM+2) > 0.1) then
         ! There is a new sink
         sstart = 0
         new_sink_info(NDIM+2) = new_sink_info(NDIM+2) - 1.0_DP
      else
         ! There is not a new sink
         sstart = 1
      end if
   end if

   allocate(sinkflag(1:ptot))
   sinkflag(1:ptot) = -1

! Find id of sink to which particle is closest and within the sink radius
! ----------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dr,drsqd,drsqd2,rads,rs_dp,s,ss)
   do p=pgasstart,pgasend
      do s=sstart,stot ! s=0 represents new sink
         if (s==0) then
            rs_dp(1:NDIM) = new_sink_info(1:NDIM)
            rads = real(new_sink_info(NDIM + 1),DP)
         else
            if (.NOT. sink(s)%accrete) cycle
            rs_dp(1:NDIM) = real(sink(s)%r(1:NDIM),DP)
            rads = real(sink(s)%radius,DP)
         end if
         call distance2_dp(rs_dp(1:NDIM),p,dr(1:NDIM),drsqd)

         if (drsqd < rads*rads) then
            if (sinkflag(p) == -1) then
               sinkflag(p) = s
            else
               ss = sinkflag(p)
               if (ss == 0) then
                  rs_dp(1:NDIM) = new_sink_info(1:NDIM)
               else
                  rs_dp(1:NDIM) = real(sink(ss)%r(1:NDIM),DP)
               end if
               call distance2_dp(rs_dp(1:NDIM),p,dr(1:NDIM),drsqd2)
               if (drsqd < drsqd2) sinkflag(p) = s
            end if
         end if

      end do
   end do
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------------

   debug_timing("SINK_TRANSFER_PARTICLES SEND")
   allocate(exportlist(1:ptot,0:lastrank))
   numsendto = 0
   numdel = 0
   pleft = ptot
   
   ! Go through each particle; add to appropriate list to send
   ! --------------------------------------------------------------------------
   do p=1,ptot

      ! Test if the particle is to be accreted
      s = sinkflag(p)
      if (s /= -1) then
         ! Prevent us from having less than pp_gather+1 particles in normal SPH
         ! or less than 2*pp_gather particles in grad-h SPH
#ifdef GRAD_H_SPH
         if (pleft <= 2*pp_gather) then
#else
         if (pleft <= pp_gather+1) then
#endif
            write (6,*) "Too few particles left in task ", rank, &
                       &" due to sink accretion!"
         end if
         pleft = pleft - 1
         if (s==0) then
            d = floor(new_sink_info(NDIM+2))
         else
            d = sink(s)%domain
         end if
         if (d==rank) cycle
         numsendto(d) = numsendto(d) + 1
         exportlist(numsendto(d),d) = p
      end if
   end do
   ! --------------------------------------------------------------------------

   ! Allocate storage to send particles from
#ifdef DEBUG_TRANSFER
   write (6,*) rank, ": allocating sendparticles from 1 to ", sum(numsendto)+1
#endif
   allocate(sendparticles(1:sum(numsendto)+1)) ! +1 so even if sending no ptcls...
   allocate(deletelist(1:sum(numsendto)))

   ! Pack the particle data into the sendparticles array
   slots(0) = 1
   do d=1,lastrank
      ! Slot domain d starts at in sendparticles (???)
       slots(d) = slots(d-1) + numsendto(d-1)  
   end do

   ! --------------------------------------------------------------------------
   do d=0,lastrank
      if (d==rank) cycle
      i = slots(d)
      do p=1,numsendto(d)
         sendparticles(i) = sph(exportlist(p,d))
         i = i + 1
         numdel = numdel + 1
         deletelist(numdel) = exportlist(p,d)
      end do
   end do
   ! --------------------------------------------------------------------------

  ! Sort deletelist of particles that have left
   if (numdel > 1) then
#ifdef HEAPSORT
      call heapsort_int(numdel,deletelist(1:numdel))
#else
      call insertion_sort_int(numdel,deletelist(1:numdel))
#endif
   end if

   ! Send those particles to the other domains
   ! --------------------------------------------------------------------------
   i = 0
   do d=0,lastrank
      if (d==rank) cycle
      call MPI_ISEND(sendparticles(slots(d)),numsendto(d),MPI_PARTICLE,&
      &d,TRANSFER_TAG,MPI_COMM_WORLD,send(i),ierr)
      i = i + 1
   end do
   ! --------------------------------------------------------------------------

   deallocate(exportlist)

   debug_timing("TRANSFER_PARTICLES REMOVE_DEAD")

   ! Now remove all dead particles
   i = sum(numsendto)
   if (i>0) then
      call remove_particles(i,deletelist,.TRUE.)
   end if
   deallocate(deletelist)

   old_ptot = ptot ! AFTER removing particles but BEFORE adding particles

   allocate(receiveparticles(0:lastrank))

   debug_timing("TRANSFER_PARTICLES RECEIVE")

   ! Wait for receives from other domains to determine total incoming number
   ! --------------------------------------------------------------------------
   do i=0,lastrank-1
      ! Test for any message - wait until one arrives (blocking probe)
      WAIT_TIME_MACRO
      call MPI_PROBE(MPI_ANY_SOURCE, TRANSFER_TAG, MPI_COMM_WORLD, stat, ierr)
      CALC_TIME_MACRO
      d = stat(MPI_SOURCE)
      call MPI_GET_COUNT(stat,MPI_PARTICLE,pcount(i),ierr)

      ! Allocate space, and post receives of incoming particles
      allocate(receiveparticles(i)%data(1:max(pcount(i),1)))

      ! Post the receive - until this is posted the receive does not actually occur
      ! and may or may not be buffered
      call MPI_IRECV(receiveparticles(i)%data(1),pcount(i),MPI_PARTICLE,&
         &d,TRANSFER_TAG,MPI_COMM_WORLD,recv(i),ierr)

   end do
   ! --------------------------------------------------------------------------

   ! If arrays are too small OR too large, expand them OR shrink them
   if ((sum(pcount) + ptot) > pmax) then
      call expand(ceiling(real(PMAXMULT,DP)*real(ptot + sum(pcount))))
   end if

   ! Wait for all sends and receives to complete
   debug_timing("WAIT TO RECEIVE PARTICLES")
   WAIT_TIME_MACRO
   call MPI_WAITALL(lastrank,recv,MPI_STATUSES_IGNORE,ierr)
   CALC_TIME_MACRO

   ! Unpack data into arrays
   call unpack_receiveparticles(pcount)

   debug_timing("WAIT TO SEND PARTICLES")

   WAIT_TIME_MACRO
   call MPI_WAITALL(lastrank,send,MPI_STATUSES_IGNORE,ierr)
   CALC_TIME_MACRO
   deallocate(sendparticles)

   ! Deallocate receiveparticles
   do i=0,lastrank-1
      deallocate(receiveparticles(i)%data)
   end do
   deallocate(receiveparticles)

#if defined(BH_TREE)
   ! Add particles to hydro and gravity trees
   if (sum(pcount) > 0) then
     
      allocate(newlist(1:sum(pcount)))
      newlist = (/(i,i=old_ptot+1,ptot)/)
      write (6,*) "add to hydro tree"
      call BH_add_particles(sum(pcount),newlist,cmax_hydro,ctot_hydro,&
                          &ltot_hydro,first_cell_hydro,last_cell_hydro,BHhydro)
#if defined(SELF_GRAVITY)
      write (6,*) "add to grav tree"
      call BH_add_particles(sum(pcount),newlist,cmax_grav,ctot_grav,&
                          &ltot_grav,first_cell_grav,last_cell_grav,BHgrav)
#endif
      deallocate(newlist)

   end if
#endif

   ! Need to restock trees after particle transfer
   call tree_update(0,nsteps)

   return
end subroutine sink_transfer_particles
