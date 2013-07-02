! FORCES.F90
! C. P. Batty & D. A. Hubber - 11/1/2007
! Calculates accelerations (hydrodynamic and gravitational) for all particles
! - controls calls to relevant subroutines
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE gravity_export_to(send,recv,acclist,acctot,chunksize)
   use type_module
   use particle_module
   use time_module
#ifdef SINKS
   use sink_module
#endif
   use mpi_communication_module
   use hydro_module
   use mpi
#ifdef DEBUG_EXPORT
   use filename_module, only : run_id, MPI_ext
#endif
   implicit none

   integer, intent(out) :: send(0:endranklist)
                  ! Array of send requests (one for each OTHER task)
   integer, intent(out) :: recv(0:endranklist)
                  ! Array of recv requests (one for each OTHER task)
   integer, intent(in) :: acctot                   ! No. of particles on acc. step
   integer, intent(in) :: acclist(1:acctot)        ! List of particles on acc. step
   integer, intent(in) :: chunksize                ! Data packet size for dynamic OpenMP
   integer :: d                     ! domain counter
   integer :: countsend(0:lastrank) ! Number of particles to each domain
   integer :: slot
   integer :: i, j                  ! loop counter
   integer :: p                     ! particle counter
#ifdef SINKS
   integer :: s                     ! sink counter
#endif
   integer :: stat(MPI_STATUS_SIZE) ! MPI Status variable
   integer :: ierr                  ! Return code
   integer :: pcount                ! Number of particles received per receive
#ifdef BH_TREE
   logical :: export                ! Whether a particle needs exporting
   real(kind=DP) :: agravp(1:NDIM)  ! Gravitational acceleration
   real(kind=DP) :: potp            ! Gravitational potential
#endif

#ifdef DEBUG_EXPORT
   character(len=100)    :: out_file ! filename extension for data output
#endif

   debug2("Calculating remote gravity or exporting particles [gravity_export_to.F90]")
   debug_timing("GRAVITY_EXPORT_TO")

  ! Zero all active sinks' 'export' gravity contribution
#if defined(SINKS)
   if (accdo_sinks) then
      do s=1,stot
         if (sink(s)%domain /= rank) cycle
         sink(s)%remote_agrav = 0._DP
         sink(s)%remote_gpot = 0._DP
      end do
   end if
#endif

! Work out which particles smoothing kernels overlap other domain's
! bounding boxes
#ifdef SINKS
   allocate(exportlist(1:ptot+stot,0:lastrank))
#else
   allocate(exportlist(1:ptot,0:lastrank))
#endif
   countsend = 0
   slot = 0

#ifdef BH_TREE
   debug_timing("REMOTE_GRAVITY")
   do i=0,lastrank-1
      WAIT_TIME_MACRO
      call MPI_WAITANY(lastrank, treerecv(0:), pcount, stat, ierr)
      CALC_TIME_MACRO
      d = stat(MPI_SOURCE)
#ifdef DEBUG_EXPORT
      write (6,*) "Rank ", rank, "has completed receiving a pruned tree from&
      & task ", d
      call MPI_GET_COUNT(stat,MPI_REMOTE_BH_GRAV,pcount,ierr)
      write (6,*) "It contains ", pcount, " cells."
#else
      call MPI_GET_COUNT(stat,MPI_REMOTE_BH_GRAV,pcount,ierr)
#endif
      ! We have finished receiving a pruned tree, we can now do gravity walks over it
      if (acctot > 0) then
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,chunksize) DEFAULT(SHARED) PRIVATE(p,agravp,potp,export)
         do j=1,acctot
            p = acclist(j)
            call BHtree_remote_grav(p,d,sph(p)%r,agravp,potp,export)
            if (export) then
               ! Throw away results, add to export list
!$OMP CRITICAL
               countsend(d) = countsend(d) + 1
               exportlist(countsend(d),d) = p
!$OMP END CRITICAL
            else
               ! Store gravitational acceleration
               sph(p)%a     = sph(p)%a + agravp(1:VDIM)
               sph(p)%gpot  = sph(p)%gpot  + potp
            end if
         end do
!$OMP END PARALLEL DO
      end if
#ifdef SINKS
      ! Now try for local sinks
      if (accdo_sinks) then
         do s=1,stot
            if (rank /= sink(s)%domain) cycle
            call BHtree_remote_grav(-s,d,sink(s)%r,agravp,potp,export)
            if (export) then
               ! Throw away results, add to export list
               countsend(d) = countsend(d) + 1
               exportlist(countsend(d),d) = -s
            else
               ! Store gravitational acceleration
               sink(s)%remote_agrav = sink(s)%remote_agrav + real(agravp(1:NDIM),PR)
               sink(s)%remote_gpot = sink(s)%remote_gpot + real(potp,PR)
            end if
         end do
      end if

   debug_timing("GRAVITY_EXPORT_TO")
#endif
   end do

#else
   do d=0,lastrank
      if (d==rank) cycle
      ! Not using the pruned trees, export all particles to all domains
      do p=pgravitystart,ptot
         if (.NOT. accdo(p)) cycle
         ! Add identifier of particle to list to add to export list later
         countsend(d) = countsend(d) + 1
         exportlist(countsend(d),d) = p
      end do
#ifdef SINKS
      ! Now try for local sinks
      if (accdo_sinks) then
         do s=1,stot
            if (rank /= sink(s)%domain) cycle
            ! Add identifier of sink to list to add to export list later
            countsend(d) = countsend(d) + 1
            exportlist(countsend(d),d) = -s
         end do
      end if
#endif
   end do
#endif

#ifdef DEBUG_EXPORT
   write (6,*) rank, ": sending this many: ", countsend
#endif
   allocate(sendgravity(1:sum(countsend)+1)) ! +1 so that last MPI send doesn't go
                                             ! out of range on zero length message

#ifdef DEBUG_EXPORT
   out_file = trim(adjustl(run_id))//".exportgravto."//trim(adjustl(MPI_ext))
   OPEN (1,file=out_file,status='unknown',form='formatted')
#endif

   do d=0,lastrank
      if (d==rank) cycle
#ifdef DEBUG_EXPORT
      write (1,*) "-------------------------------------------------------"
      write (1,*) "             From this task to task: ", d
      write (1,*) "-------------------------------------------------------"
#endif
      do i=1,countsend(d)
         p = exportlist(i,d)
         slot = slot + 1
#ifdef SINKS
         if (p > 0) then
#endif
#ifdef DEBUG_EXPORT
            write (1,*) "Particle ", p, " going in slot ", slot
#endif
            sendgravity(slot)%p = p
#ifdef DEBUG_EXPORT
            sendgravity(slot)%d = rank
#endif
            sendgravity(slot)%r = sph(p)%r
            sendgravity(slot)%h = sph(p)%h
#ifdef GRAD_H_SPH
            sendgravity(slot)%zo = sph(p)%zo
#endif
#if defined(BH_TREE) && !defined(GEOMETRIC_MAC)
            sendgravity(slot)%agravmag = sph(p)%agravmag
#endif
#ifdef DEBUG_EXPORT
            write (1,*) sendgravity(slot)
#endif
#ifdef SINKS
         else
            ! Sink particle
            s = -p
#ifdef DEBUG_EXPORT
            write (1,*) "Sink ", s, " going in slot ", slot
#endif
            sendgravity(slot)%p = p
#ifdef DEBUG_EXPORT
            sendgravity(slot)%d = rank
#endif
            sendgravity(slot)%r = sink(s)%r
            sendgravity(slot)%h = sink(s)%h
#ifdef GRAD_H_SPH
            sendgravity(slot)%zo = 0._PR
#endif
#if defined(BH_TREE) && !defined(GEOMETRIC_MAC)
            sendgravity(slot)%agravmag = sink(s)%agravmag
#endif
#if defined(BH_TREE) && defined(EIGEN_MAC)
            sendgravity(slot)%gpot = sink(s)%gpot
#endif
#ifdef DEBUG_EXPORT
            write (1,*) sendgravity(slot)
#endif
         end if
#endif
      end do
   end do

   deallocate(exportlist)

#ifdef DEBUG_EXPORT
   CLOSE (1)
#endif

! Export particles
   slot = 1
   i = 0

   do d=0,lastrank
      if (d==rank) cycle

#ifdef DEBUG_EXPORT
      write (6,*) "Rank ",rank," sending a export of gravity particles to task ", d
      write (6,*) "Message length ", countsend(d), " from slot ", slot
#endif

      call MPI_ISEND(sendgravity(slot), countsend(d), MPI_GRAVITYTYPE, &
         & d, GRAV_EXPORT_TAG, MPI_COMM_WORLD,send(i),ierr)
      slot = slot + countsend(d)
      i = i + 1
   end do

! Receive exported particles
   allocate(receivegravity(0:endranklist))
   pexportgrav = 0

   do i=0,lastrank-1
      ! Test for any message - wait until one arrives (blocking probe)
      WAIT_TIME_MACRO
      call MPI_PROBE(MPI_ANY_SOURCE, GRAV_EXPORT_TAG, MPI_COMM_WORLD, stat, ierr)
      CALC_TIME_MACRO
      d = stat(MPI_SOURCE)
      call MPI_GET_COUNT(stat,MPI_GRAVITYTYPE,pcount,ierr)

#ifdef DEBUG_EXPORT
      write (6,*) "Rank ",rank," receiving exported gravity particles from task ", d
      write (6,*) "Message length ", pcount, ", receiving to slot ", i
#endif

      allocate(receivegravity(i)%data(1:pcount+1)) ! +1 for 0 length send

      ! Record which domain and how many particles
      grav_fromeach(1,i) = d
      grav_fromeach(2,i) = pcount

! Post the receive - until this is posted the receive does not actually occur
! and may or may not be buffered
      call MPI_IRECV(receivegravity(i)%data(1),pcount,MPI_GRAVITYTYPE,&
         &d,GRAV_EXPORT_TAG,MPI_COMM_WORLD,recv(i),ierr)
      pexportgrav = pexportgrav + pcount

   end do

   return
END SUBROUTINE gravity_export_to
