! GET_GHOSTS.F90
! D. A. Hubber - 12/6/2008
! Finds potential neighbours of particles on local processor that overlap
! bounding box from other processors.  Records the position, mass and
! velocity of all potential neighbour ready to be sent back to the
! other processors.
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine get_ghosts(scatter)
   use mpi
   use mpi_communication_module
   use particle_module
   use type_module
#ifdef BH_TREE
   use tree_module
#endif
#ifdef DEBUG_GHOST_TX
   use filename_module, only : run_id, MPI_ext
#endif
   use domain_comparisons
   implicit none

! This subroutine identifies particles within the bounding boxes of other
! domains, and sends them ghost particles.
! For an SPH search, the ghost particles have position, mass and velocity.

   logical, intent(in) :: scatter   ! Do we want scatter neighbours?
   integer :: d                     ! Domain loop counter/identifier
   integer :: i, j                  ! Counter variables
   integer :: slot                  ! Slot value
   integer :: ierr                  ! Error value
   integer :: numsendto(0:lastrank)   ! Number sent to domain
   integer :: numrecvfrom(0:lastrank) ! Number received from each domain
   integer :: recv_num(0:endranklist) ! Num. of particles receive requests
   integer :: send_num(0:endranklist) ! Num. of particles send requests
   integer :: send(1:MAX_CHUNKS,0:endranklist)
                             ! Array of send requests (one for each OTHER task)
   integer :: recv(1:MAX_CHUNKS,0:endranklist)
                             ! Array of recv requests (one for each OTHER task)
   integer :: p, pp                 ! Particle counters
#ifdef DEBUG_GHOST_TX
   integer :: pcount                ! Number of particles being received
#endif
   integer :: pp_numb               ! Ghost number counter
   integer :: stat(MPI_STATUS_SIZE) ! Status variable
   integer :: tag                   ! Tag of message
   logical :: flag                  ! Logical flag
   integer :: nsend                 ! number of particles to send
   integer :: nsend_external        ! number of periodic ghosts to send
   integer :: porig                 ! Original particle of ghost
   integer :: chunk_id              ! Number of chunk to send/recv
   integer :: chunk_size_tot        ! Total size of chunk to send/recv
   integer :: chunk_size            ! Size of chunk to send/recv
   integer :: chunk_slot            ! Slot for chunk to send/recv
   integer, allocatable :: send_list(:) ! particles to send
   integer :: pghost_external       ! Number of external periodic ghosts

!    integer, allocatable       :: ghostd(:)         ! Number of ghosts in each domain

   type(ghostaccess), allocatable :: ghostsendarray(:) ! Send array of ghost particles
   type(ghostaccess), allocatable :: ghostrecvarray(:) ! Receive array of ghost particles
   type(particleaccess), allocatable :: ptclsendarray(:) ! Send array of full particles
   type(particleaccess), allocatable :: ptclrecvarray(:) ! Receive array of full particles

   real(kind=PR) :: bbmintemp(1:NDIM)
   real(kind=PR) :: bbmaxtemp(1:NDIM)

#ifdef BH_TREE
   integer :: c                         ! Cell number
   integer :: pp_pot                    ! Number of particles found in the tree walk
   integer :: treelist(1:ptot+pghost)   ! Results of the tree walk
#endif

#ifdef DEBUG_GHOST_TX
   character(len=100)    :: out_file ! filename extension for data output
#endif

   debug2("Getting ghost neighbours [get_ghosts.F90]")
   
!    allocate(ghostd(0:lastrank))
   if (scatter) then
      allocate(ptclsendarray(0:lastrank))
      allocate(ptclrecvarray(0:lastrank))
   else
      allocate(ghostsendarray(0:lastrank))
      allocate(ghostrecvarray(0:lastrank))
   end if
!    ghostd(0:lastrank) = 0
   numsendto = 0
   numrecvfrom = 0
   if (scatter) then
      tag = SCATTER_GHOST_TAG
   else
      tag = GHOST_TAG
   end if
   send = MPI_REQUEST_NULL
   recv = MPI_REQUEST_NULL
   pghost_external = 0


! Loop over all other domains
! ----------------------------------------------------------------------------
   do d=0,lastrank

      if (d == rank) cycle

#if defined(PERIODIC)
      pghost_external = size(ghost_sph_reference(d)%data)
#endif

#ifdef BH_TREE

      c = 0
      pp_pot = 0

! ----------------------------------------------------------------------------
      do
         bbmintemp(1:NDIM) = BHhydro(c)%r(1:NDIM) - BHhydro(c)%rmax
         bbmaxtemp(1:NDIM) = BHhydro(c)%r(1:NDIM) + BHhydro(c)%rmax

         if (scatter) then
            ! Include both particles within bbmin/max of other domain, and
            ! particles whose smoothing kernels overlap the activemin/max of
            ! the other domain
            flag = overlap_gather_scatter(bbmintemp,bbmaxtemp,BHhydro(c)%hrangemax,d)
         else
            ! Only include particles within bbmin
            flag = overlap_smoothing_box(bbmintemp,bbmaxtemp,d)
         end if

! Check if search sphere overlaps cell 'neighbour sphere'
! ----------------------------------------------------------------------------

         if (flag) then

! Check if c is a leaf cell or not
            if (BHhydro(c)%leaf > 0) then

! Add all particles in leaf cell to potential neighbour list
               do j=1,BHhydro(c)%leaf
                  pp_pot = pp_pot + 1
                  treelist(pp_pot) = BHhydro(c)%plist(j)
               end do
               c = BHhydro(c)%nextcell
            else if (BHhydro(c)%leaf == 0) then
               c = BHhydro(c)%ifopen
            else
               c = BHhydro(c)%nextcell
            end if

         else
            
            c = BHhydro(c)%nextcell

         end if
! ----------------------------------------------------------------------------

         if (c > ctot_hydro) exit
      end do
      
      
! ----------------------------------------------------------------------------
! #if defined(PERIODIC)
!       c = 0
! 
!       if (pghost > 0) then
!          do
!             bbmintemp(1:NDIM) = BHghost(c)%r(1:NDIM) - BHghost(c)%rmax
!             bbmaxtemp(1:NDIM) = BHghost(c)%r(1:NDIM) + BHghost(c)%rmax
! 
!             if (scatter) then
!                ! Include both particles within bbmin/max of other domain, and
!                ! particles whose smoothing kernels overlap the activemin/max of
!                ! the other domain
!                flag = overlap_gather_scatter(bbmintemp,bbmaxtemp,BHghost(c)%hrangemax,d)
!             else
!                ! Only include particles within bbmin
!                flag = overlap_smoothing_box(bbmintemp,bbmaxtemp,d)
!             end if
! 
! ! Check if search sphere overlaps cell 'neighbour sphere'
! ! ----------------------------------------------------------------------------
! 
!             if (flag) then
! 
! ! Check if c is a leaf cell or not
!                if (BHghost(c)%leaf > 0) then
! 
! ! Add all particles in leaf cell to potential neighbour list
!                   do j=1,BHghost(c)%leaf
!                      pp_pot = pp_pot + 1
!                      treelist(pp_pot) = BHghost(c)%plist(j)
!                   end do
!                   c = BHghost(c)%nextcell
!                else if (BHghost(c)%leaf == 0) then
!                   c = BHghost(c)%ifopen
!                else
!                   c = BHghost(c)%nextcell
!                end if
! 
!             else
!             
!                c = BHghost(c)%nextcell
! 
!             end if
! ! ----------------------------------------------------------------------------
! 
!             if (c > ctot_ghost) exit
!          end do
!       end if
! #endif
#endif

! Loop over all particles in current domain, check whether they are entirely
! within the domain bounding box, and then record required properties
      pp_numb = 0
#ifdef BH_TREE
      if (scatter) then
         allocate(ptclsendarray(d)%data(1:max(1,pp_pot+pghost_external)))
      else
         allocate(ghostsendarray(d)%data(1:max(1,pp_pot+pghost_external)))
      end if
      allocate(send_list(1:pp_pot+pghost_external))
      nsend = 0
      do j=1,pp_pot
         p = treelist(j)
#else
      if (scatter) then
         allocate(ptclsendarray(d)%data(1:ptot+pghost_external))
      else
         allocate(ghostsendarray(d)%data(1:ptot+pghost_external))
      end if
      allocate(send_list(1:ptot+pghost_external))
      nsend = 0
      do p=1,ptot
#endif
         if (scatter) then
            flag = inside_gather_scatter(sph(p)%r,KERNRANGE*sph(p)%h,d)
         else
            flag = inside_bounding_smoothing_box(sph(p)%r,d)
         end if
         if (flag) then
            nsend = nsend + 1
            send_list(nsend) = p
         end if
      end do
      do j=1,nsend
         p = send_list(j)
         pp_numb = pp_numb + 1
         if (scatter) then
            ptclsendarray(d)%data(pp_numb) = sph(p)
         else
            ghostsendarray(d)%data(pp_numb)%r = sph(p)%r
            ghostsendarray(d)%data(pp_numb)%m = sph(p)%m
            ghostsendarray(d)%data(pp_numb)%v = sph(p)%v
#ifdef DIV_A
            ghostsendarray(d)%data(pp_numb)%a = sph(p)%a
#endif
#if defined(SINKS) && defined(SELF_GRAVITY)
            ghostsendarray(d)%data(pp_numb)%gpot = sph(p)%gpot
#endif
#ifdef DEBUG_GHOST_TX
            ghostsendarray(d)%data(pp_numb)%porig = sph(p)%porig
#endif
            ghostsendarray(d)%data(pp_numb)%ptype = sph(p)%ptype
         end if
      end do
#if defined PERIODIC
      nsend_external = 0
      do p=1,pghost_external
         if (scatter) then
            porig = ghost_sph_reference(d)%data(p)%porig
            flag = inside_gather_scatter(ghost_sph_reference(d)%data(p)%r,&
               & KERNRANGE*sph(porig)%h,d)
         else
            flag = inside_bounding_smoothing_box(&
               & ghost_sph_reference(d)%data(p)%r,d)
         end if
         if (flag) then
            nsend_external = nsend_external + 1
            send_list(nsend + nsend_external) = p
         end if
      end do
      do j=1,nsend_external
         p = send_list(nsend + j)
         pp = ghost_sph_reference(d)%data(p)%porig
         pp_numb = pp_numb + 1
         if (scatter) then
            ptclsendarray(d)%data(pp_numb) = sph(pp)
            ptclsendarray(d)%data(pp_numb)%r = ghost_sph_reference(d)%data(p)%r
         else
            ghostsendarray(d)%data(pp_numb)%r = ghost_sph_reference(d)%data(p)%r
            ghostsendarray(d)%data(pp_numb)%m = sph(pp)%m
            ghostsendarray(d)%data(pp_numb)%v = sph(pp)%v
#ifdef DIV_A
            ghostsendarray(d)%data(pp_numb)%a = sph(pp)%a
#endif
#if defined(SINKS) && defined(SELF_GRAVITY)
            ghostsendarray(d)%data(pp_numb)%gpot = sph(pp)%gpot
#endif
#ifdef DEBUG_GHOST_TX
            ghostsendarray(d)%data(pp_numb)%porig = sph(pp)%porig
#endif
            ghostsendarray(d)%data(pp_numb)%ptype = sph(pp)%ptype
         end if
      end do
#endif
!       ghostd(d) = pp_numb
      numsendto(d) = pp_numb

#ifdef DEBUG_GHOST_TX
      write (6,*) "Rank ",rank," sending a message to task ", d
!       write (6,*) "Message length ", ghostd(d), " particles"
      write (6,*) "Message length ", numsendto(d), " particles"
#endif

      deallocate(send_list)
#if defined(PERIODIC)
      deallocate(ghost_sph_reference(d)%data)
#endif
   end do
  
  ! --------------------------------------------------------------------------
  ! Transmit/receive the number of particles to transfer to/from each domain
  i = 0
  do d=0,lastrank
     if (d==rank) cycle
     call MPI_ISEND(numsendto(d), 1, MPI_INTEGER, &
        &d, GHOST_NUM_TAG, MPI_COMM_WORLD, send_num(i), ierr)
     i = i + 1
  end do
  
  i = 0
  do d=0,lastrank
     if (d==rank) cycle
     call MPI_IRECV(numrecvfrom(d), 1, MPI_INTEGER, &
        &d, GHOST_NUM_TAG, MPI_COMM_WORLD, recv_num(i), ierr)
     i = i + 1
  end do

  ! --------------------------------------------------------------------------
  
   i = 0
   do d=0,lastrank
      if (d==rank) cycle

! Send ghost particles for this domain
      chunk_id = 1
      chunk_slot = 1
      chunk_size_tot = numsendto(d)
      do
         if (numsendto(d) == 0) exit
         chunk_size = min(chunk_size_tot, max_message_size)
         if (scatter) then
            call MPI_ISEND(ptclsendarray(d)%data(chunk_slot), chunk_size, &
               & MPI_PARTICLE, d, tag, MPI_COMM_WORLD, send(chunk_id,i), ierr)
         else
            call MPI_ISEND(ghostsendarray(d)%data(chunk_slot), chunk_size, &
               & MPI_GHOST_TYPE, d, tag, MPI_COMM_WORLD, send(chunk_id,i), ierr)
         end if
         chunk_size_tot = chunk_size_tot - chunk_size
         chunk_slot = chunk_slot + chunk_size
         if (chunk_size_tot == 0) exit
         chunk_id = chunk_id + 1
         if (chunk_id > MAX_CHUNKS) then
            write (6,*) "Too many chunks! get_ghosts"
            stop 1
         end if
      end do
      i = i + 1

   end do

#ifdef DEBUG_GHOST_TX
   out_file = trim(adjustl(run_id))//".ghostsend."//trim(adjustl(MPI_ext))
   open(1,file=out_file,status='unknown',form='formatted')
   do d=0,lastrank
      write (1,*) " ------------------ DOMAIN ", d, " ------------------"
      if (rank==d) cycle
      do p=1,numsendto(d)
         if (scatter) then
            write(1,*) p, ptclsendarray(d)%data(1:NDIM,p)
         else
            write(1,*) p, ghostsendarray(d)%data(1:NDIM,p)
         end if
      end do
   end do
   close(1)
#endif

! Post receives to receive ghost particles from other domains
! ----------------------------------------------------------------------------

!    ghostd(0:lastrank) = 0

! Receive from each domain in turn into the array ghostrecvarray

! Go through domains in the order their messages arrive
   do i=0,lastrank-1

! Test for any message - wait until one arrives (blocking probe)
      WAIT_TIME_MACRO
!       call MPI_PROBE(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, stat, ierr)
      call MPI_WAITANY(lastrank, recv_num, p, stat, ierr)
      CALC_TIME_MACRO
      d = stat(MPI_SOURCE)
!       if (scatter) then
!          call MPI_GET_COUNT(stat,MPI_PARTICLE,pcount,ierr)
!       else
!          call MPI_GET_COUNT(stat,MPI_GHOST_TYPE,pcount,ierr)
!       end if
!       ghostd(d) = pcount

#ifdef DEBUG_GHOST_TX
      write (6,*) "Rank ",rank," receiving a message from task ", d
      if (scatter) then
         call MPI_GET_COUNT(stat,MPI_PARTICLE,pcount,ierr)
      else
         call MPI_GET_COUNT(stat,MPI_GHOST_TYPE,pcount,ierr)
      end if
      write (6,*) "Message length ", pcount, " or ", numrecvfrom(d), &
                  &" particles"
      if (pcount /= numrecvfrom(d)) then
         write (6,*) "Message size does not match that expected!"
         stop 1
      end if
#endif

      if (scatter) then
         allocate(ptclrecvarray(d)%data(1:max(1,numrecvfrom(d))))
      else
         allocate(ghostrecvarray(d)%data(1:max(1,numrecvfrom(d))))
      end if

! Post the receive - until this is posted the receive does not actually occur
! and may or may not be buffered
      chunk_id = 1
      chunk_slot = 1
      chunk_size_tot = numrecvfrom(d)
      do
         if (numrecvfrom(d) == 0) exit
         chunk_size = min(chunk_size_tot, max_message_size)
#ifdef DEBUG_TRANSFER
         write (6,*) rank, ": Chunk=", chunk_id, "; chunk_slot=", chunk_slot, &
                     & "; chunk_size=", chunk_size
#endif
         if (scatter) then
            call MPI_IRECV(ptclrecvarray(d)%data(chunk_slot), chunk_size,&
               &MPI_PARTICLE, d, tag, MPI_COMM_WORLD, recv(chunk_id,i), ierr)
         else
            call MPI_IRECV(ghostrecvarray(d)%data(chunk_slot), chunk_size, &
               &MPI_GHOST_TYPE, d, tag, MPI_COMM_WORLD, recv(chunk_id,i), ierr)
         end if
         chunk_size_tot = chunk_size_tot - chunk_size
         chunk_slot = chunk_slot + chunk_size
         if (chunk_size_tot == 0) exit
         chunk_id = chunk_id + 1
      end do

   end do

  ! Check we have enough space on the end of the main array
   if ((sum(numrecvfrom)+ptot+pghost) > pmax) then
      write (6,*) "expanding main array to store ", sum(numrecvfrom)+pghost, &
         &" extra ghosts..."
      call expand(ceiling(real(PMAXMULT,DP)*real(ptot+pghost+sum(numrecvfrom))))
#if defined(BH_TREE)
      call expand_ghost_tree(pghost+sum(numrecvfrom))
#endif
   end if

! Wait for all sends and receives to be completed
! ----------------------------------------------------------------------------

! All our receives are now completed
   WAIT_TIME_MACRO
   call MPI_WAITALL(lastrank*MAX_CHUNKS,recv,MPI_STATUSES_IGNORE,ierr)
   CALC_TIME_MACRO
   slot = ptot+pghost+1
   do d=0,lastrank
      if (d==rank) cycle
      do i=1,numrecvfrom(d)
         if (scatter) then
            sph(slot) = ptclrecvarray(d)%data(i)
            sph(slot)%porig = 0
         else
            sph(slot)%r = ghostrecvarray(d)%data(i)%r
            sph(slot)%m = ghostrecvarray(d)%data(i)%m
            sph(slot)%v = ghostrecvarray(d)%data(i)%v
#ifdef DIV_A
            sph(slot)%a = ghostrecvarray(d)%data(i)%a
#endif
#if defined(SINKS) && defined(GRAVITY)
            sph(slot)%gpot = ghostrecvarray(d)%data(i)%gpot
#endif
#ifdef DEBUG_GHOST_TX
            sph(slot)%porig = ghostrecvarray(d)%data(i)%porig
#else
            sph(slot)%porig = 0
#endif
            sph(slot)%ptype = ghostrecvarray(d)%data(i)%ptype
         end if
         slot = slot + 1
      end do
   end do
   
   ! Update number of ghosts
   pghost = pghost + sum(numrecvfrom)

   if (scatter) then
      deallocate(ptclrecvarray)
   else
      deallocate(ghostrecvarray)
   end if

#ifdef DEBUG_GHOST_TX
   out_file = trim(adjustl(run_id))//".ghost."//trim(adjustl(MPI_ext))
   open(1,file=out_file,status='unknown',form='formatted')
   i = ptot+1
   do d=0,lastrank
      write (1,*) " ------------------ DOMAIN ", d, " ------------------"
      if (rank==d) cycle
      do p=i,i+numrecvfrom(d)-1
        write(1,*) p, sph(p)%r, sph(p)%porig
      end do
      i = i + 1
   end do
   close(1)
#endif

   WAIT_TIME_MACRO
   call MPI_WAITALL(lastrank,send_num,MPI_STATUSES_IGNORE,ierr)
   call MPI_WAITALL(lastrank*MAX_CHUNKS,send,MPI_STATUSES_IGNORE,ierr)
   CALC_TIME_MACRO

! Can't do this until all sends completed
   if (scatter) then
      deallocate(ptclsendarray)
   else
      deallocate(ghostsendarray)
   end if

!    deallocate(ghostd)

   return
end subroutine get_ghosts
