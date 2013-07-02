! GET_POT_LIST_BB.F90
! D. A. Hubber - 12/6/2008
! Finds potential neighbours of particles on local processor that overlap
! bounding box from other processors.  Records the position, mass and
! velocity of all potential neighbour ready to be sent back to the
! other processors.
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine get_pot_list_bb
  use mpi
  use mpi_communication_module
  use particle_module
  use type_module
  use domain_comparisons
#ifdef BH_TREE
  use tree_module
#endif
#ifdef DEBUG_GHOST_TX
  use filename_module, only : run_id, MPI_ext
#endif
  implicit none

! This subroutine identifies particles within the bounding boxes of other
! domains, and sends them ghost particles.
! For an SPH search, the ghost particles have position, mass and velocity.

  logical :: flag                  ! Logical flag
  integer :: d                     ! Domain loop counter/identifier
  integer :: i                     ! Counter variable
  integer :: slot                  ! Slot value
  integer :: ierr                  ! Error value
  integer :: send(0:endranklist)   ! Array of send requests (one for each OTHER task)
  integer :: recv(0:endranklist)   ! Array of recv requests (one for each OTHER task)
  integer :: p                     ! Particle counter
  integer :: pcount                ! Number of particles being received
  integer :: pp_numb               ! Ghost number counter
  integer :: stat(MPI_STATUS_SIZE) ! Status variable
  integer :: tag                   ! Tag of message
  integer :: nsend                 ! number of particles to send
  integer, allocatable :: ghostd(:)     ! Number of ghosts in each domain
  integer, allocatable :: send_list(:)  ! particles to send
  real(kind=PR) :: bbmintemp(1:NDIM)    ! ..
  real(kind=PR) :: bbmaxtemp(1:NDIM)    ! ..
  type(ghostaccess), allocatable :: sendarray(:)    ! Send array of ghost particles
  type(ghostaccess), allocatable :: receivearray(:) ! Receive array of ghost particles
#ifdef BH_TREE
  integer :: c                     ! Cell number
  integer :: j                     ! Loop counter
  integer :: pp_pot                ! Number of particles found in the tree walk
  integer :: treelist(1:ptot)      ! Results of the tree walk
#endif
#ifdef DEBUG_GHOST_TX
  character(len=100)    :: out_file ! filename extension for data output
#endif

  debug_timing("GET_POT_LIST_BB")
  debug2("Getting ghost neighbours [get_pot_list_bb.F90]")

  allocate(ghostd(0:lastrank))
  allocate(sendarray(0:lastrank))
  allocate(receivearray(0:lastrank))
  ghostd(0:lastrank) = 0
#ifdef PUSH_RHO
  ghost_nfrom(0:lastrank) = 0
#endif
  pghost = 0
  i = 0
  tag = GHOST_TAG
#ifdef PUSH_RHO
  allocate(ghost_send_list(0:lastrank))
#endif

! Loop over all other domains
! ----------------------------------------------------------------------------
  do d=0,lastrank

     if (d == rank) cycle

#ifdef BH_TREE

     c = 0
     pp_pot = 0

! ----------------------------------------------------------------------------
     do
        bbmintemp(1:NDIM) = BHhydro(c)%r(1:NDIM) - BHhydro(c)%rmax
        bbmaxtemp(1:NDIM) = BHhydro(c)%r(1:NDIM) + BHhydro(c)%rmax

        ! this is OK as particles here will not actually be over the periodic edge
        flag = overlap_smoothing_box(bbmintemp,bbmaxtemp,d)

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

        if (c > cmax_hydro) exit
     end do
#endif

! Loop over all particles in current domain, check whether they are entirely
! within the domain bounding box, and then record required properties
     pp_numb = 0
#ifdef BH_TREE
     allocate(sendarray(d)%data(1:max(1,pp_pot)))
     allocate(send_list(1:pp_pot))
     nsend = 0
     do j=1,pp_pot
        p = treelist(j)
#else
     allocate(sendarray(d)%data(1:ptot))
     allocate(send_list(1:ptot))
     nsend = 0
     do p=1,ptot
#endif
        flag = in_bounding_smoothing_box(sph(p)%r,d)
        if (flag) then
           nsend = nsend + 1
           send_list(nsend) = p
        end if
     end do
     do j=1,nsend
        p = send_list(j)
        pp_numb = pp_numb + 1
        sendarray(d)%data(pp_numb)%r = sph(p)%r
        sendarray(d)%data(pp_numb)%m = sph(p)%m
        sendarray(d)%data(pp_numb)%v = sph(p)%v
#ifdef DIV_A
        sendarray(d)%data(pp_numb)%a = sph(p)%a
#endif
#if defined(SINKS) && defined(GRAVITY)
        sendarray(d)%data(pp_numb)%gpot = sph(p)%gpot
#endif
#ifdef DEBUG_GHOST_TX
        sendarray(d)%data(pp_numb)%porig = sph(p)%porig
#endif
        sendarray(d)%data(pp_numb)%ptype = sph(p)%ptype
     end do
     ghostd(d) = pp_numb

#ifdef DEBUG_GHOST_TX
    write (6,*) "Rank ",rank," sending a message to task ", d
    write (6,*) "Message length ", ghostd(d), " particles"
#endif

#ifdef PUSH_RHO
     ghost_nsend(d) = nsend
     allocate(ghost_send_list(d)%data(nsend))
     ghost_send_list(d)%data(1:nsend) = send_list(1:nsend)
#endif
     deallocate(send_list)

! Send ghost particles for this domain
     call MPI_ISEND(sendarray(d)%data(1), ghostd(d), MPI_GHOST_TYPE, &
          & d, tag, MPI_COMM_WORLD, send(i), ierr)
     i = i + 1

  end do

#ifdef DEBUG_GHOST_TX
  out_file = trim(adjustl(run_id))//".ghostsend."//trim(adjustl(MPI_ext))
  open(1,file=out_file,status='unknown',form='formatted')
  do d=0,lastrank
    write (1,*) " ------------------ DOMAIN ", d, " ------------------"
    if (rank==d) cycle
    do p=1,ghostd(d)
      write(1,*) p, sendarray(d)%data(1:NDIM,p)
    end do
  end do
  close(1)
#endif

! Post receives to receive ghost particles from other domains
! ----------------------------------------------------------------------------

! Change meaning of pp_numb to a position of the first particle in the current domain
  ghostd(0:lastrank) = 0

! Receive from each domain in turn into the array receivearray

! Go through domains in the order their messages arrive
! ---------------------------------------------------------------------------
  do i=0,lastrank-1

! Test for any message - wait until one arrives (blocking probe)
    WAIT_TIME_MACRO
    call MPI_PROBE(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, stat, ierr)
    CALC_TIME_MACRO
    d = stat(MPI_SOURCE)
    call MPI_GET_COUNT(stat,MPI_GHOST_TYPE,pcount,ierr)
    ! If doing a sink search, check if we are receiving particles from
    ! this task
    ghostd(d) = pcount
    pghost = pghost + ghostd(d)

#ifdef DEBUG_GHOST_TX
    write (6,*) "Rank ",rank," receiving a message from task ", d
    write (6,*) "Message length ", pcount, ", or ", ghostd(d), " particles"
#endif

    allocate(receivearray(d)%data(1:max(1,ghostd(d))))

! Post the receive - until this is posted the receive does not actually occur
! and may or may not be buffered
    call MPI_IRECV(receivearray(d)%data(1),pcount,MPI_GHOST_TYPE,&
         &d,tag,MPI_COMM_WORLD,recv(i),ierr)

  end do
! ---------------------------------------------------------------------------

  ! Check we have enough space on the end of the main array
  if ((pghost+ptot) > pmax) then
    call expand(ceiling(real(PMAXMULT,DP)*real(ptot+pghost)))
  end if

! Wait for all sends and receives to be completed
! ----------------------------------------------------------------------------

! All our receives are now completed
  WAIT_TIME_MACRO
  call MPI_WAITALL(lastrank,recv,MPI_STATUSES_IGNORE,ierr)
  CALC_TIME_MACRO
  slot = ptot+1
  do d=0,lastrank
    if (d==rank) cycle
    do i=1,ghostd(d)
       sph(slot)%r = receivearray(d)%data(i)%r
       sph(p)%m = receivearray(d)%data(i)%m
       sph(p)%v = receivearray(d)%data(i)%v
#ifdef DIV_A
       sph(p)%a = receivearray(d)%data(i)%a
#endif
#if defined(SINKS) && defined(GRAVITY)
       sph(p)%gpot = receivearray(d)%data(i)%gpot
#endif
#ifdef DEBUG_GHOST_TX
       sph(p)%porig = receivearray(d)%data(i)%porig
#endif
       sph(p)%ptype = receivearray(d)%data(i)%ptype
       slot = slot + 1
    end do
    slot = slot + ghostd(d)
  end do

#ifdef PUSH_RHO
  ghost_nfrom = ghostd
#endif

  deallocate(receivearray)

#ifdef DEBUG_GHOST_TX
  out_file = trim(adjustl(run_id))//".ghost."//trim(adjustl(MPI_ext))
  open(1,file=out_file,status='unknown',form='formatted')
  i = ptot+1
  do d=0,lastrank
    write (1,*) " ------------------ DOMAIN ", d, " ------------------"
    if (rank==d) cycle
    do p=i,i+ghostd(d)-1
      write(1,*) p, sph(p)%r, sph(p)%porig
    end do
    i = i + 1
  end do
  close(1)
#endif

  WAIT_TIME_MACRO
  call MPI_WAITALL(lastrank,send,MPI_STATUSES_IGNORE,ierr)
  CALC_TIME_MACRO

! Can't do this until all sends completed
  deallocate(sendarray)

  deallocate(ghostd)

  return
end subroutine get_pot_list_bb

#ifdef PUSH_RHO
! ============================================================================
! PUSH_GHOST_RHO
! ..
! ============================================================================
subroutine push_ghost_rho
  use mpi
  use mpi_communication_module
  use hydro_module, only : rho
  use domain_comparisons
  implicit none

  integer :: d                     ! Domain loop counter/identifier
  integer :: i, j                  ! Counter variables
  integer :: ierr                  ! Error value
  integer :: send(0:endranklist)   ! Array of send requests (one for each OTHER task)
  integer :: recv(0:endranklist)   ! Array of recv requests (one for each OTHER task)
  integer :: sends                 ! number of actual sends
  integer :: recvs                 ! number of actual recvs
  integer :: p                     ! Particle counter
  integer :: pcount                ! Number of particles being received
  integer :: stat(MPI_STATUS_SIZE) ! Status variable
  integer :: tag                   ! Tag of message
  integer, allocatable           :: ghostd(:)       ! Number of ghosts in each domain
  type(real_access), allocatable :: sendarray(:)    ! Send array of ghost particles
  type(real_access), allocatable :: receivearray(:) ! Receive array of ghost particles

  sends = 0
  tag = GHOST_RHO_TAG
  allocate(sendarray(0:lastrank))
  allocate(receivearray(0:lastrank))
  allocate(ghostd(0:lastrank))

  ! Send new densities
  ! -------------------------------------------------------------------------
  do d=0,lastrank
     if (d==rank) cycle
     if (ghost_nsend(d)==0) cycle

     ! Load data into sendarray
     allocate(sendarray(d)%data(1:ghost_nsend(d)))
     do j=1,ghost_nsend(d)
        p = ghost_send_list(d)%data(j)
        sendarray(d)%data(j) = sph(p)%rho
     end do

     call MPI_ISEND(sendarray(d)%data(1), ghost_nsend(d), MPI_REAL_PR, &
          & d, tag, MPI_COMM_WORLD, send(sends), ierr)

     sends = sends + 1
  end do
  ! -------------------------------------------------------------------------

  ! Receive new densities
  recvs = 0
  do d=0,lastrank
     if (d==rank) cycle
     if (ghost_nfrom(d)==0) cycle
     allocate(receivearray(d)%data(1:ghost_nfrom(d)))

     call MPI_IRECV(receivearray(d)%data(1),ghost_nfrom(d),MPI_REAL_PR,&
        &d,tag,MPI_COMM_WORLD,recv(recvs),ierr)

     recvs = recvs + 1
  end do

  ! Wait for receives to finish
  WAIT_TIME_MACRO
  call MPI_WAITALL(recvs,recv(0:recvs-1),MPI_STATUSES_IGNORE,ierr)
  CALC_TIME_MACRO

  ! Load data into ghost_rho_array
  i = 1
  do d=0,lastrank
    if (d==rank) cycle
    if (ghost_nfrom(d)==0) cycle
    ghost_rho_array(i:i+ghost_nfrom(d)-1) = receivearray(d)%data(1:ghost_nfrom(d))
    i = i + ghost_nfrom(d)
  end do

  ! Wait for sends to finish
  WAIT_TIME_MACRO
  call MPI_WAITALL(sends,send(0:sends-1),MPI_STATUSES_IGNORE,ierr)
  CALC_TIME_MACRO

  ! Clean up
  deallocate(sendarray)
  deallocate(receivearray)
  deallocate(ghost_send_list)

  return
end subroutine push_ghost_rho
#endif
