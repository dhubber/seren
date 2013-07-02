! TRANSFER_PARTICLES.F90
! A. McLeod - 31/07/08
! ...
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine transfer_particles()
  use neighbour_module, only : pp_gather
  use particle_module
  use mpi_communication_module
  use time_module
  use mpi
  use domain_comparisons
  use type_module
#ifdef DEBUG_TRANSFER
  use filename_module
#endif
  implicit none

  integer :: i                             ! Loop counters
  integer :: ierr                          ! MPI error value
  integer :: d                             ! Domain counter
  integer :: p                             ! Particle counter
  integer :: pleft                         ! Number of particles left
  integer :: numsendto(0:lastrank)         ! Number sent to domain
  integer :: numrecvfrom(0:lastrank)       ! Number received from each domain
  integer :: numdel                        ! Number to delete
  integer :: recv_num(0:endranklist)       ! Num. of particles receive requests
  integer :: send_num(0:endranklist)       ! Num. of particles send requests
  integer :: pcount(0:endranklist)         ! Particles in receive_particles
  integer :: recv(1:MAX_CHUNKS,0:endranklist)   ! Array of receive requests
  integer :: send(1:MAX_CHUNKS,0:endranklist)   ! Array of send requests
  integer :: slots(0:lastrank)             ! Slots
  integer :: stat(MPI_STATUS_SIZE)         ! MPI Status variable
  integer :: chunk_id                      ! Number of chunk to send/recv
  integer :: chunk_size_tot                ! Total size of chunk to send/recv
  integer :: chunk_size                    ! Size of chunk to send/recv
  integer :: chunk_slot                    ! Slot for chunk to send/recv
  integer, allocatable :: deletelist(:)    ! List of ptcls to delete
  integer :: p_est                         ! Estimated particle number
  real(kind=DP) :: p_est_max               ! Estimated space required
#ifdef DEBUG_TRANSFER
  character(len=256) :: out_file           ! ..
#endif

  ! If there's only one MPI task, then return immediately
  if (numtasks==1) return

  debug_timing("TRANSFER_PARTICLES SEND")
  allocate(exportlist(1:ptot,0:lastrank))
  numsendto = 0
  numrecvfrom = 0
  numdel = 0
  pleft = ptot
  send = MPI_REQUEST_NULL
  recv = MPI_REQUEST_NULL

#ifdef DEBUG_TRANSFER
  out_file = trim(adjustl(run_id))//".sending."//trim(adjustl(MPI_ext))
  open(1,file=out_file,position='APPEND',form='formatted')
  write (1,*) " ==================== ", nsteps, " ===================="
  write (1,*) "    ptot = ", ptot
  write (6,*) rank, ": Domain boundary is: ",&
           &domain_bbmin(1:NDIM,rank), domain_bbmax(1:NDIM,rank)
#endif

  ! Go through each particle; test whether it has left the domain
  ! --------------------------------------------------------------------------
  ploop: do p=1,ptot

     ! Test if the particle is outside our own domain:
     if (.NOT. inside_domain_box(sph(p)%r,rank)) then
#ifdef DEBUG_TRANSFER
        write (6,*) rank, ": Particle ",p," r = ",&
        &sph(p)%r, "; outside domain boundary!"
#endif
        ! Prevent us from having less than pp_gather+1 particles in normal SPH
        ! or less than 2*pp_gather particles in grad-h SPH
#ifdef GRAD_H_SPH
        if (pleft <= 2*pp_gather) cycle
#else
        if (pleft <= pp_gather+1) cycle
#endif
        pleft = pleft - 1
        ! Find which domain the particle came from
        do d=0,lastrank
           if (d==rank) cycle
           ! Test if the particle is inside domain d
           if (inside_domain_box(sph(p)%r,d)) then
#ifdef DEBUG_TRANSFER
              write (6,*) rank, ": sending to domain ",d,": ", &
                       &domain_bbmin(1:NDIM,d), domain_bbmax(1:NDIM,d)
#endif
              ! Add to a list to send to domain d
              numsendto(d) = numsendto(d) + 1
              exportlist(numsendto(d),d) = p
              cycle ploop
           end if
        end do
        write (6,*) " ************************* STOP *************************"
        write (6,*) "In task ", rank
        write (6,*) "Particle ", p, " outside domain boundary but not sent anywhere!"
        write (6,*) "Particle is at position ", sph(p)%r
        write (6,*) "Domain boundaries:"
        do d=0,lastrank
           write (6,*) "Domain ", d
           write (6,*) "min: ", domain_bbmin(1:NDIM,d)
           write (6,*) "max: ", domain_bbmax(1:NDIM,d)
        end do
        stop
     end if
  end do ploop

  ! --------------------------------------------------------------------------
  ! Transmit/receive the number of particles to transfer to/from each domain
  i = 0
  do d=0,lastrank
     if (d==rank) cycle
     call MPI_ISEND(numsendto(d), 1, MPI_INTEGER, &
        &d, TRANSFER_NUM_TAG, MPI_COMM_WORLD, send_num(i), ierr)
     i = i + 1
  end do
  
  i = 0
  do d=0,lastrank
     if (d==rank) cycle
     call MPI_IRECV(numrecvfrom(d), 1, MPI_INTEGER, &
        &d, TRANSFER_NUM_TAG, MPI_COMM_WORLD, recv_num(i), ierr)
     i = i + 1
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
     ! Slot domain d starts at slots(d) in sendparticles
     slots(d) = slots(d-1) + numsendto(d-1)
  end do

  ! --------------------------------------------------------------------------
  do d=0,lastrank
    if (d==rank) cycle
    i = slots(d)
    do p=1,numsendto(d)
      sendparticles(i) = sph(exportlist(p,d))
#ifdef DEBUG_TRANSFER
      write (6,*) rank, ": packing particle ",p," to sendparticles slot ", i
#endif
      i = i + 1
      numdel = numdel + 1
      deletelist(numdel) = exportlist(p,d)
    end do
  end do
  ! --------------------------------------------------------------------------

  ! If we have more than one sent particle, sort list of particles to delete
  if (numdel > 1) then
#ifdef HEAPSORT
     call heapsort_int(numdel,deletelist(1:numdel))
#else
     call insertion_sort_int(numdel,deletelist(1:numdel))
#endif
  end if

#ifdef DEBUG_TRANSFER
  write (6,*) rank, ": numsendto = ", numsendto
  write (6,*) rank, ":     slots = ", slots
#endif

  ! Send those particles to the other domains
  ! --------------------------------------------------------------------------
  i = 0
  do d=0,lastrank
    if (d==rank) cycle
#ifdef DEBUG_TRANSFER
    write (6,*) rank, ": Sending ", numsendto(d), " particles to task ",&
    & d, " from slot ", slots(d)
#endif
    chunk_id = 1
    chunk_slot = slots(d)
    chunk_size_tot = numsendto(d)
    do
       if (numsendto(d) == 0) exit
       chunk_size = min(chunk_size_tot, max_message_size)
#ifdef DEBUG_TRANSFER
       write (6,*) rank, ": Chunk=", chunk_id, "; chunk_slot=", chunk_slot, &
                   & "; chunk_size=", chunk_size
#endif
       call MPI_ISEND(sendparticles(chunk_slot),chunk_size,MPI_PARTICLE,&
          &d,TRANSFER_TAG,MPI_COMM_WORLD,send(chunk_id,i),ierr)
       chunk_size_tot = chunk_size_tot - chunk_size
       chunk_slot = chunk_slot + chunk_size
       if (chunk_size_tot == 0) exit
       chunk_id = chunk_id + 1
       if (chunk_id > MAX_CHUNKS) then
          write (6,*) "Too many chunks! transfer_particles"
          stop 1
       end if
    end do
    i = i + 1
  end do
  ! --------------------------------------------------------------------------

  deallocate(exportlist)

#ifdef DEBUG_TRANSFER
  write (6,*) rank, ": Removing dead particles..."
#endif
  debug_timing("TRANSFER_PARTICLES REMOVE_DEAD")

  ! Now remove all dead particles
  i = sum(numsendto)
  if (i>0) then
    call remove_particles(i,deletelist,.FALSE.)
  end if
  deallocate(deletelist)

  allocate(receiveparticles(0:endranklist))

  debug_timing("TRANSFER_PARTICLES RECEIVE")

  ! Wait for receives from other domains to determine total incoming number
  ! --------------------------------------------------------------------------
  do i=0,lastrank-1
    !! Test for any message - wait until one arrives (blocking probe)
    ! Test for knowledge of number of incoming particles from any task
    WAIT_TIME_MACRO
    !call MPI_PROBE(MPI_ANY_SOURCE, TRANSFER_TAG, MPI_COMM_WORLD, stat, ierr)
    call MPI_WAITANY(lastrank, recv_num, p, stat, ierr)
    CALC_TIME_MACRO
    d = stat(MPI_SOURCE)
    !call MPI_GET_COUNT(stat,MPI_PARTICLE,pcount(i),ierr)
    pcount(i) = numrecvfrom(d)

#ifdef DEBUG_TRANSFER
    write (6,*) "Rank ",rank," receiving particles from task ", d
    !write (6,*) "Message length ", pcount(i)
    write (6,*) "Message length ", numrecvfrom(d)
#endif

    ! Allocate space, and post receives of incoming particles
    !allocate(receiveparticles(i)%data(1:max(pcount(i),1)))
    allocate(receiveparticles(i)%data(1:max(numrecvfrom(d),1)))

    ! Post the receive - until this is posted the receive does not actually
    ! occur and may or may not be buffered
    
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
       call MPI_IRECV(receiveparticles(i)%data(chunk_slot), chunk_size, &
                     &MPI_PARTICLE, d, TRANSFER_TAG, MPI_COMM_WORLD, &
                     recv(chunk_id,i), ierr)
       chunk_size_tot = chunk_size_tot - chunk_size
       chunk_slot = chunk_slot + chunk_size
       if (chunk_size_tot == 0) exit
       chunk_id = chunk_id + 1
    end do

  end do
  ! --------------------------------------------------------------------------

  ! If arrays are too small OR too large, expand them OR shrink them
  ! Use last value of pghost here (larger of h_find and scatter ghosts)
  !p_est = sum(pcount) + last_pghost + ptot
  p_est = sum(numrecvfrom) + last_pghost + ptot
  p_est_max = real(PMAXMULT,DP) * real(p_est,DP)
  if (p_est > pmax .OR. ceiling(real(PMAXMULT,DP)*p_est_max) < pmax) then
     call expand(ceiling(p_est_max))
  end if

  ! Wait for all sends and receives to complete
  debug_timing("WAIT TO RECEIVE PARTICLES")
#ifdef DEBUG_TRANSFER
  write (6,*) rank, ": Waiting for receives to complete..."
#endif
  WAIT_TIME_MACRO
  call MPI_WAITALL(lastrank*MAX_CHUNKS,recv,MPI_STATUSES_IGNORE,ierr)
  CALC_TIME_MACRO

  ! Unpack data into arrays
  call unpack_receiveparticles(pcount)
!   call unpack_receiveparticles(numrecvfrom)

  debug_timing("WAIT TO SEND PARTICLES")
#ifdef DEBUG_TRANSFER
    write (6,*) rank, ": Waiting for sends to complete..."
#endif

  WAIT_TIME_MACRO
  call MPI_WAITALL(lastrank,send_num,MPI_STATUSES_IGNORE,ierr)
  call MPI_WAITALL(lastrank*MAX_CHUNKS,send,MPI_STATUSES_IGNORE,ierr)
  CALC_TIME_MACRO
  deallocate(sendparticles)

  ! Deallocate receiveparticles
  do i=0,lastrank-1
     deallocate(receiveparticles(i)%data)
  end do
  deallocate(receiveparticles)

#ifdef DEBUG_TRANSFER
  write (1,*) "------------------------------"
  close(1)
#endif

  ! Need to re-build trees after particle transfer
  call tree_update(nsteps,nsteps)

  return
end subroutine transfer_particles
