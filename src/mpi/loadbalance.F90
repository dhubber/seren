! DECOMPOSITION.F90
! C. P. Batty - 11/6/2008
! Domain decomposition subroutine for MPI
! ============================================================================

#include "macros.h"

! ============================================================================
subroutine loadbalance()
   use MPI
   use mpi_communication_module
   use mpi_decomposition_module
   use particle_module
   use hydro_module
   use type_module
   use time_module
   use periodic_module
   use filename_module, only : run_id
#if defined(SINKS)
   use sink_module, only : stot, sink
#endif
   implicit none

   integer :: acctot_list(0:lastrank)   ! To avoid temporaries warning
   integer :: d                         ! domain counter
   integer :: p, i, j, k, x, y, z       ! Loop counters
   integer :: slowest_task              ! Rank of slowest task
   integer :: s(1:NDIM)                 ! Grid slot
   integer :: mingrid(1:NDIM)           ! Min positions of box in grid
   integer :: maxgrid(1:NDIM)           ! Max positions of box in grid
   integer :: dimen                     ! Dimension tree has split along
   integer :: mpicount                  ! MPI number of elements
   integer :: ierr                      ! MPI error value

   integer :: second_h_list(0:lastrank) ! To avoid temporaries warning
   integer :: ptot_list(0:lastrank)     ! To avoid temporaries warning
   integer :: predict_acctot_list(0:lastrank) ! To avoid temporaries warning
   real(kind=PR) :: overlap_list(0:lastrank)  ! To avoid temporaries warning
   real(kind=PR) :: costs_list(0:lastrank)    ! To avoid temporaries warning
   real(kind=PR) :: rowsum(0:gridsize)        ! Sum along a dimen of a part 
                                              ! of temppdensitygrid

   real(kind=PR) :: left                ! Number to left of current position
   real(kind=PR) :: ln                  ! Particle numbers (left branch)
   real(kind=PR) :: aux1, aux2, aux3    ! just some variables
   real(kind=PR) :: Vfrac               ! Volume filling fraction
   real(kind=PR) :: dx(1:NDIM)          ! Grid spacing
   real(kind=PR) :: ds                  ! Size of box for interpolation
   real(kind=PR) :: curposr(1:NDIM)     ! Temporary position variable
   real(kind=PR) :: minr(1:NDIM)        ! minimum particle positions
   real(kind=PR) :: maxr(1:NDIM)        ! maximum particle positions
   real(kind=PR) :: minrtemp(1:NDIM)    ! temp. minimum position
   real(kind=PR) :: maxrtemp(1:NDIM)    ! temp. maximum position
   real(kind=PR) :: COM(1:NDIM)         ! Average particle density position
   real(kind=PR) :: MOI(1:NDIM)         ! 'Moment of inertia' of ptcl density
   real(kind=PR) :: aux(1:NDIM)         ! Temporary variable
   real(kind=PR) :: totalpdensity       ! Total of work-weighted density
   real(kind=PR) :: local_vol           ! Volume of local non-overlapped domain
   real(kind=PR) :: d_vol               ! Volume of domain
   real(kind=PR) :: overlap             ! Fraction of domain which is overlapped
   real(kind=PR) :: p_costs(1:6)        ! ptot, acctot, predict_acctot, 
                                        ! second_h, costs, overlap
   real(kind=DP), allocatable :: receivecalctimes(:)   ! Calc times for tasks
   real(kind=DP), allocatable :: receivetottimes(:)    ! Total times for tasks
   real(kind=DP), allocatable :: receivewaittimes(:)   ! Wait times for tasks
   real(kind=DP), allocatable :: receiveexporttimes(:) ! Export times for tasks
   real(kind=PR), allocatable :: receive_costs(:,:)    ! Number of particles, acctot

#if NDIM==3
   real (kind=PR), allocatable :: temppdensity(:,:,:)
                   ! Local particle density grid
#elif NDIM==2
   real (kind=PR), allocatable :: temppdensity(:,:)
                   ! Local particle density grid
#else
   real (kind=PR), allocatable :: temppdensity(:)
                   ! Local particle density grid
#endif

   debug2("Load balancing [loadbalance.F90]")

   if (rank==0) then
#if NDIM==3
      allocate(temppdensity(0:gridsize,0:gridsize,0:gridsize))
#elif NDIM==2
      allocate(temppdensity(0:gridsize,0:gridsize))
#else
      allocate(temppdensity(0:gridsize))
#endif
   end if

! Reset load balance markers
   do_load_balance = .FALSE.
   nloadbalance = nsteps
   last_loadbalance_time = MPI_WTIME() - since_loadbalance_time

#ifdef PERIODIC
! Rationalize particle positions
   do p=1,ptot
      call unwrap_particle_position(sph(p)%r)
   end do
#endif

! Find the minimum and maximum particle positions in each dimension
   minr = BIG_NUMBER
   maxr = -BIG_NUMBER
   do p=1,ptot
      minr = min(minr,sph(p)%r(1:NDIM))
      maxr = max(maxr,sph(p)%r(1:NDIM))
   end do

#ifdef SINKS
   do i=1,stot
      minr = min(minr,sink(i)%r)
      maxr = max(maxr,sink(i)%r)
   end do
#endif

! Reduce minr/maxr amongst all tasks
   call MPI_ALLREDUCE(MPI_IN_PLACE,minr,NDIM,MPI_REAL_PR,&
      & MPI_MIN,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE,maxr,NDIM,MPI_REAL_PR,&
      & MPI_MAX,MPI_COMM_WORLD,ierr)

#ifdef PERIODIC
! Force 'root' of domain tree to take the size of the periodic box
   where (leftwall) minr = periodic_min
   where (rightwall) maxr = periodic_max
#endif
! Fudge very slightly where not a boundary
   dx = 0.0001d0 * (maxr - minr) / real(gridsize,PR)
   where (.not. leftwall) minr = minr - dx
   where (.not. rightwall) maxr = maxr + dx
! Recalculate dx on fudged values
   dx = (maxr - minr) / real(gridsize,PR)

! Calculate volume of total and non-overlapped domain
   local_vol = product(local_max-local_min)
   d_vol = product(domain_bbmax(1:NDIM,rank)-domain_bbmin(1:NDIM,rank))
   overlap = (d_vol - local_vol) / d_vol
   overlap = overlap * 100._PR

   allocate(receivecalctimes(0:lastrank))
   allocate(receivetottimes(0:lastrank))
   allocate(receivewaittimes(0:lastrank))
   allocate(receiveexporttimes(0:lastrank))
   allocate(receive_costs(1:6,0:lastrank))

! Now gather wait times
   call MPI_GATHER(calctime,1,MPI_DOUBLE_PRECISION,receivecalctimes,1,&
   &MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

   call MPI_GATHER(last_loadbalance_time,1,MPI_DOUBLE_PRECISION,&
        &receivetottimes,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

   call MPI_GATHER(waittime,1,MPI_DOUBLE_PRECISION,receivewaittimes,1,&
        &MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

   call MPI_GATHER(exporttime,1,MPI_DOUBLE_PRECISION,receiveexporttimes,1,&
        &MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

   p_costs = (/real(ptot,PR),real(sum_acctot,PR),real(predict_acctot,PR),&
        &real(sum_second_h,PR),sum_costs,overlap/)

   call MPI_GATHER(p_costs,6,MPI_REAL_PR,receive_costs,6,&
   &MPI_REAL_PR,0,MPI_COMM_WORLD,ierr)
   ! sum_acctot = 0 ! We now do this in loadbalance_time
   sum_second_h = 0
   sum_costs = 0._PR
   slowest_task = -1

   if (rank==0) then
      ptot_list = int(receive_costs(1,0:lastrank))
      acctot_list = int(receive_costs(2,0:lastrank))
      sum_acctot = sum(acctot_list)
      predict_acctot_list = int(receive_costs(3,0:lastrank))
      second_h_list = int(receive_costs(4,0:lastrank))
      costs_list = receive_costs(5,0:lastrank)
      overlap_list = receive_costs(6,0:lastrank)
      slowest_task = maxloc(receivecalctimes,1) - 1

#if defined(DEBUG_MPI) || defined(DEBUG_LOADBALANCE)
      write (6,*) "Load balancing efficiency: ", &
      &sum(receivecalctimes)/sum(receivetottimes)*100._PR, "%"
      write (6,*) "receivecalctimes = ", receivecalctimes
      write (6,*) "taskbias      =    ", taskbias
      write (6,*) "Load balancing: times and task biases"
      write (6,*) "Calculation times"
      write (6,'(9999(F14.2,1X:))') receivecalctimes
      write (6,*) "Total integration times: "
      write (6,'(9999(F14.2,1X:))') receivetottimes
      write (6,*) "MPI waiting times: "
      write (6,'(9999(F14.2,1X:))') receivewaittimes
      write (6,*) "Time calculating for exported particles: "
      write (6,'(9999(F14.2,1X:))') receiveexporttimes
      write (6,*) "Percentage of domain which is overlapped: "
      write (6,'(9999(F14.4,1X:))') overlap_list
      write (6,*) "Number of particles in each domain: "
      write (6,'(9999(I14,1X:))') ptot_list
      write (6,*) "Fraction of particles in each domain: "
      write (6,'(9999(F14.8,1X:))') real(ptot_list) / real(sum(ptot_list))
      write (6,*) "Active particle calculations in each domain: "
      write (6,'(9999(I14,1X:))') acctot_list
      write (6,*) "Predicted active particle calculations in each domain: "
      write (6,'(9999(I14,1X:))') predict_acctot_list
      write (6,*) "Second h gather calculations in each domain: "
      write (6,'(9999(I14,1X:))') second_h_list
      write (6,*) "Cost of calculations in each domain: "
      write (6,'(9999(F14.4,1X:))') costs_list
      write (6,*) "Fraction of cost in each domain: "
      write (6,'(9999(F14.8,1X:))') costs_list / sum(costs_list)
      write (6,*) "Average particle cost per calculation in each domain: "
      write (6,'(9999(F14.2,1X:))') costs_list / real(max(acctot_list,1),PR)
      write (6,*) "Slowest task is task: ", slowest_task
#endif
   end if

#if defined(DEBUG_MPI)
   if (rank==0) then
      write (6,*) "Load balance: Reallocating particles"
   end if
#endif

!    call MPI_SCATTER(receivecalctimes,1,MPI_DOUBLE_PRECISION,calctime,1,&
!    &MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

   deallocate(receivecalctimes)
   deallocate(receivetottimes)
   deallocate(receivewaittimes)
   deallocate(receiveexporttimes)
   deallocate(receive_costs)

#if NDIM==3
   allocate (MPIpdensitygrid(0:gridsize,0:gridsize,0:gridsize))
#elif NDIM==2
   allocate (MPIpdensitygrid(0:gridsize,0:gridsize))
#else
   allocate (MPIpdensitygrid(0:gridsize))
#endif
   ! Fill grid with particles
   call fill_grid(dx,minr,.FALSE.)
   
   ! -------------------------------------------------
   ! Reduce grid by sum
   mpicount = (gridsize+1)**NDIM
   debug_timing("REDUCE MPIpdensitygrid")
   if (rank==0) then
      call MPI_REDUCE(MPI_IN_PLACE,MPIpdensitygrid,mpicount,MPI_REAL_PR,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   else
      call MPI_REDUCE(MPIpdensitygrid,0,mpicount,MPI_REAL_PR,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   end if
   debug_timing("LOADBALANCE")

#if defined(DEBUG_LOADBALANCE)
   if (rank==0) then
      write (6,*) "sum(MPIpdensitygrid) = ", sum(MPIpdensitygrid)
      write (6,*) "size(MPItreemin(0)%data) = ", size(MPItreemin(0)%data)
      write (6,*) "MPItreemin(0)%data(0)%xyz = ", MPItreemin(0)%data(0)%xyz
      write (6,*) "MPItreemax(0)%data(0)%xyz = ", MPItreemax(0)%data(0)%xyz
   end if
#endif
   
#ifdef FIXED_BOUNDARIES
   if (.FALSE.) then
#else
   if (rank==0) then
#endif

      ! Create a tree from taskbias, summing upwards
      ! except that taskbias was a stupid idea, so always =1 now
      taskbias_tree(MPItreedepth)%data = 0._DP
      do d=0,lastrank
         ! Load taskbias_tree into correct places in tree
         taskbias_tree(MPItreedepth)%data(MPIgeometry(d)) = 1.0_PR
      end do
      do i=MPItreedepth-1, 0, -1
         ! Go backwards through the levels, summing up
         do d=0,(2**i)-1
            taskbias_tree(i)%data(d) = &
            & taskbias_tree(i+1)%data(d*2) + taskbias_tree(i+1)%data((d*2)+1)
         end do
      end do

!       open(unit=7,file="fillgrid.dat",status="REPLACE",form="FORMATTED",iostat=ierr)
!       do j=1,gridsize
!          do k=1,gridsize
!             write (7,*) j, k, MPIpdensitygrid(j,k,64)
!          end do
!       end do
!       close(7)

      do i=1,MPItreedepth
         d = (2**i)-1
         do j=0,d-1,2
   !           write (6,*) "Level ", i, " and slots ", j, " & ", (j+1)
            ! Load limits on min and max from higher level
            minrtemp = MPItreemin(i-1)%data(j/2)%xyz    ! Note integer division
            maxrtemp = MPItreemax(i-1)%data(j/2)%xyz    ! Note integer division
#ifdef DEBUG_LOADBALANCE
             write (6,*) "MPItreemin(i-1)%data(j/2)%xyz = ", MPItreemin(i-1)%data(j/2)%xyz
             write (6,*) "MPItreemax(i-1)%data(j/2)%xyz = ", MPItreemax(i-1)%data(j/2)%xyz
             write (6,*) "minrtemp = ", minrtemp
             write (6,*) "maxrtemp = ", maxrtemp
#endif

            ! Limit minr and maxr to min and max of grid
            minrtemp = max(minrtemp,minr)
            maxrtemp = min(maxrtemp,maxr)
#ifdef DEBUG_LOADBALANCE
!              write (6,*) "minrtemp = ", minrtemp
!              write (6,*) "maxrtemp = ", maxrtemp
#endif

            ! Use MPIpdensitygrid to calculate new min/max
            mingrid = int(((1.000001d0*minrtemp)-minr)/dx)
            maxgrid = int(((1.000001d0*maxrtemp)-minr)/dx)
#ifdef DEBUG_LOADBALANCE
             write (6,*) "mingrid = ", mingrid
             write (6,*) "maxgrid = ", maxgrid
#endif
! Start of section dealing with the grid
! ----------------------------------------------------------------------------
            temppdensity = 0._PR
            minrtemp = (minrtemp-minr)/dx
               ! Minimum in real space of overlap area, in terms of grid units
            maxrtemp = (maxrtemp-minr)/dx
               ! Maximum in real space of overlap area, in terms of grid units
#ifdef DEBUG_LOADBALANCE
             write (6,*) "minrtemp = ", minrtemp
             write (6,*) "maxrtemp = ", maxrtemp
#endif
            ! Iterate through all grid cells that are overlapped by the box
            do x=mingrid(1),maxgrid(1)
               s(1) = x
#if NDIM==2 || NDIM==3
               do y=mingrid(2),maxgrid(2)
                  s(2) = y
#endif
#if NDIM==3
                  do z=mingrid(3),maxgrid(3)
                     s(3) = z
#endif
                     if (all(s>mingrid).AND.all(s<maxgrid)) then
                        ! This box is entirely within the overlapped area
                        Vfrac = real(1,PR)
                     else
                     ! Calculate fraction of box occupied
#ifdef DEBUG_LOADBALANCE
!                      write (6,*) "Outside overlapped area"
!                      write (6,*) "Currently extracting from grid cell ", s
!                      write (6,*) "minrtemp, maxrtemp = ", minrtemp, maxrtemp
!                      write (6,*) "s=", s, " ; s+1 = ", s+1
!                      write (6,*) "max(real(s,PR),minrtemp) = ", max(real(s,PR),minrtemp)
!                      write (6,*) "min(real(s+1,PR),maxrtemp) = ", min(real(s+1,PR),maxrtemp)
#endif
                        Vfrac = product( min(real(s+1,PR),maxrtemp) - max(real(s,PR),minrtemp) )
                        Vfrac = max (Vfrac, 0._PR)
                     end if
#if NDIM==3
                     temppdensity(s(1),s(2),s(3)) = &
                     &MPIpdensitygrid(s(1),s(2),s(3)) * Vfrac
#elif NDIM==2
                     temppdensity(s(1),s(2)) = &
                     &MPIpdensitygrid(s(1),s(2)) * Vfrac
#else
                     temppdensity(s(1)) = &
                     &MPIpdensitygrid(s(1)) * Vfrac
#endif

#if NDIM==3
                  end do
#endif
#if NDIM==2 || NDIM==3
               end do
#endif
            end do
            totalpdensity = sum(real(temppdensity,PR))
#if defined(DEBUG_LOADBALANCE)
            write (6,*) "totalpdensity = ", totalpdensity
#endif

            ! Completely inherit min and max from upper level
            MPItreemax(i)%data(j)%xyz = MPItreemax(i-1)%data(j/2)%xyz
            MPItreemax(i)%data(j+1)%xyz = MPItreemax(i-1)%data(j/2)%xyz
            MPItreemin(i)%data(j)%xyz = MPItreemin(i-1)%data(j/2)%xyz
            MPItreemin(i)%data(j+1)%xyz = MPItreemin(i-1)%data(j/2)%xyz

            if (i==MPItreedepth.AND.MPIrevgeometry(j+1)>=numtasks) then
               ! This branch now only has one cell
               ! Make the right-hand branch of zero size (min = max of
               ! left hand branch)
               MPItreemin(i)%data(j+1)%xyz = MPItreemax(i)%data(j)%xyz
               cycle
            end if

            ! Calculate which dimension we SHOULD split along
            ! We cannot use exact COM and MOI, but we can do something
            ! equivalent with the (timestep-weighted) density

            ! Calculate COM and moments of inertia around COM
            COM = real(0.0,PR)
            MOI = real(0.0,PR)

            do x=mingrid(1),maxgrid(1)
               curposr(1) = (real((x-mingrid(1)),PR) + 0.5_PR) * dx(1)
#if NDIM==2 || NDIM==3
               do y=mingrid(2),maxgrid(2)
                  curposr(2) = (real((y-mingrid(2)),PR) + 0.5_PR) * dx(2)
#endif
#if NDIM==3
                  do z=mingrid(3),maxgrid(3)
                     curposr(3) = (real((z-mingrid(3)),PR) + 0.5_PR) * dx(3)
#endif
                     ! Add to COM and MOI
#if NDIM==3
                     aux(1:NDIM) = temppdensity(x,y,z) * curposr
#elif NDIM==2
                     aux(1:NDIM) = temppdensity(x,y) * curposr
#else
                     aux(1:NDIM) = temppdensity(x) * curposr
#endif
                     COM = COM + aux
                     aux(1:NDIM) = aux(1:NDIM) * curposr
                     MOI = MOI + aux
#if NDIM==3
                  end do
#endif
#if NDIM==2 || NDIM==3
               end do
#endif
            end do

            COM = COM / totalpdensity

            ! Translate MOI to COM: Icom = Iparallelaxis - m*d^2
            MOI = MOI - (totalpdensity*(COM**2))

#ifdef DEBUG_LOADBALANCE
            write (6,*) "COM = ", COM

            write (6,*) "MOI = ", MOI
#endif

            ! Axis to cut is the one with the largest MOI
            ! BUT don't change if difference is less than 10%
            k = maxloc(MOI,1)

!             if (k/=dimen) then
!             ! If axis to cut is the same anyway, don't worry
!                if (MOI(k) > 1.1*MOI(dimen)) dimen = k
!                ! If axis with largest MOI has MOI at least 10% larger than MOI of current splitting axis,
!                ! cut along axis with largest MOI.
!             end if
            dimen = k

#ifdef DEBUG_LOADBALANCE
             write (6,*) "We have decided on a dimension: ", dimen
#endif

#if NDIM==3
            ! Calculate rowsum along dimen
            x = mod(dimen,NDIM) + 1
            y = mod(dimen+1,NDIM) + 1
            rowsum = sum(sum(temppdensity,max(x,y)),min(x,y))
#ifdef DEBUG_LOADBALANCE
            write (6,*) "x & y: ", x, y
            write (6,*) "rowsum:", rowsum
#endif
#elif NDIM==2
            x = mod(dimen,NDIM) + 1
            rowsum = sum(temppdensity,x)
#else
            rowsum = temppdensity
#endif

            ! Calculate splitting point
            aux1 = taskbias_tree(i)%data(j)
            aux1 = aux1 / (aux1 + taskbias_tree(i)%data(j+1))
            ln = sum(rowsum) * aux1
            left = 0._PR
#ifdef DEBUG_LOADBALANCE
            write (6,*) "mingrid(dimen) = ", mingrid(dimen)
            write (6,*) "maxgrid(dimen) = ", maxgrid(dimen)
            write (6,*) "minrtemp(dimen) = ", minrtemp(dimen)
            write (6,*) "maxrtemp(dimen) = ", maxrtemp(dimen)
            write (6,*) "dimen:", dimen
#endif
            do k=mingrid(dimen),maxgrid(dimen)
               ! Go through rowsum until we get half-way, then interpolate splitting point
#ifdef DEBUG_LOADBALANCE
               write (6,*) "k ",k,"; left=", left
#endif
               if (left + rowsum(k) >= ln) then
#ifdef DEBUG_LOADBALANCE
                  write (6,*) "This will bring left over total, interpolate"
                  write (6,*) "minrtemp(dimen) = ", minrtemp(dimen)
                  write (6,*) "maxrtemp(dimen) = ", maxrtemp(dimen)
                  write (6,*) "real(k,PR)", real(k,PR)
                  write (6,*) "real(k+1,PR)", real(k+1,PR)
                  write (6,*) "ln = ", ln, "; left = ", left, "; rowsum(k) = ", rowsum(k)
#endif
                  ! Interpolate
                  ds = 1._PR
                  aux1 = 0._PR
                  if (real(k,PR) < minrtemp(dimen)) then
                     ds = ds - (minrtemp(dimen)-real(k,PR))
                     aux1 = minrtemp(dimen) - real(k,PR)
#ifdef DEBUG_LOADBALANCE
                     write (6,*) "real(k,PR) < minrtemp(dimen)"
                     write (6,*) "ds = ds - ", minrtemp(dimen)-real(k,PR)
                     write (6,*) "ds = ", ds
#endif
                  end if
                  if (real(k+1,PR) > maxrtemp(dimen)) then
                     ds = ds - (real(k+1,PR)-maxrtemp(dimen))
#ifdef DEBUG_LOADBALANCE
                     write (6,*) "real(k+1,PR) > maxrtemp(dimen)"
                     write (6,*) "ds = ds - ", real(k+1,PR)-maxrtemp(dimen)
                     write (6,*) "ds = ", ds
#endif
                  end if
                  ds = ds * dx(dimen)
#ifdef DEBUG_LOADBALANCE
                  write (6,*) "dx(", dimen, ") = ", dx(dimen)
                  write (6,*) "ds = ", ds
                  write (6,*) "Component due to minrtemp being in same cell as splitting point:", aux1
#endif
                  aux1 = aux1 * dx(dimen)
                  aux2 = ds * ((ln - left) / rowsum(k))
#ifdef DEBUG_LOADBALANCE
                  write (6,*) "aux1 = ", aux1
                  write (6,*) "Interpolation component: ", ((ln-left)/rowsum(k))
                  write (6,*) "Scaled: ", ((ln-left)/rowsum(k)) * ds / dx(dimen)
                  write (6,*) "aux2 = ", aux2
#endif
                  aux3 = real(k,PR) * dx(dimen)
#ifdef DEBUG_LOADBALANCE
                  write (6,*) "Grid component = ", real(k,PR)
                  write (6,*) "aux3 = ", aux3
#endif
                  aux3 = aux1 + aux2 + aux3 + minr(dimen)
                                        ! This is the splitting point
                  !aux3 = aux3 + minr(dimen)        ! This is the splitting point
#ifdef DEBUG_LOADBALANCE
                  write (6,*) "splitting point = ", aux3
                  if (aux3 < MPItreemin(i)%data(j)%xyz(dimen).OR.&
                     &aux3 > MPItreemax(i)%data(j)%xyz(dimen)) then
                     write (6,*) "Splitting point outside min/max!"
                     stop
                  end if
#endif
                  exit
               else
                  left = left + rowsum(k)
               end if
            end do

            MPItreemax(i)%data(j)%xyz(dimen) = aux3
            MPItreemin(i)%data(j+1)%xyz(dimen) = aux3

! End of section dealing with the grid
! ----------------------------------------------------------------------------

#ifdef DEBUG_LOADBALANCE
             write (6,*) "Level ", i, ", cells ", j, " & ", (j+1), "; splitting point ", aux3
             write (6,*) "**********************************************************************"
!              write (6,*) "MPItreemin = "
!              do x=0,MPItreedepth
!                write (6,*) MPItreemin(x)%data(:)%xyz
!              end do
!              write (6,*) "MPItreemax = "
!              do x=0,MPItreedepth
!                write (6,*) MPItreemax(x)%data(:)%xyz
!              end do
#endif

         end do
      end do

! Copy domain boundaries from tree
      do d=0,lastrank
         i = MPIgeometry(d)
         domain_bbmin(1:NDIM,d) = MPItreemin(MPItreedepth)%data(i)%xyz
         domain_bbmax(1:NDIM,d) = MPItreemax(MPItreedepth)%data(i)%xyz
#ifdef DEBUG_LOADBALANCE
         write (6,'(A,I0,A,3F10.7)') "domain_bbmin(",d,") = ", domain_bbmin(1:NDIM,d)
         write (6,'(A,I0,A,3F10.7)') "domain_bbmax(",d,") = ", domain_bbmax(1:NDIM,d)
#endif
      end do

   end if

   ! Broadcast domain boundaries to all tasks
   call MPI_BCAST(domain_bbmin(1,0),NDIM*numtasks,MPI_REAL_PR,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(domain_bbmax(1,0),NDIM*numtasks,MPI_REAL_PR,0,MPI_COMM_WORLD,ierr)

   local_min = domain_bbmin(1:NDIM,rank)
   local_max = domain_bbmax(1:NDIM,rank)

   deallocate (MPIpdensitygrid)
   if (rank==0) deallocate(temppdensity)

   calctime = 0._DP
   exporttime = 0._DP
   waittime = 0._DP
   since_loadbalance_time = MPI_WTIME()

#ifdef DEBUG_LOADBALANCE
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   if (rank==0) write (6,*) "LOADBALANCE COMPLETED"
#endif

! ----------------------------------------------------------------------------

  return
end subroutine loadbalance
