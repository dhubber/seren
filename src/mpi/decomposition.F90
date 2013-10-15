! DECOMPOSITION.F90
! C. P. Batty - 11/6/2008
! Domain decomposition subroutine for MPI
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE decomposition
   use filename_module, only : run_id, run_dir, fileform_ext, restart
   use mpi_communication_module
   use mpi_decomposition_module
   use particle_module
   use hydro_module
   use scaling_module
   use type_module
   use periodic_module
   use time_module
   use neighbour_module
   use sink_module
   implicit none

   character(len=20) :: data_id(1:50)   ! Char ids of arrays written
   character(len=20) :: data_id_root(1:50) ! data_id for root only
   character(len=100) :: dom_file       ! filename for domain data output
   character(len=7) :: file_ext         ! filename extension for data output
   character(len=20) :: format_id       ! File format (for verification)
   character(len=20) :: unit_data(1:50) ! Unit data
   integer :: curpos                    ! Current position in tree counter
   integer :: dimen                     ! Dimension tree has split along
   integer :: d                         ! domain counter
   integer :: i, j, k, x, y, z          ! Loop counters
   integer :: idata(1:50)               ! Integer data
   integer :: mingrid(1:NDIM)           ! Min positions of box in grid
   integer :: maxgrid(1:NDIM)           ! Max positions of box in grid
   integer :: ndata                     ! Number of arrays written
   integer :: ndata_root                ! Number of arrays written (root only)
   integer :: nunit                     ! Number of units
   integer :: p                         ! particle counter
   integer :: pdom(0:lastrank)          ! ptot in each domain
   integer :: pboundary_dom(0:lastrank) ! pboundary in each domain
   integer :: picm_dom(0:lastrank)      ! picm in each domain
   integer :: pgas_dom(0:lastrank)      ! pgas in each domain
   integer :: pcdm_dom(0:lastrank)      ! pcdm in each domain
   integer :: pdust_dom(0:lastrank)     ! pdust in each domain
   integer :: pion_dom(0:lastrank)      ! pion in each domain
   integer :: skip                      ! Skip counter for going through tree
   integer :: sl(1:NDIM)                ! Grid slot
   integer :: typedata(1:5,1:50)        ! type data header array
   integer :: typedata_root(1:5,1:50)   ! type data header array (root only)
   integer, allocatable :: domain(:)    ! domain indicator for particle
   integer, allocatable :: idummy(:)    ! ..
   integer(kind=ILP) :: ilpdata(1:50)   ! Long integer data
   real(kind=PR) :: aux(1:NDIM)         ! Temporary variable
   real(kind=PR) :: aux1, aux2, aux3    ! just some variables
   real(kind=PR) :: COM(1:NDIM)         ! Average particle density position
   real(kind=PR) :: curposr(1:NDIM)     ! Temporary position variable
   real(kind=PR) :: ds                  ! Size of box for interpolation
   real(kind=PR) :: dx(1:NDIM)          ! Grid spacing
   real(kind=DP) :: dpdata(1:50)        ! ..
   real(kind=PR) :: left                ! Number to left of current position
   real(kind=PR) :: ln                  ! Particle numbers (left branch)
   real(kind=PR) :: minr(1:NDIM)        ! minimum particle positions
   real(kind=PR) :: maxr(1:NDIM)        ! maximum particle positions
   real(kind=PR) :: minrtemp(1:NDIM)    ! temp. minimum position
   real(kind=PR) :: maxrtemp(1:NDIM)    ! temp. maximum position
   real(kind=PR) :: MOI(1:NDIM)         ! 'Moment of inertia' of ptcl density
   real(kind=PR) :: rowsum(0:gridsize)  ! Sum along a dimen of a part of pdensitygrid
   real(kind=PR) :: rdata(1:50)         ! Real data array
   real(kind=PR) :: totalpdensity       ! Total of work-weighted density
   real(kind=PR) :: Vfrac               ! Volume filling fraction
   real(kind=PR), allocatable :: rdummy1(:)     ! real dummy array
   real(kind=PR), allocatable :: rdummy3(:,:)   ! real vector dummy array
#if NDIM==3
   real(kind=PR), allocatable :: temppdensity(:,:,:)
                   ! Local particle density grid
#elif NDIM==2
   real(kind=PR), allocatable :: temppdensity(:,:)
                   ! Local particle density grid
#else
   real(kind=PR), allocatable :: temppdensity(:)
                   ! Local particle density grid

#endif
#if defined(SINKS)
   integer, parameter :: sink_data_length=12+NDIM+VDIM+2*DMDT_RANGE
   logical :: ldummy(1:2)                       ! Logical dummy array
   integer :: idummy2(1:2)                      ! Integer dummy array 2
   integer :: s                                 ! Sink counter
   real(kind=PR) :: raux(1:sink_data_length)    ! Aux. variable
#endif

   ! Write output that is suppressed in read routines for MPI code
   write(6,*) "SPH Particles  = ", ptot,"    Sink Particles = ", stot
   write(6,*) "Gas            = ", pgas
   write(6,*) "Boundary       = ", pboundary
   write(6,*) "Intercloud     = ", picm
   write(6,*) "Dark matter    = ", pcdm
   write(6,*) "Dust           = ", pdust
   write(6,*) "Ions           = ", pion
  
! Scale saved variables
   time = time / tscale
   lastsnap = lastsnap / tscale
   mgas_orig  = mgas_orig / mscale

#if NDIM==3
   allocate(temppdensity(0:gridsize,0:gridsize,0:gridsize))
#elif NDIM==2
   allocate(temppdensity(0:gridsize,0:gridsize))
#else
   allocate(temppdensity(0:gridsize))
#endif

! Set totalptot now
   totalptot = ptot

   allocate(domain(1:ptot))
   pdom(0:lastrank)          = 0
   pboundary_dom(0:lastrank) = 0
   picm_dom(0:lastrank)      = 0
   pgas_dom(0:lastrank)      = 0
   pcdm_dom(0:lastrank)      = 0
   pdust_dom(0:lastrank)     = 0
   pion_dom(0:lastrank)      = 0

   MPItreedepth = 0
   d = 1

! Determine the depth of the binary decomposition tree
   do
      if (d >= numtasks) exit
      MPItreedepth = MPItreedepth+1
      d = d * 2
   end do

! Allocate variables
   allocate(MPItreeoccupation(0:MPItreedepth))
   allocate(MPItreemin(0:MPItreedepth))
   allocate(MPItreemax(0:MPItreedepth))
   allocate(MPIgeometry(0:lastrank))
   allocate(MPIrevgeometry(0:(2**MPItreedepth)-1))
   allocate(taskbias_tree(0:MPItreedepth))

! Allocate and zero tree components
   do d=0,MPItreedepth
      allocate(MPItreeoccupation(d)%data(0:(2**d)-1))
      MPItreeoccupation(d)%data = 0
      allocate(MPItreemin(d)%data(0:(2**d)-1))
      allocate(MPItreemax(d)%data(0:(2**d)-1))
      do i=0,2**d-1
        MPItreemin(d)%data(i)%xyz = -BIG_NUMBER
        MPItreemax(d)%data(i)%xyz = BIG_NUMBER
      end do
      allocate(taskbias_tree(d)%data(0:(2**d)-1))
   end do

! Find the minimum and maximum particle positions in each dimension
   minr = BIG_NUMBER
   maxr = -BIG_NUMBER
   do p=1,ptot
      minr = min(minr,minimal_sph(p)%r(1:NDIM))
      maxr = max(maxr,minimal_sph(p)%r(1:NDIM))
   end do

   where (leftwall)
      minr = periodic_min
      MPItreemin(0)%data(0)%xyz = minr
   end where
   where (rightwall)
      maxr = periodic_max
      MPItreemax(0)%data(0)%xyz = maxr
   end where

   ! Fudge very slightly where boundary is not periodic
   dx = 0.0001d0 * (maxr - minr) / real(gridsize,PR)
   where (.not. leftwall) minr = minr - dx
   where (.not. rightwall) maxr = maxr + dx

   ! Recalculate dx on fudged values
   dx = (maxr - minr) / real(gridsize,PR)

   MPIgeometry = 0
   MPIrevgeometry = 0

   ! Check tree geometry is representable as an integer - if not, something is very wrong...
   if ( bit_size(MPIgeometry(1)) <= MPItreedepth ) then
      write (6,*) "Tree too deep to represent as an integer!"
      stop
   end if

   ! Define mapping of task numbers to bottom level tree cell numbers
   do d=0,lastrank
      do i=0,MPItreedepth-1
         call mvbits(d,i,1,MPIgeometry(d),(MPItreedepth-1)-i)
      end do
   end do

   ! Create MPIrevgeometry
   do d=0,(2**MPItreedepth)-1
      do i=0,MPItreedepth-1
         call mvbits(d,i,1,MPIrevgeometry(d),(MPItreedepth-1)-i)
      end do
   end do

! For each thread, determine its position in the geometry and calculate the tree occupation numbers
   MPItreeoccupation(0)%data(0) = numtasks
   do d=0,lastrank
      skip = 0
      ! Update tree occupation numbers
      do i=(MPItreedepth-1),0,-1 ! Reverse order do
         curpos = MPItreedepth - i
         if (btest(MPIgeometry(d),i)) then
            skip = skip+1
         end if
         MPItreeoccupation(curpos)%data(skip) = MPItreeoccupation(curpos)%data(skip) + 1
         skip = skip * 2
      end do
   end do

#if NDIM==3
   allocate (MPIpdensitygrid(0:gridsize,0:gridsize,0:gridsize))
#elif NDIM==2
   allocate (MPIpdensitygrid(0:gridsize,0:gridsize))
#else
   allocate (MPIpdensitygrid(0:gridsize))
#endif

   ! Calculate (for root only, replicate for other tasks) the grid
   phydrostart = pboundary + 1
   nlevel_sinks = 1_ILP
   call fill_grid(dx,minr,.TRUE.)

   do i=1,MPItreedepth
      d = (2**i)-1
      do j=0,d-1,2
!           write (6,*) "Level ", i, " and slots ", j, " & ", (j+1)
         ! Load limits on min and max from higher level
         minrtemp = MPItreemin(i-1)%data(j/2)%xyz    ! Note integer division
         maxrtemp = MPItreemax(i-1)%data(j/2)%xyz    ! Note integer division
!           write (6,*) "MPItreemin(i-1)%data(j/2)%xyz = ", MPItreemin(i-1)%data(j/2)%xyz
!           write (6,*) "MPItreemax(i-1)%data(j/2)%xyz = ", MPItreemax(i-1)%data(j/2)%xyz
!           write (6,*) "minrtemp = ", minrtemp
!           write (6,*) "maxrtemp = ", maxrtemp

         ! Limit minr and maxr to min and max of grid
         minrtemp = max(minrtemp,minr)
         maxrtemp = min(maxrtemp,maxr)
!           write (6,*) "minrtemp = ", minrtemp
!           write (6,*) "maxrtemp = ", maxrtemp

         ! Use MPIpdensitygrid to calculate new min/max
         mingrid = int(((1.000001d0*minrtemp)-minr)/dx)
         maxgrid = int(((1.000001d0*maxrtemp)-minr)/dx)
!           write (6,*) "mingrid = ", mingrid
!           write (6,*) "maxgrid = ", maxgrid
         temppdensity = real(0,PR)
         minrtemp = (minrtemp-minr)/dx     ! Minimum in real space of overlap area, in terms of grid units
         maxrtemp = (maxrtemp-minr)/dx     ! Maximum in real space of overlap area, in terms of grid units
!           write (6,*) "minrtemp = ", minrtemp
!           write (6,*) "maxrtemp = ", maxrtemp

         ! Iterate through all grid cells that are overlapped by the box
         do x=mingrid(1),maxgrid(1)
            sl(1) = x
#if NDIM==2 || NDIM==3
            do y=mingrid(2),maxgrid(2)
               sl(2) = y
#endif
#if NDIM==3
               do z=mingrid(3),maxgrid(3)
                  sl(3) = z
#endif
                  if (all(sl>mingrid).AND.all(sl<maxgrid)) then
                     ! This box is entirely within the overlapped area
                     Vfrac = real(1,PR)
                  else
                  ! Calculate fraction of box occupied
!                   write (6,*) "Outside overlapped area"
!                   write (6,*) "Currently extracting from grid cell ", sl
!                   write (6,*) "minrtemp, maxrtemp = ", minrtemp, maxrtemp
!                   write (6,*) "sl=", sl, " ; sl+1 = ", sl+1
!                   write (6,*) "max(real(sl,PR),minrtemp) = ", max(real(sl,PR),minrtemp)
!                   write (6,*) "min(real(sl+1,PR),maxrtemp) = ", min(real(sl+1,PR),maxrtemp)
                     Vfrac = product( min(real(sl+1,PR),maxrtemp) - max(real(sl,PR),minrtemp) )
                     Vfrac = max (Vfrac, 0._PR)
                  end if
#if NDIM==3
                  temppdensity(sl(1),sl(2),sl(3)) = &
                     & MPIpdensitygrid(sl(1),sl(2),sl(3)) * Vfrac
#elif NDIM==2
                  temppdensity(sl(1),sl(2)) = &
                     & MPIpdensitygrid(sl(1),sl(2)) * Vfrac
#else
                  temppdensity(sl(1)) = MPIpdensitygrid(sl(1)) * Vfrac
#endif

#if NDIM==3
               end do
#endif
#if NDIM==2 || NDIM==3
            end do
#endif
         end do

         ! Completely inherit min and max from upper level
         MPItreemax(i)%data(j)%xyz = MPItreemax(i-1)%data(j/2)%xyz
         MPItreemax(i)%data(j+1)%xyz = MPItreemax(i-1)%data(j/2)%xyz
         MPItreemin(i)%data(j)%xyz = MPItreemin(i-1)%data(j/2)%xyz
         MPItreemin(i)%data(j+1)%xyz = MPItreemin(i-1)%data(j/2)%xyz

         if (i==MPItreedepth.AND.MPIrevgeometry(j+1)>=numtasks) then
            ! This branch now only has one cell
            !MPItreemin(dimen,i,j+1) = MPItreemax(dimen,i,j)
            MPItreemin(i)%data(j+1)%xyz = MPItreemax(i)%data(j)%xyz
            cycle
         end if

         ! Calculate which dimension we should split along
         ! We cannot use exact COM and MOI, but we can do something equivalent
         ! with the work-weighted density

         ! Calculate COM and moments of inertia around COM
         totalpdensity = sum(temppdensity)
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

!          write (6,*) "COM = ", COM
! 
!          write (6,*) "MOI = ", MOI

         ! Axis to cut is the one with the largest MOI
         ! BUT don't change if difference is less than 10%
         k = maxloc(MOI,1)

!          if (k/=dimen) then
!          ! If axis to cut is the same anyway, don't worry
!             if (MOI(k) > 1.1*MOI(dimen)) dimen = k
!             ! If axis with largest MOI has MOI at least 10% larger than MOI of current splitting axis,
!             ! cut along axis with largest MOI.
!          end if
         dimen = k

!        write (6,*) "We have decided on a dimension: ", dimen

#if NDIM==3
         ! Calculate rowsum along dimen
         x = mod(dimen,NDIM) + 1
         y = mod(dimen+1,NDIM) + 1
         rowsum = sum(sum(temppdensity,max(x,y)),min(x,y))
!           write (6,*) "x & y: ", x, y
!           write (6,*) "rowsum:", rowsum
#elif NDIM==2
         x = mod(dimen,NDIM) + 1
         rowsum = sum(temppdensity,x)
#else
         rowsum = temppdensity
#endif

         ! Calculate splitting point
         aux1 = real(MPItreeoccupation(i)%data(j),PR)
         aux1 = aux1 / (aux1 + real(MPItreeoccupation(i)%data(j+1),PR))
         ln = sum(rowsum) * aux1
         left = 0._PR
         do k=0,gridsize-1
            ! Go through rowsum until we get half-way, then interpolate splitting point
!             write (6,*) "k ",k,"; left=", left
            if (left + rowsum(k) > ln) then
!               write (6,*) "This will bring left over total, interpolate"
               ! Interpolate
               ds = 1._PR
               aux1 = 0._PR
               if (real(k,PR) < minrtemp(dimen)) then
                  ds = ds - (minrtemp(dimen)-real(k,PR))
                  aux1 = minrtemp(dimen) - real(k,PR)
               end if
               if (real(k+1,PR) > maxrtemp(dimen)) then
                  ds = ds - (real(k+1,PR)-maxrtemp(dimen))
               end if
               ds = ds * dx(dimen)
!                write (6,*) "ln = ", ln, "; left = ", left, "; rowsum(k) = ", rowsum(k)
!                write (6,*) "dx(", dimen, ") = ", dx(dimen)
               aux1 = aux1 * dx(dimen)
               aux2 = ds * ((ln - left) / rowsum(k))
               aux3 = real(k,PR) * dx(dimen)
               aux3 = aux1 + aux2 + aux3 + minr(dimen)
               exit
            else
               left = left + rowsum(k)
            end if
         end do

         MPItreemax(i)%data(j)%xyz(dimen) = aux3
         MPItreemin(i)%data(j+1)%xyz(dimen) = aux3

!           write (6,*) "Level ", i, ", cells ", j, " & ", (j+1), "; splitting point ", aux2
!           write (6,*) "MPItreemin = ", MPItreemin(i)%data(j)%xyz
!           write (6,*) "MPItreemax = ", MPItreemax(i)%data(j)%xyz

      end do
   end do

   deallocate (MPIpdensitygrid)

! Allocate particles to their task
   if (pboundary > 0) then
      do p=pboundarystart,pboundaryend
         call find_task(minimal_sph(p)%r,d)
         domain(p) = d
         pdom(d) = pdom(d) + 1
         pboundary_dom(d) = pboundary_dom(d) + 1
      end do
   end if
   if (picm > 0) then
      do p=picmstart,picmend
         call find_task(minimal_sph(p)%r,d)
         domain(p) = d
         pdom(d) = pdom(d) + 1
         picm_dom(d) = picm_dom(d) + 1
      end do
   end if
   if (pgas > 0) then
      do p=pgasstart,pgasend
         call find_task(minimal_sph(p)%r,d)
         domain(p) = d
         pdom(d) = pdom(d) + 1
         pgas_dom(d) = pgas_dom(d) + 1
      end do
   end if
   if (pcdm > 0) then
      do p=pcdmstart,pcdmend
         call find_task(minimal_sph(p)%r,d)
         domain(p) = d
         pdom(d) = pdom(d) + 1
         pcdm_dom(d) = pcdm_dom(d) + 1
      end do
   end if
   if (pdust > 0) then
      do p=pduststart,pdustend
         call find_task(minimal_sph(p)%r,d)
         domain(p) = d
         pdom(d) = pdom(d) + 1
         pdust_dom(d) = pdust_dom(d) + 1
      end do
   end if
   if (pion > 0) then
      do p=pionstart,pionend
         call find_task(minimal_sph(p)%r,d)
         domain(p) = d
         pdom(d) = pdom(d) + 1
         pion_dom(d) = pion_dom(d) + 1
      end do
   end if

! Copy domain boundaries from tree
   do d=0,lastrank
      i = MPIgeometry(d)
      domain_bbmin(1:NDIM,d) = MPItreemin(MPItreedepth)%data(i)%xyz
      domain_bbmax(1:NDIM,d) = MPItreemax(MPItreedepth)%data(i)%xyz
   end do

#ifdef DEBUG_GHOST
   write (6,*) "Initial domain boundaries"
   do d=0,lastrank
      write (6,'(A,I0,A,3F10.7)') "domain_bbmin(1:NDIM,",d,") = ", domain_bbmin(1:NDIM,d)
      write (6,'(A,I0,A,3F10.7)') "domain_bbmax(1:NDIM,",d,") = ", domain_bbmax(1:NDIM,d)
   end do
#endif

! ----------------------------------------------------------------------------

! Diagnostics
  do p=0,lastrank
     write(6,*) "Domain", p, "contains", pdom(p), "particles"
  end do
  if (sum(pdom) /= ptot) stop "Particles in domains do not add up to ptot!"

! ----------------------------------------------------------------------------

! Write to files in seren_unform format

! Allocate main dummy arrays and zero them
  allocate(idummy(1:ptot))
  allocate(rdummy1(1:ptot))
  idata(1:50)     = 0
  ilpdata(1:50)   = 0_ILP
  dpdata(1:50)    = 0.0_DP
  rdata(1:50)     = 0.0_PR
  unit_data(1:50) = ''
  data_id(1:50)   = ''
  ndata           = 0
  
! Set unit character ids
  unit_data(1)  = runit
  unit_data(2)  = munit
  unit_data(3)  = tunit
  unit_data(4)  = vunit
  unit_data(5)  = aunit
  unit_data(6)  = rhounit
  unit_data(7)  = sigmaunit
  unit_data(8)  = Punit
  unit_data(9)  = funit
  unit_data(10) = Eunit
  unit_data(11) = momunit
  unit_data(12) = angmomunit
  unit_data(13) = angvelunit
  unit_data(14) = dmdtunit
  unit_data(15) = Lunit
  unit_data(16) = kappaunit
  unit_data(17) = Bunit
  unit_data(18) = Qunit
  unit_data(19) = Junit
  unit_data(20) = uunit
  unit_data(21) = tempunit
  nunit = 21

! Set array ids and array information data if there are any SPH particles
! ----------------------------------------------------------------------------
  if (ptot > 0) then
     ndata = ndata + 1;    data_id(ndata) = 'porig'
     typedata(1:5,ndata) = (/1,1,ptot,2,0/)
     
     ndata = ndata + 1;    data_id(ndata) = 'r'
     typedata(1:5,ndata) = (/NDIM,1,ptot,4,1/)
     
     ndata = ndata + 1;    data_id(ndata) = 'm'
     typedata(1:5,ndata) = (/1,1,ptot,4,2/)
     
     ndata = ndata + 1;    data_id(ndata) = 'h'
     typedata(1:5,ndata) = (/1,1,ptot,4,1/)
     
     ndata = ndata + 1;    data_id(ndata) = 'v' 
     typedata(1:5,ndata) = (/VDIM,1,ptot,4,4/)
     
#if defined(MHD)
     ndata = ndata + 1;    data_id(ndata) = 'B' 
     typedata(1:5,ndata) = (/BDIM,1,ptot,4,4/)
#endif
     
     ndata = ndata + 1;    data_id(ndata) = 'rho'
     typedata(1:5,ndata) = (/1,1,ptot,4,6/)
#if defined(HYDRO)
     ndata = ndata + 1;    data_id(ndata) = 'temp'
     typedata(1:5,ndata) = (/1,1,ptot,4,21/)
     
     ndata = ndata + 1;    data_id(ndata) = 'u'
     typedata(1:5,ndata) = (/1,1,ptot,4,20/)
#endif
  end if
  
  ndata_root = ndata
  typedata_root = typedata
  data_id_root = data_id

#if defined(SINKS)
  if (stot > 0) then
     ndata_root = ndata_root + 1;    data_id_root(ndata_root) = 'sink_v1'
        typedata_root(1:5,ndata_root) = (/0,1,stot,7,0/)
  end if
#endif

! Set important header information
  idata(2)    = stot
  idata(20)   = nunit
  idata(30)   = DMDT_RANGE
  idata(31)   = pgas_orig
  idata(32)   = pp_gather
  idata(41)   = numtasks
  ilpdata(1)  = snapshot
  ilpdata(2)  = nsteps
  ilpdata(3)  = ntempnext
  ilpdata(4)  = ndiagnext
  ilpdata(5)  = nsnapnext
  ilpdata(6)  = nsinknext
  rdata(1)    = real(h_fac,PR)
  rdata(2)    = real(gamma,PR)
  rdata(3)    = real(mu_bar,PR)
  dpdata(1)   = time*tscale
  dpdata(2)   = lastsnap*tscale
  dpdata(3)   = mgas_orig*mscale

  format_id = 'SERENBINARYDUMPV2'

! Write domains to separate files --------------------------------------------
  do d=0,lastrank

     write(file_ext, "(I0)") d
     dom_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//&
     &trim(adjustl(fileform_ext))//".tmp.MPI."//trim(adjustl(file_ext))

     open(1, file=dom_file, status="unknown", form="unformatted")

     idata(1)    = pdom(d)
     idata(3)    = pboundary_dom(d)
     idata(4)    = picm_dom(d)
     idata(5)    = pgas_dom(d)
     idata(6)    = 0 !pcdm
     idata(7)    = 0 !pdust
     idata(8)    = 0 !pion
     idata(40)   = d
     if (d == 0) then
        idata(21) = ndata_root
        typedata_root(3,1:ndata) = pdom(d) ! Set correct number of particles
        ! note that we only go to ndata as if we have sinks we don't want to
        ! change plast in the typedata header for them
     else
        idata(21) = ndata
        typedata(3,1:ndata) = pdom(d) ! Set correct number of particles
     end if

! Write information identifying format and precision of file
     write(1) format_id
     write(1) PR
     write(1) NDIM
     write(1) VDIM
     write(1) BDIM

! Write header information to file
     write(1) idata
     write(1) ilpdata
     write(1) rdata
     write(1) dpdata
     if (d==0) then
        if (nunit > 0) write(1) unit_data(1:nunit)
        if (ndata_root > 0) write(1) data_id_root(1:ndata_root)
        if (ndata_root > 0) write(1) typedata_root(1:5,1:ndata_root)
     else
        if (nunit > 0) write(1) unit_data(1:nunit)
        if (ndata > 0) write(1) data_id(1:ndata)
        if (ndata > 0) write(1) typedata(1:5,1:ndata)
     end if

! Original ids
! ----------------------------------------------------------------------------
     i = 0
     do p=1,ptot
        if (domain(p) == d) then
           i = i + 1
           idummy(i) = minimal_sph(p)%porig
        end if
     end do
     write(1) idummy(1:i)
     
     
! Positions
! ----------------------------------------------------------------------------
     allocate(rdummy3(1:NDIM,1:ptot))
     i = 0
     do p=1,ptot
        if (domain(p) == d) then
           i = i + 1
           rdummy3(1:NDIM,i) = minimal_sph(p)%r*real(rscale,PR)
        end if
     end do
     write(1) rdummy3(1:NDIM,1:i)
     deallocate(rdummy3)
     
! Mass
! ----------------------------------------------------------------------------
     
     i = 0
     do p=1,ptot
        if (domain(p) == d) then
           i = i + 1
           rdummy1(i) = minimal_sph(p)%m*real(mscale,PR)
        end if
     end do
     write(1) rdummy1(1:i)
     
! Smoothing length
! ----------------------------------------------------------------------------
     
     i = 0
     do p=1,ptot
        if (domain(p) == d) then
           i = i + 1
           rdummy1(i) = minimal_sph(p)%h*real(rscale,PR)
        end if
     end do
     write(1) rdummy1(1:i)

! Velocities
! ----------------------------------------------------------------------------
     allocate(rdummy3(1:VDIM,1:ptot))
     i = 0
     do p=1,ptot
        if (domain(p) == d) then
           i = i + 1
           rdummy3(1:VDIM,i) = minimal_sph(p)%v(1:VDIM)*vscale
        end if
     end do
     write(1) rdummy3(1:VDIM,1:i)
     deallocate(rdummy3)

#if defined(MHD)
! Magnetic fields
! ----------------------------------------------------------------------------
     allocate(rdummy3(1:BDIM,1:ptot))
     i = 0
     do p=1,ptot
        if (domain(p) == d) then
           i = i + 1
           rdummy3(1:BDIM,i) = minimal_sph(p)%B(1:BDIM)*Bscale
        end if
     end do
     write(1) rdummy3(1:BDIM,1:i)
     deallocate(rdummy3)
#endif

! Density
! ----------------------------------------------------------------------------
     i = 0
     do p=1,ptot
        if (domain(p) == d) then
           i = i + 1
           rdummy1(i) = minimal_sph(p)%rho*rhoscale
        end if
     end do
     write(1) rdummy1(1:i)

! Temperature
! ----------------------------------------------------------------------------
#if defined(HYDRO)
     i = 0
     do p=1,ptot
        if (domain(p) == d) then
           i = i + 1
           rdummy1(i) = minimal_sph(p)%temp
        end if
     end do
     write(1) rdummy1(1:i)
#endif

! Internal energy
! ----------------------------------------------------------------------------
#if defined(HYDRO)
     i = 0
     do p=1,ptot
        if (domain(p) == d) then
           i = i + 1
           rdummy1(i) = minimal_sph(p)%u*Escale
        end if
     end do
     write(1) rdummy1(1:i)
#endif

! Sinks
! ----------------------------------------------------------------------------
#if defined(SINKS)
     if (stot > 0 .AND. d==0) then
        write(1) 2,2,0,sink_data_length,0,0
        do s=1,stot
           ldummy(1)                     = sink(s)%accrete
           ldummy(2)                     = sink(s)%static
           idummy2(1)                    = sink(s)%id
           idummy2(2)                    = sink(s)%ncreate
           raux(1:sink_data_length)      = 0.0_PR
           raux(1)                       = real(sink(s)%tcreate*tscale,PR)
           raux(2:NDIM+1)                = sink(s)%r(1:NDIM)*real(rscale,PR)
           raux(NDIM+2:NDIM+VDIM+1)      = sink(s)%v(1:NDIM)*real(vscale,PR)
           raux(NDIM+VDIM+2)             = sink(s)%m*real(mscale,PR)
           raux(NDIM+VDIM+3)             = sink(s)%h*real(rscale,PR)
           raux(NDIM+VDIM+4)             = sink(s)%radius*real(rscale,PR)
           raux(NDIM+VDIM+5:NDIM+VDIM+7) = real(sink(s)%angmom(1:3)*angmomscale,PR)
           raux(NDIM+VDIM+8)             = real(sink(s)%dmdt*dmdtscale,PR)
           raux(NDIM+VDIM+9)             = real(sink(s)%star_radius*rscale,PR)
           raux(NDIM+VDIM+10)            = real(sink(s)%luminosity*Lscale,PR)
           raux(NDIM+VDIM+11)            = real(sink(s)%temperature,PR)
           if (DMDT_RANGE > 0) then
              raux(NDIM+VDIM+12:NDIM+VDIM+11+DMDT_RANGE) = &
                 & real(sink(s)%macc(1:DMDT_RANGE)*mscale,PR)
              raux(NDIM+VDIM+12+DMDT_RANGE:NDIM+VDIM+11+2*DMDT_RANGE) = &
                 & real(sink(s)%tacc(1:DMDT_RANGE)*tscale,PR)
           end if
           raux(NDIM+VDIM+12+2*DMDT_RANGE) = real(sink(s)%mmax*mscale,PR)
           write(1) ldummy
           write(1) idummy2
           write(1) raux
        end do
     end if
#endif

     close(1)

  end do
  deallocate(domain)
  deallocate(temppdensity)
! ----------------------------------------------------------------------------

  return
END SUBROUTINE decomposition
