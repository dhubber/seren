! DEBUG_TESTS.F90
! A. McLeod - 30/07/2008
! Subroutines to test stuff for debugging purposes
! ============================================================================

#include "macros.h"

! ============================================================================

module debug_tests

   contains

! ----------------------------------------------------------------------------
! Subroutines to check that particles have not been duplicated or overlap

   subroutine check_particle_uniqueness(local)
      use particle_module
#if defined(USE_MPI)
      use mpi_communication_module
      use mpi
#endif
      implicit none

      logical, intent(in)  :: local             ! Is this purely a local check?
                                                ! Irrelevant for non-MPI
      integer, allocatable :: porigcount(:)     ! Histogram of porig numbers
                                                ! (should be zero/one per porig!)
      integer              :: maxporig          ! Maximum porig
      integer              :: ierr              ! MPI error value
      integer              :: p                 ! Particle counter
      logical              :: OK                ! Failure = FALSE

      debug_timing("CHECK PARTICLE UNIQUENESS")

      OK = .TRUE.

      maxporig = maxval(sph(1:ptot)%porig)

#if defined(USE_MPI)
      if (.NOT. local) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, maxporig, 1, MPI_INTEGER, MPI_MAX, &
            &MPI_COMM_WORLD, ierr)
      end if
#endif
     
      allocate (porigcount(1:maxporig))
      porigcount = 0

      do p=1,ptot
         porigcount(sph(p)%porig) = porigcount(sph(p)%porig) + 1
      end do

#if defined(USE_MPI)
      if (.NOT. local) then
         if (rank==0) then
            call MPI_REDUCE(MPI_IN_PLACE, porigcount(1), maxporig, MPI_INTEGER, &
               &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         else
            call MPI_REDUCE(porigcount(1), 0, maxporig, MPI_INTEGER, &
               &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         end if
     
         if (rank/=0) return
      end if
#endif

      do p=1,maxporig
         if (porigcount(p)>1) then
            write (6,*) "Particle porig==",p," exists ",porigcount(p)," times!"
            OK = .FALSE.
         end if
      end do

      if (.NOT. OK) stop 1

      write (6,*) "Particle uniqueness verified"

      return
   end subroutine check_particle_uniqueness

   subroutine doublecheck_particle_uniqueness(local)
      use particle_module
      use scaling_module
#if defined(USE_MPI)
      use mpi_communication_module
      use mpi
#endif
      implicit none
     
      logical, intent(in)  :: local             ! Is this purely a local check?
                                                ! Irrelevant for non-MPI
      integer, allocatable :: porigcount(:)     ! Histogram of porig numbers
                                                ! (should be zero/one per porig!)
      real(kind=PR), allocatable :: r(:,:)      ! Particle locations
      real (kind=PR)       :: dr(1:NDIM)        ! Particle-particle vector
      real (kind=PR)       :: drsqd             ! Particle-particle dist. sqd.
      integer              :: maxporig          ! Maximum porig
      integer              :: ierr              ! MPI error value
      integer              :: p                 ! Particle counter
      integer              :: i                 ! Loop counter
      integer              :: j                 ! Loop counter
      logical              :: OK                ! Failure = FALSE

      debug_timing("DOUBLECHECK PARTICLE UNIQUENESS")

      OK = .TRUE.

      maxporig = maxval(sph(1:ptot)%porig)
     
#if defined(USE_MPI)
      if (.NOT. local) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, maxporig, 1, MPI_INTEGER, MPI_MAX, &
            &MPI_COMM_WORLD, ierr)
      end if
#endif

      allocate (porigcount(1:maxporig))
      porigcount = 0
      allocate (r(1:NDIM,1:maxporig))
      r = 0._PR

      do p=1,ptot
         porigcount(sph(p)%porig) = porigcount(sph(p)%porig) + 1
      end do
      do p=1,ptot
         r(1:NDIM,sph(p)%porig) = sph(p)%r
      end do

#if defined(USE_MPI)
      if (.NOT. local) then
         if (rank==0) then
            call MPI_REDUCE(MPI_IN_PLACE, porigcount(1), maxporig, MPI_INTEGER, &
               &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         else
            call MPI_REDUCE(porigcount(1), 0, maxporig, MPI_INTEGER, &
               &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         end if
         if (rank==0) then
            call MPI_REDUCE(MPI_IN_PLACE, r(1,1), NDIM*maxporig, MPI_REAL_PR,&
            &MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         else
            call MPI_REDUCE(r(1,1), 0, NDIM*maxporig, MPI_REAL_PR, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         end if
      
         if (rank/=0) return
        
      end if
#endif

      do p=1,maxporig
         if (porigcount(p)>1) then
            write (6,*) "Particle porig==",p," exists ",porigcount(p)," times!"
            OK = .FALSE.
         end if
      end do

      write (6,*) "About to check particle-particle spacing."
      write (6,*) "Warning - this will take an extremely long time for larger"
      write (6,*) "numbers of particles as it is an n^2 operation."

      do i=1,maxporig-1
         if (porigcount(i) == 0) cycle
         write (6,*) i
         do j=i+1,maxporig
            if (porigcount(j) == 0) cycle
            dr = r(1:NDIM,i) - r(1:NDIM,j)
            drsqd = sum(dr**2)
#if defined(BHTREE)
            if (drsqd <= min(maxval(spacing(drsqd)),2._PR**(-LMAX))) then
#else
            if (drsqd <= maxval(spacing(r))) then
#endif
               write (6,*) "Particle porig==",i," too close to porig==",j,&
                  &", separation ",sqrt(drsqd),"!"
               write (6,*) "r(:,",i,") = ", r(:,i)
               write (6,*) "r(:,",j,") = ", r(:,j)
               write (6,*) "|dr| (non-code units) = ", sqrt(drsqd)*real(rscale,PR)
               OK = .FALSE.
            end if
         end do
      end do

      if (.NOT. OK) stop 1

      write (6,*) "Particle position uniqueness verified"

      return
   end subroutine doublecheck_particle_uniqueness

! ----------------------------------------------------------------------------
! Subroutines to check the BH tree

#ifdef BH_TREE
   SUBROUTINE check_tree_integrity
      use time_module
      use type_module
      use particle_module
      use tree_module
      implicit none

      integer  :: c, cc, p, i, leaf
      integer  :: countrefs(1:ptot)
      integer, allocatable :: count_tree_refs(:)
      logical  :: OK

      debug_timing("CHECK TREE INTEGRITY")

   ! Check tree integrity

#if defined(HYDRO)
      allocate(count_tree_refs(0:cmax_hydro))
      countrefs = 0
      count_tree_refs = 0
      c = 0
      do
         count_tree_refs(c) = count_tree_refs(c) + 1
         if (BHhydro(c)%leaf > 0) then
            do i=1,BHhydro(c)%leaf
               p = BHhydro(c)%plist(i)
               if (p > ptot) then
                  write (6,*) "Hydro cell ", c, " references a particle ", p
                  write (6,*) "but ptot = ", ptot, "!"
                  if (p < pmax) then
                     write (6,*) "porig(p) = ", sph(p)%porig
                  end if
                  write (6,*) "Cell properties:"
                  write (6,*) "c = ", c
                  write (6,*) "ifopen = ", BHhydro(c)%ifopen
                  write (6,*) "nextcell = ", BHhydro(c)%nextcell
                  write (6,*) "leaf = ", BHhydro(c)%leaf
                  write (6,*) "plist = ", BHhydro(c)%plist
                  stop
               end if
               countrefs(p) = countrefs(p) + 1
            end do
            c = BHhydro(c)%nextcell
         else if (BHhydro(c)%leaf == 0) then
            c = BHhydro(c)%ifopen
         else
            c = BHhydro(c)%nextcell
         end if

         if (c > cmax_hydro) exit
      end do

      do p=1,ptot
         if (countrefs(p) > 0 .AND. (p<1.OR.p>phydroend)) then
            write (6,*) "Particle ",p," outside of 1:phydroend! ("&
               &,1,":",phydroend,")"
            stop 1
         else if (countrefs(p) /= 1) then
            write (6,*) "Particle ",p," referenced ", countrefs(p), " times by the hydro tree!"
            stop 1
         end if
      end do
      
      if (any(count_tree_refs(ctot_hydro+1:cmax_hydro) > 0)) then
         write (6,*) "Tree cell > ctot_hydro is referenced!"
         do c=ctot_hydro+1,cmax_hydro
            if (count_tree_refs(c) == 0) cycle
            write (6,*) "Cell properties:"
            write (6,*) "c = ", c
            write (6,*) "ifopen = ", BHhydro(c)%ifopen
            write (6,*) "nextcell = ", BHhydro(c)%nextcell
            write (6,*) "leaf = ", BHhydro(c)%leaf
            write (6,*) "plist = ", BHhydro(c)%plist
         end do
         stop 1
      end if
      
      do c=0,ctot_hydro
         if (count_tree_refs(c) == 0) then
            if (BHhydro(c)%leaf == -1) cycle ! dead cells can be dereferenced
            write (6,*) "Hydro tree cell ", c, " is not referenced!"
            write (6,*) "Cell properties:"
            write (6,*) "c = ", c
            write (6,*) "ifopen = ", BHhydro(c)%ifopen
            write (6,*) "nextcell = ", BHhydro(c)%nextcell
            write (6,*) "leaf = ", BHhydro(c)%leaf
            write (6,*) "plist = ", BHhydro(c)%plist
            do cc=0,ctot_hydro
               if (BHhydro(cc)%ifopen == c) then
                  write (6,*) "Cell ", cc, " ifopen => ", c
                  write (6,*) "Cell properties:"
                  write (6,*) "c = ", cc
                  write (6,*) "ifopen = ", BHhydro(cc)%ifopen
                  write (6,*) "nextcell = ", BHhydro(cc)%nextcell
                  write (6,*) "leaf = ", BHhydro(cc)%leaf
                  write (6,*) "plist = ", BHhydro(cc)%plist
               end if
               if (BHhydro(cc)%nextcell == c) then
                  write (6,*) "Cell ", cc, " nextcell => ", c
                  write (6,*) "Cell properties:"
                  write (6,*) "c = ", cc
                  write (6,*) "ifopen = ", BHhydro(cc)%ifopen
                  write (6,*) "nextcell = ", BHhydro(cc)%nextcell
                  write (6,*) "leaf = ", BHhydro(cc)%leaf
                  write (6,*) "plist = ", BHhydro(cc)%plist
               end if
            end do
            stop 1
         else if (count_tree_refs(c) > 1) then
            if (BHhydro(c)%leaf /= 0) then
               write (6,*) "Hydro tree non-node cell (leaf=",BHhydro(c)%leaf,")&
                  &referenced ", count_tree_refs(c), "times!"
               write (6,*) "Cell properties:"
               write (6,*) "c = ", c
               write (6,*) "ifopen = ", BHhydro(c)%ifopen
               write (6,*) "nextcell = ", BHhydro(c)%nextcell
               write (6,*) "leaf = ", BHhydro(c)%leaf
               write (6,*) "plist = ", BHhydro(c)%plist
               stop 1
            end if
         end if
      end do
      
      deallocate(count_tree_refs)
#endif

#if defined(SELF_GRAVITY)
      allocate(count_tree_refs(0:cmax_grav))
      countrefs = 0
      count_tree_refs = 0
      c = 0
      do
         count_tree_refs(c) = count_tree_refs(c) + 1
         if (BHgrav(c)%leaf > 0) then
            do i=1,BHgrav(c)%leaf
               p = BHgrav(c)%plist(i)
               if (p > ptot) then
                  write (6,*) "Grav cell ", c, " references a particle ", p
                  write (6,*) "but ptot = ", ptot, "!"
                  if (p < pmax) then
                     write (6,*) "porig(p) = ", sph(p)%porig
                  end if
                  write (6,*) "Cell properties:"
                  write (6,*) "c = ", c
                  write (6,*) "ifopen = ", BHgrav(c)%ifopen
                  write (6,*) "nextcell = ", BHgrav(c)%nextcell
                  write (6,*) "leaf = ", BHgrav(c)%leaf
                  write (6,*) "plist = ", BHgrav(c)%plist
                  stop
               end if
               countrefs(p) = countrefs(p) + 1
            end do
            c = BHgrav(c)%nextcell
         else if (BHgrav(c)%leaf == 0) then
            c = BHgrav(c)%ifopen
         else
            c = BHgrav(c)%nextcell
         end if

         if (c > cmax_grav) exit
      end do

      do p=1,ptot
         if (countrefs(p) > 0 .AND. (p<pgravitystart.OR.p>pgravityend)) then
            write (6,*) "Particle ",p," outside of pgravitystart:pgravityend! ("&
               &,pgravitystart,":",pgravityend,")"
            stop
         else if (countrefs(p) /= 1 .AND. (p>=pgravitystart.AND.p<=pgravityend)) then
            write (6,*) "Particle ",p," referenced ", countrefs(p), " times by the gravity tree!"
            stop
         end if
      end do
      
      if (any(count_tree_refs(ctot_grav+1:cmax_grav) > 0)) then
         write (6,*) "Tree cell > ctot_grav is referenced!"
         do c=ctot_grav+1,cmax_grav
            if (count_tree_refs(c) == 0) cycle
            write (6,*) "Cell properties:"
            write (6,*) "c = ", c
            write (6,*) "ifopen = ", BHgrav(c)%ifopen
            write (6,*) "nextcell = ", BHgrav(c)%nextcell
            write (6,*) "leaf = ", BHgrav(c)%leaf
            write (6,*) "plist = ", BHgrav(c)%plist
         end do
         stop 1
      end if
      
      do c=0,ctot_grav
         if (count_tree_refs(c) == 0) then
            if (BHgrav(c)%leaf == -1) cycle ! dead cells can be dereferenced
            write (6,*) "Grav tree cell ", c, " is not referenced!"
            write (6,*) "Cell properties:"
            write (6,*) "c = ", c
            write (6,*) "ifopen = ", BHgrav(c)%ifopen
            write (6,*) "nextcell = ", BHgrav(c)%nextcell
            write (6,*) "leaf = ", BHgrav(c)%leaf
            write (6,*) "plist = ", BHgrav(c)%plist
            do cc=0,ctot_grav
               if (BHgrav(cc)%ifopen == c) then
                  write (6,*) "Cell ", cc, " ifopen => ", c
                  write (6,*) "Cell properties:"
                  write (6,*) "c = ", cc
                  write (6,*) "ifopen = ", BHhydro(cc)%ifopen
                  write (6,*) "nextcell = ", BHhydro(cc)%nextcell
                  write (6,*) "leaf = ", BHhydro(cc)%leaf
                  write (6,*) "plist = ", BHhydro(cc)%plist
               end if
               if (BHgrav(cc)%nextcell == c) then
                  write (6,*) "Cell ", cc, " nextcell => ", c
                  write (6,*) "Cell properties:"
                  write (6,*) "c = ", cc
                  write (6,*) "ifopen = ", BHhydro(cc)%ifopen
                  write (6,*) "nextcell = ", BHhydro(cc)%nextcell
                  write (6,*) "leaf = ", BHhydro(cc)%leaf
                  write (6,*) "plist = ", BHhydro(cc)%plist
               end if
            end do
            stop 1
         else if (count_tree_refs(c) > 1) then
            if (BHgrav(c)%leaf /= 0) then
               write (6,*) "Grav tree non-node cell (leaf=",BHgrav(c)%leaf,")&
                  &referenced ", count_tree_refs(c), "times!"
               write (6,*) "Cell properties:"
               write (6,*) "c = ", c
               write (6,*) "ifopen = ", BHgrav(c)%ifopen
               write (6,*) "nextcell = ", BHgrav(c)%nextcell
               write (6,*) "leaf = ", BHgrav(c)%leaf
               write (6,*) "plist = ", BHgrav(c)%plist
               stop 1
            end if
         end if
      end do
      
      deallocate(count_tree_refs)
#endif

      return
   end subroutine check_tree_integrity

   subroutine check_overlap_hydro_tree(BHtree,cmax_tree)
      use tree_module
      use particle_module
      implicit none

      integer, intent(in)            :: cmax_tree ! Total size of tree
      type(BHhydro_node), intent(in) :: BHtree(0:cmax_tree) ! Tree to build
      integer  :: c, p, pp, i, j
      logical  :: OK
      real(kind=PR) :: rp(1:NDIM), rpp(1:NDIM), drsqd, space_p, space_pp, space_sqd

      OK = .TRUE.
      c = 0
      do
         if (BHtree(c)%leaf > 0) then
            do i=1,BHtree(c)%leaf
               p = BHtree(c)%plist(i)
               rp = sph(p)%r
               space_p = maxval(spacing(rp))
               do j=i+1,BHtree(c)%leaf
                  pp = BHtree(c)%plist(j)
                  rpp = sph(pp)%r
                  space_pp = maxval(spacing(rp))
                  drsqd = sum((rp-rpp)**2)
                  space_sqd = (max(space_p,space_pp))**2
                  if (drsqd < 2._PR*space_sqd) then
                     write (6,*) "Possible collision between particle ", p, " and ", pp, "!"
                     write (6,*) "rp = ", rp
                     write (6,*) "rpp = ", rpp
                     write (6,*) "drsqd = ", drsqd, "; drmag = ", sqrt(drsqd)
                     write (6,*) "space_p = ", space_p, "; space_pp = ", space_pp
                     OK = .FALSE.
                  end if
               end do
            end do
            c = BHtree(c)%nextcell
         else if (BHtree(c)%leaf == 0) then
            c = BHtree(c)%ifopen
         else
            c = BHtree(c)%nextcell
         end if

         if (c > cmax_tree) exit
      end do

      if (.NOT. OK) stop

      write (6,*) "Overlap of particles in tree verified (within leaves)"

      return
   end subroutine check_overlap_hydro_tree
#endif

! ----------------------------------------------------------------------------
! MPI-only subroutines

#if defined(USE_MPI)
   subroutine check_particles_in_boxes
      use particle_module
      use domain_comparisons
      use mpi_communication_module
      implicit none
      integer         :: p, d
      logical         :: OK
      real(kind=PR)   :: rp(1:NDIM)

      debug_timing("CHECK PARTICLES IN BOXES")

      OK = .TRUE.

      if (any(domain_bbmin(1:NDIM,rank) == domain_bbmax(1:NDIM,rank))) then
         write (6,*) "This domain's boundaries are invalid!"
         write (6,'(A, 3F10.7)') "Domain min: ", domain_bbmin(1:NDIM,rank)
         write (6,'(A, 3F10.7)') "Domain max: ", domain_bbmax(1:NDIM,rank)
      end if

      do p=1,ptot
         rp = sph(p)%r
#if defined(PERIODIC)
         call unwrap_particle_position(rp)
#endif
            if (.NOT. inside_domain_box(rp,rank)) then
            write (6,*) "Rank: ",rank," - Particle ", p, " outside domain box at ", rp, "!"
            OK = .FALSE.
         end if
      end do

      if (.NOT. OK) then
         write (6,*) "Particles outside box!"
         write (6,'(A, 3F10.7)') "Domain min: ", domain_bbmin(1:NDIM,rank)
         write (6,'(A, 3F10.7)') "Domain max: ", domain_bbmax(1:NDIM,rank)
         stop
      else
         if (rank==0) write (6,*) "Checked particles in boxes"
      end if

      return
   end subroutine check_particles_in_boxes
#endif

end module debug_tests