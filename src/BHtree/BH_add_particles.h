SUBROUTINE BH_ADD_PART_SUB_NAME(n_new,new_list,cmax,ctot,ltot,&
                            &first_cell,last_cell,BHtree)
  use interface_module, only : distance3
  use definitions
  use particle_module
  use tree_module
  implicit none

  integer, intent(in) :: n_new                         ! No. of new particles
  integer, intent(in) :: new_list(1:n_new)             ! List of new particles
  integer, intent(in) :: cmax                          ! Max. no of cells
  integer, intent(inout) :: ctot                       ! Total no of cells
  integer, intent(in)    :: ltot                       ! Total no. of levels
  integer, intent(inout) :: first_cell(0:LMAX)         ! 1st cell in level
  integer, intent(inout) :: last_cell(0:LMAX)          ! Last cell in level
  type(BH_ADD_PART_TYPE), intent(inout) :: BHtree(0:cmax)  ! Tree array

  integer :: p                                ! Particle we are adding
  integer :: i                                ! Auxilary counter
  integer :: l                                ! Level counter
  integer :: c                                ! Cell counter
  integer :: cc                               ! Auxiliary cell counter
  integer :: pp                               ! Auxiliary leaf counter
  integer :: leaf                             ! Leaf number of cell
  integer :: start_level                      ! First cell on level
  integer :: end_level                        ! Last cell on level
  real(kind=PR) :: dr(1:NDIM)                 ! Relative displacement vector
  real(kind=PR) :: drsqd                      ! Distance squared
  real(kind=PR) :: maxdistsqd                 ! Maximum gather dist. squared
  real(kind=PR) :: rp(1:NDIM)                 ! Position of particle p
  real(kind=PR) :: least_dist_sqd             ! Smallest distance to a cell sqd
  integer :: least_dist_c                     ! Cell with smallest distance

#if defined(DEBUG_BH_ADD_PARTICLES)
  write (6,*) "Adding ", n_new, " new particles to tree..."
  write (6,*) new_list
#endif

! Main particle loop
! ============================================================================
  p_loop: do i=1,n_new
     p = new_list(i)
#if defined(DEBUG_BH_ADD_PARTICLES)
     write (6,*) "adding new particle ", p
#endif
     rp = sph(p)%r
#if defined(DEBUG_BH_ADD_PARTICLES)
     write (6,*) "rp = ", rp
#endif

     ! Start on first level
     l = 0
     c = BHtree(0)%ifopen
     
#if defined(DEBUG_BH_ADD_PARTICLES)
     write (6,*) "first cell = 0"
     write (6,*) "then ifopen = ", c
#endif

   ! Main level loop
   ! =========================================================================
     level_loop: do
        l = l + 1
#if defined(DEBUG_BH_ADD_PARTICLES)
        write (6,*) "level = ", l
#endif
        start_level = first_cell(l)
        end_level = last_cell(l)
#if defined(DEBUG_BH_ADD_PARTICLES)
        write (6,*) "start_level = ", start_level
        write (6,*) "end_level = ", end_level
#endif
        least_dist_sqd = BIG_NUMBER
        least_dist_c = -1
        
        do ! child cell loop
#if defined(DEBUG_BH_ADD_PARTICLES)
           write (6,*) "cell ", c
           if (c > ctot) then
              write (6,*) "Tree cell outside of ctot!"
              call flush(6)
              stop
           end if
#endif
           maxdistsqd = BH_ADD_PART_DIST
           leaf = BHtree(c)%leaf
#if defined(DEBUG_BH_ADD_PARTICLES)
           write (6,*) "leaf = ", leaf
#endif
           ! Compute distance from particle to cell centre, plus the maximum 
           ! radial extent (squared) of all particles within the cell
#if defined(PERIODIC) && !defined(GHOST_PARTICLES)
           call distance3(rp(1:NDIM),BHtree(c)%r(1:NDIM),dr(1:NDIM),drsqd)
#else
           dr(1:NDIM) = BHtree(c)%r(1:NDIM) - rp(1:NDIM)
           drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif
#if defined(DEBUG_BH_ADD_PARTICLES)
           write (6,*) "maxdistsqd = ", maxdistsqd
           write (6,*) "drsqd = ", drsqd
#endif
#if (BH_ADD_PART_GRAV == 0)
           ! If this is a gravity tree, there is no maxdistsqd
           ! so none of this can ever happen
           if (drsqd < maxdistsqd) then
#if defined(DEBUG_BH_ADD_PARTICLES)
              write (6,*) "Let's use this cell"
#endif
              ! If this is a non-full leaf cell or a dead cell, insert here
              ! If this is a node cell, step into now
              if (leaf == 0) then
                 ! Step into
#if defined(DEBUG_BH_ADD_PARTICLES)
                 if (BHtree(c)%ifopen > ctot) then
                    write (6,*) "ifopen tree cell outside of ctot!"
                    write (6,*) "c = ", c
                    write (6,*) "ifopen = ", BHtree(c)%ifopen
                    write (6,*) "nextcell = ", BHtree(c)%nextcell
                    write (6,*) "leaf = ", BHtree(c)%leaf
                    write (6,*) "plist = ", BHtree(c)%plist
                    call flush(6)
                    stop
                 end if
#endif
                 c = BHtree(c)%ifopen
#if defined(DEBUG_BH_ADD_PARTICLES)
                 write (6,*) "Node cell, stepping into, ifopen = ", c
#endif
                 cycle level_loop
              elseif (leaf < LEAFMAX) then
                 ! Insert particle!
                 if (leaf == -1) then
#if defined(DEBUG_BH_ADD_PARTICLES)
                    write (6,*) "Dead cell - bingo!"
#endif
                    ! Dead cell
                    leaf = 1
                 else
#if defined(DEBUG_BH_ADD_PARTICLES)
                    write (6,*) "Non-full leaf cell - bingo!"
#endif
                    ! Non-dead cell
                    leaf = leaf + 1
                 end if
                 BHtree(c)%leaf = leaf
                 BHtree(c)%plist(leaf) = p
#if defined(DEBUG_BH_ADD_PARTICLES)
                 write (6,*) "Cell properties:"
                 write (6,*) "c = ", c
                 write (6,*) "ifopen = ", BHtree(c)%ifopen
                 write (6,*) "nextcell = ", BHtree(c)%nextcell
                 write (6,*) "leaf = ", BHtree(c)%leaf
                 write (6,*) "plist = ", BHtree(c)%plist
#endif
                 cycle p_loop
#if defined(DEBUG_BH_ADD_PARTICLES)
              else
                 write (6,*) "We found a great cell, but it was full. Continue..."
#endif
              end if
           end if
#endif
           
           if (drsqd < least_dist_sqd .AND. leaf < LEAFMAX) then
#if defined(DEBUG_BH_ADD_PARTICLES)
              write (6,*) "Currently the closest cell that is not a full &
                          &leaf cell. Record..."
#endif
              ! If this is the closest node or non-full leaf cell
              ! we have found so far, record
              least_dist_sqd = drsqd
              least_dist_c = c
#if defined(DEBUG_BH_ADD_PARTICLES)
           else if (drsqd < least_dist_sqd) then
               write (6,*) "Reject this cell as it is full!"
           else
               write (6,*) "Do not record as drsqd (", drsqd, &
                  &") >= least_dist_sqd (", least_dist_sqd, ")"
#endif
           end if
           
#if defined(DEBUG_BH_ADD_PARTICLES)
           if (BHtree(c)%nextcell == ctot) write (6,*) "this seems a bad idea"
#endif
           
           if (BHtree(c)%nextcell < start_level .OR. &
              &BHtree(c)%nextcell > ctot) then
#if defined(DEBUG_BH_ADD_PARTICLES)
              write (6,*) "next cell takes us somewhere else: ", BHtree(c)%nextcell
#endif
              ! We have run over all child cells of our entry point
              ! Decide on a winner
              if (least_dist_c == -1) then
                 ! No non-full leaf cells on this level
                 ! Wholescale tree organization required :(
#if defined(DEBUG_BH_ADD_PARTICLES)
                 write (6,*) "Having to add new cell..."
#endif
                 if (ctot + 1 == cmax) stop "Out of cells!"
                 ! Add another cell. New cell id = c + 1
#if defined(DEBUG_BH_ADD_PARTICLES)
                 write (6,*) "Old cell properties:"
                 write (6,*) "c - 1 = ", c - 1
                 write (6,*) "ifopen = ", BHtree(c-1)%ifopen
                 write (6,*) "nextcell = ", BHtree(c-1)%nextcell
                 write (6,*) "leaf = ", BHtree(c-1)%leaf
                 write (6,*) "plist = ", BHtree(c-1)%plist
                 write (6,*) "c = ", c
                 write (6,*) "ifopen = ", BHtree(c)%ifopen
                 write (6,*) "nextcell = ", BHtree(c)%nextcell
                 write (6,*) "leaf = ", BHtree(c)%leaf
                 write (6,*) "plist = ", BHtree(c)%plist
                 write (6,*) "c + 1 = ", c + 1
                 write (6,*) "ifopen = ", BHtree(c+1)%ifopen
                 write (6,*) "nextcell = ", BHtree(c+1)%nextcell
                 write (6,*) "leaf = ", BHtree(c+1)%leaf
                 write (6,*) "plist = ", BHtree(c+1)%plist
                 write (6,*) "New cell id = c + 1 = ", c + 1
                 write (6,*) "old first_cell: ", first_cell(0:ltot)
                 write (6,*) "old last_cell: ", last_cell(0:ltot)
#endif
                 ! Increase first_cell/last_cell appropriately
                 where (first_cell > c) first_cell = first_cell + 1
                 where (last_cell > c) last_cell = last_cell + 1
#if defined(DEBUG_BH_ADD_PARTICLES)
                 write (6,*) "new first_cell: ", first_cell(0:ltot)
                 write (6,*) "new last_cell: ", last_cell(0:ltot)
#endif
                 ! Adjust ifopen and nextcell
                 do cc=0,ctot
                    if (BHtree(cc)%ifopen > c) then
                       if (BHtree(cc)%ifopen /= LARGEST_INT) then
                          BHtree(cc)%ifopen = BHtree(cc)%ifopen + 1
                       end if
                    end if
                    if (BHtree(cc)%nextcell > c) then
                       if (BHtree(cc)%nextcell /= LARGEST_INT) then
                          BHtree(cc)%nextcell = BHtree(cc)%nextcell + 1
                       end if
                    end if
                 end do
                 ! Shuffle along cells.
                 do cc=ctot,c+1,-1
                    BHtree(cc + 1) = BHtree(cc)
                 end do
                 BHtree(c+1)%nextcell = BHtree(c)%nextcell
                 BHtree(c+1)%ifopen = LARGEST_INT ! DELETE ME (shouldn't matter)
                 BHtree(c)%nextcell = c + 1
                 leaf = 1
                 BHtree(c+1)%leaf = leaf
                 BHtree(c+1)%plist(leaf) = p
                 ctot = ctot + 1
#if defined(DEBUG_BH_ADD_PARTICLES)
                 write (6,*) "New cell properties:"
                 write (6,*) "c - 1 = ", c - 1
                 write (6,*) "ifopen = ", BHtree(c-1)%ifopen
                 write (6,*) "nextcell = ", BHtree(c-1)%nextcell
                 write (6,*) "leaf = ", BHtree(c-1)%leaf
                 write (6,*) "plist = ", BHtree(c-1)%plist
                 write (6,*) "c = ", c
                 write (6,*) "ifopen = ", BHtree(c)%ifopen
                 write (6,*) "nextcell = ", BHtree(c)%nextcell
                 write (6,*) "leaf = ", BHtree(c)%leaf
                 write (6,*) "plist = ", BHtree(c)%plist
                 write (6,*) "c + 1 = ", c + 1
                 write (6,*) "ifopen = ", BHtree(c+1)%ifopen
                 write (6,*) "nextcell = ", BHtree(c+1)%nextcell
                 write (6,*) "leaf = ", BHtree(c+1)%leaf
                 write (6,*) "plist = ", BHtree(c+1)%plist
#endif
                 cycle p_loop
              else
                 ! Use cell least_dist_c
                 c = least_dist_c
#if defined(DEBUG_BH_ADD_PARTICLES)
                 write (6,*) "Having to use non-ideal cell ", c
#endif
                 leaf = BHtree(c)%leaf
                 if (leaf == 0) then
#if defined(DEBUG_BH_ADD_PARTICLES)
                 write (6,*) "Node cell, stepping into, ifopen = ", c
#endif
                    ! Step into
#if defined(DEBUG_BH_ADD_PARTICLES)
                    if (BHtree(c)%ifopen > ctot) then
                       write (6,*) "ifopen tree cell outside of ctot!"
                       write (6,*) "c = ", c
                       write (6,*) "ifopen = ", BHtree(c)%ifopen
                       write (6,*) "nextcell = ", BHtree(c)%nextcell
                       write (6,*) "leaf = ", BHtree(c)%leaf
                       write (6,*) "plist = ", BHtree(c)%plist
                       call flush(6)
                       stop
                    end if
#endif
                    c = BHtree(c)%ifopen
                    cycle level_loop
                 elseif (leaf < LEAFMAX) then
                    ! Insert particle!
                    if (leaf == -1) then
#if defined(DEBUG_BH_ADD_PARTICLES)
                       write (6,*) "Dead cell - bingo!"
#endif
                       ! Dead cell
                       leaf = 1
                    else
#if defined(DEBUG_BH_ADD_PARTICLES)
                       write (6,*) "Non-full leaf cell - bingo!"
#endif
                       ! Non-dead cell
                       leaf = leaf + 1
                    end if
                    BHtree(c)%leaf = leaf
                    BHtree(c)%plist(leaf) = p
#if defined(DEBUG_BH_ADD_PARTICLES)
                    write (6,*) "Cell properties:"
                    write (6,*) "c = ", c
                    write (6,*) "ifopen = ", BHtree(c)%ifopen
                    write (6,*) "nextcell = ", BHtree(c)%nextcell
                    write (6,*) "leaf = ", BHtree(c)%leaf
                    write (6,*) "plist = ", BHtree(c)%plist
#endif
                    cycle p_loop
                 else
                    stop "If cell was full, shouldn't have added it!"
                 end if
              end if
           end if
           
           ! If we haven't come to the end of this set of child cells, continue
#if defined(DEBUG_BH_ADD_PARTICLES)
           if (BHtree(c)%nextcell > ctot) then
              write (6,*) "nextcell tree cell outside of ctot!"
              write (6,*) "c = ", c
              write (6,*) "ifopen = ", BHtree(c)%ifopen
              write (6,*) "nextcell = ", BHtree(c)%nextcell
              write (6,*) "leaf = ", BHtree(c)%leaf
              write (6,*) "plist = ", BHtree(c)%plist
              call flush(6)
              stop
           end if
#endif
           c = BHtree(c)%nextcell
           
        end do
     
     end do level_loop
! ============================================================================
  end do p_loop

#if defined(DEBUG_BH_ADD_PARTICLES)
  write (6,*) "we are done; BHtree(0)%ifopen = ", BHtree(0)%ifopen
#endif
  
  return
END SUBROUTINE BH_ADD_PART_SUB_NAME
