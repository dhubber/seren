! H_GATHER.F90
! C. P. Batty & D. A. Hubber - 11/12/2006
! Calculates smoothing length of particle p using pure gather method by 
! setting h so the kernel contains exactly pp_gather neighbours.  Uses old 
! value of h as first guess for new h.  Then walks the tree to find all 
! potential neighbours within search limits. Finally performs a limited 
! sort of the potential neighbours to obtain exact value of h.  
! Returns new value of h as subroutine parameter. 
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE h_gather(p,hguess,rp,typemask)
  use interface_module, only : distance2,gather_neib_on_fly,&
       &heapsort_real,insertion_sort_real
  use particle_module
  use neighbour_module
  use kernel_module
  use time_module
  use type_module
  use tree_module
  implicit none

  integer, intent(in) :: p                ! chosen particle
  real(kind=PR), intent(inout) :: hguess  ! guess of smoothing length
  real(kind=PR), intent(in) :: rp(1:NDIM) ! position of particle p
  logical, optional, intent(in) :: typemask(1:nmaxtypes) ! part. types to include?

  integer :: i                            ! counter in neighbour search
  integer :: pp                           ! neighbouring particles (p')
  integer :: pp_high                      ! no. of neighbours within 2 h_high
  integer :: pp_low                       ! no. of neighbours within 2 h_low
  integer :: pp_max                       ! length of pot list array
  integer :: pp_pot                       ! no. of potential neighbours
  integer :: pp_tot                       ! no. of neighbours 
  integer, allocatable :: pp_high_list(:) ! list of neighbours < 2 h_high
  integer, allocatable :: pp_potlist(:)   ! potential neighbour list
  real(kind=PR) :: dr(1:NDIM)             ! vector displacements (p'-p)
  real(kind=PR) :: drsqd                  ! particle separation squared
  real(kind=PR) :: hrangesqd_high         ! particle radius (2*h_high) squared
  real(kind=PR) :: hrangesqd_low          ! particle radius (2*h_low) squared
  real(kind=PR) :: h_high                 ! high value estimate of h
  real(kind=PR) :: h_low                  ! low value of h
  real(kind=PR) :: h_old                  ! old value of h
  real(kind=PR) :: h_search               ! value of h for tree search  
  real(kind=PR), allocatable :: rsqd_high(:)  ! drsqd values for high list

  debug3("Calculating h [h_gather.F90] for particle ", p)

  pp_pot   = 0
  pp_max   = LISTSIZE
  h_old    = hguess
  h_high   = HMULT*h_old
  h_low    = h_old/HMULT
  h_search = 0.0_PR
  hrangesqd_high = KERNRANGESQD*h_high*h_high
  allocate(pp_high_list(1:pp_max))
  allocate(rsqd_high(1:pp_max))


! Main searching loop
! ============================================================================
  do

     ! Obtain potential neighbour list either by direct sum or the tree
     ! -----------------------------------------------------------------------
     do 

        if (h_high > h_search) then
           h_search = HMULT*h_high
           pp_pot = 0
           if (allocated(rsqd_high)) deallocate(rsqd_high)
           if (allocated(pp_high_list)) deallocate(pp_high_list)
           if (allocated(pp_potlist)) deallocate(pp_potlist)
           
           ! Obtain list only from certain particle types if mask defined
           if (present(typemask)) then
              call gather_neib_on_fly(p,pp_max,pp_pot,pp_potlist,&
                   &rp(1:NDIM),KERNRANGE*h_search,typemask)
           else
              call gather_neib_on_fly(p,pp_max,pp_pot,pp_potlist,&
                   &rp(1:NDIM),KERNRANGE*h_search)
           end if

           ! Increase pp_max to ptot in spurious cases which return 
           ! ridiculously large potential neighbour lists.
           if (pp_pot < 0) then
              h_search = 0.0_PR
              pp_pot   = 0
              pp_max   = ptot
           else if (pp_pot <= pp_gather) then
              pp_pot   = 0
              h_search = 0.0_PR
              h_low    = h_high
              h_high   = HMULT*h_high
              hrangesqd_high = KERNRANGESQD*h_high*h_high
           end if
        end if

        if (pp_pot > 0) exit
     end do
     ! -----------------------------------------------------------------------

     pp_tot  = 0
     pp_high = 0
     pp_low  = 0
     hrangesqd_high = KERNRANGESQD*h_high*h_high
     hrangesqd_low  = KERNRANGESQD*h_low*h_low
     if (.not. allocated(pp_high_list)) allocate(pp_high_list(1:pp_max))
     if (.not. allocated(rsqd_high)) allocate(rsqd_high(1:pp_max))


     ! Loop over potential neighbour list and find the number of particles 
     ! with r < 2*h_low and the number and list of particles with 
     ! 2*h_low < r < 2*h_high
     ! -----------------------------------------------------------------------
     do i=1,pp_pot
        pp = pp_potlist(i)
        if (pp == p) cycle
        call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
        
        ! Is potential neighbour within KERNRANGE*h_high of particle p
        if (drsqd < hrangesqd_high) then   
           
           ! Now find particles within KERNRANGE*h_low and 
           ! KERNRANGE*h_low < r < 2*h_high
           if (drsqd < hrangesqd_low) then  
              pp_low = pp_low + 1
           else                         
              pp_high = pp_high + 1
              pp_high_list(pp_high) = pp
              rsqd_high(pp_high) = drsqd
           end if
        end if
        if (pp_high == pp_max) exit
     end do
     ! -----------------------------------------------------------------------
     
     pp_tot = pp_low + pp_high   
     
     ! If potential list is too small or too large, adjust h limits
     ! -----------------------------------------------------------------------
     if (pp_tot < pp_gather) then
        h_low  = h_high
        h_high = HMULT*h_high
        hrangesqd_high = KERNRANGESQD*h_high*h_high
        cycle
        
     ! else if h_low is too large, decrease limits
     else if (pp_low >= pp_gather) then
        h_high = h_low
        h_low  = h_low/HMULT
        cycle
        
     ! If the limits contain enough neighbours, exit loop
     else if (pp_low < pp_gather .and. pp_tot >= pp_gather) then 
        exit

     end if
     ! -----------------------------------------------------------------------
     
     ! If minimum h flag is set and it is clear that h is less than hmin,
     ! set h to hmin here and return
#if defined(MINIMUM_H)
     if (pp_tot >= pp_gather .and. h_high < hmin) then
        hguess = hmin
        return
     end if
#endif

  end do
! ============================================================================


! Sort h_high list to determine h
#if defined(HEAPSORT)
  call heapsort_real(pp_high,pp_high_list,rsqd_high)
#else
  call insertion_sort_real(pp_high,&
       &pp_high_list(1:pp_high),rsqd_high(1:pp_high))
#endif
  drsqd  = rsqd_high(pp_gather - pp_low)
  hguess = INVKERNRANGE*sqrt(drsqd)

! If minimum h flag is set, ensure h >= hmin 
#if defined(MINIMUM_H)
  if (hguess < hmin) hguess = hmin
#endif

! Record potential neighbour list in order to speed up SPH calculations
#if defined(NEIGHBOUR_LISTS)
  pptot(p) = pp_pot
  if (pp_pot <= pp_limit) then
     pplist(1:pp_pot,p) = pp_potlist(1:pp_pot)
  end if
#endif

#if defined(DEBUG_H_GATHER)
  write(6,*) "Number of neighbours for particle", p, "=", pp_tot
  write(6,*) "pp_low : ",pp_low,"  pp_high : ",pp_high
  write(6,*) "Smoothing length =", hguess,h_low,h_high,h_old
  if ((abs(hguess - h_old)/h_old) > 0.1_PR) &
       & write(6,*) "Smoothing length changed by more than 10%"
#endif

  deallocate(rsqd_high)
  deallocate(pp_potlist)
  deallocate(pp_high_list)

  return
END SUBROUTINE h_gather
