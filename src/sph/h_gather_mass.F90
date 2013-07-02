! H_GATHER_MASS.F90
! C. P. Batty & D. A. Hubber - 11/12/2006
! Calculates smoothing length of particle p using pure gather method by 
! setting h so kernel contains at least a specified amount of mass 
! (p_gather * mp).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE h_gather(p,hguess,rp)
  use interface_module, only : insertion_sort_real
  use particle_module
  use neighbour_module
  use kernel_module
  use time_module
  implicit none

  integer, intent(in) :: p                ! chosen particle
  real(kind=PR), intent(inout) :: hguess  ! guess of smoothing length
  real(kind=PR), intent(in) :: rp(1:NDIM) ! position of particle p

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
  real(kind=PR) :: hold                   ! old value of h
  real(kind=PR) :: h_search               ! value of h for tree search  
  real(kind=PR) :: mgather                ! ..
  real(kind=PR) :: mhigh                  ! ..
  real(kind=PR) :: mlow                   ! ..
  real(kind=PR) :: mp                     ! ..
  real(kind=PR) :: mpp                    ! ..
  real(kind=PR) :: mtemp                  ! ..
  real(kind=PR) :: mtot                   ! ..
  real(kind=PR), allocatable :: mlist(:)      ! ..
  real(kind=PR), allocatable :: mlist_high(:) ! ..
  real(kind=PR), allocatable :: rsqd_high(:)  ! drsqd values for high list

  debug3("Calculating h [h_gather.F90] for particle ", p)

  mp = sph(p)%m
  mgather = mp*real(pp_gather,PR)

  pp_pot   = 0
  pp_max   = LISTSIZE
  hold     = hguess
  h_high   = HMULT*hold
  h_low    = hold / HMULT
  h_search = 0.0_PR
  hrangesqd_high = KERNRANGESQD*h_high*h_high
  allocate(mlist(1:pp_max))
  allocate(mlist_high(1:pp_max))
  allocate(pp_potlist(1:pp_max))
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
#if defined(BH_TREE)
           call BHhydrowalk_hgather(rp(1:NDIM),KERNRANGE*h_search,&
                &pp_pot,pp_max,pp_potlist)
#elif defined(BINARY_TREE)
           call binary_neibfind(p,h_search,pp_pot,pp_max,pp_potlist)
#else
           do pp=1,ptot
              if (pp == p) cycle
              call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
              if (drsqd < hrangesqd_high) then
                 pp_pot = pp_pot + 1
                 pp_potlist(pp_pot) = pp
              endif
           end do
#endif
        end if
           
        ! Precaution against ridiculously high neighbour numbers
        if (pp_pot < 0) then
           h_search = 0.0_PR
           pp_pot   = 0
           pp_max   = ptot
           deallocate(rsqd_high)
           deallocate(pp_high_list)
           deallocate(pp_potlist)
           deallocate(mlist_high)
           deallocate(mlist)
           allocate(mlist(1:pp_max))
           allocate(mlist_high(1:pp_max))
           allocate(pp_potlist(1:pp_max))
           allocate(pp_high_list(1:pp_max))
           allocate(rsqd_high(1:pp_max))
        else if (pp_pot <= pp_gather) then
           h_search = 0.0_PR
           pp_pot   = 0
           h_low    = h_high
           h_high   = HMULT*h_high
           hrangesqd_high = KERNRANGE*KERNRANGE*h_high*h_high
        end if
        
        if (pp_pot > 0) exit
     end do
     ! -----------------------------------------------------------------------

     mhigh   = 0.0_PR
     mlow    = 0.0_PR
     mtot    = 0.0_PR
     pp_tot  = 0
     pp_high = 0
     pp_low  = 0
     hrangesqd_high = KERNRANGE*KERNRANGE*h_high*h_high
     hrangesqd_low  = KERNRANGE*KERNRANGE*h_low*h_low

     ! Loop over potential neighbour list and find the number of particles 
     ! with r < 2*h_low and the number and list of particles with 
     ! 2*h_low < r < 2*h_high
     ! -----------------------------------------------------------------------
     do i=1,pp_pot
        pp = pp_potlist(i)
        if (pp == p) cycle
        call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
        mpp = sph(pp)%m
        
        ! Is potential neighbour within 2*h_high of particle p
        if (drsqd < hrangesqd_high) then   
           
           ! Now find particles within 2*h_low and 2*h_low < r < 2*h_high
           if (drsqd < hrangesqd_low) then  
              mlow = mlow + mpp
              pp_low = pp_low + 1
           else
              pp_high = pp_high + 1
              pp_high_list(pp_high) = pp
              rsqd_high(pp_high) = drsqd        
              mhigh = mhigh + mpp
              mlist_high(pp_high) = mpp
           end if
        end if
        if (pp_high == pp_max) exit
     end do
     ! -----------------------------------------------------------------------
     
     pp_tot = pp_low + pp_high   
     mtot = mlow + mhigh
     
     ! If potential list is too small or too large, adjust h limits
     ! -----------------------------------------------------------------------
     if (mtot < mgather) then
        h_low = h_high
        h_high = HMULT * h_high
        hrangesqd_high = KERNRANGE*KERNRANGE*h_high*h_high
        cycle
        
     ! else if h_low is too large, decrease limits
     else if (mlow >= mgather) then
        h_high = h_low
        h_low = h_low / HMULT
        cycle
        
     ! If the limits contain enough neighbours, exit loop
     else if (mlow < mgather .and. mtot >= mgather) then 
        exit
     end if
     ! -----------------------------------------------------------------------
     
     ! If minimum h flag is set and it is clear that h is less than hmin,
     ! set h to hmin here and return
#if defined(MINIMUM_H)
     if (mtot >= mgather .and. h_high < hmin) then
        hguess = hmin
        return
     end if
#endif

  end do
! ============================================================================


! Sort potential neighbours to determine correct h
  call insertion_sort_real(pp_high,pp_high_list(1:pp_high),rsqd_high(1:pp_high))
  mtemp = mlow

! Loop through high list and find appropriate h
  do i=1,pp_high
     mtemp = mtemp + mlist_high(i)
     if (mtemp < mgather) cycle
     drsqd = rsqd_high(i)
     hguess = 0.5_PR*sqrt(drsqd)
     exit
  end do

! If minimum h flag is set, ensure h >= hmin 
#if defined(MINIMUM_H)
  if (hguess < hmin) hguess = hmin
#endif

#if defined(DEBUG_H_GATHER)
  write(6,*) "Number of neighbours for particle", p, "=", pp_gather
  write(6,*) "Smoothing length =", hguess,h_low,h_high,hold
  if ( (abs(hguess - hold)/hold) > 0.1) &
  	& write(6,*) "Smoothing length changed by more than 10%"
#endif

  deallocate(rsqd_high)
  deallocate(pp_high_list)
  deallocate(pp_potlist)
  deallocate(mlist_high)
  deallocate(mlist)

  return
END SUBROUTINE h_gather
