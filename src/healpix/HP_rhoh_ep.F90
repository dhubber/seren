! HP_RHOH_EP.F90
! C. P. Batty & D. A. Hubber - 11/12/2006
! Calculates smoothing length of particle p  using pure gather method by 
! setting h so kernel contains exactly pp_gather neighbours.  Uses old value 
! of h as first guess for new h.  Then walks the tree to find all potential 
! neighbours within search limits. Finally performs a limited sort of 
! particles to obtain exact value of h.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE HP_rhoh_ep(rep,hep,rhoep)
  use interface_module, only : binary_gather_walk,BHhydrowalk_hgather,&
       &distance2,heapsort_real,insertion_sort_real,w0,w2
  use particle_module
  use neighbour_module
  use kernel_module
  use time_module
  use tree_module
  implicit none

  real(kind=PR), intent(in) :: rep(1:NDIM)        ! position of e.p.  
  real(kind=PR), intent(inout) :: hep             ! smoothing length
  real(kind=PR), intent(inout) :: rhoep           ! density at position rep

  integer :: i                            ! counter in neighbour search
  integer :: pp                           ! neighbouring particles (p')
  integer :: pp_high                      ! no. of neighbours within 2 h_high
  integer :: pp_low                       ! no. of neighbours within 2 h_low
  integer :: pp_max                       ! length of pot list array
  integer :: pp_pot                       ! no. of potential neighbours
  integer :: pp_tot                       ! no. of neighbours 
  integer, allocatable :: pp_low_list(:)  ! list of neighbours < 2 h_low
  integer, allocatable :: pp_high_list(:) ! list of neighbours < 2 h_high
  integer, allocatable :: pp_potlist(:)   ! potential neighbour list
  real(kind=PR) :: dr(1:NDIM)             ! vector displacements (p'-p)
  real(kind=PR) :: drmag                  ! distance
  real(kind=PR) :: drsqd                  ! particle separation squared
  real(kind=PR) :: hrangesqd_high         ! particle radius (2*h_high) squared
  real(kind=PR) :: hrangesqd_low          ! particle radius (2*h_low) squared
  real(kind=PR) :: hfactor                ! invh^NDIM
  real(kind=PR) :: h_high                 ! high value estimate of h
  real(kind=PR) :: h_low                  ! low value of h
  real(kind=PR) :: h_old                  ! old value of h
  real(kind=PR) :: h_search               ! value of h for tree search  
  real(kind=PR) :: invhp                  ! 1 / h
  real(kind=PR) :: mpp                    ! Mass of neighbour
  real(kind=PR), allocatable :: rsqd_high(:)  ! drsqd values for high list
  real(kind=PR), allocatable :: rsqd_low(:)   ! drsqd values for low list

  pp_pot   = 0
  pp_max   = LISTSIZE
  h_old    = hep
  h_high   = HMULT*h_old
  h_low    = h_old / HMULT
  h_search = 0.0_PR
  hrangesqd_high = KERNRANGESQD*h_high*h_high
  allocate(pp_potlist(1:pp_max))
  allocate(pp_high_list(1:pp_max))
  allocate(pp_low_list(1:pp_max))
  allocate(rsqd_high(1:pp_max))
  allocate(rsqd_low(1:pp_max))

  
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
           call BHhydrowalk_hgather(rep(1:NDIM),KERNRANGE*h_search,&
                &pp_pot,pp_max,pp_potlist,ctot_hydro,BHhydro)
#else
           do pp=1,ptot
              if (pp==p) cycle
              call distance2(rep(1:NDIM),pp,dr(1:NDIM),drsqd)
              if (drsqd < hrangesqd_high) then
                 pp_pot = pp_pot + 1
                 pp_potlist(pp_pot) = pp
              end if
           end do
#endif
        end if
           
        ! Precaution against ridiculously high neighbour numbers
        if (pp_pot < 0) then
           h_search = 0.0_PR
           pp_pot = 0
           pp_max = ptot
           deallocate(rsqd_low)
           deallocate(rsqd_high)
           deallocate(pp_low_list)
           deallocate(pp_high_list)
           deallocate(pp_potlist)
           allocate(pp_potlist(1:pp_max))
           allocate(pp_high_list(1:pp_max))
           allocate(pp_low_list(1:pp_max))
           allocate(rsqd_high(1:pp_max))
           allocate(rsqd_low(1:pp_max))
        else if (pp_pot <= pp_gather) then
           h_search = 0.0_PR
           pp_pot = 0
           h_low = h_high
           h_high = HMULT*h_high
           hrangesqd_high = KERNRANGESQD*h_high*h_high
        end if
        if (pp_pot > 0) exit
     end do
     ! -----------------------------------------------------------------------

     pp_tot  = 0
     pp_high = 0
     pp_low  = 0
     hrangesqd_high = KERNRANGESQD*h_high*h_high
     hrangesqd_low  = KERNRANGESQD*h_low*h_low

     ! Loop over potential neighbour list and find the number of particles 
     ! with r < 2*h_low and the number and list of particles with 
     ! 2*h_low < r < 2*h_high
     ! -----------------------------------------------------------------------
     do i=1,pp_pot
        pp = pp_potlist(i)
        call distance2(rep(1:NDIM),pp,dr(1:NDIM),drsqd)
        
        ! Is potential neighbour within 2*h_high of particle p
        if (drsqd < hrangesqd_high) then   
           
           ! Now find particles within 2*h_low and 2*h_low < r < 2*h_high
           if (drsqd < hrangesqd_low) then  
              pp_low = pp_low + 1
              pp_low_list(pp_low) = pp
              rsqd_low(pp_low) = drsqd
           else                         
              pp_high = pp_high + 1
              pp_high_list(pp_high) = pp
              rsqd_high(pp_high) = drsqd        
           end if
        end if

        if (pp_high == LISTSIZE) exit
     end do
     ! -----------------------------------------------------------------------
     
     pp_tot = pp_low + pp_high   
     
     ! If potential list is too small or too large, adjust h limits
     ! -----------------------------------------------------------------------
     if (pp_tot < pp_gather) then
        h_low = h_high
        h_high = HMULT*h_high
        hrangesqd_high = KERNRANGESQD*h_high*h_high
        cycle
        
     ! else if h_low is too large, decrease limits
     else if (pp_low >= pp_gather) then
        h_high = h_low
        h_low = h_low / HMULT
        cycle
        
     ! If the limits contain enough neighbours, exit loop
     else if (pp_low < pp_gather .and. pp_tot >= pp_gather) then 
        exit
     end if
     ! -----------------------------------------------------------------------
     
  end do
! ============================================================================


! Sort h_high list to determine h
#if defined(HEAPSORT)
  call heapsort_real(pp_high,pp_high_list,rsqd_high)
#else
  call insertion_sort_real(pp_high,pp_high_list(1:pp_high),rsqd_high(1:pp_high))
#endif
  drsqd  = rsqd_high(pp_gather - pp_low)
  hep = INVKERNRANGE*sqrt(drsqd)


! Combine near and far sorted lists
  if (pp_gather - pp_low >= 1) then
     do i=1,pp_gather - pp_low
        pp_low_list(i+pp_low) = pp_high_list(i)
        rsqd_low(i+pp_low) = rsqd_high(i)
     end do
  end if


! Now calculate the density using the gather list
! ============================================================================
  rhoep = 0.0_PR
  invhp = 1.0_PR / hep
  hfactor = invhp**(NDIMPR)

! Now loop over all neighbours to find contributions
! ----------------------------------------------------------------------------
  do i=1,pp_gather
     pp      = pp_low_list(i)
     call distance2(rep(1:NDIM),pp,dr(1:NDIM),drsqd)
     mpp     = sph(pp)%m
     drmag   = sqrt(drsqd)
     rhoep   = rhoep + mpp*hfactor*w0(drmag*invhp)
  end do
! ----------------------------------------------------------------------------


#ifdef DEBUG_H_GATHER
  write(6,*) "Number of neighbours for particle", p, "=", pp_gather
  write(6,*) "Smoothing length =", hep,h_low,h_high,h_old
  if ( (abs(hep - h_old)/h_old) > 0.1) &
  	& write(6,*) "Smoothing length changed by more than 10%"
#endif
  
  deallocate(rsqd_low)
  deallocate(rsqd_high)
  deallocate(pp_low_list)
  deallocate(pp_high_list)
  deallocate(pp_potlist)
  
  return
END SUBROUTINE HP_rhoh_ep
