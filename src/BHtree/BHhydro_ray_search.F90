! BHHYDRO_RAY_SEARCH.F90
! D. A. Hubber - 24/01/2008
! Walks hydro tree to determine the id of all particles that are intersected 
! by a ray of source position rsource and direction runit.  Also returns 
! the impact parameters of all particles intersected.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE BHhydro_ray_search(rsource,rdir,pp_tot,pp_max,pp_list,blist)
  use interface_module, only : distance2,distance3
  use definitions
  use particle_module, only : ptot,sph
  use tree_module, only : BHhydro,ctot_hydro
  implicit none

  integer, intent(in)  :: pp_max               ! Max. no. of pot. neighbours
  integer, intent(out) :: pp_tot               ! No. of particles in ray
  integer, intent(out) :: pp_list(1:pp_max)    ! Particle ray list
  real(kind=PR),intent(in) :: rsource(1:NDIM)  ! Position of ray source
  real(kind=PR),intent(in) :: rdir(1:NDIM)     ! Direction of ray
  real(kind=PR),intent(out) :: blist(1:pp_max) ! Impact parameters

  integer :: i                                 ! Auxilary leaf counter
  integer :: c                                 ! Cell counter
  integer :: pp                                ! Particle id
  real(kind=PR) :: bsqd                        ! Impact parameter squared
  real(kind=PR) :: dr(1:NDIM)                  ! Relative displacement vector
  real(kind=PR) :: drsqd                       ! Distance squared
  real(kind=PR) :: maxdistsqd                  ! Max. distance squared
  real(kind=PR) :: sdotr                       ! Distance (sqd) along ray

! Initialize variables.  Start on root cell (c=0) and work our way down.
  c = 0
  pp_tot = 0

! Loop through all cells in the tree
! ============================================================================
  do
     call distance3(rsource(1:NDIM),BHhydro(c)%r(1:NDIM),dr(1:NDIM),drsqd)
     sdotr = dot_product(dr(1:NDIM),rdir(1:NDIM))
     bsqd = drsqd - sdotr*sdotr
     maxdistsqd = (BHhydro(c)%rmax + BHhydro(c)%hrangemax)**2

     ! Check if cell contains particles that are all behind the ray
     !if (sdotr < 0.) then
     !   if (drsqd >= maxdistsqd) then
     !      c = BHhydro(c)%nextcell
     !      cycle
     !   end if
     !end if


     ! Check if search sphere overlaps cell 'neighbour sphere'
     ! -----------------------------------------------------------------------
     if (bsqd < maxdistsqd) then 

        ! Check if c is a leaf cell or not
        if (BHhydro(c)%leaf > 0) then

           ! If templist has become too long, flag and return
           if (pp_tot + BHhydro(c)%leaf > pp_max) then
              pp_tot = -1
              return
           end if
              
           ! Now loop through all particles in leaf cell and check if 
           ! the ray intercepts their smoothing kernels
           do i=1,BHhydro(c)%leaf
              pp = BHhydro(c)%plist(i)
              call distance2(rsource(1:NDIM),pp,dr(1:NDIM),drsqd)
              sdotr = dot_product(dr(1:NDIM),rdir(1:NDIM))
              bsqd = drsqd - sdotr*sdotr
              if (bsqd < KERNRANGESQD*sph(pp)%h**2 .and. sdotr > 0) then
                 pp_tot = pp_tot + 1
                 pp_list(pp_tot) = pp
                 blist(pp_tot) = sqrt(bsqd)
              end if
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
     ! -----------------------------------------------------------------------

     if (c > ctot_hydro) exit
  end do
! ============================================================================

#ifdef DEBUG_BHTREEWALK
  write(6,*) hp,pp_pot,pp_list(1:pp_pot)
#endif

  return
END SUBROUTINE BHhydro_ray_search
