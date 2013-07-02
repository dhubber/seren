! H_RHO_ITERATION.F90
! C. P. Batty & D. A. Hubber - 11/12/2006
! Calculates density and smoothing length of particle p using the 
! grad-h method of Price & Monaghan (2004).  Also calculates all other 
! gather-only SPH properties (e.g. velocity divergence, velocity curl, etc).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE h_rho_iteration(p,hp,typemask)
  use interface_module, only : copy_neighbour_list,distance2,&
       &gather_neib_on_fly,w0,womega
  use particle_module
  use hydro_module
  use neighbour_module
  use kernel_module
  use time_module
  use sink_module
  use type_module
  implicit none

  integer, intent(in) :: p                            ! particle id
  real(kind=PR), intent(inout) :: hp                  ! h of particle p
  logical, optional, intent(in) :: typemask(1:ntypes) ! part. types to include?

  integer :: i                          ! counter in neighbour search
  integer :: iteration                  ! Number of iterations
  integer :: pp                         ! neighbouring particles (p')
  integer :: pp_max                     ! length of pot list array
  integer :: pp_pot                     ! number of potential neighbours
  integer, allocatable :: pp_potlist(:) ! list of potential neighbours
  real(kind=PR) :: dr(1:NDIM)           ! vector displacements (p'-p)
  real(kind=PR) :: drmag                ! Distance
  real(kind=PR) :: drsqd                ! p'-p separation squared
  real(kind=PR) :: dt                   ! Time since last force update 
  real(kind=PR) :: dt_new               ! Latest timestep
  real(kind=PR) :: dt_old               ! Previous timestep
  real(kind=PR) :: fbisection           ! Aux variable for bisection method
  real(kind=PR) :: frho                 ! Newton-Rhapson function
  real(kind=PR) :: frhoprime            ! Newton-Rhapson function gradient
  real(kind=PR) :: hrangesqd            ! particle radius (2*h_new) squared
  real(kind=PR) :: hfactor              ! invhp**(NDIM)
  real(kind=PR) :: h_high               ! high value estimate of h
  real(kind=PR) :: h_lower              ! Lower bound on h for bisection method
  real(kind=PR) :: h_new                ! New value of h
  real(kind=PR) :: h_old                ! old value of smoothing length
  real(kind=PR) :: h_upper              ! Upper bound on h for bisection method
  real(kind=PR) :: invhp                ! (1 / hp)
  real(kind=PR) :: mp                   ! Mass of particle p
  real(kind=PR) :: mpp                  ! Mass of particle pp
  real(kind=PR) :: rhotemp              ! Auxilary summation variable for rho
  real(kind=PR) :: rp(1:NDIM)           ! position of particle p
!#if defined(GRAD_H_SPH)
  real(kind=PR) :: omega_p              ! Omega correction factor for p
!#endif
#if defined(SINKS)
  integer :: s                          ! Sink counter
#endif

  debug3("Calculating h [h_rho_iteration.F90] for particle ", p)

! Make local copy of old smoothing length and particle position
  rp(1:NDIM) = sph(p)%r(1:NDIM)
  mp         = sph(p)%m
  hp         = sph(p)%h
  rhotemp    = sph(p)%rho
#if defined(GRAD_H_SPH)
  omega_p    = sph(p)%omega
#else
  omega_p    = 1.0_PR
#endif

! Gather neibs over a larger volume
  iteration = 0
  h_old     = hp
  h_new     = hp
  h_high    = 0.0_PR
  h_upper   = rextent
  h_lower   = 0.0_PR
  pp_pot    = 0

! Estimate the new smoothing length h
! Determine time interval since last acceleration calculation
#if defined(EULER) || defined(LEAPFROG_KDK)
  dt = real(sph(p)%laststep,PR)
#elif defined(RUNGE_KUTTA2)
  dt = 0.5_DP*real(n - sph(p)%nlast,PR)*real(timestep,PR)
#elif defined(LEAPFROG_DKD)
  dt = 0.5_DP*(real(sph(p)%laststep,PR) + &
       &real(n - sph(p)%nlast,PR)*real(timestep,PR))
#endif
  dt = 0.0_PR
#if defined(GRAD_H_SPH)
  hp = hp * (1.0_PR - ((sph(p)%div_v / (NDIMPR * omega_p))*dt))
#else
  hp = hp * (1.0_PR - ((sph(p)%div_v / (NDIMPR))*dt))
#endif
  invhp = 1.0_PR / hp
  hfactor = invhp**(NDIM)
  hrangesqd = KERNRANGESQD*hp*hp


! Main iteration loop
! ============================================================================
  do

     ! Obtain new potential neighbour list either by direct sum or the tree
     ! -----------------------------------------------------------------------
     do 

#if defined(DEBUG_H_RHO_ITERATION)
        write(6,*) "Iterating h/rho for :",p,iteration,pp_pot,pp_gather
        write(6,*) "h : ",hp,h_high,rextent,h_upper,h_lower
#endif

        ! Compute new potential neighbour list if new search radius is
        ! greater than that used in last tree walk
        if (hp > h_high .or. pp_pot <= 0) then
           pp_pot = 0        
           h_high = hp*HMULT
           if (allocated(pp_potlist)) deallocate(pp_potlist)

           ! Obtain list only from certain particle types if mask defined
           if (present(typemask)) then
              call gather_neib_on_fly(p,pp_max,pp_pot,pp_potlist,&
                   &rp(1:NDIM),KERNRANGE*h_high,typemask)
           else
              call gather_neib_on_fly(p,pp_max,pp_pot,pp_potlist,&
                   &rp(1:NDIM),KERNRANGE*h_high)
           end if
        end if
        
        ! Precaution against ridiculously high neighbour numbers
        if (pp_pot < 0) then
           h_high = 0.0_PR
           pp_pot = 0
           pp_max = ptot
        else if (pp_pot < pp_gather) then
           h_high = 0.0_PR
           pp_pot = 0
           hp     = hp*HMULT
        end if

        ! If we have a reasonable potential neighbour list, exit loop
        if (pp_pot > 0) exit

     end do
     ! -----------------------------------------------------------------------

     h_new   = hp
     invhp   = 1.0 / hp
     hfactor = invhp**(NDIM)
     hrangesqd = KERNRANGESQD*hp*hp

     ! First, include self contributions
#if defined(H_NUMBER)
     rhotemp = hfactor*w0(0.0_PR)
#else
     rhotemp = mp*hfactor*w0(0.0_PR)
#endif
#if defined(GRAD_H_SPH)
     omega_p = mp*hfactor*invhp*womega(0.0_PR)
#endif

     ! Now loop over all potential neighbours and calculate SPH quantities
     ! -----------------------------------------------------------------------
     do i=1,pp_pot
        pp = pp_potlist(i)
        if (pp == p) cycle
#if defined(PERIODIC) && !defined(GHOST_PARTICLES)
        call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
#else
        dr(1:NDIM) = sph(pp)%r(1:NDIM) - rp(1:NDIM)
        drsqd = dot_product(dr(1:NDIM),dr(1:NDIM))
#endif
        if (drsqd >= hrangesqd) cycle
        drmag = sqrt(drsqd) !+ SMALL_NUMBER
#if defined(H_NUMBER)
        rhotemp = rhotemp + hfactor*w0(drmag*invhp)
#else
        mpp = sph(pp)%m
        rhotemp = rhotemp + mpp*hfactor*w0(drmag*invhp)
#endif
#if defined(GRAD_H_SPH)
        omega_p = omega_p + mpp*hfactor*invhp*womega(drmag*invhp)
#endif
     end do
     ! -----------------------------------------------------------------------

#if defined(H_NUMBER)
     rhotemp = mp*rhotemp
#endif
#if defined(GRAD_H_SPH)
     omega_p = 1.0_PR + ((hp * omega_p) / (NDIMPR * rhotemp))
#endif

     ! Safety check to make sure fixed-point iteration has not gone crazy
     ! If so, go straight to bisection method
     if (.not.(h_new > 0.0_PR .and. h_new < rextent)) then
        iteration = GRADH_ITERATION_MAX
        
     ! If iteration has not gone crazy, first convergence check here
     else if ((abs(hp - h_fac*(mp/rhotemp)**(INVNDIM))/hp < H_CON) &
          & .and. iteration >= 1) then
        exit
     end if


     ! Calculate new value of h_new depending on no. of iterations performed
     ! -----------------------------------------------------------------------
     if (iteration < GRADH_ITERATION_MAX) then

        ! If the change is too big (or invalid), use a smaller step
        frho = mp*(h_fac/hp)**(NDIM) - rhotemp
#if defined(GRAD_H_SPH)
        frhoprime = -NDIMPR*rhotemp*omega_p/hp
#else
        frhoprime = -NDIMPR*rhotemp/hp
#endif
        if (abs(frhoprime) > SMALL_NUMBER) then
           h_new = hp - frho/frhoprime
        else
           h_new = h_fac*(mp/rhotemp)**(INVNDIM)
        end if
        
        ! If values are unacceptable, use a simple fixed-point iteration
        if (.not. (h_new >= 0.0_PR .and. h_new <= rextent)) &
             &h_new = h_fac*(mp/rhotemp)**(INVNDIM)
     
     ! If too many iterations have been attempted, prepare for bisection
     else if (iteration == GRADH_ITERATION_MAX) then
        h_upper = rextent
        h_lower = 0.0_PR
        h_new   = 0.5_PR*(h_lower + h_upper)
 
     ! Bisection method (after GRADH_ITERATION_MAX iterations)
     else if (iteration < 2*GRADH_ITERATION_MAX) then
        if (rhotemp > SMALL_NUMBER) then
           fbisection = hp**(NDIM) - (h_fac**(NDIM))*mp/rhotemp
        else
           fbisection = 1.0_PR
        end if
        if (fbisection > 0.0_PR) then
           h_upper = hp
        else
           h_lower = hp
        end if
        h_new = 0.5_PR*(h_lower + h_upper)

     ! If bisection is not converging, simply use the old value
     ! (Not correct for grad-h, but allows continuation of simulation)
     else
        h_new = h_old
        hp = h_new
        exit

     end if
     ! -----------------------------------------------------------------------

     hp = h_new
     iteration = iteration + 1

#if defined(MINIMUM_H)
     if (h_new < hp .and. h_new <= hmin) then
        h_new = hmin
        hp = hmin
        exit
     end if
#endif
     
  end do
! ============================================================================


! If h is less than allowed minimum, reaffirm it again here for safety
#if defined(MINIMUM_H)
  if (hp < hmin) hp = hmin
#endif

! Record potential neighbour list in order to speed up SPH calculations
#if defined(NEIGHBOUR_LISTS)
  pptot(p) = pp_pot
  if (pp_pot <= pp_limit) pplist(1:pp_pot,p) = pp_potlist(1:pp_pot)
#endif

  if (allocated(pp_potlist)) deallocate(pp_potlist)

#if defined(MINIMUM_H)
  if (.not. (sph(p)%h >= hmin)) then
     write(6,*) "Invalid smoothing lengths : ",p,sph(p)%h,hmin,rhotemp,&
          &h_fac*(mp/rhotemp)**(INVNDIM)/hp,sph(p)%r(1:NDIM)
     stop
  end if
#endif

  return
END SUBROUTINE h_rho_iteration

