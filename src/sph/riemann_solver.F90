! RIEMANN_SOLVER.F90
! D. A. Hubber - 09/11/2010
! Solves the Riemann problem between two particles, p and pp, assuming the 
! particles as the LHS and RHS states of the shocktube.
! Isothermal Riemann solver (Balsara 19??; Cha & Whitworth 2003)
! Iterative Riemann solver  (Van Leer 1979; Cha & Whitworth 2003)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE riemann_solver(gamma_eff,vp,vpp,dr_unit,&
     &rho_p,rho_pp,press_p,press_pp,Pstar,vstar)
  use definitions
  use hydro_module, only :gamma
  implicit none

  real(kind=PR), intent(in) :: dr_unit(1:NDIM) ! Unit vector joining particles
  real(kind=PR), intent(in) :: gamma_eff       ! Effective value of gamma
  real(kind=PR), intent(in) :: press_p         ! SPH pressure of particle p
  real(kind=PR), intent(in) :: press_pp        !  "      "     "     "    pp
  real(kind=PR), intent(in) :: rho_p           ! SPH density of particle p
  real(kind=PR), intent(in) :: rho_pp          !  "     "    "     "     pp
  real(kind=PR), intent(in) :: vp(1:NDIM)      ! Velocity of particle p
  real(kind=PR), intent(in) :: vpp(1:NDIM)     !     "    "     "     pp
  real(kind=PR), intent(out) :: Pstar          ! Pressure solution
  real(kind=PR), intent(out) :: vstar          ! Velocity solution

  logical :: doit                      ! ..
  integer :: niterations               ! No. of iterations of Riemann solver
  real(kind=PR) :: Ca                  ! Sound speed of state a
  real(kind=PR) :: Cb                  ! Sound speed of state b
  real(kind=PR) :: Ci                  ! Intermediate state sound speed
  real(kind=PR) :: gamma_aux1          ! Aux. gamma variable
  real(kind=PR) :: gamma_aux2          ! Aux. gamma variable
  real(kind=PR) :: Pold                ! Old pressure
  real(kind=PR) :: va                  ! LHS velocity
  real(kind=PR) :: vb                  ! RHS velocity
  real(kind=PR) :: Wa                  ! ..
  real(kind=PR) :: Wb                  ! ..
  real(kind=PR), parameter :: pfloor = 1.0E-19_PR   ! Pressure floor

! Initialise quantities before iterating
  gamma_aux1  = (gamma_eff - 1.0_PR)/(2.0_PR*gamma_eff)
  gamma_aux2  = (gamma_eff + 1.0_PR)/(2.0_PR*gamma_eff)
  va          = -dot_product(vp(1:NDIM),dr_unit(1:NDIM))
  vb          = -dot_product(vpp(1:NDIM),dr_unit(1:NDIM))
  Ca          = sqrt(gamma_eff*press_p*rho_p)
  Cb          = sqrt(gamma_eff*press_pp*rho_pp)
  Pstar       = (Ca*press_pp + Cb*press_p - Ca*Cb*(va - vb))/(Ca + Cb)
  vstar       = (Ca*va + Cb*vb - (press_p - press_pp))/(Ca + Cb)
  
! If initial guess using linear Riemann solver is good enough, do not 
! use the iterative solver
  if (abs((Pstar - press_p)/press_p) < 0.1_PR .and. &
       abs((Pstar - press_pp)/press_pp) < 0.1_PR) goto 100

 ! write(6,*) "Linear Riemann solver? : ",&
 !      &(Pstar - press_p)/Pstar,(Pstar - press_pp)/Pstar,Pstar,press_p,press_pp

! Non-iterative solution for two rarefactions
  if (Pstar < press_p .and. Pstar < press_pp) then
     Ci = Ca / (rho_p*press_p**(gamma_aux1))
     Pstar = (0.5_PR*(gamma_eff - 1.0_PR)*(vb - va) + Ca/rho_p + Cb/rho_pp)&
          &/(Ci + Cb/(rho_pp*press_pp**(gamma_aux1)))
     if (Pstar < pfloor) then
        Pstar = pfloor
     else
        Pstar = Pstar**(1.0_PR/gamma_aux1)
     end if
     !write(6,*) "Rarefaction : ",Pstar,press_p,press_pp
     goto 100
  end if

  Pold = Pstar;
  Wa   = Ca*wave(Pstar,press_p,gamma_aux1,gamma_aux2)
  Wb   = -Cb*wave(Pstar,press_pp,gamma_aux1,gamma_aux2)
  doit = .true.
  niterations = 0
     
  
! Iterate towards consistent solution
! ----------------------------------------------------------------------------
  do while(doit)
     Pstar = (Wa*press_pp - Wb*press_p + Wa*Wb*(va - vb))/(Wa - Wb)

     if (abs(Pstar - Pold) > 0.01_PR*Pstar .and. Pstar > 0.0_PR) then
        Pold = Pstar
        Wa = Ca*wave(Pstar,press_p,gamma_aux1,gamma_aux2)
        Wb = -Cb*wave(Pstar,press_pp,gamma_aux1,gamma_aux2)
        niterations = niterations + 1

        if (niterations > 50) then
           write(6,*) "Convergence failure in Riemann solver"

           ! Use linear guess
           Pstar = (Ca*press_pp + Cb*press_p - Ca*Cb*(va - vb))/(Ca + Cb)
           Pstar = max(Pstar,pfloor)
           Wb = -Cb
           Wa = Ca
           doit = .false.
        else if (niterations > 20) then
           write(6,*) "Slow convergence in Riemann solver; nit : ",niterations
        end if
     else
        doit = .false.
     end if
  end do

! Calculate resolved velocity
  vstar = (Wa*va - Wb*vb - (press_p - press_pp))/(Ca - Cb)

    if (Pstar < 0.0_PR) then
       write(6,*) "Riemann problem solved!! : ",press_p,press_pp,Pstar,Ca,Cb,Wa,Wb,vstar,va - vb
       write(6,*) "Riemann solver failed!!"
     stop
  end if

! ----------------------------------------------------------------------------

100 return

contains

FUNCTION wave(Pstar,Pi,gamma_aux1,gamma_aux2)
  use definitions
  implicit none

  real(kind=PR), intent(in) :: Pstar       ! Pressure from Riemann solver
  real(kind=PR), intent(in) :: Pi          ! Pressure of particle i
  real(kind=PR), intent(in) :: gamma_aux1  ! ..
  real(kind=PR), intent(in) :: gamma_aux2  ! ..
  
  real(kind=PR) :: wave                    ! ..
  real(kind=PR) :: xaux                    ! ..
  
  xaux = Pstar/Pi

! Use linear expression for wave
  if (abs(xaux - 1.0_PR) < 1.0E-3_PR) then
     wave = 1.0_PR + 0.5_PR*gamma_aux2*(xaux - 1.0_PR)

! Else use non-linear expression
  else
     ! Shock
     if (xaux >= 1.0_PR) then
        wave = sqrt(1.0_PR + gamma_aux2*(xaux - 1.0_PR))
     else
        wave = gamma_aux2*(1.0_PR - xaux)/(1.0_PR - xaux**(gamma_aux2))
     end if
  end if

  return
END FUNCTION wave

END SUBROUTINE riemann_solver



! ============================================================================
! ISOTHERMAL_RIEMANN_SOLVER
! Simpler isothermal Riemann solver (Balsara 19??)
! ============================================================================
SUBROUTINE isothermal_riemann_solver(vp,vpp,dr_unit,&
     &rho_p,rho_pp,sound_p,Pstar,vstar)
  use definitions
  implicit none

  real(kind=PR), intent(in) :: dr_unit(1:NDIM) ! Unit displacement vector
  real(kind=PR), intent(in) :: rho_p           ! Density of p
  real(kind=PR), intent(in) :: rho_pp          ! Density of pp
  real(kind=PR), intent(in) :: sound_p         ! Sound speed of p (and pp)
  real(kind=PR), intent(in) :: vp(1:NDIM)      ! Velocity of p
  real(kind=PR), intent(in) :: vpp(1:NDIM)     ! Velocity of pp
  real(kind=PR), intent(out) :: Pstar          ! Pressure solution
  real(kind=PR), intent(out) :: vstar          ! Velocity solution

  real(kind=PR) :: va                          ! LHS velocity
  real(kind=PR) :: vb                          ! RHS velocity
  real(kind=PR) :: X                           ! Aux. var. for Riemann solver

  va = dot_product(vp(1:NDIM),dr_unit(1:NDIM))
  vb = dot_product(vpp(1:NDIM),dr_unit(1:NDIM))
  X = sqrt(rho_p)*sqrt(rho_pp) / (sqrt(rho_p) + sqrt(rho_pp))
  Pstar = 0.25_PR*(X*(va - vb) + sqrt(X*X*(va - vb)**2 + &
       &4*sound_p*sound_p*X*(sqrt(rho_p) + sqrt(rho_pp))))**2
  !vstar = vl - (Pstar - sound_p*sound_p*rho_p)/(sqrt(Pstar*rho_p))
  !vstar = vr + (Pstar - sound_p*sound_p*rho_pp)/(sqrt(Pstar*rho_pp))

  return
END SUBROUTINE isothermal_riemann_solver

