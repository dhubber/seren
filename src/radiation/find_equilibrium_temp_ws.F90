! FIND_EQUILIBRIUM_TEMP_WS.F90
! D. Stamatellos & D. A. Hubber - 17/8/2009
! Contains subroutines and functions for use in new energy equation 
! and radiative transfer method (Stamatellos et al. 2007).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE find_equilibrium_temp_ws(p)
  use interface_module, only : ambient_temp,ebalance,eosenergy,eosmu,&
       &find_idens,find_itemp,find_temp_from_energy,getkappa
  use particle_module
  use scaling_module
  use constant_module
  use hydro_module
  use Eos_module
  use Tprof_module
  use type_module
  implicit none

  integer, intent(in) :: p            ! Particle id
  logical       :: ZERO               ! Logical flag
  integer       :: itempEq            ! ...
  integer       :: itempnext=0        ! ...
  integer       :: itempindex=0       ! ..
  real(kind=PR) :: balanceHigh        ! Rate of change of internal energy 
  real(kind=PR) :: balanceLow         !    for bisection method
  real(kind=PR) :: balance            ! dudt for equilibirum temp(?)
  real(kind=PR) :: minbalance         ! dudt for equilibirum temp(?)
  real(kind=PR) :: column2_p          ! Local copy of column2 for p
  real(kind=PR) :: dtemp              ! temperature difference for bisection
  real(kind=PR) :: dt_thermal         ! thermalization timescale
  real(kind=PR) :: dudt_Eq            ! equilibrium heating rate
  real(kind=PR) :: dudt_hydro         ! compressive heating rate
  real(kind=PR) :: dudt_p             ! ..
  real(kind=PR) :: dudt_hc_p          ! ..
  real(kind=PR) :: dudt_rad_p         ! radiative heating rate
  real(kind=PR) :: dudt_tot           ! net heating rate
  real(kind=PR) :: kappaEq            ! Rosseland pseudo-mean opacities
  real(kind=PR) :: kappaT             !    for current and equilibrium temp.
  real(kind=PR) :: kappaHigh          ! Rosseland pseudo-mean opacities 
  real(kind=PR) :: kappaLow           !    for bisection method
  real(kind=PR) :: kapparEq           ! Rosseland  opacities
  real(kind=PR) :: kapparT            !    for current and equilibrium temp.
  real(kind=PR) :: kapparHigh         ! Rosseland  opacities 
  real(kind=PR) :: kapparLow          !    for bisection method
  real(kind=PR) :: kappapEq           ! Planck mean opacities for current
  real(kind=PR) :: kappapT            !    and equilibrium temp
  real(kind=PR) :: kappapHigh         ! high Planck opacity on grid
  real(kind=PR) :: kappapLow          ! low Planck opacity on grid
  real(kind=PR) :: rho_p              ! Density of particle p
  real(kind=PR) :: Teq                ! equilibrium temperature
  real(kind=PR) :: To                 ! ambient temperature
  real(kind=PR) :: Thigh              ! High temp on grid
  real(kind=PR) :: Tlow               ! Low temp on grid
  real(kind=PR) :: Thigh_temp         ! High temp on grid
  real(kind=PR) :: Tlow_temp          ! Low temp on grid
  real(kind=PR) :: kappaCool          ! ..
  real(kind=PR) :: kappapCool         ! ..
  real(kind=PR) :: kapparCool         ! ..
  real(kind=PR) :: kappaHeat          ! ..
  real(kind=PR) :: kappapHeat         ! ..
  real(kind=PR) :: kapparHeat         ! ..
  real(kind=PR) :: balanceCool        ! ..
  real(kind=PR) :: balanceHeat        ! ..
  real(kind=PR) :: minkappaT          ! ..
  real(kind=PR) :: minkappapT         ! ..
  real(kind=PR) :: minkapparT         ! .. 
  real(kind=PR) :: ueq_p              ! equilibrium specific energy
#if defined(DIFFUSION)
  real(kind=PR) :: dudt_diff          ! diffusion rate in dt
#endif
#if defined(STAR_HEATING)
  real(kind=PR) :: dr(1:NDIM)         ! relative displacement vector
  real(kind=PR) :: dudt_star          ! diffusion rate in dt
  real(kind=PR) :: dudt_star_lost     ! ..
#endif
#if defined(SPIEGEL_TEST)
  real(kind=PR), parameter :: EPSIL = 0.000001_PR
#else
  real(kind=PR), parameter :: EPSIL = 0.001_PR
#endif 
#if defined(DEBUG_RAD)
  real(kind=PR) :: drsqd               ! radius squared
  real(kind=PR) :: dr(1:NDIM)          ! relative displacement vector
#endif

  debug3("Calculating equation of state [eos.F90] for particle ",p)


#if defined(DEBUG_RAD)
  dr(1:NDIM) = sph(p)%r(1:NDIM)
  drsqd= (dr(1)*dr(1) + dr(2)*dr(2))*(6.684e-14*rscale*rcgs)**2
#endif
#if defined(DIFFUSION)
  dudt_diff = 0.0_PR
#endif
#if defined(STAR_HEATING)
  dudt_star = 0.0_PR
  dudt_star_lost = 0.0_PR
#endif
  rho_p = sph(p)%rho
  ZERO = .false.

! Calculate ambient temperature
  call ambient_temp(p,To)

  call find_idens(rho_p,sph(p)%idens)
  call find_itemp(sph(p)%temp,sph(p)%itemp)
  if (sph(p)%idens <= 1) sph(p)%idens = 2
  if (sph(p)%idens >= dim_dens) sph(p)%idens = dim_dens - 1
  if (sph(p)%itemp <= 1) sph(p)%itemp = 2
  if (sph(p)%itemp >= dim_temp) sph(p)%itemp = dim_temp - 1

#if defined(RAD_WS_SINK_POT)
  sph(p)%column2 = (fcolumn**2)*sph(p)%gpot*rho_p
#else
  sph(p)%column2 = (fcolumn**2)*sph(p)%sphgpot*rho_p
#endif

  column2_p = sph(p)%column2

! Rosseland and Planck opacities at rho_p and temperature temp(p)
  call getkappa(rho_p,sph(p)%temp,sph(p)%idens,kappaT,kapparT,kappapT)

! calculate radiative heating/cooling rate
  call ebalance(dudt_rad_p,0.0_PR,To,sph(p)%temp,kappapT,kappaT,column2_p)

#if defined(DIFFUSION)
  dudt_diff = sph(p)%du_dt_diff
  dudt_p    = sph(p)%dudt

! make sure there is no double counting in hybrid method
  if (dudt_diff*dudt_rad_p <= 0.0_PR) then
     dudt_hydro =  dudt_p !+ dudt_diff
     dudt_hc_p  =  dudt_rad_p !+dudt_diff
  else
     dudt_hydro = dudt_p + dudt_diff
     dudt_hc_p  = dudt_diff + dudt_rad_p
  endif

! both positive
!  if (dudt_diff>0 .and. dudt_rad_p>0 .and. dudt_rad_p>dudt_diff) then 
!     dudt_hydro = dudt_p
!     dudt_hc_p  = dudt_rad_p
!  endif

!  if (dudt_diff>0 .and. dudt_rad_p>0 .and. dudt_rad_p<dudt_diff) then
!     dudt_hydro = dudt_p+(dudt_diff-dudt_rad_p)
!     dudt_hc_p  = (dudt_diff-dudt_rad_p)
!  endif

!both negative
!  if (dudt_diff<0 .and. dudt_rad_p<0 .and. dudt_rad_p>dudt_diff) then
!     dudt_hydro = dudt_p+(dudt_diff-dudt_rad_p)
!     dudt_hc_p  = (dudt_diff-dudt_rad_p)
!  endif

!  if (dudt_diff<0 .and. dudt_rad_p<0 .and. dudt_rad_p<dudt_diff) then
!      dudt_hydro = dudt_p
!      dudt_hc_p  = dudt_rad_p
!  endif
  
#else
  dudt_hydro = sph(p)%dudt
#endif

#if defined(STAR_HEATING)
  dr(1:NDIM) = sph(p)%r(1:NDIM)
  if (sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)) < 5.0e-3) then
     dudt_star = 4.0_PR*PI*((3*5.0d-3)**2/(6.684d-14*rscale*rcgs)**2)&
          &*rad_const*(6000d0)**4/sph(p)%m/119.
     dudt_hydro = dudt_hydro + dudt_star      
  endif
#endif 

#if defined(DEBUG_RAD)
  rad_info(1,p) = kappapT
  rad_info(2,p) = kappaT
  rad_info(3,p) = kappapT*sqrt(column2_p)
  rad_info(4,p) = kappaT*sqrt(column2_p)
  rad_info(5,p) = sqrt(column2_p)
  rad_info(6,p) = sph(p)%m*sph(p)%gpot
#endif

! Find equilibrium temperature
! ----------------------------------------------------------------------------

! Calculate energy balance assuming current temperature as starting point
  call ebalance(balance,0.0_PR,To,sph(p)%temp,kappapT,kappaT,column2_p)

!  call getkappa(rho_p,eos_temp(itemp(p)-1),idens(p),kappaCool,kapparCool,kappapCool)
!  call ebalance(balanceCool,0.0_PR,To,eos_temp(itemp(p)-1),kappapCool,kappaCool,column2_p)

!  call getkappa(rho_p,eos_temp(itemp(p)+1),idens(p),kappaHeat,kapparHeat,kappapHeat)
!  call ebalance(balanceHeat,0.0_PR,To,eos_temp(itemp(p)+1),kappapHeat,kappaHeat,column2_p)

  call getkappa(rho_p,5.0_PR,sph(p)%idens,minkappaT,minkapparT,minkappapT)
  call ebalance(minbalance,0.0_PR,To,5.0_PR,minkappapT,minkappaT,column2_p)
 
  itempindex = 0

! Search for solutions at temperatures < temp(p)  (i.e. when the compressional
! heating is less than the radiative cooling resulting in net cooling)
! ============================================================================
!  if (dudt_hydro <= -balance .and. balanceCool>=balance) then 
  if (dudt_hydro <= -balance) then

     ! Ensure that temperature is not too low due to expansion 
     ! (5K is set as minimum due to external radiation field such as CMB)
     if (dudt_hydro <= -minbalance)  then
        dudt_hydro = 0.0_PR
        Teq = 5.0_PR
        goto 40
     end if

     ! If temperature drops too low, set to minimum value for look-up tables
     if (sph(p)%itemp <= 2) then     
        Teq = eos_temp(2)
        goto 40
     end if

     ! Low and high temperatures for search
     Tlow  = real(eos_temp(sph(p)%itemp - 1),PR) 
     Thigh = real(eos_temp(sph(p)%itemp + 1),PR)
    
     itempindex = 2 
     itempnext = sph(p)%itemp - itempindex
     
     call getkappa(rho_p,Tlow,sph(p)%idens,kappaLow,kapparLow,kappapLow)
     call eBalance(balanceLow,dudt_hydro,To,Tlow,kappapLow,kappaLow,column2_p)
     call getkappa(rho_p,Thigh,sph(p)%idens,kappaHigh,kapparHigh,kappapHigh)
     call eBalance(balanceHigh,dudt_hydro,To,Thigh,&
          &kappapHigh,kappaHigh,column2_p)     

     ! -----------------------------------------------------------------------
    do while (balanceLow*balanceHigh > 0.0_PR )

        if (itempnext <= 2) then     
           Teq = eos_temp(2)
           goto 40
        end if

        Thigh       = Tlow     
        kappaHigh   = kappaLow
        kappapHigh  = kappapLow
        balanceHigh = balanceLow
        Tlow        = real(eos_temp(itempnext),PR)
        itempindex  = itempindex + 1
        itempnext   = sph(p)%itemp - itempindex 
        
        call getkappa(rho_p,Tlow,sph(p)%idens,kappaLow,kapparLow,kappapLow)
        call eBalance(balanceLow,dudt_hydro,To,Tlow,&
             &kappapLow,kappaLow,column2_p)

     end do
    
     dtemp = Thigh - Tlow

! Else search for solutions where the compressional heating is greater than 
! the radiative cooling (resulting in temperature increase)
! ============================================================================
!  else if (dudt_hydro > -balance .and. balanceHeat<=balance) then
!  else if (dudt_hydro > -balance) then
  else

    ! Search for solutions at temperatures > temp(p)
     if (sph(p)%itemp >= dim_temp - 1 ) then     
        Teq = eos_temp(dim_temp - 1)
        goto 40
     endif

     Tlow  = real(eos_temp(sph(p)%itemp - 1),PR)
     Thigh = real(eos_temp(sph(p)%itemp + 1),PR)
     itempindex = 2
     itempnext  = sph(p)%itemp + itempindex
     
     call getkappa(rho_p,Tlow,sph(p)%idens,kappaLow,kapparLow,kappapLow)
     call eBalance(balanceLow,dudt_hydro,To,Tlow,kappapLow,kappaLow,column2_p)
     call getkappa(rho_p,Thigh,sph(p)%idens,kappaHigh,kapparHigh,kappapHigh)
     call eBalance(balanceHigh,dudt_hydro,To,Thigh,&
          &kappapHigh,kappaHigh,column2_p)     
     
     ! -----------------------------------------------------------------------
    do while (balanceLow*balanceHigh > 0.0_PR)       

       if (itempnext >= dim_temp - 1) then     
          Teq = eos_temp(dim_temp - 1)
          goto 40
       endif
       
       Tlow       = Thigh
       kappaLow   = kappaHigh
       kappapLow  = kappapHigh
       balanceLow = balanceHigh
       !if (itempnext == 0) write(6,*) "Something wrong here : ",&
       !     &p,itempindex,itempnext,dim_temp
       Thigh      = real(eos_temp(itempnext),PR)
       itempindex = itempindex + 1 
       
       itempnext  = sph(p)%itemp + itempindex
       
       call getkappa(rho_p,Thigh,sph(p)%idens,kappaHigh,kapparHigh,kappapHigh)
       call eBalance(balanceHigh,dudt_hydro,To,Thigh,&
            &kappapHigh,kappaHigh,column2_p)
       
     end do 

     dtemp = Thigh - Tlow

  end if

! ============================================================================

  Teq = 0.5_PR*(Tlow + Thigh)
  dtemp = Thigh - Tlow


! ----------------------------------------------------------------------------
! do while ((.not. ZERO) .and. (dtemp>1))  
  do while ((.not. ZERO) .and. abs(2.0_PR*dtemp/(Thigh + Tlow)) > EPSIL) 
     
     call getkappa(rho_p,Teq,sph(p)%idens,kappaEq,kapparEq,kappapEq)
     call eBalance(balance,dudt_hydro,To,Teq,kappapEq,kappaEq,column2_p)       
     
     Thigh_temp = Thigh
     Tlow_temp = Tlow
     
     ! If we've found the equilibirum temperature exactly, then discontinue 
     ! iteration.  Else, use bisection method to refine guess of temperature.
     if (balance == 0.0_PR) then
        ZERO = .true.
     else 
        if (balanceLow*balance < 0.0_PR) then
           Thigh       = Teq
           balanceHigh = balance
           kappaHigh   = kappaEq
           kappapHigh  = kappapEq
        else 
           Tlow       = Teq
           balanceLow = balance
           kappaLow   = kappaEq
           kappapLow  = kappapEq
        end if
     end if
     Teq = 0.5_PR*(Tlow + Thigh)
     dtemp = (Thigh - Tlow)     
  end do
  
  Teq = 0.5_PR*(Tlow + Thigh)  
! ----------------------------------------------------------------------------

! Make sure equilibrium temperature is not too low
40 if (Teq < 5.0_PR) Teq = 5.0_PR

! Now we have calculated the equilibrium temperature, calculate all other 
! properties at this temperature
  call getkappa(rho_p,Teq,sph(p)%idens,kappaEq,kapparEq,kappapEq)
  call find_itemp(Teq,itempEq)  
  ueq_p = eosenergy(rho_p,Teq,sph(p)%idens,itempEq)

  call ebalance(dudt_Eq,0.0_PR,To,Teq,kappapEq,kappaEq,column2_p)
  
! This is Teq^4/column2kappaEq-temp(p)^4/column2 kappaT  
  dudt_tot = -(dudt_Eq - dudt_rad_p) 

! Calculate thermalization timescale
  if (dudt_tot == 0.0_PR) then 
     dt_thermal = 1.E20_PR
  else
#if defined(DIFFUSION) && defined(STAR_HEATING)
     dt_thermal = (ueq_p - sph(p)%u) / &
          & (sph(p)%dudt + dudt_rad_p + dudt_diff + dudt_star)
#elif defined(DIFFUSION)
     !dt_thermal = (ueq_p - u(p))/(du_dt(p) + dudt_rad_p + dudt_diff)
     dt_thermal = (ueq_p - sph(p)%u)/(sph(p)%dudt + dudt_hc_p)
#else
     dt_thermal = (ueq_p - sph(p)%u)/(sph(p)%dudt + dudt_rad_p)
#endif
  end if
    
! Record equilibrium internal energy and thermalization timescale
  sph(p)%dt_therm = dt_thermal
  sph(p)%ueq = ueq_p
#if defined(DEBUG_DUDTRAD)
  sph(p)%dudt_rad = dudt_rad_p
#endif
#if defined(DEBUG_RAD)
  rad_info(7,p) = sph(p)%dudt
  rad_info(8,p) = dudt_rad_p
  rad_info(9,p) = dudt_tot
  rad_info(10,p) = dudt_Eq
#endif

  return
END SUBROUTINE find_equilibrium_temp_ws



! ============================================================================
! EBALANCE
! D. Stamatellos - 3/1/2008
! Calculates net heating rate due to hydro (i.e. expansion/contraction) 
! and radiative effects (i.e heating/cooling)
! ============================================================================
SUBROUTINE ebalance(energybalance,du_dt,To,T,kappapT,kappaT,column2)
  use definitions
  use Eos_module, only : rad_const
  implicit none

  real(kind=PR), intent(in)  :: du_dt          ! Compressional heating rate
  real(kind=PR), intent(in)  :: To             ! 'Background' temperature
  real(kind=PR), intent(in)  :: T              ! Current gas temperature
  real(kind=PR), intent(in)  :: kappapT        ! Planck mean opacity
  real(kind=PR), intent(in)  :: kappaT         ! Rosseland mean opacity
  real(kind=PR), intent(in)  :: column2        ! Column density squared
  real(kind=PR), intent(out) :: energybalance  ! Net heating/cooling rate
  
  energybalance = du_dt - 4.0_PR*rad_const*(T**4 - To**4)/ &
       &((column2*kappaT) + (1.0_PR/kappapT)) 
  
  return
END SUBROUTINE ebalance



! ============================================================================
! GETKAPPA
! D. Stamatellos - 3/1/2008
! Obtains value of Rosseland mean opacity from inputted table.  
! Interpolates between grid points for more accurate value.  
! ============================================================================
SUBROUTINE getkappa(dens,temp,idens,kappa_mean,kappa_rosseland,kappa_planck)
  use interface_module, only : find_itemp
  use definitions
  use Eos_module, only : bdens,btemp,densmin,kappa,kappap,kappar,tempmin
  implicit none

  integer, intent(in) :: idens                   ! ..
  real(kind=PR), intent(in) :: temp              ! ..
  real(kind=PR), intent(in) :: dens              ! ..
  real(kind=PR), intent(out) :: kappa_mean       ! ..
  real(kind=PR), intent(out) :: kappa_planck     ! ..
  real(kind=PR), intent(out) :: kappa_rosseland  ! ..

  integer :: itemp                               ! ..
  real(kind=PR) :: deltadens                     ! ..
  real(kind=PR) :: deltatemp                     ! ..
  
! Find table location
  call find_itemp(temp,itemp)
  
! Find linear interpolation weightings inbetween table points
  deltadens = bdens*log10(dens/densmin) + 1.0_PR - real(idens,PR)
  deltatemp = btemp*log10(temp/tempmin) + 1.0_PR - real(itemp,PR)

! Linearly interpolate between logarithmic table points
 
 kappa_mean = kappa(idens,itemp)*deltadens*deltatemp + &
       & kappa(idens+1,itemp+1)*(1.0_PR - deltadens)*(1.0_PR - deltatemp) + &
       & kappa(idens,itemp+1)*deltadens*(1.0_PR - deltatemp) + &
       & kappa(idens+1,itemp)*(1.0_PR - deltadens)*deltatemp

 kappa_rosseland = kappar(idens,itemp)*deltadens*deltatemp + &
       & kappar(idens+1,itemp+1)*(1.0_PR - deltadens)*(1.0_PR - deltatemp) + &
       & kappar(idens,itemp+1)*deltadens*(1.0_PR - deltatemp) + &
       & kappar(idens+1,itemp)*(1.0_PR - deltadens)*deltatemp
 
  kappa_planck = kappap(idens,itemp)*deltadens*deltatemp + &
       & kappap(idens+1,itemp+1)*(1.0_PR - deltadens)*(1.0_PR - deltatemp) + &
       & kappap(idens,itemp+1)*deltadens*(1.0_PR - deltatemp) + &
       & kappap(idens+1,itemp)*(1.0_PR - deltadens)*deltatemp

  return
END SUBROUTINE getkappa



! ============================================================================
! FIND_IDENS
! D. Stamatellos - 3/1/2008
! Finds grid index corresponding to the real inputted density variable
! ============================================================================
SUBROUTINE find_idens(dens,idens_index)
  use definitions
  use Eos_module, only: eos_dens,dim_dens,bdens,densmin
  implicit none

  real(kind=PR), intent(in)   :: dens       ! Particle density
  integer, intent(out)  :: idens_index      ! Density table index

  idens_index = int(bdens*log10(dens/densmin)) + 1
  if (idens_index < 1) idens_index = 1
  if (idens_index >= dim_dens) idens_index = dim_dens - 1

  return
END SUBROUTINE find_idens



! ============================================================================
! FIND_ITEMP
! D. Stamatellos - 3/1/2008
! Finds grid index corresponding to real inputted temperature
! ============================================================================
SUBROUTINE find_itemp(temp,itemp_index)
  use definitions
  use Eos_module, only: eos_temp,dim_temp,btemp,tempmin
  implicit none

  real(kind=PR), intent(in)   :: temp       ! Particle temperature
  integer, intent(out)  :: itemp_index      ! Temperature table index

  itemp_index = int(btemp*log10(temp/tempmin)) + 1
  if (itemp_index < 1) itemp_index = 1
  if (itemp_index >= dim_temp) itemp_index = dim_temp - 1
  
  return
END SUBROUTINE find_itemp



! ============================================================================
! FIND_TEMP_FROM_ENERGY
! D. Stamatellos - 3/1/2008
! Finds temp and itemp from energy
! ============================================================================
SUBROUTINE find_temp_from_energy(idens,energy,itemp,temp)
  use definitions
  use Eos_module, only: eos_energy,eos_temp,dim_temp
  implicit none

  integer, intent(in)  :: idens             ! ..
  real(kind=PR), intent(in)     :: energy   ! ..
  integer, intent(out) :: itemp             ! ..
  real(kind=PR), intent(out)    ::temp      ! ..
  
  integer :: mid                            ! ..
  integer :: lo                             ! ..
  integer :: hi                             ! ..
  real(kind=PR) slope                       ! ..

! Find energy_index
  lo = 1; hi = dim_temp; itemp = 0
    
! ----------------------------------------------------------------------------
  if (energy < eos_energy(idens,1)) then
     temp = eos_temp(1) 
     itemp = 1
  else if (energy > eos_energy(idens,dim_temp)) then
     temp = eos_temp(dim_temp)
     itemp = dim_temp
  else

     do 
        mid = (lo + hi)/2
        if (hi - lo <= 1) then
           itemp = lo
           exit
        else if (eos_energy(idens,mid).gt.energy) then
           hi = mid
        else
           lo = mid
        end if
     end do
     
     if (itemp >= dim_temp) itemp = dim_temp - 1

     slope = (eos_energy(idens,itemp+1) - eos_energy(idens,itemp))/ &
          &(eos_temp(itemp+1) - eos_temp(itemp))
     temp = real(eos_temp(itemp) + (energy - eos_energy(idens,itemp))/slope,PR)
     
  end if
! ----------------------------------------------------------------------------
  
  return
END SUBROUTINE find_temp_from_energy



! ============================================================================
! EOSENERGY
! D. Stamatellos - 3/1/2008
! Calculates energy ....
! ============================================================================
FUNCTION eosenergy(dens,temp,idens,itemp)
  use definitions 
  use Eos_module, only: eos_dens,eos_temp,eos_energy,&
       &bdens,btemp,densmin,tempmin
  implicit none

  integer, intent(in) :: idens              ! ..
  integer, intent(in) :: itemp              ! ..
  real(kind=PR), intent(in)  :: dens,temp   ! ..
  
  real(kind=PR) :: eosenergy                ! ..
  real(kind=PR) :: deltadens                ! ..
  real(kind=PR) :: deltatemp                ! ..

  deltadens = bdens*log10(dens/densmin) + 1.0_PR - real(idens,PR)
  deltatemp = btemp*log10(temp/tempmin) + 1.0_PR - real(itemp,PR)

  eosenergy = eos_energy(idens,itemp)*&
       &(1.0_PR - deltadens)*(1.0_PR - deltatemp) + &
       &eos_energy(idens+1,itemp)*deltadens*(1.0_PR - deltatemp) + &
       &eos_energy(idens,itemp+1)*(1.0_PR - deltadens)*deltatemp + &
       &eos_energy(idens+1,itemp+1)*deltadens*deltatemp

END FUNCTION eosenergy



! ============================================================================
! EOSMU
! D. Stamatellos - 3/1/2008
! Calculates value of mu (mean molecular mass) for a given temperature 
! and density.  
! ============================================================================
FUNCTION eosmu(dens,temp,idens,itemp)
  use definitions
  use Eos_module, only: eos_dens,eos_temp,eos_mu,bdens,btemp,densmin,tempmin
  implicit none

  real(kind=PR), intent(in)   :: dens,temp
  integer, intent(in):: idens,itemp

  real(kind=PR) :: eosmu
  real(kind=PR) :: deltadens
  real(kind=PR) :: deltatemp

  deltadens = bdens*log10(dens/densmin) + 1.0_PR - real(idens,PR)
  deltatemp = btemp*log10(temp/tempmin) + 1.0_PR - real(itemp,PR)

  eosmu = eos_mu(idens,itemp)*(1.0_PR - deltadens)*(1.0_PR - deltatemp) + &
       &eos_mu(idens+1,itemp)*deltadens*(1.0_PR - deltatemp) + &
       &eos_mu(idens,itemp+1)*(1.0_PR - deltadens)*deltatemp + &
       &eos_mu(idens+1,itemp+1)*deltadens*deltatemp

END FUNCTION eosmu
