! UNITS.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Calculates all units and scaling variables.  
! To convert code units to specified unit (e.g. pc, myr, msun), 
! multiply by rscale, mscale, tscale etc.
! To convert code units to S.I. units, multiply by rscale*r_SI, 
! mscale*m_SI, tscale*t_SI etc.
! To convert code units to cgs units, multiply by rscale*rcgs, 
! mscale*mcgs, tscale*tcgs etc.
! (NOTE: rscale and mscale are chosen as parameters, and all other 
! scales are calculated based upon these)
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE units
  use interface_module, only : paramerror
  use HP_module, only : intmax, N_LyC, a_star, Xfrac
  use Eos_module, only : rad_const
  use scaling_module
  use constant_module
  use hydro_module
  implicit none

  real(kind=DP) :: denom         ! Auxilary variable to store denominator
  real(kind=DP) :: invdenom2     ! 1 / denom**2
  real(kind=DP) :: num           ! Auxilary variable to store numerator

  debug1("Setting up scaling units for simulation [units.F90]")


! Set all scaling variables to unity if we are using the dimensionless flag
! ----------------------------------------------------------------------------
  if (dimensionless) then
     rscale      = 1.0_DP; r_SI      = 1.0_DP; rcgs      = 1.0_DP; runit      = ""
     mscale      = 1.0_DP; m_SI      = 1.0_DP; mcgs      = 1.0_DP; munit      = ""
     tscale      = 1.0_DP; t_SI      = 1.0_DP; tcgs      = 1.0_DP; tunit      = ""
     vscale      = 1.0_DP; v_SI      = 1.0_DP; vcgs      = 1.0_DP; vunit      = ""
     ascale      = 1.0_DP; a_SI      = 1.0_DP; acgs      = 1.0_DP; aunit      = ""
     rhoscale    = 1.0_DP; rho_SI    = 1.0_DP; rhocgs    = 1.0_DP; rhounit    = ""
     sigmascale  = 1.0_DP; sigma_SI  = 1.0_DP; sigmacgs  = 1.0_DP; sigmaunit  = ""
     Pscale      = 1.0_DP; P_SI      = 1.0_DP; Pcgs      = 1.0_DP; Punit      = ""
     fscale      = 1.0_DP; f_SI      = 1.0_DP; fcgs      = 1.0_DP; funit      = ""
     Escale      = 1.0_DP; E_SI      = 1.0_DP; Ecgs      = 1.0_DP; Eunit      = ""
     momscale    = 1.0_DP; mom_SI    = 1.0_DP; momcgs    = 1.0_DP; momunit    = ""
     angmomscale = 1.0_DP; angmom_SI = 1.0_DP; angmomcgs = 1.0_DP; angmomunit = ""
     angvelscale = 1.0_DP; angvel_SI = 1.0_DP; angvelcgs = 1.0_DP; angvelunit = ""
     dmdtscale   = 1.0_DP; dmdt_SI   = 1.0_DP; dmdtcgs   = 1.0_DP; dmdtunit   = ""
     Lscale      = 1.0_DP; L_SI      = 1.0_DP; Lcgs      = 1.0_DP; Lunit      = ""
     kappascale  = 1.0_DP; kappa_SI  = 1.0_DP; kappacgs  = 1.0_DP; kappaunit  = ""
     Bscale      = 1.0_DP; B_SI      = 1.0_DP; Bcgs      = 1.0_DP; Bunit      = ""
     Qscale      = 1.0_DP; Q_SI      = 1.0_DP; Qcgs      = 1.0_DP; Qunit      = ""
     Jscale      = 1.0_DP; J_SI      = 1.0_DP; Jcgs      = 1.0_DP; Junit      = ""
     uscale      = 1.0_DP; u_SI      = 1.0_DP; ucgs      = 1.0_DP; uunit      = ""
     tempscale   = 1.0_DP; temp_SI   = 1.0_DP; tempcgs   = 1.0_DP; tempunit   = ""
     dudtscale   = 1.0_DP; dudt_SI   = 1.0_DP; dudtcgs   = 1.0_DP; dudtunit   = ""
     Pconst      = 1.0_PR
     Pconst2     = 1.0_PR
     rad_const   = 1.0_PR
     sound_const = 1.0_PR
     newsound_const = 1.0_PR
     
! ----------------------------------------------------------------------------
else

   ! Length unit in S.I. units (m)
   if (runit=="mpc") then
      r_SI = 1.0E6_DP*r_pc
   else if (runit=="kpc") then
      r_SI = 1.0E3_DP*r_pc
   else if (runit=="pc") then
      r_SI = r_pc
   else if (runit=="au") then
      r_SI = r_au
   else if (runit=="r_sun") then
      r_SI = r_sun
   else if (runit=="r_earth") then
      r_SI = r_earth
   else if (runit=="km") then
      r_SI = 1000.0_DP
   else if (runit=="m") then
      r_SI = 1.0_DP
   else if (runit=="cm") then
      r_SI = 0.01_DP
   else
      call paramerror("Unrecognised runit")
   end if
   rcgs = 100.0_DP*r_SI
   
   
   ! Mass units in S.I. units (kg)
   if (munit=="m_sun") then
      m_SI = m_sun
   else if (munit=="m_jup") then
      m_SI = m_jup
   else if (munit=="m_earth") then
      m_SI = m_earth
   else if (munit=="kg") then
      m_SI = 1.0_DP
   else if (munit=="g") then
      m_SI = 1.0E-3_DP
   else
      call paramerror("Unrecognised munit")
   end if
   mcgs = 1.0E3_DP*m_SI
   
   
   ! Calculate time units in S.I. units (s).  Scale the time unit so that 
   ! G = 1 (for gravity switched both on and off).  
   tscale = ((rscale*r_SI)**1.50_DP) / sqrt(mscale*m_SI*G_const)
   if (tunit=="gyr") then
      t_SI = 1.0D3*myr
   else if (tunit=="myr") then
      t_SI = myr
   else if (tunit=="yr") then
      t_SI = yr
   else if (tunit=="day") then
      t_SI = day
   else if (tunit=="sec") then
      t_SI = 1.0_DP
   else
      call paramerror("Unrecognised tunit")
   end if
   tcgs = t_SI
   tscale = tscale / t_SI
   
   
   ! Velocity units in S.I. units (i.e. m/s)
   vscale = (rscale*r_SI) / (tscale*t_SI)
   if (vunit=="km_s") then
      v_SI = 1000.0_DP
   else if (vunit=="pc_myr") then
      v_SI = r_pc / myr
   else if (vunit=="au_yr") then
      v_SI = r_au / yr
   else if (vunit=="m_s") then
      v_SI = 1.0_DP
   else if (vunit=="cm_s") then
      v_SI = 0.01_DP
   else
      call paramerror("Unrecognised vunit")
   end if
   vcgs = 100.0_DP*v_SI
   vscale = vscale / v_SI
   
   
   ! Acceleration units in S.I. units (i.e. m/s^2)
   ascale = (rscale*r_SI) / ((tscale*t_SI)*(tscale*t_SI))
   if (aunit=="km_s2") then
      a_SI = 1000.0_DP
   else if (aunit=="km_s_myr") then
      a_SI = 1000.0_DP / myr
   else if (aunit=="pc_myr2") then
      a_SI = r_pc / (myr * myr)
   else if (aunit=="au_yr2") then
      a_SI = r_au / (yr*yr)
   else if (aunit=="m_s2") then
      a_SI = 1.0_DP
   else if (aunit=="cm_s2") then
      a_SI = 0.01_DP
   else
      call paramerror("Unrecognised aunit")
   end if
   acgs = 100.0_DP*a_SI
   ascale = ascale / a_SI
   
   
   ! Density units in S.I. units (i.e. kg/m^3)
   denom = (rscale*r_SI)
   rhoscale = (mscale*m_SI) / (denom)**3
   if (rhounit=="m_sun_pc3") then
      denom = r_pc*r_pc*r_pc
      rho_SI = m_sun / denom
   else if (rhounit=="kg_m3") then
      rho_SI = 1.0_DP
   else if (rhounit=="g_cm3") then
      rho_SI = 1.0E3_DP
   else
      call paramerror("Unrecognised rhounit")
   end if
   rhocgs = 1.0E-3_DP*rho_SI
   rhoscale = rhoscale / rho_SI
   
   
   ! Column density units in S.I. units (i.e. kg/m^2)
   denom = (rscale*r_SI)
   sigmascale = (mscale*m_SI) / (denom)**2
   if (sigmaunit=="m_sun_pc2") then
      denom = r_pc*r_pc
      sigma_SI = m_sun / denom
   else if (sigmaunit=="kg_m2") then
      sigma_SI = 1.0_DP
   else if (sigmaunit=="g_cm2") then
      sigma_SI = 10.0_DP
   else
      call paramerror("Unrecognised sigmaunit")
   end if
   sigmacgs = 1.0E-1_DP*sigma_SI
   sigmascale = sigmascale / sigma_SI
   
   
   ! Pressure units in S.I. units (i.e. kg/m/s^2)
   denom = tscale*t_SI 
   invdenom2 = (1.0_DP / denom)**2
   Pscale = invdenom2 * (mscale*m_SI) / (rscale*r_SI)
   if (Punit=="Pa") then
      P_SI = 1.0_DP
   else if (Punit=="bar") then
      P_SI = 1.0D5
   else if (Punit=="g_cms2") then
      P_SI = 10.0_DP
   else
      call paramerror("Unrecognised Punit")
   end if
   Pcgs = 0.1_DP*P_SI
   Pscale = Pscale / P_SI
   
   
   ! Force units in S.I. units (i.e. N)
   fscale = mscale*m_SI*rscale*r_SI / (tscale*tscale*t_SI*t_SI)
   if (funit=="N") then
      f_SI = 1.0_DP
   else if (funit=="dyn") then
      f_SI = 1.0E-5_DP
   else
      call paramerror("Unrecognised funit")
   end if
   fcgs = 1.0E5_DP*f_SI
   fscale = fscale / f_SI
   
   
   ! Energy units in S.I. units (i.e. J)
   Escale = (mscale*m_SI*rscale*r_SI*rscale*r_SI) / (tscale*t_SI*tscale*t_SI)
   if (Eunit=="J") then
      E_SI = 1.0_DP
   else if (Eunit=="erg") then
      E_SI = 1.0E-7_DP
   else if (Eunit=="GJ") then
      E_SI = 1.0E12_DP
   else if (Eunit=="10^40erg") then
      E_SI = 1.0E33_DP
   else
      call paramerror("Unrecognised Eunit")
   end if
   Ecgs = 1.0E7_DP*E_SI
   Escale = Escale / E_SI
   
   
   ! Momentum units in S.I. units (i.e. kg m/s)
   momscale = mscale*m_SI*rscale*r_SI / (tscale*t_SI)
   if (momunit=="m_sunkm_s") then
      mom_SI = m_sun * 1000.0_DP
   else if (momunit=="m_sunau_yr") then
      mom_SI = (m_sun * r_au) / yr
   else if (momunit=="kgm_s") then
      mom_SI = 1.0_DP
   else if (momunit=="gcm_s") then
      mom_SI = 1.0E-5_DP
   else
      call paramerror("Unrecognised momunit")
   end if
   momcgs = 1.0E5_DP*mom_SI
   momscale = momscale / mom_SI
   
   
   ! Angular momentum units in S.I. units (i.e. kg m^2/s)
   angmomscale = (mscale*m_SI*rscale*r_SI*rscale*r_SI) / (tscale*t_SI)
   if (angmomunit=="m_sunkm2_s") then
      angmom_SI = m_sun * 1000.0_DP * 1000.0_DP
   else if (angmomunit=="m_sunau2_yr") then
      angmom_SI = m_sun * r_au * r_au / yr
   else if (angmomunit=="kgm2_s") then
      angmom_SI = 1.
   else if (angmomunit=="gcm2_s") then
      angmom_SI = 1.0E-7_DP
   else
      call paramerror("Unrecognised angmomunit")
   end if
   angmomcgs = 1.0E7_DP*angmom_SI
   angmomscale = angmomscale / angmom_SI
   
   
   ! Angular velocity units in S.I. units (i.e. rad/s)
   angvelscale = 1.0_DP / (tscale*t_SI)
   if (angvelunit=="rad_s") then
      angvel_SI = 1.0_DP
   else
      call paramerror("Unrecognised angvelunit")
   end if
   angvelcgs = angvel_SI
   angvelscale = angvelscale / angvel_SI
   
   
   ! Accretion rate units in S.I. units (i.e. kg/s)
   dmdtscale = mscale*m_SI/(tscale*t_SI)
   if (dmdtunit=="m_sun_myr") then
      dmdt_SI = m_sun / myr
   else if (dmdtunit=="m_sun_yr") then
      dmdt_SI = m_sun / yr
   else if (dmdtunit=="kg_s") then
      dmdt_SI = 1.0_DP
   else if (dmdtunit=="g_s") then
      dmdt_SI = 1.0E-3_DP
   else
      call paramerror("Unrecognised dmdtunit")
   end if
   dmdtcgs = 1.0E3_DP*dmdt_SI
   dmdtscale = dmdtscale / dmdt_SI
   
   
   ! Luminosity units in S.I. units (i.e. J/s)
   Lscale = Escale*E_SI/(tscale*t_SI)
   if (Lunit=="L_sun") then
      L_SI = L_sun
   else if (Lunit=="J_s") then
      L_SI = 1.0_DP
   else if (Lunit=="ergs_s") then
      L_SI = 1.0E-7_DP
   else
      call paramerror("Unrecognised Lunit")
   end if
   Lcgs = 1.0E7_DP*L_SI
   Lscale = Lscale / L_SI
   
   
   ! Opacity units in S.I. units (HARDWIRED TO CM2_G FOR NOW)
   kappascale = (rscale*rscale*r_SI*r_SI)/(mscale*m_SI)
   if (kappaunit=="m2_kg") then
      kappa_SI = 1.0_DP 
   else if (kappaunit=="cm2_g") then
      kappa_SI = 0.10_DP
   else
      call paramerror("Unrecognised kappaunit")
   end if
   kappacgs = 10.0_DP*kappa_SI
   kappascale = kappascale / kappa_SI
   
   
   ! Charge units in S.I. units (i.e. coulombs)
   Qscale = sqrt(mscale*m_SI*rscale*r_SI / mu_0)
   if (Qunit=="C") then
      Q_SI = 1.0_DP
   else if (Qunit=="e") then
      Q_SI = e_charge
   else
      call paramerror("Unrecognised Qunit")
   end if
   Qcgs = Q_SI
   Qscale = Qscale / Q_SI
   
   
   ! Magnetic field in S.I. units (i.e. tesla)
   Bscale = mscale*m_SI / (tscale*t_SI*Qscale*Q_SI)
   if (Bunit=="tesla") then
      B_SI = 1.0_DP
   else if (Bunit=="gauss") then
      B_SI = 0.00010_DP
   else
      call paramerror("Unrecognised Bunit")
   end if
   Bcgs = B_SI
   Bscale = Bscale / B_SI
   
   
   ! Current density units in S.I. units (i.e. C /s/m^2)
   Jscale = Qscale*Q_SI / (tscale*t_SI*rscale*r_SI*rscale*r_SI)
   if (Junit=="C_s_m2") then
      J_SI = 1.0_DP
   else
      call paramerror("Unrecognised Junit")
   end if
   Jcgs = J_SI
   Jscale = Jscale / J_SI
   
   
   ! Specific internal energy units in S.I. units (i.e. J/kg)
   uscale = (rscale*r_SI*rscale*r_SI) / (tscale*t_SI*tscale*t_SI)
   if (uunit=="J_kg") then
      u_SI = 1.0_DP
   else if (uunit=="erg_g") then
      u_SI = 1.0E-4_DP
   else
      call paramerror("Unrecognised uunit")
   end if
   ucgs = 1.0E4_DP*u_SI
   uscale = uscale / u_SI
   
   
   ! Temperature units in S.I. units (i.e. K)
   if (tempunit=="K") then
      temp_SI = 1.0_DP
   else
      call paramerror("Unrecognised tempunit")
   end if
   tempcgs = temp_SI
   tempscale = (m_hydrogen*uscale*u_SI)/(k_boltzmann*temp_SI)
   
   
   ! Rate of change of specific internal energy units in S.I. units (i.e. J/kg/s)
   dudtscale = uscale*u_SI/(tscale*t_SI)
   if (dudtunit=="J_kg_s") then
      dudt_SI = 1.0_DP
   else if (dudtunit=="erg_g_s") then
      dudt_SI = 1.0E-4_DP
   else
      call paramerror("Unrecognised dudtunit")
   end if
   dudtcgs = 1.0E4_DP*dudt_SI
   dudtscale = dudtscale / dudt_SI
   
   
   ! Sound constant in S.I. units (i.e. m/s/K^0.5)
   sound_const    = real(sqrt(k_boltzmann / (mu_bar * m_hydrogen)),PR)
   sound_const    = sound_const / real(vscale*v_SI,PR) ! (/K^0.5)
   newsound_const = real(sqrt(k_boltzmann/m_hydrogen),PR)
   newsound_const = newsound_const / real(vscale*v_SI,PR) ! (/K^0.5)
   num            = rscale*r_SI
   denom          = Escale*E_SI
   rad_const      = real(stefboltz*num*num*tscale*t_SI/denom,PR)
   
   
   ! Pressure constant in S.I. units
   Pconst  = real(k_boltzmann/(mu_bar*m_hydrogen),PR)
   Pconst  = Pconst * real((rho_SI*rhoscale)/(Pscale*P_SI),PR)
   Pconst2 = real(k_boltzmann/m_hydrogen,PR)
   Pconst2 = Pconst2 * real((rho_SI*rhoscale)/(Pscale*P_SI),PR)
   
   
   ! Ionization front integral limit
   !  auxscale = (rhoscale*rho_SI)**2
   !  auxscale = auxscale*(rscale*r_SI)**3
   !  intmax = N_LyC*(m_hydrogen**2)/(4.*PI*a_star)
   a_star = a_star*tscale*t_SI / (rscale*rcgs)**3
   N_LyC  = N_LyC*tscale*t_SI
   intmax = N_LyC*(m_hydrogen/(mscale*m_SI*Xfrac))**2/(4.*PI*a_star) 
   ! NEED TO DIVIDE BY X^2 WHICH IS FRACTION OF MASS OF HYDROGEN (ASK TB)
   !  intmax = intmax / auxscale
   
end if
! ----------------------------------------------------------------------------


! Write a summary of the scaling factors to the screen
#if defined(DEBUG1)
MPI_ROOT
  write(6,*) "Length units           = ", rscale,      trim(runit)
  write(6,*) "Mass units             = ", mscale,      trim(munit)
  write(6,*) "Time units             = ", tscale,      trim(tunit)
  write(6,*) "Velocity units         = ", vscale,      trim(vunit)
  write(6,*) "Acceleration units     = ", ascale,      trim(aunit)
  write(6,*) "Density units          = ", rhoscale,    trim(rhounit)
  write(6,*) "Surface density units  = ", sigmascale,  trim(sigmaunit)
  write(6,*) "Pressure units         = ", Pscale,      trim(Punit)
  write(6,*) "Force units            = ", fscale,      trim(funit)
  write(6,*) "Energy units           = ", Escale,      trim(Eunit)
  write(6,*) "Momentum units         = ", momscale,    trim(momunit)
  write(6,*) "Angular momentum units = ", angmomscale, trim(angmomunit)
  write(6,*) "Angular velocity units = ", angvelscale, trim(angvelunit)
  write(6,*) "Opacity units          = ", kappascale,  trim(kappaunit)
  write(6,*) "B-field units          = ", Bscale,      trim(Bunit)
  write(6,*) "Charge units           = ", Qscale,      trim(Qunit)
  write(6,*) "Current density units  = ", Jscale,      trim(Junit)
  write(6,*) "Internal energy units  = ", uscale,      trim(uunit)
  write(6,*) "Temperature units      = ", tempscale,   trim(tempunit)
  write(6,*) "dudt units             = ", dudtscale,   trim(dudtunit)
MPI_END
#endif

  return
END SUBROUTINE units
