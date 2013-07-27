module chemistry_constants

implicit none

double precision :: abundance_C, abundance_Si
parameter ( abundance_C = 1.4D-4 )
parameter ( abundance_Si = 1.5D-5 )

double precision :: rel_tol_H2, rel_tol_Hp
parameter ( rel_tol_H2 = 1D-4 )
parameter ( rel_tol_Hp = 1D-4 )

double precision :: abs_tol_H2, abs_tol_Hp
parameter ( abs_tol_H2 = 1D-7 )
parameter ( abs_tol_Hp = 1D-12 )

double precision :: energy_rate_limit, rel_tol_eng
parameter ( energy_rate_limit = 0.1 ) 
parameter ( rel_tol_eng = 1D-4 )

double precision :: ABHE, PROTONMASS, KBOLTZ, EV, GAMMA_GAS, YEAR
parameter ( ABHE = 0.1 )
parameter ( PROTONMASS = 1.6726D-24 )
parameter ( KBOLTZ = 1.3806D-16 )
parameter ( EV = 1.602176565D-12 )
parameter ( GAMMA_GAS = 5./3. )
parameter ( YEAR = 31557600.0D0 )

double precision :: temp_1ev
parameter (temp_1ev = EV / KBOLTZ)

double precision :: cIIa10, cIIe10
parameter (cIIa10 = 2.291d-6)
parameter (cIIe10 = 9.125d1 * KBOLTZ)

!
! These things might change
!

double precision :: dust_to_gas_ratio, dust_to_gas_ratio_units_solar, deff, phi_pah
double precision :: isrf_chi, shielding_length
double precision :: AV_conversion_factor
double precision :: cr_ion_rate

parameter ( phi_pah = 0.5)
parameter ( dust_to_gas_ratio = 0.01 )
parameter ( dust_to_gas_ratio_units_solar = 1.0 ) ! the ratio above is the solar value
parameter ( deff = 1.0 )
parameter ( AV_conversion_factor = 5.348D-22 )

parameter ( cr_ion_rate = 3.0D-17 )
parameter ( isrf_chi = 1.7 )
parameter ( shielding_length = 3.85136e+19 ) ! 6 pc -- could also set it to the cloud radius?

double precision h2_form_kin
parameter ( h2_form_kin = 0.12 ) ! fraction of H2 formation heating that goes into kinetic energy

double precision :: opratio, fortho, fpara
parameter ( opratio = 2.4d0 )
parameter ( fortho  = opratio / (1.d0 + opratio) )
parameter ( fpara   = 1.d0 / (1.d0 + opratio) )

double precision :: CMB_temp
parameter ( CMB_temp = 2.7 )

end module chemistry_constants

