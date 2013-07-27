!
! All the rate coefficients are stored here, in an attempt to keep the 
! iterative loop as clean as possible. Note that many of the rate coeff are
! extremely costly to compute (and they can be done often, for each particle,
! each timestep). You might want to think about makeing them into look-up
! tables, that are computed once at the beginning of the simulation. 
!


function k_form_dust(temp, tdust)
!
! The formation of H2 on dust grains, taken from Hollenbach and McKee 1979 
!
use chemistry_constants 
implicit none
double precision :: k_form_dust
double precision :: temp, tdust
double precision :: temp2, tdust2
double precision :: numerator, denominator

temp2 = temp/100.0D0
tdust2 = tdust/100.0D0
numerator = dust_to_gas_ratio_units_solar * deff * 3d-17 * temp2**.5d0
denominator = 1.0 + 0.4*sqrt(temp2 + tdust2) + 0.2*temp2 + 0.08*temp2*temp2
k_form_dust = numerator / denominator / (1d0 + 1d4 * dexp(-6d2 / tdust))
end function k_form_dust




function k_diss_h2h(temp)
!
! Collisional dissociation by H2. Note that this rate is only valid in the
! low density limit (i.e. n < n_crit). Should be fine for our ISM purposes.
!
implicit none
double precision :: k_diss_h2h, temp

k_diss_h2h = 6.67d-12 * sqrt(temp) * exp(-(1.0d0 + 63590.0D0/temp))
end function k_diss_h2h




function k_photo_diss(numdens, numdens_H2, temp)
!
! The phtoto-dissociation rate for H2. 
!
use chemistry_constants
implicit none
double precision :: k_photo_diss
double precision :: numdens, numdens_H2, temp
double precision :: column_H2, total_column
double precision :: k_ph_0, tau_d_1000, x, b5, fshield

total_column = numdens * shielding_length
column_H2 = numdens_H2 * shielding_length
k_ph_0 = 3.3d-11 * isrf_chi
tau_d_1000 = 2.0d-21 * total_column   ! optical depth of dust at 1000 Angstroms.
x = column_H2 / 5.0d14
b5  = dsqrt(KBOLTZ * temp / PROTONMASS) / 1d5
fshield = 0.965d0/(1.0d0 + x/b5)**2 + 0.035d0/sqrt(1.0d0 + x) * exp(-8.5d-4 * sqrt(1.0d0 + x))
k_photo_diss = fshield * exp(-tau_d_1000) * k_ph_0
end function k_photo_diss



function k_recomb_b(temp)
!
! The case B recombination coeff. From Ferland et al, 1992, ApJ, 387, 95
!
implicit none
double precision :: k_recomb_b, temp, tinv

tinv = 1.0D0 / temp
k_recomb_b = 2.753d-14 * (315614d0 * tinv)**1.500d0 / ((1d0 + (115188d0 * tinv)**0.407d0)**2.242d0)
end function k_recomb_b




function k_diss_h2h2(temp)
!
! The collisional dissociation coeff for H2-H2 collisions. Using the low-density
! case given in Glover et al. 2010
!
implicit none
double precision :: k_diss_h2h2, temp

k_diss_h2h2 = 5.996D-30 * temp**(4.1881) / (1.0D0 + 6.761D-6 * temp)**5.6881
k_diss_h2h2 = k_diss_h2h2 * exp(-54657.4D0 / temp)
end function k_diss_h2h2
   



function k_colldiss_H_with_e(temp)
!
! Collisional ionization of HI (H + e => H+ + 2e)
! From A97; based on data from J87
!
use chemistry_constants
implicit none
double precision :: k_colldiss_H_with_e
double precision :: temp, te, lnte

te    = temp / temp_1ev
lnte  = dlog(te)
k_colldiss_H_with_e = dexp( -32.71396786d0 + 13.5365560d0  * lnte &
                            - 5.73932875d0  * (lnte**2) &
                            + 1.56315498d0  * (lnte**3) &
                            - 0.28770560d0  * (lnte**4) &
                            + 3.48255977d-2 * (lnte**5) &
                            - 2.63197617d-3 * (lnte**6) &
                            + 1.11954395d-4 * (lnte**7) &
                            - 2.03914985d-6 * (lnte**8) ) 
end function k_colldiss_H_with_e





function k_hp_recomb_dust(temp, yn, abe)
!
! Rate co-eff for the H+ recombination on dist grains 
!
use chemistry_constants
implicit none
double precision :: k_hp_recomb_dust
double precision :: shielding_column, AV, G_dust
double precision :: temp, yn, abe
double precision :: phi, hgrvar1, hgrvar2
double precision :: temp_dep_1, temp_dep_2

shielding_column = yn * shielding_length
AV = AV_conversion_factor * dust_to_gas_ratio_units_solar * shielding_column
G_dust = isrf_chi * dexp(-2.5d0 * AV)
if (abe .eq. 0d0) then
! If the fractional ionization is zero, then there won't be any recombination,
! so the value we use for phi doesn't matter too much -- 1d20 is simply an
! arbitrary large number
!
   phi = 1d20
else
   phi = G_dust * sqrt(temp) / (yn * abe)
endif

if (phi .lt. 1d-6) then
   k_hp_recomb_dust = 1.225d-13 * dust_to_gas_ratio_units_solar
else
   temp_dep_1 = 5.087d2 * temp**1.586d-2
   temp_dep_2 = - 0.4723d0 - 1.102d-5 * log(temp)
   hgrvar1  = 8.074d-6 * phi**1.378d0
   hgrvar2  = (1d0 + temp_dep_1 * phi**temp_dep_2)
   k_hp_recomb_dust = 1.225d-13 * dust_to_gas_ratio_units_solar / (1d0 + hgrvar1 * hgrvar2)
endif
end function k_hp_recomb_dust

