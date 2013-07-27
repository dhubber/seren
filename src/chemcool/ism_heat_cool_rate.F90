function ism_heat_cool_rate(numdens, numdens_H2, numdens_Hp, numdens_H, numdens_e, temp, tdust)
!
! The various non-chemical heating and cooling rates are collected
! in this function. 
!
use chemistry_constants

implicit none

double precision :: ism_heat_cool_rate
double precision :: heating_rate, cooling_rate
double precision :: numdens, numdens_H2, numdens_Hp, numdens_H, numdens_e
double precision :: temp, temp2, tdust
double precision :: cplus_cooling, cplus_colldeex_h, cplus_colldeex_h2, cplus_colldeex_e
double precision :: cIIc10, cIIc01
double precision :: cIIb10, cIIb01
double precision :: cIIvar0, cIIvar1, cIIvar2, cIIvar3
double precision :: cIIn0, cIIn1, abcII 
double precision :: f, hh, coeff
double precision :: lyman_alpha_cooling
double precision :: gas_to_dust_rate
double precision :: PEvar0, PEvar1, PEvar2, eps, G_dust
double precision :: pe_heating
double precision :: shielding_column, AV


heating_rate = 0
cooling_rate = 0
temp2 = temp * 1.0D-02

!
! C+ cooling 
!

cplus_cooling = 0
! CII - HI (HM89 below 2000K; K86 above 2000K)
if (temp .le. 2d3) then
   cplus_colldeex_h =  8d-10 * temp2**0.07d0
else
   cplus_colldeex_h = 3.113619d-10 * temp2**0.385d0
endif

! CII - H2 (Below 250K, we use the fit from WBV96 to the data from FL77
! Above 250K, we assume that the rate scales the same as the low-temp. HI rate)
if (temp .lt. 250d0) then
   f  = (4.7d-10 + 4.6d-13 * temp) * fortho
   hh = (2.5d-10 * temp**0.12d0) * fpara
else
   f  = (5.85d-10 * temp**0.07d0) * fortho
   hh = (4.85d-10 * temp**0.07d0) * fpara
endif
!print *, "fortho, fpara", fortho, fpara
cplus_colldeex_h2 = f + hh

!  CII - electron (WB02)
if (temp .le. 2d3) then
   cplus_colldeex_e  = 3.86d-7 / dsqrt(temp2)
else
   cplus_colldeex_e = 2.426206d-7 / temp2**0.345d0
endif

cIIc10 = cplus_colldeex_h*numdens_H + cplus_colldeex_h2*numdens_H2 + cplus_colldeex_e*numdens_e 
!print *, "cIIc10   ", cIIc10

if ( cIIc10.eq.0 ) then
   cplus_cooling = 0
else
   !print *, "temp, exp", temp, dexp(-91.25d0 / temp)
   cIIc01 = cIIc10 * 2d0 * dexp(-91.25d0 / temp)
   call compute_stim(cIIa10, cIIe10, CMB_temp, cIIb10)
   !print *, "cIIa10, cIIe10, CMB_temp, cIIb10", cIIa10, cIIe10, CMB_temp, cIIb10 
   cIIb01 = 2d0 * cIIb10
   cIIvar0 = (cIIc10 + cIIa10 + cIIb10)
   cIIvar1 = (cIIc01 + cIIb01)
   cIIn0 = cIIvar0 / (cIIvar0 + cIIvar1)
   cIIn1 = cIIvar1 / (cIIvar0 + cIIvar1)
   cIIvar2 = (cIIa10 + cIIb10) * cIIe10 * cIIn1
   cIIvar3 = cIIb01 * cIIe10 * cIIn0
   !print *, "var0, var1", cIIvar0, cIIvar1
   !print *, "var2, var3", cIIvar2, cIIvar3
   abcII = abundance_C ! we assume all C is in C+   
   !print *, "abcII * numdens ", abcII * numdens
   cplus_cooling = (cIIvar2 - cIIvar3) * abcII * numdens
end if

!
! Lyman-alpha cooling
!

lyman_alpha_cooling = 0
if ( temp .gt. 1d3 ) then
   lyman_alpha_cooling = numdens_H * numdens_e * 7.5d-19 * dexp(-1.18348d5 / temp) / (1d0 + dsqrt(temp / 1d5))
endif


!
! Dust-gas transfer cooling/heating from HM89, eqn 2.15

coeff = 3.8d-33 * sqrt(temp) * ( 1d0 - 0.8 * dexp(-75d0 / temp)  )
gas_to_dust_rate = coeff * (temp - tdust) * numdens * numdens * dust_to_gas_ratio_units_solar

!
! Photoelectric heating
!

! If there's no UV field, or if the electron density is very low (in which
! case the photoheating efficiency will also be very low), then we set the
! rates to zero. Otherwise, we compute the heating rate using the
! Bakes & Tielens (1994) formula (as modified by Wolfire et al, 2003).
!

shielding_column = numdens * shielding_length
AV = AV_conversion_factor * dust_to_gas_ratio_units_solar * shielding_column

G_dust = isrf_chi * dexp(-2.5d0 * AV)
if (G_dust .eq. 0d0 .or. dust_to_gas_ratio_units_solar .eq. 0d0 ) then
   pe_heating = 0d0
elseif (numdens_e .lt. 1d-9 * G_dust * dsqrt(temp) / phi_pah) then
   pe_heating = 0d0
else
   PEvar0 = G_dust * dsqrt(temp) / phi_pah / numdens_e
   PEvar1 = (1d0 + 4d-3 * PEvar0**0.73)
   PEvar2 = (1d0 + 2d-4 * PEvar0)

   eps = (4.9d-2 / PEvar1) + ( (3.7d-2 * (temp / 1d4)**0.7d0) / PEvar2)
   pe_heating = 1.3d-24 * eps * G_dust * numdens * dust_to_gas_ratio_units_solar
endif

!
! Add up all the processes
!


!print *, "Heating/cooling rates:"
!print *, "pe_heating          ", pe_heating
!print *, "cplus_cooling       ", cplus_cooling
!print *, "lyman_alpha_cooling ", lyman_alpha_cooling
!print *, "gas_to_dust_rate    ", gas_to_dust_rate

!print *, "Lyman alpha term ", (7.5d-19 * dexp(-1.18348d5 / temp) / (1d0 + dsqrt(temp / 1d5)))
!print *, "nh2,nh+,nh,ne", numdens, numdens_H2, numdens_Hp, numdens_H, numdens_e

!stop

heating_rate = heating_rate + pe_heating
cooling_rate = cooling_rate + cplus_cooling 
cooling_rate = cooling_rate + lyman_alpha_cooling 
cooling_rate = cooling_rate + gas_to_dust_rate
ism_heat_cool_rate = heating_rate - cooling_rate

end function ism_heat_cool_rate
