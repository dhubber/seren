subroutine solve_chem_timestep(dt, iteration_success, density, energy_density_init, pdv_and_shocks, abundances)

use chemistry_constants

!
! pcc - 4th July 2013
!
! Solves the reaction network outlined in Glover & Mac Low 2007a
! using a 1st order backwards-differencing scheme, with Gauss-Seidel
! iteration. Basically a stripped down version of Simon Glover's chemistry
! cooling module that is used in ZEUS, FLASH, GADGET2, and Arepo. More 
! details can be found in Glover et al. 2010, and the references therein. 
! 
! Note that this subroutine may need to be called several times within the code's
! timestep, to "sub-cycle" the chemistry and cooling. This is taken care of in the 
! driver subroutine "do_chemcool_step" 
!

implicit none

! iteration control
integer :: iiter, iiter_max
parameter (iiter_max = 100)
integer :: iteration_success, num_converged
double precision :: rel_diff_H2, rel_diff_Hp, rel_diff_eng

! abundances, etc
double precision :: abundances(2)    ! the ones we track (H2 and H+)
double precision :: numdens_H2_init, numdens_Hp_init, numdens_H_init, numdens_e_init
double precision :: numdens_H2, numdens_Hp, numdens_H, numdens_e
double precision :: abundance_H2, abundance_Hp, abundance_H, abundance_e
double precision :: abundance_H2_old, abundance_Hp_old, energy_density_old
double precision :: temp, tdust
double precision :: temp_init, tdust_init
double precision :: energy_density_init, energy_density
double precision :: density, numdens, numdens_tot
double precision :: pdv_and_shocks


! rate equation
double precision :: dt
double precision :: C_H2, D_H2
double precision :: C_Hp, D_Hp
double precision :: chem_heat_rate, chem_cool_rate

! rate co-efficients
double precision :: k_form_dust, k_diss_h2h, k_diss_h2h2, k_photo_diss
double precision :: k_colldiss_H_with_e, k_recomb_b, k_hp_recomb_dust

! the reaction rates typically enter the main rate equation AND the 
! heating/cooling functions associated with that rate. To save calling
! the costly rate functions twice, we store them in the following variable
double precision :: reaction_rate

! energy equation stuff
double precision :: energy_diff, ism_dteng 
double precision :: ism_heat_cool_rate
double precision :: heating_h2_form_dust


!
! Chemical rate equations work in number density of colliders, but it
! is more convenient to store the abundances of the colliders, so we
! need to convert from abunance to number density first. We'll convert
! back again at the end of the iteration.
! Unpack the initial abundances and energy/temp, etc...
! Also work out the e and H abundance from conservation equations
!

tdust_init = 15. ! Doesn't change in this implementation!
numdens = density /  ((1.0 + 4.0 * ABHE) * PROTONMASS)
numdens_tot = (1.0D0 + ABHE - abundances(1) + abundances(2)) * numdens
temp_init = energy_density_init * (GAMMA_GAS - 1.0D0) / ( numdens_tot * KBOLTZ )

numdens_H2_init = numdens * abundances(1)

numdens_Hp_init = numdens * abundances(2)

abundance_H = max(1d0 - 2d0 * abundances(1) - abundances(2), 0d0)
numdens_H_init = numdens * abundance_H

abundance_e = abundances(2) + abundance_C + abundance_Si
numdens_e_init = numdens * abundance_e

!print *, "START:", real(dt), abundances(1), abundances(2), abundance_H, abundance_e

!
! Need to update a few things before entering the interation. These
! same variables will be updated again at the END of the loop, ready
! ready for the next iteration (assuming that there is one)
!

numdens_H2 = numdens_H2_init
numdens_Hp = numdens_Hp_init
numdens_H = numdens_H_init
numdens_e = numdens_e_init
temp = temp_init
tdust = tdust_init

!print *, "1) In chemistry solver: "
!print *, "1) Abund: ", abundances(1), abundances(2)
!print *, "1) num dens:", numdens, numdens_H2, numdens_Hp
!print *, "1) Temps :", temp, tdust
!print *, "1) energy density: ", energy_density
!
! Entering iterative loop
!

iteration_success = 0
do iiter = 1, iiter_max

   !
   ! Initialise the heating cooling
   !

   chem_heat_rate = 0
   chem_cool_rate = 0

   !
   ! The H2 formation / destruction reactions come first
   !

   C_H2 = 0
   D_H2 = 0

   ! Reaction 1: H2 formation on dust grains
   reaction_rate = k_form_dust(temp, tdust) * numdens_H * numdens_H
   C_H2 = C_H2 + reaction_rate
   chem_heat_rate = chem_heat_rate + 4.48d0 * EV * h2_form_kin * reaction_rate

   ! Reaction 2: H2 dissociation (H2 + H -> 3H) Taken from Glover et al. (2010)
   reaction_rate = k_diss_h2h(temp)
   D_H2 = D_H2 + reaction_rate * numdens_H
   chem_cool_rate = chem_cool_rate + 4.48d0 * EV * reaction_rate * numdens_H * numdens_H2 

   ! Reaction 3: H2 dissociation (H2 + H2 -> 2H + H2) 
   reaction_rate = k_diss_h2h2(temp)
   D_H2 = D_H2 + reaction_rate * numdens_H2 
   chem_cool_rate = chem_cool_rate + 4.48d0 * EV * reaction_rate * numdens_H2 * numdens_H2

   ! Reaction 4: H2 photo-dissociation (H2 + UV-photon -> 2H)
   reaction_rate = k_photo_diss(numdens, numdens_H2, temp)
   D_H2 = D_H2 + reaction_rate
   chem_heat_rate = chem_heat_rate + 6.4D-13 * reaction_rate * numdens_H2

   !
   ! Do the H2 reaction....
   !

   numdens_H2 = (numdens_H2_init + C_H2 * dt) / (1 + D_H2 * dt)

   !
   ! update the H abundance from the conservation laws
   !

   abundance_H2 =  numdens_H2/numdens
   abundance_H = max(1d0 - 2d0 * abundance_H2 - abundance_Hp, 0d0)
   numdens_H = numdens * abundance_H

   !
   ! Now the H+ formation / destruction reactions
   !

   C_Hp = 0
   D_Hp = 0

   ! Reaction 5: Cosmic-ray ionisation rate (H + cr -> H+ + e)
   C_Hp = C_Hp + cr_ion_rate * numdens_H
   chem_heat_rate = chem_heat_rate + 3.2D-28*(cr_ion_rate/1.0D-17) * numdens_H 

   ! Reaction 6: Collisional ionisation with electrons (H + e -> H+ + 2e)
   reaction_rate = k_colldiss_H_with_e(temp) 
   C_Hp = C_Hp + reaction_rate * numdens_H * numdens_e
   chem_cool_rate = chem_cool_rate + 13.6d0 * eV * reaction_rate * numdens_H * numdens_e 

   ! Reaction 7: Recombination (H+ + e -> H + photon)
   reaction_rate = k_recomb_b(temp)
   D_Hp = D_Hp + reaction_rate * numdens_e 
   chem_cool_rate = chem_cool_rate + KBOLTZ * temp * reaction_rate * numdens_Hp * numdens_e 

   ! Reaction 8: Recombination on dust grains (H+ + e + grain -> H + grain)
   reaction_rate = k_hp_recomb_dust(temp, numdens, abundance_e)
   D_Hp = D_Hp + reaction_rate * numdens_e 
   chem_cool_rate = chem_cool_rate + KBOLTZ * temp * reaction_rate * numdens_Hp * numdens_e

   !
   ! Do the H+ reaction 
   !

   numdens_Hp = (numdens_Hp_init + C_Hp * dt) / (1 + D_Hp * dt)

   !
   ! The energy from these reactions.
   !
 
   ism_dteng = ism_heat_cool_rate(numdens, numdens_H2, numdens_Hp, numdens_H, numdens_e, temp, tdust)
   energy_diff = abs(chem_heat_rate - chem_cool_rate + ism_dteng + pdv_and_shocks) * dt / energy_density_init
   if ( energy_diff.gt.energy_rate_limit ) then
      !
      ! The energy wants to change by a greater amount than we allow. Need to reduce the
      ! timestep. We now send an appropriate value back to the driving subrountine, and
      ! exit this attempt without updating the internal energy (or temp) or the abundances.
      !
      iteration_success = -1
      print *, "Dt too big for d_eng/dt! Aborting to try again with a smaller dt"
      print *, energy_density_init, energy_diff
      exit
   else
      energy_density = energy_density_init + (chem_heat_rate - chem_cool_rate + ism_dteng + pdv_and_shocks) * dt
   end if


   !
   ! update the H/e abundance from the conservation laws
   !

   abundance_H2 = numdens_H2 / numdens
   abundance_Hp = numdens_Hp / numdens
   abundance_H = max(1d0 - 2d0 * abundance_H2 - abundance_Hp, 0d0)
   abundance_e = abundance_Hp + abundance_C + abundance_Si

   numdens_H = numdens * abundance_H
   numdens_e = numdens * abundance_e

   !
   ! End of the iteration? Check convergence. If the abundances of the species
   ! are below some minimum value, 'abs_tol_X', we don't care, and can leave the
   ! loop (essentially equivalent to taking a Euler step). Otherwise we need to 
   ! check for convergence for that species. The convergence is then decided by
   ! checking if the difference between the old and new abundance is within some
   ! tolerance parameter (set in chemistry_constants).
   !

   num_converged = 0
   if ( iiter.gt.1 ) then 
      ! check the H2 convergence
      if ( abundance_H2.gt.abs_tol_H2 ) then
         rel_diff_H2 = abs(abundance_H2 -  abundance_H2_old) / abundance_H2
         if ( rel_diff_H2.lt.rel_tol_H2 ) num_converged = num_converged + 1
      else
         num_converged = num_converged + 1
      end if 
      ! check the H+ convergence
      if ( abundance_Hp.gt.abs_tol_Hp ) then
         rel_diff_Hp = abs(abundance_Hp -  abundance_Hp_old) / abundance_Hp
         if ( rel_diff_Hp.lt.rel_tol_Hp ) num_converged = num_converged + 1
      else
         num_converged = num_converged + 1
      end if 
      ! check the energy convergence -- no absolute check here, just relative
      rel_diff_eng = abs(energy_density -  energy_density_old) / energy_density
      if ( rel_diff_eng.lt.rel_tol_eng ) num_converged = num_converged + 1
      !
      ! If we're converged, then we can exit the loop
      !
      if ( num_converged.eq.3 ) then
         iteration_success = iiter
         exit
      end if
   end if 

   !
   ! Still not converged. Need to go around again! Get things ready...
   !

   abundance_H2_old = abundance_H2
   abundance_Hp_old = abundance_Hp
   energy_density_old = energy_density
   numdens_tot = (1.0D0 + ABHE - abundance_H2 + abundance_Hp) * numdens
   temp = energy_density * (GAMMA_GAS - 1.0D0) / numdens_tot / KBOLTZ
end do

!
! convert back from number densities of species to abundance
!

!print *, "2) In chemistry solver: "
!print *, "2) Abund: ", abundance_H2, abundance_Hp
!print *, "2) num dens:", numdens, numdens_H2, numdens_Hp
!print *, "2) Temps :", temp, tdust
!print *, "2) energy density: ", energy_density
!print *, "2) iteration_success", iteration_success

if ( iteration_success.gt.0 ) then
   abundance_H2 = numdens_H2 / numdens
   abundance_Hp = numdens_Hp / numdens
   abundances(1) = abundance_H2
   abundances(2) = abundance_Hp
   energy_density_init = energy_density
end if

end subroutine solve_chem_timestep


