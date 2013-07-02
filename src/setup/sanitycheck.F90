! SANITYCHECK.F90
! C. P. Batty & D. A. Hubber - 1/4/2007
! Checks all compiler flags and all input parameters for any conflicts, 
! bad or illegal values etc.  Will cause the program to stop during 
! runtime and print an error message to screen, rather than compile time. 
! ============================================================================ 

#include "macros.h"

! ============================================================================ 
SUBROUTINE sanitycheck
  use interface_module, only : comperror,paramerror
  use particle_module
  use periodic_module
  use scaling_module
  use neighbour_module
  use hydro_module
  use time_module
  use type_module
  use sink_module
  use tree_module
  use Nbody_module
  use filename_module
  use HP_module
  use turbulence_module
  implicit none

  debug2("Checking compiler flags and parameter values [sanitycheck.F90]")

! Check main compiler flags 
! ----------------------------------------------------------------------------
#if NDIM < 1 || NDIM > 3
  call comperror("NDIM out of range")
#endif

#if defined(PERIODIC_X) || defined(PERIODIC_Y) || defined(PERIODIC_Z)
#if !defined(PERIODIC)
  call comperror("PERIODIC flag not on with X/Y/Z")
#endif
#endif

#if !defined(SPH_SIMULATION) && !defined(NBODY_SPH_SIMULATION) && !defined(NBODY_SIMULATION)
  call comperror("No simulation mode activated")
#endif

#if !defined(M4_KERNEL) && !defined(QUINTIC_KERNEL) && !defined(GAUSSIAN_3H_KERNEL) && !defined(LINEAR_KERNEL)
  call comperror("No valid SPH kernel function selected")
#endif

!#if defined(LINEAR_KERNEL)
!  call comperror("WARNING : The linear kernel is a 'toy' kernel used only &
!       &for testing certain elements of the SPH code, and shouldn't be &
!       &used for actual simulations.  If you wish to use this option, then &
!       &comment out this error message in setup/sanitycheck.F90.")
!#endif


! SPH-only options
! ----------------------------------------------------------------------------
#if defined(SPH_SIMULATION) || defined(NBODY_SPH_SIMULATION)

#if !defined(EULER) && !defined(RUNGE_KUTTA2) && !defined(LEAPFROG_KDK) && !defined(LEAPFROG_DKD)
  call comperror("No valid SPH integration scheme selected")
#endif

#if defined(EULER)
  call comperror("WARNING : The Euler integration scheme is a RUBBISH &
       &integration scheme and should not be used unless you want to find &
       &out how awful it is, or you are having a laugh.  If you really &
       &wish to use this option, then comment out this error message in &
       &setup/sanitycheck.F90.")
#endif

#if !defined(ENERGY_EQN)
  if (typeinfo(boundaryid)%eos == "energy_eqn" .or. &
       &typeinfo(icmid)%eos == "energy_eqn" .or. &
       &typeinfo(gasid)%eos == "energy_eqn") call comperror(&
       &"Energy equation selected but not activated in Makefile")
#endif

#if !defined(ENTROPY_EQN)
  if (typeinfo(boundaryid)%eos == "entropy_eqn" .or. &
       &typeinfo(icmid)%eos == "entropy_eqn" .or. &
       &typeinfo(gasid)%eos == "entropy_eqn") call comperror(&
       &"Entropy equation selected but not activated in Makefile")
#endif

#if defined(IONIZING_UV_RADIATION)
#if !defined(ISOTHERMAL) && !defined(BAROTROPIC)
!  call comperror("IONIZING_UV_RADIATION flag on without an allowed EOS")
#endif
#endif

#if defined(COOLING_HEATING) && !defined(ENERGY_EQN)
  call comperror("COOLING_HEATING flag on and ENERGY_EQN flag off")
#endif

#if defined(RAD_WS) && !defined(SELF_GRAVITY)
  call comperror("RAD_WS flag on and SELF_GRAVITY flag off")
#endif

#if defined(RAD_WS)
#if !defined(AMBIENT_HEATING) && !defined(HDISC_HEATING) && !defined(HDISC_HEATING_3D_SINGLE) && !defined(STAR_HEATING) && !defined(STAR_SIMPLE_HEATING)
  call comperror("No heating mechanism selected for RAD_WS")
#endif
#endif

#if defined(FLUX_LIMITED_DIFFUSION) && !defined(RAD_WS)
  call comperror("FLUX_LIMITED_DIFFUSION flag on and RAD_WS flag off")
#endif

#if defined(ARTIFICIAL_CONDUCTIVITY) && !defined(ENERGY_EQN)
  call comperror("ARTIFICIAL_CONDUCTIVITY activated but not ENERGY_EQN")
#endif

#if !defined(GRAVITY) && !defined(HYDRO)
  call comperror("GRAVITY and HYDRO flags both off")
#endif

#endif
! ----------------------------------------------------------------------------


#if defined(GRAVITY) && NDIM == 1
  call comperror("NDIM == 1 for GRAVITY flag")
#endif

#if defined(EWALD) && !defined(PERIODIC)
  call comperror("PERIODIC flag off, EWALD flag on")
#endif

#if defined(BH_TREE) && defined(SELF_GRAVITY)
#if !defined(GEOMETRIC_MAC) && !defined(GADGET_MAC) && !defined(GADGET2_MAC) && !defined(EIGEN_MAC)
  call comperror("No recognised MAC for BH_TREE and SELF_GRAVITY")
#endif
#endif

#if defined(BINARY_TREE)
  call comperror("Binary tree not completely implemented; currently disabled")
#endif

#if defined(BINARY_TREE) && LEAFMAX == 1
  call comperror("LEAFMAX == 1 does not work for the binary tree")
#endif

#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
!#if !defined(SINKS)
!  call comperror("SINKS flag off and NBODY flag on")
!#endif
#endif

#if defined(SINK_GRAVITY_ONLY) && !defined(SINKS)
  call comperror("SINKS flag off and SINK_GRAVITY_ONLY flag on")
#endif

#if defined(SINK_GRAVITY_ONLY) && defined(SELF_GRAVITY)
  call comperror("SINK_GRAVITY_ONLY and SELF_GRAVITY flags on")
#endif

!#if defined(BH_TREE) && defined(REORDER_TREE)
!  call comperror("REORDER not working correctly with TREE or ALL options")
!#endif

#if !defined(ADAPTIVE_TIMESTEP_LEVELS) && !defined(RESTRICTED_TIMESTEP_LEVELS) && !defined(FIXED_TIMESTEP_LEVELS)
  call comperror("No TIMESTEP option selected")
#endif

#if defined(LEAPFROG_DKD) && defined(CHECK_NEIGHBOUR_TIMESTEPS) && defined(IMMEDIATE_TIMESTEP_REDUCTION)
  call comperror("Combination of SPH_INTEGRATION=LFDKD and &
       &CHECK_NEIB_TIMESTEP=2 potentially dangerous due to secular &
       &error increase (particularly prominent in disk simulations).  &
       &If you wish to use this option, then comment out this error &
       &message in setup/sanity_check.F90.")
#endif

#if defined(SINKS) && defined(SINK_REMOVE_ANGMOM)
!  call comperror("SINK_REMOVE_ANGMOM under development - currently disabled")
#endif

#if defined(HEALPIX) && defined(MULTIPLE_SINK_SOURCES)
#if defined(SINGLE_STATIC_SOURCE) || defined(SINGLE_SINK_SOURCE)
  call comperror("Both single and multiple HEALPix sources selected in Makefile")
#endif
#endif

#if defined(TURBULENT_FORCING)
#if !(NDIM==3)
  call comperror("Turbulent forcing only available in 3D")
#endif
#if !(VDIM==3)
  call comperror("Turbulent forcing only available in 3D")
#endif
#if !defined(RESTRICTED_TIMESTEP_LEVELS)
  call compwarning("Turbulent forcing should be used with &
                   &RESTRICTED_TIMESTEP_LEVELS")
#endif
  call compwarning("Warning - when SEREN is used with turbulent forcing, &
                   &large amounts of stack space are required. If SEREN &
                   &segfaults during startup, ensure the available &
                   &stacksize is unlimited.")
#endif

! Things that don't work with MPI
! ----------------------------------------------------------------------------

#if defined(USE_MPI)
  call compwarning("Use of MPI is experimental!")
#if defined(PERIODIC)
#if defined(PERIODIC_X) || defined(PERIODIC_Y) || defined(PERIODIC_Z)
  call compwarning("Use of MPI and periodic boundary conditions is &
     &even more experimental!")
#endif
#endif

#if defined(NEIGHBOUR_LISTS)
  call compwarning("NEIGHBOUR_LISTS may or may not work with MPI (untested)")
#endif

#if !(NDIM==3)
  call compwarning("MPI only tested in 3 dimensions")
#endif

#if defined(REMOVE_OUTLIERS)
  call compwarning("REMOVE_OUTLIERS may or may not work with MPI (untested)")
#endif

#if !defined(GHOST_PARTICLES)
  call comperror("Must use ghost particles with MPI")
#endif

#if defined(OSPH) || defined(SM2012_SPH) || defined(RTSPH) || defined(RPSPH)
  call comperror("Only standard SPH and grad-h SPH supported with MPI")
#endif

#if defined(ENTROPIC_FUNCTION)
  call comperror("ENTROPIC_FUNCTION not supported with MPI")
#endif

#if defined(SELF_GRAVITY) && defined(PERIODIC)
  call comperror("Use of self-gravity and periodic boundary conditions &
     &is not just a bad idea - with MPI it probably won't work")
#endif

#if defined(N_BODY)
  call comperror("NBODY gravity not supported with MPI")
#endif

#if defined(RAD_WS)
  call comperror("RAD_WS not supported with MPI")
#endif

#if defined(NBODY_SPH_SIMULATION) || defined(NBODY_SIMULATION)
  call comperror("NBODY_SIMULATION / NBODY_SPH_SIMULATION is not supported &
     &with MPI")
#endif

#if defined(CHECK_NEIGHBOUR_TIMESTEPS)
  call comperror("CHECK_NEIB_TIMESTEP is not supported with MPI")
#endif

#if defined(REORDER_TREE)
  call comperror("REORDER_TREE is not supported with MPI")
#endif

#if defined(BINARY_TREE)
  call comperror("BINARY_TREE is not supported with MPI")
#endif

#if defined(CELL_WALK)
  call comperror("CELL_WALK is not supported with MPI")
#endif

#if defined(HEALPIX)
  call comperror("HEALPIX is not MPI-parallelized")
#endif
#endif


! Now check all input parameters 
! ----------------------------------------------------------------------------
  if (rseed < 1) call paramerror("rseed < 1")
  if (sph_endtime < 0.0_DP) call paramerror("sph_endtime < 0")
  if (firstsnap < 0.0_DP) call paramerror("firstsnap < 0")
  if (snaptime <= 0.0_DP) call paramerror("snaptime <= 0")
  if (noutputstep < 1) call paramerror("noutputstep < 1")
  if (ntempstep < 1) call paramerror("ntempstep < 1")
  if (ndiagstep < 1) call paramerror("ndiagstep < 1")
  if (nsinkstep < 1) call paramerror("nsinkstep < 1")
  if (nsnapstep < 1) call paramerror("nsnapstep < 1")
  if (courant_mult < 0.0_DP) call paramerror("courant_mult < 0")
  if (courant_mult > 1.0_DP) call paramerror("courant_mult too big (> 1)")
  if (accel_mult < 0.0_DP) call paramerror("accel_mult < 0")
  if (accel_mult > 1.0_DP) call paramerror("accel_mult too big (> 1)")
  if (sink_mult < 0.0_DP) call paramerror("sink_mult < 0")
  if (sink_mult > 1.0_DP) call paramerror("sink_mult too big (> 1)")
  if (nlevels < 1) call paramerror("nlevels < 1")
  if (nlevels > 20) call paramerror("nlevels too large (> 20)")
  if (dt_fixed < 0.0_DP) call paramerror("dt_fixed < 0")
#if defined(NBODY_SIMULATION)
  if (nbody_endtime < 0.0_DP) call paramerror("nbody_endtime < 0")
  if (nbody_frac < 0.0_DP .or. nbody_frac > 1.0_DP) &
       & call paramerror("nbody_frac out of range")
  if (npec < 1) call paramerror("npec < 1")
  if (gammapertmax < 0.0_PR .or. gammapertmax > 1.0_PR) &
       &call paramerror("gammapertmax < 0 or gammapertmax > 1")
#endif
#if defined(NBODY_SIMULATION) || defined(NBODY_SPH_SIMULATION)
  if (nbody_timemult < 0.0_DP) call paramerror("nbody_timemult < 0")
  if (nbody_timemult > 1.0_DP) call paramerror("nbody_timemult too big (> 1)")
#endif
#if defined(SPH_SIMULATION) && defined(NBODY_SIMULATION)
  if (nbody_endtime < sph_endtime) call paramerror("nbody_endtime < endtime")
#endif
  if (rscale <= 0.0_DP) call paramerror("rscale <= 0")
  if (mscale <= 0.0_DP) call paramerror("mscale <= 0")
#if defined(PERIODIC_X)
  if (periodic_min(1) >= periodic_max(1)) call paramerror("periodic_x")
#endif
#if defined(PERIODIC_Y)
  if (periodic_min(2) >= periodic_max(2)) call paramerror("periodic_y")
#endif
#if defined(PERIODIC_Z)
  if (periodic_min(3) >= periodic_max(3)) call paramerror("periodic_z")
#endif
#if !defined(H_RHO)
  if (pp_gather < 1) call paramerror("pp_gather < 1")
#endif
#if defined(MINIMUM_H)
  if (hmin < 0.0_PR) call paramerror("hmin < 0")
#endif
#if defined(H_RHO)
  if (h_fac <= 0.0_PR) call paramerror("h_fac < 0")
  if (h_fac >= 2.0_PR) call paramerror("h_fac > 0 (clumping instability)")
#endif
#if defined(HYDRO)
  if (isotemp <= 0.0_PR) call paramerror("isotemp <= 0")
  if (rhobary <= 0.0_PR) call paramerror("rhobary <= 0")
  if (gamma <= 1.0_PR .or. gamma > 5.0_PR/3.0_PR) &
       &call paramerror("gamma <= 1 or gamma > 5/3")
  if (mu_bar <= 0.0_PR) call paramerror("mu_bar <= 0")
  if (Kpoly < 0.0_PR) call paramerror("Kpoly < 0")
#if defined(EXTERNAL_PRESSURE)
  if (Pext < 0.0_PR) call paramerror("Pext < 0")
#endif
#if defined(ARTIFICIAL_VISCOSITY)
  if (alpha < 0.0_PR) call paramerror("alpha < 0")
  if (beta < 0.0_PR) call paramerror("beta < 0")
#if defined(VISC_TD)
  if (alpha_min < 0.0_PR) call paramerror("alpha_min < 0")
  if (alpha_min > alpha) call paramerror("alpha_min > alpha")
#endif
#endif
#endif
#if defined(BH_TREE)
  if (abserror < 0.0_PR) call paramerror("abserror < 0")
  if (abserror > 1.0_PR) call paramerror("abserror > 1")
  if (thetamaxsqd < 0.0_PR) call paramerror("thetamaxsqd < 0")
  if (thetamaxsqd > 1.0_PR) call paramerror("thetamaxsqd too big (> 1)")
  if (nbuildstep <= 0) call paramerror("nbuildstep <= 0")
#endif
#if defined(SINKS)
  if (rhosink <= 0.0_PR) call paramerror("rhosink <= 0")
  if (sinkrad <= 0.0_PR) call paramerror("sinkrad <= 0")
  if (nsearchstep < 1) call paramerror("nsearchstep < 1")
#if !defined(FIXED_ABSOLUTE_SINKRAD)
  if (sinkrad < KERNRANGE) call paramerror("sinkrad < KERNRANGE")
#endif
#if defined(SMOOTH_ACCRETION) || defined(SINK_REMOVE_ANGMOM)
  if (alpha_ss < 0.0_PR) call paramerror("alpha_ss < 0")
#endif
#if defined(SMOOTH_ACCRETION)
  if (smooth_accrete_frac < 0.0_PR .or. smooth_accrete_frac > 1.0_PR) &
       call comperror("smooth_accrete_frac < 0 or > 1")
  if (smooth_accrete_dt < 0.0_PR .or. smooth_accrete_dt > 1.0_PR) &
       call comperror("smooth_accrete_dt < 0 or > 1")
  if (in_file_form /= "sf" .and. in_file_form /= "seren_form" .and. &
       &in_file_form /= "su" .and. in_file_form /= "seren_unform") &
       call paramerror("Cannot use file formats other than native SEREN &
       &format with smooth-accretion sinks")
#endif
#endif
#if defined(REMOVE_OUTLIERS)
  if (rholost < 0.0_PR) call paramerror("rholost < 0.0")
  if (rad_lost < 0.0_PR) call paramerror("rad_lost < 0.0")
#endif
#if defined(IONIZING_UV_RADIATION)
  if (nionallstep <= 0) call paramerror("nionallstep <= 0")
  if (f1 < 0.0_PR) call paramerror("f1 < 0.0")
  if (f2 < 0.0_PR) call paramerror("f1 < 0.0")
  if (f2 > 2.0_PR) call paramerror("f2 > 2.0")
  if (f3 < 0.0_PR) call paramerror("f1 < 0.0")
  if (f4 < 0.0_PR) call paramerror("f1 < 0.0")
  if (Tneut < 0.0_PR) call paramerror("Tneut < 0.0")
  if (Tneut > Tion) call paramerror("Tneut > Tion")
  if (lmax_hp > HP_LEVELS) call paramerror("lmax_hp too big : Too many HEALPix levels")
#endif
#if defined(STELLAR_WIND)
  if (M_loss < 0.0_PR) call paramerror("M_loss < 0.0")
  if (v_wind < 0.0_PR) call paramerror("v_wind < 0.0")
#endif
#if defined(TURBULENT_FORCING)
  if (comp_frac < 0.0_PR) call paramerror("comp_frac < 0.0")
  if (comp_frac > 1.0_PR) call paramerror("comp_frac > 1.0")
  if (turb_T <= 0.0_PR) call paramerror("turb_T <= 0.0")
  if (turb_Ndt <= 0) call paramerror("turb_Ndt <= 0")
  if (mod(turb_T / (turb_Ndt * dt_fixed),1.0) > min(turb_T,dt_fixed)*0.001) then
    call paramwarn("turb_T / turbNdt should be a multiple of dt_fixed")
  end if
  if (turb_rms <= 0.0) call paramerror("turb_rms <= 0.0")
  if (turb_min(1) >= turb_max(1)) call paramerror("turb_x")
  if (turb_min(2) >= turb_max(2)) call paramerror("turb_y")
  if (turb_min(3) >= turb_max(3)) call paramerror("turb_z")
#endif

  return
END SUBROUTINE sanitycheck


! ============================================================================ 
SUBROUTINE comperror(errmsg)
  implicit none

  character(len=*), intent(in) :: errmsg

  write(6,*) "Compiler flag error : ", errmsg
  stop

  return
END SUBROUTINE comperror


! ============================================================================ 
SUBROUTINE compwarning(warnmsg)
  use definitions
  implicit none

  character(len=*), intent(in) :: warnmsg

#if defined(USE_MPI)
  if (rank /= 0) return
#endif
  write(6,*) "Compiler flag warning : ", warnmsg

  return
END SUBROUTINE compwarning


! ============================================================================ 
SUBROUTINE paramerror(errmsg)
  implicit none

  character(len=*), intent(in) :: errmsg

  write(6,*) "Parameter error : ", errmsg
  stop

  return
END SUBROUTINE paramerror


! ============================================================================ 
SUBROUTINE paramwarn(warnmsg)
  implicit none

  character(len=*), intent(in) :: warnmsg

  write(6,*) "Parameter warning : ", warnmsg

  return
END SUBROUTINE paramwarn
