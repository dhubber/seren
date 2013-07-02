! DEFAULT_PARAMETERS.F90
! D. A. Hubber & K. Rawiraswattana - 22/08/2010
! Set default values for parameters in case parameters file is not read in.
! (Mainly used for initial conditions and test programs where certain 
! unused parameters still need to tbe set to some value).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE default_parameters
  use interface_module
  use periodic_module
  use scaling_module
  use filename_module
  use time_module
  use sink_module
  use Tprof_module
  use particle_module
  use type_module
  use neighbour_module, only : pp_gather,hmin,h_fac
  use hydro_module, only : alpha,beta,alpha_min,cooling_law,&
       &gamma,isotemp,Kpoly,mu_bar,mu_ion,rhobary,Pext
  use Nbody_module, only : gammapertmax,nbody_endtime,&
       &nbody_frac,nbody_timemult,npec
  use tree_module, only : thetamaxsqd, abserror
  use HP_module, only : a_star,f1,f2,f3,f4,lmax_hp,M_loss,N_Lyc,&
       &rstatic,Tion,Tneut,v_wind,Xfrac
  use turbulence_module
  implicit none

  debug1("Setting default parameter values [default_parameters.F90]")

  nparams = 0

! File variables
  call add_string_parameter("run_dir",run_dir,".")
  call add_string_parameter("run_id",run_id,"test1")
  call add_string_parameter("in_file",in_file,"")
  call add_string_parameter("in_file_form",in_file_form,"sf")
  call add_string_parameter("out_file_form",out_file_form,"sf")

! Misc. variables  
  call add_logical_parameter("restart",restart,.false.)
  call add_logical_parameter("com_frame",com_frame,.false.)
  call add_integer_parameter("rseed",rseed,1)
  call add_integer_parameter("ptrack",ptrack,1)
  
! Time variables
  call add_double_parameter("sph_endtime",sph_endtime,1.0_DP)
  call add_double_parameter("nbody_endtime",nbody_endtime,1.0_DP)
  call add_double_parameter("nbody_sph_endtime",nbody_sph_endtime,1.0_DP)
  call add_double_parameter("firstsnap",firstsnap,0.0_DP)
  call add_double_parameter("snaptime",snaptime,1.0_DP)
  call add_long_integer_parameter("noutputstep",noutputstep,100_ILP)
  call add_long_integer_parameter("ntempstep",ntempstep,1000_ILP)
  call add_long_integer_parameter("ndiagstep",ndiagstep,32_ILP)
  call add_long_integer_parameter("nsinkstep",nsinkstep,32_ILP)
  call add_long_integer_parameter("nsnapstep",nsnapstep,10000_ILP)
  call add_double_parameter("courant_mult",courant_mult,0.2_DP)
  call add_double_parameter("accel_mult",accel_mult,0.5_DP)
  call add_double_parameter("sink_mult",sink_mult,0.5_DP)
  call add_double_parameter("nbody_timemult",nbody_timemult,0.95_DP)
  call add_long_integer_parameter("nlevels",nlevels,1_ILP)
  call add_double_parameter("dt_fixed",dt_fixed,1.0_DP)

! Scaling preferences
  call add_unit_parameter("runit",runit,"pc")
  call add_unit_parameter("munit",munit,"m_sun")
  call add_unit_parameter("tunit",tunit,"myr")
  call add_unit_parameter("vunit",vunit,"km_s")
  call add_unit_parameter("aunit",aunit,"km_s2")
  call add_unit_parameter("rhounit",rhounit,"g_cm3")
  call add_unit_parameter("sigmaunit",sigmaunit,"g_cm2")
  call add_unit_parameter("Punit",Punit,"Pa")
  call add_unit_parameter("funit",funit,"N")
  call add_unit_parameter("Eunit",Eunit,"J")
  call add_unit_parameter("momunit",momunit,"m_sunkm_s")
  call add_unit_parameter("angmomunit",angmomunit,"m_sunau2_yr")
  call add_unit_parameter("angvelunit",angvelunit,"rad_s")
  call add_unit_parameter("dmdtunit",dmdtunit,"m_sun_yr")
  call add_unit_parameter("Lunit",Lunit,"L_sun")
  call add_unit_parameter("kappaunit",kappaunit,"cm2_g")
  call add_unit_parameter("Bunit",Bunit,"tesla")
  call add_unit_parameter("Qunit",Qunit,"C")
  call add_unit_parameter("Junit",Junit,"C_s_m2")
  call add_unit_parameter("uunit",uunit,"J_kg")
  call add_unit_parameter("tempunit",tempunit,"K")
  call add_unit_parameter("dudtunit",dudtunit,"J_kg_s")
  call add_double_parameter("rscale",rscale,1.0_DP)
  call add_double_parameter("mscale",mscale,1.0_DP)

! Periodic boundary conditions (PERIODIC_X/Y/Z)
  call add_real_parameter("periodic_min(1)",periodic_min(1),0.0_PR)
  call add_real_parameter("periodic_max(1)",periodic_max(1),1.0_PR)
  call add_real_parameter("periodic_min(2)",periodic_min(2),0.0_PR)
  call add_real_parameter("periodic_max(2)",periodic_max(2),1.0_PR)
  call add_real_parameter("periodic_min(3)",periodic_min(3),0.0_PR)
  call add_real_parameter("periodic_max(3)",periodic_max(3),1.0_PR)
  call add_real_parameter("rspheremax",rspheremax,1.0_PR)
  call add_integer_parameter("psphere",psphere,0)

! Neighbour and smoothing length parameters
  call add_integer_parameter("pp_gather",pp_gather,50)
  call add_real_parameter("hmin",hmin,0.0_PR)
  call add_real_parameter("h_fac",h_fac,1.2_PR)

! Equation of state parameters
  call add_string_parameter("boundaryeos",&
       &typeinfo(boundaryid)%eos,"isothermal")
  call add_string_parameter("icmeos",typeinfo(icmid)%eos,"isothermal")
  call add_string_parameter("gaseos",typeinfo(gasid)%eos,"isothermal")
  call add_real_parameter("isotemp",isotemp,10.0_PR)
  call add_real_parameter("rhobary",rhobary,1.0E-14_PR)
  call add_real_parameter("gamma",gamma,1.4_PR)
  call add_real_parameter("mu_bar",mu_bar,2.35_PR)
  call add_real_parameter("Kpoly",Kpoly,0.4246_PR)
  call add_real_parameter("Pext",Pext,0.0_PR)

! Cooling law parameters
  call add_string_parameter("cooling_law",cooling_law,"none")

! Viscosity factors: alpha, beta and alpha_min
  call add_real_parameter("alpha",alpha,1.0_PR)
  call add_real_parameter("beta",beta,2.0_PR)
  call add_real_parameter("alpha_min",alpha_min,0.1_PR)

! Tree parameters
  call add_real_parameter("abserror",abserror,0.1_PR)
  call add_real_parameter("thetamaxsqd",thetamaxsqd,0.5_PR)
  call add_long_integer_parameter("nbuildstep",nbuildstep,8_ILP)

! Sink particle creation parameters
  call add_real_parameter("rhosink",rhosink,1.0E-10_PR)
  call add_real_parameter("sinkrad",sinkrad,3.0_PR)
  call add_long_integer_parameter("nsearchstep",nsearchstep,16_ILP)
  call add_logical_parameter("rho_search",rho_search,.true.)
  call add_logical_parameter("potmin_search",potmin_search,.true.)
  call add_logical_parameter("hill_sphere_search",hill_sphere_search,.true.)
  call add_logical_parameter("energy_search",energy_search,.false.)
  call add_logical_parameter("thermal_search",thermal_search,.false.)
  call add_logical_parameter("div_v_search",div_v_search,.true.)
  call add_logical_parameter("div_a_search",div_a_search,.true.)
  call add_logical_parameter("timescale_search",timescale_search,.false.)

! Sink accretion and feedback parameters parameters
  call add_logical_parameter("energy_accrete",energy_accrete,.false.)
  call add_real_parameter("alpha_ss",alpha_ss,0.01_PR)
  call add_double_parameter("smooth_accrete_frac",&
       &smooth_accrete_frac,5.0E-3_DP)
  call add_double_parameter("smooth_accrete_dt",smooth_accrete_dt,5.0E-3_DP)
  call add_double_parameter("smooth_accrete_dt",smooth_accrete_rho,10.0_DP)
  call add_real_parameter("f_accretion",f_accretion,0.75_PR)
  call add_real_parameter("feedback_tdelay",feedback_tdelay,0.0_PR)
  call add_real_parameter("feedback_minmass",feedback_minmass,10.0_PR)
  call add_real_parameter("star_radius",star_radius,3.0_PR)
  call add_real_parameter("alpha_EA",alpha_EA,0.1_PR)
  call add_double_parameter("dmdt_regular",dmdt_regular,1.0E-07_DP)

! SPH particle removal criteria
  call add_logical_parameter("rho_remove",rho_remove,.false.)
  call add_logical_parameter("energy_remove",energy_remove,.false.)
  call add_logical_parameter("rad_remove",rad_remove,.false.)
  call add_real_parameter("rholost",rholost,1.0E-21_PR)
  call add_real_parameter("rad_lost",rad_lost,1.0E10_PR)

! N-body simulation parameters
  call add_integer_parameter("npec",npec,1)
  call add_double_parameter("nbody_frac",nbody_frac,0.95_DP)
  call add_double_parameter("gammapertmax",gammapertmax,1.0E-6_DP)

! Radiative cooling (Stamatellos et al. 2007) parameters
  call add_string_parameter("eos_opa_file",eos_opa_file,"eos.dat")
  call add_real_parameter("z_factor",z_factor,1.0_PR)
  call add_real_parameter("ptemp0",ptemp0,250.0_PR)
  call add_real_parameter("temp_inf",temp_inf,10.0_PR)
  call add_real_parameter("ptemp_r0",ptemp_r0,0.25_PR)
  call add_real_parameter("ptemp_q",ptemp_q,0.75_PR)

! HEALPix ionizing radiation and stellar wind parameters
  call add_long_integer_parameter("nionallstep",nionallstep,8_ILP)
  call add_real_parameter("f1",f1,0.5_PR)
  call add_real_parameter("f2",f2,1.0_PR)
  call add_real_parameter("f3",f3,1.0_PR)
  call add_real_parameter("f4",f4,0.8_PR)
  call add_real_parameter("Tneut",Tneut,10.0_PR)
  call add_real_parameter("Tion",Tion,10000.0_PR)
  call add_real_parameter("Xfrac",Xfrac,0.7_PR)
  call add_real_parameter("mu_ion",mu_ion,1.0_PR)
  call add_double_parameter("a_star",a_star,2.7E-13_DP)
  call add_double_parameter("N_LyC",N_LyC,1.0E49_DP)
  call add_real_parameter("rstatic(1)",rstatic(1),0.0_PR)
  call add_real_parameter("rstatic(2)",rstatic(2),0.0_PR)
  call add_real_parameter("rstatic(3)",rstatic(3),0.0_PR)
  call add_integer_parameter("lmax_hp",lmax_hp,11)
  call add_real_parameter("M_loss",M_loss,1.0E-8_PR)
  call add_real_parameter("v_wind",v_wind,2000.0_PR)

! Turbulent forcing parameters
  call add_real_parameter("comp_frac",comp_frac, 0.5_PR)
  call add_double_parameter("turb_T", turb_T, 0.1_DP)
  call add_integer_parameter("turb_Ndt", turb_Ndt, 100)
  call add_real_parameter("turb_rms", turb_rms, 1.0_PR)
  call add_real_parameter("turb_min(1)",turb_min(1),0.0_PR)
  call add_real_parameter("turb_max(1)",turb_max(1),1.0_PR)
  call add_real_parameter("turb_min(2)",turb_min(2),0.0_PR)
  call add_real_parameter("turb_max(2)",turb_max(2),1.0_PR)
  call add_real_parameter("turb_min(3)",turb_min(3),0.0_PR)
  call add_real_parameter("turb_max(3)",turb_max(3),1.0_PR)

  return
END SUBROUTINE default_parameters



! ============================================================================
! ADD_INTEGER_PARAMETER
! D. A. Hubber - 21/08/2010
! Add integer variable parameter to input parameter list
! ============================================================================
SUBROUTINE add_integer_parameter(param_name,ipointer,idefault)
  use definitions
  use filename_module
  implicit none

  character(len=*), intent(in) :: param_name     ! name of parameter
  integer, target, intent(in) :: ipointer        ! integer variable pointer
  integer, intent(in) :: idefault                ! default value of parameter

  nparams = nparams + 1
  params(nparams)%var_name = param_name
  params(nparams)%var_type = "i"
  params(nparams)%var_i    => ipointer
  params(nparams)%var_i    = idefault

  return
END SUBROUTINE add_integer_parameter



! ============================================================================
! ADD_LONG_INTEGER_PARAMETER
! D. A. Hubber - 21/08/2010
! Add long-integer variable parameter to input parameter list
! ============================================================================
SUBROUTINE add_long_integer_parameter(param_name,ipointer,idefault)
  use definitions
  use filename_module
  implicit none

  character(len=*), intent(in) :: param_name        ! name of parameter
  integer(kind=ILP), target, intent(in) :: ipointer ! integer variable pointer
  integer(kind=ILP), intent(in) :: idefault         ! default value of parameter

  nparams = nparams + 1
  params(nparams)%var_name = param_name
  params(nparams)%var_type = "j"
  params(nparams)%var_j    => ipointer
  params(nparams)%var_j    = idefault

  return
END SUBROUTINE add_long_integer_parameter



! ============================================================================
! ADD_REAL_PARAMETER
! D. A. Hubber - 21/08/2010
! Add real (PR) variable parameter to input parameter list
! ============================================================================
SUBROUTINE add_real_parameter(param_name,rpointer,rdefault)
  use definitions
  use filename_module
  implicit none

  character(len=*), intent(in) :: param_name     ! name of parameter
  real(kind=PR), target, intent(in) :: rpointer  ! real variable pointer
  real(kind=PR), intent(in) :: rdefault          ! default value of parameter

  nparams = nparams + 1
  params(nparams)%var_name = param_name
  params(nparams)%var_type = "r"
  params(nparams)%var_r    => rpointer
  params(nparams)%var_r    = rdefault

  return
END SUBROUTINE add_real_parameter



! ============================================================================
! ADD_DOUBLE_PARAMETER
! D. A. Hubber - 21/08/2010
! Add double-precision real variable parameter to input parameter list
! ============================================================================
SUBROUTINE add_double_parameter(param_name,dpointer,ddefault)
  use definitions
  use filename_module
  implicit none

  character(len=*), intent(in) :: param_name     ! name of parameter
  real(kind=DP), target, intent(in) :: dpointer  ! DP variable pointer
  real(kind=DP), intent(in) :: ddefault          ! default value of parameter

  nparams = nparams + 1
  params(nparams)%var_name = param_name
  params(nparams)%var_type = "d"
  params(nparams)%var_d    => dpointer
  params(nparams)%var_d    = ddefault

  return
END SUBROUTINE add_double_parameter



! ============================================================================
! ADD_STRING_PARAMETER
! D. A. Hubber - 21/08/2010
! Add string variable parameter to input parameter list
! ============================================================================
SUBROUTINE add_string_parameter(param_name,cpointer,cdefault)
  use definitions
  use filename_module
  implicit none

  character(len=*), intent(in) :: param_name        ! name of parameter
  character(len=*), target, intent(in) :: cpointer  ! DP variable pointer
  character(len=*), intent(in) :: cdefault          ! default value of string

  nparams = nparams + 1
  params(nparams)%var_name = param_name
  params(nparams)%var_type = "c"
  params(nparams)%var_c    => cpointer
  params(nparams)%var_c    = cdefault

  return
END SUBROUTINE add_string_parameter



! ============================================================================
! ADD_UNIT_PARAMETER
! D. A. Hubber - 21/08/2010
! Add unit-string variable parameter to input parameter list
! ============================================================================
SUBROUTINE add_unit_parameter(param_name,upointer,udefault)
  use definitions
  use filename_module
  implicit none

  character(len=*), intent(in) :: param_name        ! name of parameter
  character(len=*), target, intent(in) :: upointer  ! DP variable pointer
  character(len=*), intent(in) :: udefault          ! default value of string

  nparams = nparams + 1
  params(nparams)%var_name = param_name
  params(nparams)%var_type = "u"
  params(nparams)%var_u    => upointer
  params(nparams)%var_u    = udefault

  return
END SUBROUTINE add_unit_parameter



! ============================================================================
! ADD_LOGICAL_PARAMETER
! D. A. Hubber - 21/08/2010
! Add logical variable parameter to input parameter list
! ============================================================================
SUBROUTINE add_logical_parameter(param_name,lpointer,ldefault)
  use definitions
  use filename_module
  implicit none

  character(len=*), intent(in) :: param_name        ! name of parameter
  logical, target, intent(in) :: lpointer           ! logical variable pointer
  logical, intent(in) :: ldefault                   ! default value of param

  nparams = nparams + 1
  params(nparams)%var_name = param_name
  params(nparams)%var_type = "l"
  params(nparams)%var_l    => lpointer
  params(nparams)%var_l    = ldefault

  return
END SUBROUTINE add_logical_parameter
