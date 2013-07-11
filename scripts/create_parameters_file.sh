#! /bin/bash -
# ============================================================================
# CREATE_PARAMETERS.SH
# D. A. Hubber - 10/2/2009
# Creates the parameters file based on set values defined for the selected 
# test, and default values for undefined parameters.
# ============================================================================
echo '=====================' > params.dat
echo 'SEREN parameters file' >> params.dat
echo '=====================' >> params.dat
echo 'Contains all input parameters to perform a simulation using Seren.' >> params.dat  
echo '=====================' >> params.dat
echo 'File parameters' >> params.dat
echo '=====================' >> params.dat
echo 'run_id = '${run_id:=TEST} >> params.dat
echo 'run_dir = '${run_dir:=\".\"} >> params.dat
echo 'in_file = '${in_file:=test.dat} >> params.dat
echo 'in_file_form = '${in_file_form:=dragon_form} >> params.dat
echo 'out_file_form = '${out_file_form:=out_file_form} >> params.dat
echo '=====================' >> params.dat
echo 'Misc. parameters' >> params.dat
echo '=====================' >> params.dat
echo 'restart = '${restart:=.FALSE.} >> params.dat
echo 'com_frame = '${com_frame:=.FALSE.} >> params.dat
echo 'rseed =  '${rseed:=1} >> params.dat
echo 'ptrack = '${ptrack:=1} >> params.dat
echo '=====================' >> params.dat
echo 'Time parameters' >> params.dat
echo '=====================' >> params.dat
echo 'sph_endtime = '${sph_endtime:=1.0} >> params.dat
echo 'nbody_sph_endtime = '${nbody_sph_endtime:=1.0} >> params.dat
echo 'nbody_endtime = '${nbody_endtime:=1.0} >> params.dat
echo 'firstsnap = '${firstsnap:=0.1} >> params.dat
echo 'snaptime = '${snaptime:=0.1} >> params.dat
echo 'noutputstep ='${noutputstep:=16} >> params.dat
echo 'ntempstep = '${ntempstep:=2000} >> params.dat
echo 'ndiagstep = '${ndiagstep:=4} >> params.dat
echo 'nsinkstep = '${nsinkstep:=2} >> params.dat
echo 'nsnapstep = '${nsnapstep:=20000} >> params.dat
echo 'courant_mult = '${courant_mult:=0.2} >> params.dat
echo 'accel_mult = '${accel_mult:=0.5} >> params.dat
echo 'sink_mult = '${sink_mult:=0.5} >> params.dat
echo 'nbody_timemult = '${nbody_timemult:=0.1} >> params.dat
echo 'nlevels = '${nlevels:=1} >> params.dat
echo 'dt_fixed = '${dt_fixed:=1.0} >> params.dat
echo '=====================' >> params.dat
echo 'Unit parameters' >> params.dat
echo '=====================' >> params.dat
echo 'dimensionless = '${dimensionless:=0} >> params.dat
echo 'runit = '${runit:=pc} >> params.dat
echo 'munit = '${munit:=m_sun} >> params.dat
echo 'tunit = '${tunit:=myr} >> params.dat
echo 'vunit = '${vunit:=km_s} >> params.dat
echo 'aunit = '${aunit:=km_s2} >> params.dat
echo 'rhounit = '${rhounit:=g_cm3} >> params.dat
echo 'sigmaunit = '${sigmaunit:=g_cm2} >> params.dat
echo 'Punit = '${Punit:=Pa} >> params.dat
echo 'funit = '${funit:=N} >> params.dat
echo 'Eunit = '${Eunit:=J} >> params.dat
echo 'momunit = '${momunit:=m_sunkm_s} >> params.dat
echo 'angmomunit = '${angmomunit:=m_sunau2_yr} >> params.dat
echo 'angvelunit = '${angvelunit:=rad_s} >> params.dat
echo 'dmdtunit = '${dmdtunit:=m_sun_yr} >> params.dat
echo 'Lunit = '${Lunit:=L_sun} >> params.dat
echo 'kappaunit = '${kappaunit:=cm2_g} >> params.dat
echo 'Bunit = '${Bunit:=tesla} >> params.dat
echo 'Qunit = '${Qunit:=C} >> params.dat
echo 'Junit = '${Junit:=C_s_m2} >> params.dat
echo 'uunit = '${uunit:=J_kg} >> params.dat
echo 'tempunit = '${tempunit:=K} >> params.dat
echo 'dudtunit = '${dudtunit:=J_kg_s} >> params.dat
echo 'rscale = '${rscale:=1.0} >> params.dat
echo 'mscale = '${mscale:=1.0} >> params.dat
echo '=====================' >> params.dat
echo 'Periodic parameters' >> params.dat
echo '=====================' >> params.dat
echo 'periodic_min(1) = '${periodic_min_x:=0.0} >> params.dat
echo 'periodic_max(1) = '${periodic_max_x:=1.0} >> params.dat
echo 'periodic_min(2) = '${periodic_min_y:=0.0} >> params.dat
echo 'periodic_max(2) = '${periodic_max_y:=1.0} >> params.dat
echo 'periodic_min(3) = '${periodic_min_z:=0.0} >> params.dat
echo 'periodic_max(3) = '${periodic_max_z:=1.0} >> params.dat
echo 'rspheremax = '${rspheremax:=0.0} >> params.dat
echo 'psphere = '${psphere:=0} >> params.dat
echo '=====================' >> params.dat
echo 'SPH parameters' >> params.dat
echo '=====================' >> params.dat
echo 'pp_gather = '${pp_gather:=50} >> params.dat
echo 'hmin = '${hmin:=0.0} >> params.dat
echo 'h_fac = '${h_fac:=1.2} >> params.dat
echo '=====================' >> params.dat
echo 'Thermal parameters' >> params.dat
echo '=====================' >> params.dat
echo 'boundaryeos = '${boundaryeos:=isothermal} >> params.dat
echo 'icmeos = '${icmeos:=isothermal} >> params.dat
echo 'gaseos = '${gaseos:=isothermal} >> params.dat
echo 'isotemp = '${isotemp:=10.0} >> params.dat
echo 'rhobary = '${rhobary:=1.e-14} >> params.dat
echo 'gamma = '${gamma:=1.4} >> params.dat
echo 'mu_bar = '${mu_bar:=2.35} >> params.dat
echo 'Kpoly = '${Kpoly:=1.0} >> params.dat
echo 'Pext = '${Pext:=0.0} >> params.dat
echo '=====================' >> params.dat
echo 'Cooling parameters' >> params.dat
echo '=====================' >> params.dat
echo 'cooling_law = '${cooling_law:=none} >> params.dat
echo '=====================' >> params.dat
echo 'Viscosity parameters' >> params.dat
echo '=====================' >> params.dat
echo 'alpha = '${alpha:=1.0} >> params.dat
echo 'beta = '${beta:=2.0} >> params.dat
echo 'alpha_min = '${alpha_min:=0.1} >> params.dat
echo '=====================' >> params.dat
echo 'Tree parameters' >> params.dat
echo '=====================' >> params.dat
echo 'abserror = '${abserror:=0.1} >> params.dat
echo 'thetamaxsqd = '${thetamaxsqd:=0.2} >> params.dat
echo 'nbuildstep = '${nbuildstep:=8} >> params.dat
echo '========================' >> params.dat
echo 'Sink creation parameters' >> params.dat
echo '========================' >> params.dat
echo 'rhosink = '${rhosink:=1.0e-11} >> params.dat
echo 'sinkrad = '${sinkrad:=2.0} >> params.dat
echo 'nsearchstep = '${nsearchstep:=16} >> params.dat
echo 'rho_search = '${rho_search:=.FALSE.} >> params.dat
echo 'potmin_search = '${potmin_search:=.FALSE.} >> params.dat
echo 'hill_sphere_search = '${hill_sphere_search:=.FALSE.} >> params.dat
echo 'energy_search = '${energy_search:=.FALSE.} >> params.dat
echo 'thermal_search = '${thermal_search:=.FALSE.} >> params.dat
echo 'div_v_search = '${div_v_search:=.FALSE.} >> params.dat
echo 'div_a_search = '${div_a_search:=.FALSE.} >> params.dat
echo 'timescale_search = '${timescale_search:=.FALSE.} >> params.dat
echo '======================================' >> params.dat
echo 'Sink accretion and feedback parameters' >> params.dat
echo '======================================' >> params.dat
echo 'energy_accrete = '${energy_accrete:=.FALSE.} >> params.dat
echo 'alpha_ss = '${alpha_ss:=0.01} >> params.dat
echo 'smooth_accrete_frac = '${smooth_accrete_frac:=0.005} >> params.dat
echo 'smooth_accrete_dt = '${smooth_accrete_dt:=0.005} >> params.dat
echo 'f_accretion = '${f_accretion:=0.75} >> params.dat
echo 'feedback_tdelay = '${feedback_tdelay:=0.0} >> params.dat
echo 'feedback_minmass = '${feedback_minmass:=10.0} >> params.dat
echo 'star_radius = '${star_radius:=3.0} >> params.dat
echo 'alpha_EA = '${alpha_EA:=0.1} >> params.dat
echo 'dmdt_regular = '${dmdt_regular:=1.0e-7} >> params.dat
echo 'z_factor = '${z_factor:=1.0} >> params.dat
echo '===============================' >> params.dat
echo 'SPH particle removal parameters' >> params.dat
echo '===============================' >> params.dat
echo 'rho_remove = '${rho_remove:=.FALSE.} >> params.dat
echo 'energy_remove = '${energy_remove:=.FALSE.} >> params.dat
echo 'rad_remove = '${rad_remove:=.FALSE.} >> params.dat
echo 'rholost = '${rholost:=0.0} >> params.dat
echo 'rad_lost = '${rad_lost:=0.0} >> params.dat
echo '=================' >> params.dat
echo 'N-body parameters' >> params.dat
echo '=================' >> params.dat
echo 'npec = '${npec:=1} >> params.dat
echo 'nbody_frac = '${nbody_frac:=0.95} >> params.dat
echo 'gammapertmax = '${gammapertmax:=1.0E-4} >> params.dat
echo '=============================' >> params.dat
echo 'Polytropic cooling parameters' >> params.dat
echo '=============================' >> params.dat
echo 'eos_opa_file = '${eos_opa_file:=eos.bell.cc.dat} >> params.dat
echo 'ptemp0 = '${ptemp0:=250.0} >> params.dat
echo 'temp_inf = '${temp_inf:=10.0} >> params.dat
echo 'ptemp_r0 = '${ptemp_r0:=0.25} >> params.dat
echo 'ptemp_q = '${ptemp_q:=0.75} >> params.dat
echo 'fcolumn = '${fcolumn:=0.104} >> params.dat
echo '=====================' >> params.dat
echo 'HEALPix parameters' >> params.dat
echo '=====================' >> params.dat
echo 'nionallstep = '${nionallstep:=8} >> params.dat
echo 'f1 = '${f1:=0.25} >> params.dat
echo 'f2 = '${f2:=1.0} >> params.dat
echo 'f3 = '${f3:=1.0} >> params.dat
echo 'f4 = '${f4:=0.5} >> params.dat
echo 'Tneut = '${Tneut:=10.0} >> params.dat
echo 'Tion = '${Tion:=10000.0} >> params.dat
echo 'Xfrac = '${Xfrac:=0.7} >> params.dat
echo 'mu_ion = '${mu_ion:=0.678} >> params.dat
echo 'a_star = '${a_star:=2.0D-13} >> params.dat
echo 'N_LyC = '${N_LyC:=1.D+48} >> params.dat
echo 'rstatic(1) = '${rstatic1:="0.0"} >> params.dat
echo 'rstatic(2) = '${rstatic2:="0.0"} >> params.dat
echo 'rstatic(3) = '${rstatic3:="0.0"} >> params.dat
echo 'lmax_hp = '${lmax_hp:=7} >> params.dat
echo 'M_loss = '${M_loss:=1.e-6} >> params.dat
echo 'v_wind = '${v_wind:=1000.0} >> params.dat
echo '=====================' >> params.dat
echo 'Turbulence parameters' >> params.dat
echo '=====================' >> params.dat
echo 'comp_frac = '${comp_frac:=0.5} >> params.dat
echo 'turb_T = '${turb_T:=1.0} >> params.dat
echo 'turb_Ndt = '${turb_Ndt:=100} >> params.dat
echo 'turb_rms = '${turb_rms:=1.0} >> params.dat
echo 'turb_min(1) = '${turb_min1:="0.0"} >> params.dat
echo 'turb_max(1) = '${turb_max1:="0.0"} >> params.dat
echo 'turb_min(2) = '${turb_min2:="0.0"} >> params.dat
echo 'turb_max(2) = '${turb_max2:="0.0"} >> params.dat
echo 'turb_min(3) = '${turb_min3:="0.0"} >> params.dat
echo 'turb_max(3) = '${turb_max3:="0.0"} >> params.dat
