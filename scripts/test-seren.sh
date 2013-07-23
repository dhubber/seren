#! /bin/bash -
# ============================================================================
# TEST-SEREN.SH
# D. A. Hubber - 10/2/2009
#
# Automated script to run a complete suite of tests (or just individual tests) 
# of the Seren SPH code.
#
# Usage:
# test-seren.sh [-comp]                  : Fortran compiler
#               [-dp]                    : double precision
#               [-openmp]                : parallelise with OpenMP
#               [-mpi]                   : parallelise with MPI
#               [-n NPROC]               : Set no. of MPI processes
#               [-fast]                  : fast compiler flags
#               [-debugX]                : set DEBUG mode to X
#               [-clean]                 : Clean-up files once sim. has ended
#               [-all]                   : perform all tests
#               [-listX]                 : perform test suite X
#               [-test test1 test2 ..]   : perform specified tests 1, 2 ..
#
# Perform all tests in suite             : test-seren.sh -all
# Perform pre-selected list of tests, X  : test-seren.sh -listX
# Perform tests 'test1', 'test2' etc..   : test-seren.sh -test test1 test2 ..
#
# List of current tests
# ADSOD-3D-AB                 : Adiabatic Sod test with standard AV.
# ADSOD-3D-AB-COND            : As ADSOD-3D-AB but with art. conductivity.
# ADSOD-3D-GRADH-AB-COND      : As ADSOD-3D-AB-COND but with 'grad-h' SPH.
# ADSOD-3D-GRADH-MON-COND     : As ADSOD-3D-GRADH-AB-COND, but with MON97 AV.
# BURRAU1                     : Burrau 3-body test
# COL-3D-AB                   : Colliding flows in 3D with AB art. visc.
# COL-3D-AB-TD                : As COL-3D-AB but with time-dependent visc.
# COL-3D-AB-BAL               : As COL-3D-AB but with the Balsara switch.
# COL-3D-MON                  : Colliding flows in 3D with AB art. visc.
# EIGEN1-MONO-KSGRAV          : BH Eigen MAC test with KS gravity (4000).
# EIGEN1-QUAD-KSGRAV          : As above, but with quadrupole moment terms.
# EIGEN1-OCT-KSGRAV           : As above, but with octupole moment terms.
# EIGEN1-MONO-NBODY           : BH Eigen MAC test with N-body gravity.
# EIGEN1-QUAD-NBODY           : As above, but with quadrupole moment terms.
# EIGEN1-OCT-NBODY            : As above, but with octupole moment terms.
# FIGURE8                     : Figure-8 3-body test
# FREEFALL1-GRADH-BH          : Freefall collapse test with grad-h gravity
# FREEFALL1-NBODY-BH          : As above, but with no kernel-softening.
# GEO1-MONO-KSGRAV            : BH gemoetric MAC test with KS gravity (4000).
# GEO1-QUAD-KSGRAV            : As above, but with quadrupole moment terms.
# GEO1-OCT-KSGRAV             : As above, but with octupole moment terms.
# GEO1-MONO-NBODY             : BH geometric MAC test with N-body gravity.
# GEO1-QUAD-NBODY             : As above, but with quadrupole moment terms.
# GEO1-OCT-NBODY              : As above, but with octupole moment terms.
# ISOFREEFALL1-GRADH-BH       : As FREEFALL1-GRADH-BH, but with isothermal gas
# NTSI1-2D-MON-COND-CONSTH    : Non-linear thin shell instability
# POLYRAD1-AB                 : Polytropic-cooling collapse test with AB AV.
# POLYRAD1-AB-FLD             : As above, but with flux-limited diffusion
# SEDOV1-3D-GRADH             : Sedov blast wave test with global timesteps.
# SEDOV2-3D-GRADH             : As above, but with individual timesteps.
# SEDOV3-3D-GRADH             : As above, but with neighbour-checking.
# SHEAR-2D-GRADH-AB           : 2D shear-flow
# SIT1-AB-BH-SINK             : Boss-Bodenheimer test (standard AV, sinks).
# SIT1-GRADH-AB-BH-SINK       : As above, but with 'grad-h' SPH.
# SIT1-AB-BH-SMOOTH_SINK      : ..
# STATPOLY1-AB-CONSTH         : Relax polytrope to hydrostatic balance.
# STATPOLY1-AB-GRADH          : As above but with grad-h SPH (1000 particles)
# STATPOLY2-AB-GRADH          : As above but with 100,000 particles
# STATPOLY3-AB-GRADH          : As above but with 114 particles
# ============================================================================


# Set-up lists of tests
LISTADSOD="ADSOD-3D-AB ADSOD-3D-AB-COND ADSOD-3D-GRADH-AB-COND ADSOD-3D-GRADH-MON-COND"
LISTCOL="COL-3D-AB COL-3D-AB-BAL COL-3D-AB-TD COL-3D-MON"
LISTEIGEN1SPH="EIGEN1-QUAD-KSGRAV EIGEN1-OCT-KSGRAV"
LISTEIGEN1NBODY="EIGEN1-QUAD-NBODY EIGEN1-OCT-NBODY"
LISTFREEFALL="FREEFALL1-GRADH-BH FREEFALL-NBODY-BH"
LISTGEO1SPH="GEO1-MONO-KSGRAV GEO1-QUAD-KSGRAV GEO1-OCT-KSGRAV"
LISTGEO1NBODY="GEO1-MONO-NBODY GEO1-QUAD-NBODY GEO1-OCT-NBODY"
LISTGRAV="$LISTEIGEN1SPH $LISTEIGEN1NBODY $LISTGEO1SPH $LISTGEO1NBODY"
LISTKH="KH-2D-GRADH-COND"
LISTNTSI="NTSI1-2D-MON-COND-CONSTH"
LISTPOLYRAD1="POLYRAD1-AB"
LISTSEDOV="SEDOV1-3D-GRADH SEDOV2-3D-GRADH SEDOV3-3D-GRADH"
LISTSHEAR2D="SHEAR-2D-GRADH-AB SHEAR-2D-GRADH-AB-TD"
LISTSTATPOLY="STATPOLY1-AB-CONSTH STATPOLY1-AB-GRADH STATPOLY2-AB-GRADH STATPOLY3-AB-GRADH"
ALL_TESTS="$LISTADSOD $LISTCOL $LISTGRAV $LISTPOLYRAD1 $LISTSEODV $LISTSHEAR2D $LISTSTATPOLY"
LIST1="STATPOLY1-AB-GRADH ADSOD-3D-GRADH-MON-COND COL-3D-AB-TD SEDOV3-3D-GRADH NTSI1-2D-MON-COND-CONSTH KH-2D-GRADH-COND BURRAU1 FIGURE8 SIT1-GRADH-AB-BH-SINK POLYRAD1-AB"
LISTPAPER="ADSOD-3D-GRADH-MON-COND COL-3D-AB SEDOV1-3D-GRADH SEDOV2-3D-GRADH SEDOV3-3D-GRADH FREEFALL1-GRADH-BH STATPOLY1-AB-GRADH STATPOLY2-AB-GRADH STATPOLY3-AB-GRADH FIGURE8 BURRAU""LISTGEO1SPH"
#ALL_TESTS="COL-3D-AB COL-3D-AB-BAL COL-3D-AB-TD COL-3D-MON97 ADSOD-3D-AB ADSOD-3D-AB-COND ADSOD-3D-GRADH-AB-COND SIT1-AB-SINK POLYRAD1-AB GEO1-MONO-KSGRAV GEO1-QUAD-KSGRAV GEO1-OCT-KSGRAV GEO1-MONO-NBODY GEO1-QUAD-NBODY GEO1-OCT-NBODY NTSI1-2D-MON-COND-CONSTH STATPOLY1-AB-CONST"

echo '-----------------'
echo 'Seren test script'
echo '-----------------'
echo 'No. of arguments : '$#

# Unset key variables for safety
unset F90 PRECISION OPENMP DEBUG TESTLIST CLEAN


# Process arguments for script
# ----------------------------------------------------------------------------
if test $# -gt 0
then
    while [ $# -gt 0 ]
      do 
      case $1 in
          -gfortran)
             F90=gfortran
             ;;
          -f90)
             F90=f90
             ;;
          -ifort)
             F90=ifort
             ;;
          -g95)
             F90=g95
             ;;
          -f95)
             F90=f95
             ;;
          -dp)
             PRECISION=DOUBLE
             ;;
          -openmp)
             OPENMP=1
             ;;
          -mpi)
             MPI=1
             MPIF90=mpif90
             GHOST_PARTICLES=1
             ;;
          -n1)
             NPROC=1
             ;;
          -n2)
             NPROC=2
             ;;
          -n4)
             NPROC=4
             ;;
          -n8)
             NPROC=8
             ;;
          -fast)
             COMPILER_MODE=FAST
             ;;
          -debug)
             COMPILER_MODE=DEBUG
             ;;
          -debug0)
             DEBUG=0
             ;;
          -debug1)
             DEBUG=1
             ;;
          -debug2)
             DEBUG=2
             ;;
          -debug3)
             DEBUG=3
             ;;
          -clean)
             CLEAN=1
             ;;
	  -all)
             TESTLIST=$ALL_TESTS
	     echo "Running all known tests: "$TESTLIST
	     break;
	     ;;
	  -list1)
             TESTLIST=$LIST1
	     echo "Running test list 1: "$TESTLIST
	     break;
	     ;;
	  -listPAPER)
             TESTLIST=$LISTPAPER
	     echo "Running test list 1: "$TESTLIST
	     break;
	     ;;
	  -listADSOD)
             TESTLIST=$LISTADSOD
	     echo "Running adiabatic Sod test list: "$TESTLIST
	     break;
	     ;;
          -listCOL)
             TESTLIST=$LISTCOL
             echo "Running colliding flows test list: "$TESTLIST
             break;
             ;;
          -listRAD)
             TESTLIST=$LISTRAD
             echo "Running Polytropic cooling test list: "$TESTLIST
             break;
             ;;
          -listEIGEN1NBODY)
             TESTLIST=$LISTEIGEN1NBODY
             echo "Running tree grav Eigen MAC N-body test list: "$TESTLIST
             break;
             ;;
          -listEIGEN1SPH)
             TESTLIST=$LISTEIGEN1SPH
             echo "Running tree grav Eigen MAC KS test list: "$TESTLIST
             break;
             ;;
          -listFREEFALL)
             TESTLIST=$LISTFREEFALL
             echo "Running Freefall collapse test list: "$TESTLIST
             break;
             ;;
          -listGEO1SPH)
             TESTLIST=$LISTGEO1SPH
             echo "Running tree grav geometric MAC KS test list: "$TESTLIST
             break;
             ;;
          -listGEO1NBODY)
             TESTLIST=$LISTGEO1NBODY
             echo "Running tree grav geometric MAC N-body test list: "$TESTLIST
             break;
             ;;
          -listGRAV)
             TESTLIST=$LISTGRAV
             echo "Running tree gravity test list: "$TESTLIST
             break;
             ;;
          -listSEDOV)
             TESTLIST=$LISTSEDOV
             echo "Running Sedov blast wave test list: "$TESTLIST
             break;
             ;;
          -listSHEAR2D)
             TESTLIST=$LISTSHEAR2D
             echo "Running Shear-flow test list: "$TESTLIST
             break;
             ;;
	  -test)
	     shift
	     if test $# -gt 0
		 then
		 TESTLIST="$@"        
		 echo "Running selected test via argument list: "$TESTLIST
	     else
		 exit 1
	     fi
	     break
	     ;;
          -*)
	     error "Unrecognised option: $1"
	     exit 0
	     ;;
	  *) 
	  exit 1
	  ;;
    esac
    shift
    done
else
    echo 'No arguments - Exiting script'
    exit 1
fi



# Initialise variables for test.
# ----------------------------------------------------------------------------
cd ..
SEREN_DIR=$(pwd)
SCRIPT_DIR=$SEREN_DIR/scripts
#IC_DIR=$SCRIPT_DIR/IC
IC_DIR=./icfiles
TESTFILE_DIR=$SCRIPT_DIR/test-files
PLOTFILE_DIR=$SCRIPT_DIR/plot-files
cd $SCRIPT_DIR

if [ ! -e results ]
then
   mkdir results
fi

echo 'SCRIPT_DIR : '$SCRIPT_DIR
echo 'SEREN_DIR : '$SEREN_DIR



# Now loop over all tests
# ============================================================================
for i in $TESTLIST
do

  # Read in Makefile and parameter options for current test.
  # --------------------------------------------------------------------------
  cd $SCRIPT_DIR
  testname=$i
  TEST_DIR=$SCRIPT_DIR/$testname
  echo 'Searching for '$testname

  if test -e $TESTFILE_DIR/$testname.test
  then
      rm -rf $TEST_DIR
      mkdir $TEST_DIR
      echo 'Reading '$testname.test
      source $TESTFILE_DIR/$testname.test
      
      
  # Construct Makefile for current test.  Set default values for variables 
  # that don't yet exist (i.e. for options not required by current test).
  # Finally, add 'tail' containing the rest of the Makefile.
      cd $SCRIPT_DIR
      rm -f Makefile
      source ./create_makefile_header.sh
      echo '' >> Makefile
      echo 'DFLAGS += '${DFLAGS:="-DDEBUG_DIAGNOSTICS -DDEBUG_PLOT_DATA -DDEBUG_TRACK_ENERGY"} >> Makefile
#  echo 'DFLAGS += '${DFLAGS:="-DDEBUG_DIAGNOSTICS -DDEBUG_PLOT_DATA -DTRACK_ENERGY -Wall -Wconversion -ffpe-trap=invalid,zero,overflow,underflow,denormal"} >> Makefile
      echo 'CFLAGS += $(DFLAGS)' >> Makefile
      echo '' >> Makefile
      cat $SEREN_DIR/'makefiletail.mk' >> Makefile
      
      
  # Next, construct 'params.dat' file for current test.
      cd $SCRIPT_DIR
      rm -f params.dat
      source ./create_parameters_file.sh
      
      
  # Now place the Makefile, parameters file, initial conditions file and 
  # any other necessary files inside the test directory.
      mv Makefile $TEST_DIR/.
      mv params.dat $TEST_DIR/.
      cp $IC_DIR/$in_file.gz $TEST_DIR/.
      gunzip $TEST_DIR/$in_file.gz
      
      
  # Now make seren (or other program) in test directory and clean up any 
  # '.o' and '.mod' files.
      echo 'Compiling SEREN'
      cd $TEST_DIR
      make -s -j 12 $PROG_NAME
      rm -f *.o *.mod
      
      
  # Run seren in test directory.
      #if "$MPI" -e "1"
      #  then
      #mpirun -n $NPROC ./seren-mpi
      #else
        time ./$PROG_EX
      #fi
      
      
  # Prepare figures for current test.
  # --------------------------------------------------------------------------
      cd $PLOTFILE_DIR
      if test -e "$PLOT_FILE"
	  then
	  cp $PLOT_FILE $TEST_DIR/.
	  if test -e "$ANALYTIC"
	      then
	      cp $PLOTFILE_DIR/$ANALYTIC $TEST_DIR/.
	  fi
	  cd $TEST_DIR
	  if test "$PLOT_PROG" == "gnuplot"
	      then
	      gnuplot < $PLOT_FILE
	  fi
	  mv $testname.ps $SCRIPT_DIR/results/.
      fi
      cd $TEST_DIR
  
  # Clean-up files if required
      if [ "$CLEAN" == "1" ]
	  then
	  rm -rf $TEST_DIR
      fi


  else
      echo $testname 'does not exist'
  fi


  # Now unset (almost) all variables in preparation for next test.
  # --------------------------------------------------------------------------
  echo 'Unsetting all bash variables for next test'
  unset -v SRCDIR EXEDIR OPTIMISE COMPILER_MODE OUTPUT_LEVEL NDIM 
  unset -v INFILE_FORMAT OUTFILE_FORMAT PERIODIC X_BOUNDARY Y_BOUNDARY
  unset -v Z_BOUNDARY SPHERICAL_WALL CYLINDRICAL_WALL
  unset -v SPH_SIMULATION NBODY_SPH_SIMULATION NBODY_SIMULATION SPH 
  unset -v SPH_INTEGRATION KERNEL HFIND MINIMUM_H HYDRO ENERGY_EQN ENTROPY_EQN
  unset -v ARTIFICIAL_VISCOSITY BALSARA VISC_TD PATTERN_REC 
  unset -v ARTIFICIAL_CONDUCTIVITY EXTERNAL_PRESSURE
  unset -v RAD_WS SINK_POTENTIAL_WS AMBIENT_HEATING_WS SINK_HEATING_WS
  unset -v FLUX_LIMITED_DIFFUSION COOLING_HEATING IONIZING_RADIATION
  unset -v STELLAR_WIND PARTICLE_INJECTION_WINDS STELLAR_LUMINOSITY
  unset -v EXTERNAL_FORCE SELF_GRAVITY MEAN_H_GRAVITY EWALD 
  unset -v SINKS SINK_RADIUS SINK_REMOVE_ANGMOM SINK_GRAVITY_ONLY 
  unset -v NBODY_INTEGRATION BINARY_STATS BINARY_COM_MOTION FORCE_SPLITTING
  unset -v TREE MULTIPOLE MAC REORDER CELL_WALK SORT TIMESTEP
  unset -v CHECK_NEIB_TIMESTEP SIGNAL_VELOCITY_DT NEIGHBOURLISTS
  unset -v KERNEL_TABLES REMOVE_OUTLIERS TIMING_CODE DIMENSIONLESS TEST
  unset -v run_id run_dir in_file_form out_file_form 
  unset -v restart com_frame rseed ptrack 
  unset -v sph_endtime nbody_sph_endtime nbody_endtime firstsnap snaptime 
  unset -v noutputstep ntempstep ndiagstep nsinkstep nsnapstep
  unset -v courant_mult accel_mult sink_mult nbody_timemult nlevels dt_fixed 
  unset -v runit munit tunit vunit aunit rhounit sigmaunit Punit funit Eunit 
  unset -v momunit angmomunit angvelunit dmdtunit Lunit kappaunit Bunit Qunit 
  unset -v Junit uunit tempunit dudtunit rscale mscale 
  unset -v periodic_min_x periodic_max_x periodic_min_y periodic_max_y 
  unset -v periodic_min_z periodic_max_z rspheremax psphere
  unset -v pp_gather hmin h_fac boundaryeos icmeos gaseos
  unset -v isotemp rhobary gamma mu_bar Kpoly Pext cooling_law
  unset -v alpha beta alpha_min 
  unset -v abserror thetamaxsqd nbuildstep 
  unset -v rhosink sinkrad nsearchstep rho_search potmin_search 
  unset -v hill_sphere_search energy_search div_v_search div_a_search 
  unset -v timescale_search energy_accrete alpha_ss smooth_accrete_frac
  unset -v smooth_accrete_dt f_accretion feedback_tdelay feedback_minmass
  unset -v star_radius alpha_EA dmdt_regular z_factor
  unset -v rho_remove energy_remove rad_remove rholost rad_lost
  unset -v npec nbody_frac gammapertmax
  unset -v eos_opa_file ptemp0 temp_inf ptemp_r0 ptemp_q fcolumn
  unset -v nionallstep f1 f2 f3 f4 Tneut Tion Xfrac a_star N_LyC 
  unset -v rstatic1 rstatic2 rstatic3 lmax_hp M_loss v_wind
  unset -v ANALYTIC PLOT_PROG PLOT_FILE PROG_EX


done
# ============================================================================

echo 'Finished all tests'
exit 1
