# makefiletail.mk
# ============================================================================

# Other directories (for analysis routines)
# ----------------------------------------------------------------------------
#PGPLOTLIBS = -L$(PGPLOT_DIR) -lpgplot -L/sw/lib -lpng -laquaterm -Wl,-framework -Wl,Foundation -lSystemStubs #-lz
PGPLOTLIBS = -L$(PGPLOT_DIR) -lpgplot -lpng
X11LIBS = -L/usr/X11R6/lib -lX11


# Sub-directory linkage
# ----------------------------------------------------------------------------
MODULEDIR = $(SRCDIR)
HEADERS = $(SRCDIR)/headers
INCLUDE_DIR += $(MODULEDIR) $(HEADERS)
VPATH = $(SRCDIR)/advance:$(SRCDIR)/analyse:$(SRCDIR)/BHtree:$(SRCDIR)/binarytree:$(SRCDIR)/devel:$(SRCDIR)/dobs:$(SRCDIR)/ghosts:$(SRCDIR)/gravity:$(SRCDIR)/headers:$(SRCDIR)/ic:$(SRCDIR)/io:$(SRCDIR)/healpix:$(SRCDIR)/main:$(SRCDIR)/mhd:$(SRCDIR)/mpi:$(SRCDIR)/nbody:$(SRCDIR)/nbody_sim:$(SRCDIR)/radiation:$(SRCDIR)/setup:$(SRCDIR)/sinks:$(SRCDIR)/sorts:$(SRCDIR)/sph:$(SRCDIR)/sph_sim:$(SRCDIR)/tests:$(SRCDIR)/timestep:$(SRCDIR)/turbulence:$(SRCDIR)/user:$(SRCDIR)/newstuff


# Remove trailing whitespace from user options
# ----------------------------------------------------------------------------
F90                      := $(strip $(F90))
MPIF90                   := $(strip $(MPIF90))
OPENMP                   := $(strip $(OPENMP))
MPI                      := $(strip $(MPI))
MPI_LIBRARY              := $(strip $(MPI_LIBRARY))
COMPILER_MODE            := $(strip $(COMPILER_MODE))
OUTPUT_LEVEL             := $(strip $(OUTPUT_LEVEL))
DIAGNOSTIC_OUTPUT        := $(strip $(DIAGNOSTIC_OUTPUT))
NDIM                     := $(strip $(NDIM))
PRECISION                := $(strip $(PRECISION))
PERIODIC                 := $(strip $(PERIODIC))
X_BOUNDARY               := $(strip $(X_BOUNDARY))
Y_BOUNDARY               := $(strip $(Y_BOUNDARY))
Z_BOUNDARY               := $(strip $(Z_BOUNDARY))
GHOST_PARTICLES          := $(strip $(GHOST_PARTICLES))
SPH_SIMULATION           := $(strip $(SPH_SIMULATION))
NBODY_SIMULATION         := $(strip $(NBODY_SIMULATION))
SPH                      := $(strip $(SPH))
SPH_INTEGRATION          := $(strip $(SPH_INTEGRATION))
KERNEL                   := $(strip $(KERNEL))
HFIND                    := $(strip $(HFIND))
MINIMUM_H                := $(strip $(MINIMUM_H))
HYDRO                    := $(strip $(HYDRO))
THERMAL                  := $(strip $(THERMAL))
RAD_WS                   := $(strip $(RAD_WS))
SINK_POTENTIAL_WS        := $(strip $(SINK_POTENTIAL_WS))
AMBIENT_HEATING_WS       := $(strip $(AMBIENT_HEATING_WS))
SINK_HEATING_WS          := $(strip $(SINK_HEATING_WS))
FLUX_LIMITED_DIFFUSION   := $(strip $(FLUX_LIMITED_DIFFUSION))
IONIZING_RADIATION       := $(strip $(IONIZING_RADIATION))
PDR_CHEMISTRY            := $(strip $(PDR_CHEMISTRY))
STELLAR_WINDS            := $(strip $(STELLAR_WINDS))
ARTIFICIAL_VISCOSITY     := $(strip $(ARTIFICIAL_VISCOSITY))
VISC_TD                  := $(strip $(VISC_TD))
BALSARA                  := $(strip $(BALSARA))
ARTIFICIAL_CONDUCTIVITY  := $(strip $(ARTIFICIAL_CONDUCTIVITY))
EXTERNAL_FORCE           := $(strip $(EXTERNAL_FORCE))
SELF_GRAVITY             := $(strip $(SELF_GRAVITY))
MEANH_GRAVITY            := $(strip $(MEANH_GRAVITY))
EWALD                    := $(strip $(EWALD))
REMOVE_OUTLIERS          := $(strip $(REMOVE_OUTLIERS))
SINKS                    := $(strip $(SINKS))
SINK_RADIUS              := $(strip $(SINK_RADIUS))
SINK_REMOVE_ANGMOM       := $(strip $(SINK_REMOVE_ANGMOM))
NBODY_INTEGRATION        := $(strip $(NBODY_INTEGRATION))
BINARY_STATS             := $(strip $(BINARY_STATS))
SINK_GRAVITY_ONLY        := $(strip $(SINK_GRAVITY_ONLY))
FORCE_SPLITTING          := $(strip $(FORCE_SPLITTING))
TREE                     := $(strip $(TREE))
MULTIPOLE                := $(strip $(MULTIPOLE))
MAC                      := $(strip $(MAC))
REORDER                  := $(strip $(REORDER))
SORT                     := $(strip $(SORT))
TIMESTEP                 := $(strip $(TIMESTEP))
CHECK_NEIB_TIMESTEP      := $(strip $(CHECK_NEIB_TIMESTEP))
NEIGHBOURLISTS           := $(strip $(NEIGHBOURLISTS))
KERNEL_TABLES            := $(strip $(KERNEL_TABLES))
REMOVE_OUTLIERS          := $(strip $(REMOVE_OUTLIERS))
TURBULENT_FORCING        := $(strip $(TURBULENT_FORCING))
TIMING_CODE              := $(strip $(TIMING_CODE))
TEST                     := $(strip $(TEST))


# Object files always included in compilation list
# ----------------------------------------------------------------------------
MODULE_OBJ += definitions.o HP_types.o modules.o interface.o
SETUP_OBJ += convert_to_code_units_1.o convert_to_code_units_2.o 
SETUP_OBJ += default_parameters.o initialize_seren_variables_1.o
SETUP_OBJ += initialize_seren_variables_2.o paramstore.o read_parameters.o
SETUP_OBJ += read_arguments.o sanitycheck.o seren_setup.o
SETUP_OBJ += set_default_particle_types.o types.o units.o
SETUP_OBJ += write_makefile_options.o write_column_info.o
IO_OBJ += read_data.o write_data.o
IO_OBJ += write_rad_ws_test_data.o
IC_OBJ += ic_subroutines.o
GENERIC_OBJ += allocate_memory.o active_particle_list.o clean_up.o COM.o
GENERIC_OBJ += create_particle_list.o debug_tests.o distance2.o
GENERIC_OBJ += distance2_dp.o distance3.o distance3_dp.o errors.o heapsort.o
GENERIC_OBJ += insertion_sort.o remove_from_list.o reorder_array.o
GENERIC_OBJ += create_new_sph_particle.o sort_particle_types.o diagnostics.o

MPI_OBJ += mpi_modules.o mpi_start.o mpi_finish.o mpi_setup_types.o
MPI_OBJ += decomposition.o find_task.o fill_grid.o mpi_share_data.o
MPI_OBJ += mpi_setup_decomposition.o synchronise.o domain_comparisons.o
MPI_OBJ += loadbalance.o expanddomainboxes.o nonoverlapbox.o
MPI_OBJ += unpack_receiveparticles.o unwrap_particle_position.o
MPI_OBJ += transfer_particles.o remove_particles.o
MPI_OBJ += mpi_loadbalance_step.o loadbalance_time.o broadcastboundingboxes.o
MPI_OBJ += mpi_sph_update.o mpi_scatter_ghosts.o get_ghosts.o
MPI_OBJ += mpi_h_gather.o #BHghost_build.o


# User object files (Include your own object files here)
# ----------------------------------------------------------------------------
USER_MOD +=
USER_SUB +=


# Compiler flags
# ----------------------------------------------------------------------------
ifeq ($(F90),f95)
CFLAGS += -DCOMPILER_F95
#OPT += -I $(MODULEDIR) -I $(HEADERS)
OPT += -C=all
ifeq ($(OPENMP),1)
OPT += -fopenmp -DOPENMP
endif

else ifeq ($(F90),g95)
CFLAGS += -DCOMPILER_G95
#OPT += -I $(MODULEDIR) -I $(HEADERS)
ifeq ($(OPENMP),1)
OPT += -openmp -DOPENMP
endif
ifeq ($(COMPILER_MODE),DEBUG)
OPT += -Wall -pedantic
endif

else ifeq ($(F90),pgf90)
CFLAGS += -DCOMPILER_PGF90
#OPT += -I $(MODULEDIR) -I $(HEADERS)
OPT += -Minline=distance,distance2,distance3,gravity_sph,gravity_gradh
#OPT += -fast -fastsse -Mnontemporal
ifeq ($(OPENMP),1)
OPT += -mp=allcores -DOPENMP
endif

else ifeq ($(F90),pgf95)
CFLAGS += -DCOMPILER_PGF95
#OPT += -I $(MODULEDIR) -I $(HEADERS)
OPT += -Minline=distance,distance2,distance3,gravity_sph,gravity_gradh
#OPT += -fast -fastsse -Mnontemporal
ifeq ($(OPENMP),1)
OPT += -mp=allcores -DOPENMP
endif

else ifeq ($(F90),ifort)
CFLAGS += -DCOMPILER_IFORT -DF2003
#OPT += -I $(MODULEDIR) -I $(HEADERS)
ifeq ($(OPENMP),1)
OPT += -openmp -DOPENMP
endif
ifeq ($(COMPILER_MODE),DEBUG)
OPT += -g -warn all -debug all -check all -traceback -fp-stack-check -fpe-all=0 -ftrapuv -no-ftz
else ifeq ($(COMPILER_MODE),FAST)
OPT += -finline-functions -inline-forceinline -no-inline-min-size -inline-max-size=100 -ip -ipo -fp-model fast=2 -xHOST -O3 -no-prec-div -static-intel #-static
else ifeq ($(COMPILER_MODE),FASTDEBUG)
OPT += -g -warn all -debug all -traceback -check bounds -check format -check uninit
OPT += -finline-functions -inline-forceinline -no-inline-min-size -inline-max-size=100 -ip -ipo -fp-model fast=2 -xHOST -O3 -no-prec-div -static-intel #-static
else ifeq ($(COMPILER_MODE),STANDARD)
OPT += -finline-functions -inline-forceinline -no-inline-min-size -inline-max-size=100 -ip -ipo
endif

else ifeq ($(F90),gfortran)
CFLAGS += -DCOMPILER_GFORTRAN -DF2003
#OPT += -I $(MODULEDIR) -I $(HEADERS)
ifeq ($(OPENMP),1)
OPT += -fopenmp -DOPENMP
endif
ifeq ($(COMPILER_MODE),DEBUG)
OPT += -g -Wall -Wuninitialized -pedantic -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow#,denormal
else ifeq ($(COMPILER_MODE),FAST)
OPT += -Winline -fexpensive-optimizations -finline-functions -finline-limit=200 -funroll-loops -ftree-vectorize -ffast-math #-Ofast
#OPT += -mtune=core2 -march=core2 -mfpmath=sse -msse4.2 
else ifeq ($(COMPILER_MODE),STANDARD)
OPT += -Winline -fexpensive-optimizations -finline-functions -finline-limit=200 -funroll-loops -ftree-vectorize
endif

else ifeq ($(F90),sunf95)
CFLAGS += -DCOMPILER_SUNF95
#OPT += -I$(MODULEDIR) -I$(HEADERS)
OPT += -ansi
OPT += -inline=%auto,distance2,distance3,gravity_sph
ifeq ($(OPENMP),1)
OPT += -openmp=parallel -autopar -stackvar -DOPENMP
endif
ifneq ($(COMPILER_MODE),DEBUG)
OPT += -C -w4 -u #-xcheck=%all
endif

else
ERROR += "Invalid Fortran compiler selected : "$(F90)
endif


# Version no.
# ----------------------------------------------------------------------------
CFLAGS += -DVERSION_NO="$(VERSION_NO)"


# MPI
# ----------------------------------------------------------------------------
ifeq ($(MPI),1)
CFLAGS += -DUSE_MPI
GHOST_PARTICLES = 1
ifeq ($(MPI_LIBRARY),intel)
OPT += -mt_mpi
endif
else ifneq ($(MPI),0)
ERROR += "Invalid value for MPI : "$(MPI)
endif


# Profiling options
# ----------------------------------------------------------------------------
ifeq ($(COMPILER_MODE),DEBUG)
CFLAGS += -DCOMPILER_MODE_DEBUG
OPT += -pg -g
else ifeq ($(COMPILER_MODE),STANDARD)
CFLAGS += -DCOMPILER_MODE_STANDARD
OPT += -O3
else ifeq ($(COMPILER_MODE),FAST)
CFLAGS += -DCOMPILER_MODE_FAST
OPT += -O3
else ifeq ($(COMPILER_MODE),FASTDEBUG)
CFLAGS += -DCOMPILER_MODE_FASTDEBUG
OPT += -O3
else ifeq ($(COMPILER_MODE),0)
CFLAGS += -DCOMPILER_MODE_0
else
ERROR += "Invalid value for COMPILER_MODE : "$(COMPILER_MODE)
endif


# Diagnostic output
# ----------------------------------------------------------------------------
ifeq ($(DIAGNOSTIC_OUTPUT),1)
CFLAGS += -DDIAGNOSTIC_OUTPUT=1 -DDEBUG_DIAGNOSTICS -DDEBUG_TRACK_ENERGY
else ifeq ($(DIAGNOSTIC_OUTPUT),0)
CFLAGS += -DDIAGNOSTIC_OUTPUT=0
else
ERROR += "Invalid value for DIAGNOSTIC_OUTPUT : "$(DIAGNOSTIC_OUTPUT)
endif


# Dimensionality of the code
# ----------------------------------------------------------------------------
ifeq ($(NDIM),1)
CFLAGS += -DNDIM=1
else ifeq ($(NDIM),2)
CFLAGS += -DNDIM=2
else ifeq ($(NDIM),3)
CFLAGS += -DNDIM=3
else
ERROR += "Invalid value for NDIM : "$(NDIM)
endif


# Precision of real variables in code
# ----------------------------------------------------------------------------
ifeq ($(PRECISION),DOUBLE)
CFLAGS += -DDOUBLE_PRECISION
else ifeq ($(PRECISION),QUADRUPLE)
CFLAGS += -DQUADRUPLE_PRECISION
else ifneq ($(PRECISION),SINGLE)
ERROR += "Invalid PRECISION option : "$(PRECISION)
endif


# Input file format
# ----------------------------------------------------------------------------
CFLAGS += -DDRAGON_INPUT -DSEREN_INPUT -DASCII_INPUT
IO_OBJ += read_data_ascii.o read_data_dragon_form.o read_data_dragon_unform.o
IO_OBJ += read_data_seren_form.o read_data_seren_unform.o


# Output file format
# ----------------------------------------------------------------------------
CFLAGS += -DSEREN_OUTPUT -DASCII_OUTPUT
IO_OBJ += write_data_ascii.o write_data_seren_form.o write_data_seren_unform.o


# Ghost particles
# ----------------------------------------------------------------------------
ifeq ($(GHOST_PARTICLES),1)
CFLAGS += -DGHOST_PARTICLES
SPH_OBJ += expand.o
#ERROR += "Ghost particle algorithm not activated yet"
ifeq ($(PERIODIC),1)
SPH_OBJ += copy_data_to_ghosts.o create_ghost_particle.o search_ghost_particles.o
endif
endif

# Periodic boundary conditions
# ----------------------------------------------------------------------------

ifeq ($(PERIODIC),1)
CFLAGS += -DPERIODIC -DBOUNDARY_CONDITIONS
SPH_OBJ += check_boundary_conditions.o
MPI_OBJ += mpi_search_ghost_particles.o search_external_ghost_particles.o
MPI_OBJ += create_external_ghost_particle.o
else ifneq ($(PERIODIC),0)
ERROR += "Invalid PERIODIC option selected : "$(PERIODIC)
endif

ifeq ($(X_BOUNDARY),PERIODIC)
CFLAGS += -DPERIODIC_X
else ifeq ($(X_BOUNDARY),WALL)
CFLAGS += -DWALL_X_LHS -DWALL_X_RHS
else ifeq ($(X_BOUNDARY),WALL_LHS)
CFLAGS += -DWALL_X_LHS
else ifeq ($(X_BOUNDARY),WALL_RHS)
CFLAGS += -DWALL_X_RHS
else ifneq ($(X_BOUNDARY),0)
ERROR += "Invalid X_BOUNDARY option : "$(X_BOUNDARY)
endif

ifeq ($(Y_BOUNDARY),PERIODIC)
CFLAGS += -DPERIODIC_Y
else ifeq ($(Y_BOUNDARY),WALL)
CFLAGS += -DWALL_Y_LHS -DWALL_Y_RHS
else ifeq ($(Y_BOUNDARY),WALL_LHS)
CFLAGS += -DWALL_Y_LHS
else ifeq ($(Y_BOUNDARY),WALL_RHS)
CFLAGS += -DWALL_Y_RHS
else ifneq ($(Y_BOUNDARY),0)
ERROR += "Invalid Y_BOUNDARY option : "$(Y_BOUNDARY)
endif

ifeq ($(Z_BOUNDARY),PERIODIC)
CFLAGS += -DPERIODIC_Z
else ifeq ($(Z_BOUNDARY),WALL)
CFLAGS += -DWALL_Z_LHS -DWALL_Z_RHS
else ifeq ($(Z_BOUNDARY),WALL_LHS)
CFLAGS += -DWALL_Z_LHS
else ifeq ($(Z_BOUNDARY),WALL_RHS)
CFLAGS += -DWALL_Z_RHS
else ifneq ($(Z_BOUNDARY),0)
ERROR += "Invalid Z_BOUNDARY option : "$(Z_BOUNDARY)
endif


# Simulation-mode flags
# ----------------------------------------------------------------------------
ifeq ($(SPH_SIMULATION),1)
CFLAGS += -DSPH_SIMULATION
OBJ += sph_simulation.o sph_setup.o sph_integrate.o
OBJ += sph_output.o sph_timesteps.o
OBJ += sph_grav_forces.o sph_sink_forces.o
SETUP_OBJ += initialize_sph_variables_1.o initialize_sph_variables_2.o
INCLUDE_SPH_OBJS = 1
else ifneq ($(SPH_SIMULATION),0)
ERROR += "Invalid SPH_SIMULATION option : "$(SPH_SIMULATION)
endif

ifeq ($(NBODY_SIMULATION),1)
CFLAGS += -DNBODY_SIMULATION -DSINKS
NBODY_OBJ += nbody_simulation.o nbody_setup.o
NBODY_OBJ += nbody_grav_forces.o
NBODY_OBJ += nbody_integrate.o nbody_output.o nbody_timesteps.o
INCLUDE_NBODY_OBJS = 1
GRAVITY = 1
ifeq ($(MPI),1)
ERROR += "NBODY simulation not MPI parallelized"
endif
else ifneq ($(NBODY_SIMULATION),0)
ERROR += "Invalid NBODY_SIMULATION option : "$(NBODY_SIMULATION)
endif


# SPH object files
# ============================================================================
ifeq ($(INCLUDE_SPH_OBJS),1)
SPH_OBJ += all_sph.o bounding_box.o \
gather_neib_on_fly.o get_neib.o get_neib_on_fly.o \
h_guess.o reduce_particle_timestep.o reduce_timesteps.o sph_update.o \
timestep_size.o track_particles.o tree_update.o
IO_OBJ += record_particle_data.o write_data_debug.o write_data_grid_results.o
SETUP_OBJ += initialize_thermal_properties.o


# SPH mode
# ----------------------------------------------------------------------------
ifeq ($(SPH),GRAD_H_SPH)
CFLAGS += -DGRAD_H_SPH -DH_RHO
SPH_OBJ += h_rho_iteration.o
else ifeq ($(SPH),SM2012_SPH)
CFLAGS += -DSM2012_SPH
else ifeq ($(SPH),RPSPH)
CFLAGS += -DRPSPH
else ifeq ($(SPH),RTSPH)
CFLAGS += -DRTSPH
else ifneq ($(SPH),STANDARD)
ERROR += "Invalid SPH option selected : "$(SPH)
endif


# Integration scheme used in code
# ----------------------------------------------------------------------------
SPH_OBJ += sph_advance.o advance_boundary_particle.o

ifeq ($(SPH_INTEGRATION),EULER)
CFLAGS += -DEULER
SPH_OBJ += advance_euler.o
else ifeq ($(SPH_INTEGRATION),RK2)
CFLAGS += -DRUNGE_KUTTA -DRUNGE_KUTTA2
SPH_OBJ += advance_runge_kutta.o
else ifeq ($(SPH_INTEGRATION),LFKDK)
CFLAGS += -DLEAPFROG_KDK
SPH_OBJ += advance_leapfrog_kdk.o leapfrog_kdk_correction_terms.o
else ifeq ($(SPH_INTEGRATION),LFDKD)
CFLAGS += -DLEAPFROG_DKD
SPH_OBJ += advance_leapfrog_dkd.o
else 
ERROR += "Invalid SPH_INTEGRATION option selected : "$(SPH_INTEGRATION)
endif


# Method used to calculate h if not using 'grad-h' SPH
# ----------------------------------------------------------------------------
ifeq ($(HFIND),NUMBER)
CFLAGS += -DHGATHER
SPH_OBJ += h_gather.o
else ifeq ($(HFIND),MASS)
CFLAGS += -DHMASS
SPH_OBJ += h_gather_mass.o
else ifeq ($(HFIND),CONSTANT)
CFLAGS += -DCONSTANT_H
else ifeq ($(HFIND),H_RHO)
CFLAGS += -DH_RHO
ifneq ($(SPH),GRAD_H_SPH)
SPH_OBJ += h_rho_iteration.o
endif
else ifeq ($(HFIND),H_NUMBER)
CFLAGS += -DH_RHO -DH_NUMBER
ifneq ($(SPH),GRAD_H_SPH)
SPH_OBJ += h_rho_iteration.o
endif
else 
ERROR += "Invalid HFIND option selected : "$(HFIND)
endif


# Use a minimum smoothing length
# ----------------------------------------------------------------------------
ifeq ($(MINIMUM_H),1)
CFLAGS += -DMINIMUM_H
else ifneq ($(MINIMUM_H),0)
ERROR += "Invalid MINIMUM_H option selected : "$(MINIMUM_H)
endif


# Hydrodynamic forces
# ----------------------------------------------------------------------------
ifeq ($(HYDRO),1)
CFLAGS += -DHYDRO
OBJ += sph_hydro_forces.o
ifeq ($(SPH),GRAD_H_SPH)
SPH_OBJ += hydro_gradh.o
else
SPH_OBJ += hydro.o
endif
else ifeq ($(HYDRO),0)
ERROR += "Invalid HYDRO option selected : "$(HYDRO)
endif


# Thermal physics 
# ----------------------------------------------------------------------------
SPH_OBJ += update_thermal_properties.o thermal.o thermal_properties.o

ifeq ($(ENERGY_EQN),1)
CFLAGS += -DINTERNAL_ENERGY -DENERGY_EQN
else ifneq ($(ENERGY_EQN),0)
ERROR += "Invalid ENERGY_EQN option selected : "$(ENERGY_EQN)
endif

ifeq ($(ENTROPY_EQN),1)
CFLAGS += -DENTROPIC_FUNCTION -DINTERNAL_ENERGY -DENTROPY_EQN
else ifneq ($(ENTROPY_EQN),0)
ERROR += "Invalid ENTROPY_EQN option selected : "$(ENTROPY_EQN)
endif

ifeq ($(RAD_WS),1)
CFLAGS += -DRAD_WS -DRAD -DINTERNAL_ENERGY -DENERGY_EQN -DU_IMPLICIT_SOLVER
SPH_OBJ += rad_ws_update.o find_equilibrium_temp_ws.o
SPH_OBJ += read_cooling_table_ws.o ambienttemp.o
ifeq ($(MPI),1)
ERROR += "RAD_WS not MPI parallelized"
endif
ifeq ($(FLUX_LIMITED_DIFFUSION),1)
CFLAGS += -DFLUX_LIMITED_DIFFUSION -DDIFFUSION
SPH_OBJ += conductivity.o diffusion.o
endif
ifeq ($(AMBIENT_HEATING_WS),1)
CFLAGS += -DAMBIENT_HEATING -DCONST_HEATING
endif
ifeq ($(SINK_HEATING_WS),HDISC_HEATING)
CFLAGS += -DHDISC_HEATING
endif
ifeq ($(SINK_HEATING_WS),HDISC_HEATING_3D_SINGLE)
CFLAGS += -DHDISC_HEATING_3D_SINGLE
endif
ifeq ($(SINK_HEATING_WS),STAR_HEATING)
CFLAGS += -DSTAR_HEATING
endif
ifeq ($(SINK_HEATING_WS),STAR_SIMPLE_HEATING)
CFLAGS += -DSTAR_SIMPLE_HEATING
endif
ifeq ($(SINK_POTENTIAL_WS),1)
CFLAGS += -DRAD_WS_SINK_POT
endif
else ifneq ($(RAD_WS),0)
ERROR += "Invalid RAD_WS option selected : "$(RAD_WS)
endif


# Cooling/heating
# ----------------------------------------------------------------------------
ifeq ($(COOLING_HEATING),EXPLICIT)
CFLAGS += -DCOOLING_HEATING -DEXPLICIT_COOLING_HEATING
OBJ += cooling_heating_rate.o
else ifeq ($(COOLING_HEATING),EXPONENTIAL)
CFLAGS += -DCOOLING_HEATING -DEXPONENTIAL_COOLING_HEATING
OBJ += cooling_heating_rate.o
else ifneq ($(COOLING_HEATING),0)
ERROR += "Invalid COOLING_HEATING option selected : "$(COOLING_HEATING)
endif


# HEALPix routines
# ----------------------------------------------------------------------------
ifeq ($(IONIZING_RADIATION),SINGLE_STATIC_SOURCE)
HEALPIX = 1
CFLAGS += -DIONIZING_UV_RADIATION -DSINGLE_STATIC_SOURCE
SPH_OBJ += HP_ionizing_radiation.o HP_walk_ray.o write_ionization_data.o
else ifeq ($(IONIZING_RADIATION),SINGLE_SINK_SOURCE)
HEALPIX = 1
CFLAGS += -DIONIZING_UV_RADIATION -DSINGLE_SINK_SOURCE
SPH_OBJ += HP_ionizing_radiation.o HP_walk_ray.o write_ionization_data.o
else ifneq ($(IONIZING_RADIATION),0)
ERROR += "Invalid IONIZING_RADIATION option selected : "$(IONIZING_RADIATION)
endif
ifneq ($(IONIZING_RADIATION),0)
ifeq ($(MPI),1)
ERROR += "Ionizing radiation not MPI parallelized"
endif
endif

ifeq ($(STELLAR_WIND),SINGLE_STATIC_SOURCE)
HEALPIX = 1
CFLAGS += -DSTELLAR_WIND -DMOMENTUM_WIND -DSINGLE_STATIC_SOURCE
SPH_OBJ += HP_stellar_feedback.o
else ifeq ($(STELLAR_WIND),SINGLE_SINK_SOURCE)
HEALPIX = 1
CFLAGS += -DSTELLAR_WIND -DMOMENTUM_WIND -DSINGLE_SINK_SOURCE
SPH_OBJ += HP_stellar_feedback.o
else ifneq ($(STELLAR_WIND),0)
ERROR += "Invalid STELLAR_WIND option selected : "$(STELLAR_WIND)
endif
ifneq ($(STELLAR_WIND),0)
ifeq ($(MPI),1)
ERROR += "STELLAR_WIND not MPI parallelized"
endif
endif


ifeq ($(HEALPIX),1)
CFLAGS += -DHEALPIX -DTRAPEZOIDAL_RULE
SPH_OBJ += HP_update.o HP_calculate_basis_vector.o HP_evaluation_point.o
SPH_OBJ += HP_initialize_source.o HP_inverse_positions.o HP_reorder_lists.o
SPH_OBJ += HP_rhoh_ep.o HP_split_active_rays.o HP_walk_all_rays.o
SPH_OBJ += create_HP_source.o healpix.o initialize_HP_sources.o
ifeq ($(STELLAR_MODEL),1)
CFLAGS += -DSTELLAR_FEEDBACK
SETUP_OBJ += read_stellar_model_table.o
SPH_OBJ += calculate_sink_properties.o
endif
endif


# Artificial viscosity
# ----------------------------------------------------------------------------
ifeq ($(ARTIFICIAL_VISCOSITY),AB)
CFLAGS += -DARTIFICIAL_VISCOSITY -DVISC_AB
else ifeq ($(ARTIFICIAL_VISCOSITY),MON97)
CFLAGS += -DARTIFICIAL_VISCOSITY -DVISC_MON97
else ifneq ($(ARTIFICIAL_VISCOSITY),0)
ERROR += "Invalid ARTIFICIAL_VISCOSITY option selected : "$(ARTIFICIAL_VISCOSITY)
endif

ifneq ($(ARTIFICIAL_VISCOSITY),0)
ifeq ($(VISC_TD),1)
CFLAGS += -DVISC_TD
else ifneq ($(VISC_TD),0)
ERROR += "Invalid VISC_TD option selected : "$(VISC_TD)
endif
ifeq ($(BALSARA),1)
CFLAGS += -DVISC_BALSARA
else ifneq ($(BALSARA),0)
ERROR += "Invalid BALSARA option selected : "$(BALSARA)
endif
endif


# Artificial conductivity
# ----------------------------------------------------------------------------
ifeq ($(ARTIFICIAL_CONDUCTIVITY),PRICE2008)
CFLAGS += -DARTIFICIAL_CONDUCTIVITY -DCOND_PRICE2008
else ifeq ($(ARTIFICIAL_CONDUCTIVITY),WADSLEY2008)
CFLAGS += -DARTIFICIAL_CONDUCTIVITY -DCOND_WADSLEY2008
else ifneq ($(ARTIFICIAL_CONDUCTIVITY),0)
ERROR += "Invalid ARTIFICIAL_CONDUCTIVITY option selected"
endif


# External pressure flag
# ----------------------------------------------------------------------------
ifeq ($(EXTERNAL_PRESSURE),1)
CFLAGS += -DEXTERNAL_PRESSURE
else ifneq ($(EXTERNAL_PRESSURE),0)
ERROR += "Invalid EXTERNAL_PRESSURE option selected : "$(EXTERNAL_PRESSURE)
endif


# Gravity forces
# ----------------------------------------------------------------------------
ifeq ($(SELF_GRAVITY),KS)
CFLAGS += -DGRAVITY -DSELF_GRAVITY
GRAVITY = 1
else ifeq ($(SELF_GRAVITY),NBODY)
CFLAGS += -DGRAVITY -DSELF_GRAVITY -DN_BODY
GRAVITY = 1
else ifneq ($(SELF_GRAVITY),0)
ERROR += "Invalid SELF_GRAVITY option selected"
endif
ifneq ($(SELF_GRAVITY),0)
MPI_OBJ += gravity_export_to.o do_export_gravity.o gravity_export_return.o
endif


# Sink-gravity only
# ----------------------------------------------------------------------------
ifneq ($(SINKS),0)
#SPH_OBJ += gravity_gradh.o
ifeq ($(SINK_GRAVITY_ONLY),KS)
CFLAGS += -DGRAVITY -DSINK_GRAVITY_ONLY
GRAVITY = 1
else ifeq ($(SINK_GRAVITY_ONLY),NBODY)
CFLAGS += -DGRAVITY -DSINK_GRAVITY_ONLY -DN_BODY
GRAVITY = 1
else ifneq ($(SINK_GRAVITY_ONLY),0)
ERROR += "Invalid SINK_GRAVITY_ONLY option selected : "$(SINK_GRAVITY_ONLY)
endif
endif


ifeq ($(GRAVITY),1)
SPH_OBJ += gravity_sph.o gravity_meanh.o
SPH_OBJ += gravity_nbody.o
ifeq ($(SPH),GRAD_H_SPH)
SPH_OBJ += gravity_gradh.o gravity_gradh_meanh.o
endif
endif


#ifneq ($(SELF_GRAVITY),0)
#ifneq ($(SINK_GRAVITY_ONLY),1)
#CFLAGS += -DSELF_GRAVITY
##else ifneq ($(SINK_GRAVITY_ONLY),0)
##ERROR += "Invalid value for SINK_GRAVITY_ONLY : "$(SINK_GRAVITY_ONLY)
#endif
#endif


# Tree code
# ----------------------------------------------------------------------------
ifeq ($(TREE),BH)
CFLAGS += -DBH_TREE
SPH_OBJ += BH_build_skeleton.o
SPH_OBJ += BHhydro_build.o BHhydro_stock.o BHhydro_update_hmax.o
SPH_OBJ += BHhydro_walk.o BHhydrowalk_hgather.o BHhydro_hguess.o
SPH_OBJ += BH_remove_particles.o BHhydro_ray_search.o
SPH_OBJ += BH_add_particles.o
ifeq ($(GHOST_PARTICLES),1)
SPH_OBJ += BHghost_build.o #BHghost_walk.o BHghostwalk_hgather.o
endif
ifneq ($(SELF_GRAVITY),0)
SPH_OBJ += BHgrav_build.o BHgrav_stock.o BHgrav_accel.o BHgrav_node_accel.o
SPH_OBJ += copy_BHhydro_to_BHgrav.o
MPI_OBJ += BHgrav_slice.o share_pruned_gravtrees.o BHtree_remote_grav.o
endif
ifeq ($(REORDER),TREE)
CFLAGS += -DREORDER_TREE
SPH_OBJ += BH_reorder_tree.o
endif
ifeq ($(REORDER),PARTICLES)
CFLAGS += -DREORDER_PARTICLES
SPH_OBJ += BH_reorder_particles.o
endif
ifeq ($(REORDER),ALL)
CFLAGS += -DREORDER_TREE -DREORDER_PARTICLES
SPH_OBJ += BH_reorder_particles.o BH_reorder_tree.o
endif

else ifeq ($(TREE),0)
ifneq ($(SELF_GRAVITY),0)
SPH_OBJ += direct_sph_gravity.o
endif

else 
ERROR += "Invalid TREE option selected : "$(TREE)
endif


# Multipole moment expansion in gravity tree
# ----------------------------------------------------------------------------
ifeq ($(MULTIPOLE),QUADRUPOLE)
CFLAGS += -DQUADRUPOLE
else ifeq ($(MULTIPOLE),OCTUPOLE)
CFLAGS += -DQUADRUPOLE -DOCTUPOLE
else ifneq ($(MULTIPOLE),0)
ERROR += "Invalid TREE option selected : "$(TREE)
endif


# Tree gravity MAC
# ----------------------------------------------------------------------------
ifeq ($(MAC),GEOMETRIC)
CFLAGS += -DGEOMETRIC_MAC
else ifeq ($(MAC),GADGET)
CFLAGS += -DGADGET_MAC
else ifeq ($(MAC),GADGET2)
CFLAGS += -DGADGET2_MAC
else ifeq ($(MAC),NEW)
CFLAGS += -DNEW_MAC
else ifeq ($(MAC),EIGEN)
CFLAGS += -DEIGEN_MAC
SPH_OBJ += eigenvalue_mac.o
else
ERROR += "Invalid value for MAC : "$(MAC)
endif


# SPH kernel-softeneing MAC
# ----------------------------------------------------------------------------
#ifeq ($(SPH_KS_MAC),1)
CFLAGS += -DSPH_KS_MAC
#else
#ifneq ($(SELF_GRAVITY),NBODY)
#SPH_OBJ += gravity_nbody.o
#endif
#endif


# Sink particles
# ----------------------------------------------------------------------------
ifneq ($(SINKS),0)
CFLAGS += -DSINKS -DACCRETION_RATE -DDEBUG_FORCES
SPH_OBJ += sink_advance.o sink_timestep.o sink_update.o
SPH_OBJ += write_accreted_particles.o write_sink_data.o
MPI_OBJ += sink_transfer_particles.o sink_share.o mpi_find_densest_particle.o

ifeq ($(SINKS),SMOOTH_ACC)
CFLAGS += -DSMOOTH_ACCRETION -DMINIMUM_H
SPH_OBJ += create_sink.o smooth_accrete_particles.o
SPH_OBJ += sink_accretion_properties.o
else ifeq ($(SINKS),SMOOTH_ACC2)
CFLAGS += -DSMOOTH_ACCRETION -DMINIMUM_H -DACCRETION_RADIUS
SPH_OBJ += create_sink.o smooth_accrete_particles.o
SPH_OBJ += sink_accretion_properties.o
else ifeq ($(SINKS),SIMPLE)
SPH_OBJ += create_sink.o accrete_particles.o sink_accretion_properties.o
else ifeq ($(SINKS),NO_ACC)
CFLAGS += -DNO_ACCRETION
SPH_OBJ += create_sink.o accrete_particles.o sink_accretion_properties.o
else ifneq ($(SINKS),0)
ERROR += "Invalid SINKS option selected : "$(SINKS)
endif

ifeq ($(SPH_INTEGRATION),EULER)
SPH_OBJ += advance_sink_euler.o
else ifeq ($(SPH_INTEGRATION),RK2)
SPH_OBJ += advance_sink_RK.o
else ifeq ($(SPH_INTEGRATION),LFKDK)
SPH_OBJ += advance_sink_LFKDK.o
else ifeq ($(SPH_INTEGRATION),LFDKD)
SPH_OBJ += advance_sink_LFDKD.o
else ifneq ($(SPH_INTEGRATION),0)
ERROR += "Invalid SPH_INTEGRATION option selected : "$(SPH_INTEGRATION)
endif

ifeq ($(GRAVITY),1)
SPH_OBJ += direct_sink_gravity.o
endif
ifneq ($(SELF_GRAVITY),0)
SPH_OBJ += sink_search.o
endif

ifeq ($(SINK_RADIUS),FIXED_ABSOLUTE)
CFLAGS += -DFIXED_ABSOLUTE_SINKRAD
else ifeq ($(SINK_RADIUS),FIXED_HMULT)
CFLAGS += -DFIXED_HMULT_SINKRAD
else ifeq ($(SINK_RADIUS),HMULT)
CFLAGS += -DHMULT_SINKRAD
else ifneq ($(SINK_RADIUS),0)
ERROR += "Invalid SINK_RADIUS option selected : "$(SINK_RADIUS)
endif

ifeq ($(SINK_REMOVE_ANGMOM),1)
CFLAGS += -DSINK_REMOVE_ANGMOM #-DCHECK_NEIGHBOUR_TIMESTEPS -DIMMEDIATE_TIMESTEP_REDUCTION
SPH_OBJ += redistribute_sink_angmom.o
endif

endif


# Check neighbour's timesteps
# ----------------------------------------------------------------------------
ifeq ($(CHECK_NEIB_TIMESTEP),1)
CFLAGS += -DCHECK_NEIGHBOUR_TIMESTEPS
OBJ += check_neighbour_timesteps.o
else ifeq ($(CHECK_NEIB_TIMESTEP),2)
CFLAGS += -DCHECK_NEIGHBOUR_TIMESTEPS -DIMMEDIATE_TIMESTEP_REDUCTION
OBJ += check_neighbour_timesteps.o
else ifneq ($(CHECK_NEIB_TIMESTEP),0)
ERROR += "Invalid CHECK_NEIB_TIMESTEP option selected : "$(CHECK_NEIB_TIMESTEP)
endif
ifneq ($(CHECK_NEIB_TIMESTEP),0)
ifeq ($(MPI),1)
ERROR += "CHECK_NEIB_TIMESTEP not MPI parallelized"
endif
endif


# Neighbour lists
# ----------------------------------------------------------------------------
ifeq ($(NEIGHBOURLISTS),PARTICLES)
CFLAGS += -DNEIGHBOUR_LISTS
else ifeq ($(NEIGHBOURLISTS),1)
CFLAGS += -DNEIGHBOUR_LISTS
else ifneq ($(NEIGHBOURLISTS),0)
ERROR += "Invalid NEIGHBOURLISTS option selected : "$(NEIGHBOURLISTS)
endif
ifneq ($(NEIGHBOURLISTS),0)
ifeq ($(MPI),1)
ERROR += "NEIGHBOURLISTS have never been tested for MPI, so probably don't work"
endif
endif


endif
# ============================================================================



# Include N-body source files
# ============================================================================
ifeq ($(INCLUDE_NBODY_OBJS),1)
CFLAGS += -DSINKS
NBODY_OBJ += copy_sinks_to_stars.o copy_stars_to_sinks.o
NBODY_OBJ += nbody_advance.o nbody_correction_terms.o nbody_end_step.o
NBODY_OBJ += nbody_timestep_size.o write_star_data.o nbody_diagnostics.o
SETUP_OBJ += nbody_accrete_bound_particles.o nbody_hermite4_extra_terms.o

ifeq ($(NBODY_INTEGRATION),HERMITE4)
CFLAGS += -DNBODY_HERMITE4
NBODY_OBJ += gravity_hermite4_meanh.o gravity_hermite4_gradh_meanh.o
NBODY_OBJ += gravity_hermite4.o nbody_hermite4_direct_gravity.o
else ifeq ($(NBODY_INTEGRATION),HERMITE4_TS)
CFLAGS += -DNBODY_HERMITE4 -DNBODY_HERMITE4_TS
NBODY_OBJ += gravity_hermite4_meanh.o gravity_hermite4_gradh_meanh.o
NBODY_OBJ += gravity_hermite4.o nbody_hermite4_direct_gravity.o
else ifeq ($(NBODY_INTEGRATION),LFKDK)
CFLAGS += -DNBODY_HERMITE4 -DNBODY_LEAPFROG_KDK
NBODY_OBJ += gravity_hermite4_meanh.o gravity_hermite4_gradh_meanh.o
NBODY_OBJ += gravity_hermite4.o nbody_hermite4_direct_gravity.o
else
ERROR += "Invalid value for NBODY_INTEGRATION : "$(NBODY_INTEGRATION)
endif

ifeq ($(NBODY_SIMULATION),1)
ifneq ($(SPH_SIMULATION),1)
SPH_OBJ += gravity_sph.o
endif
endif

ifeq ($(BINARY_STATS),1)
CFLAGS += -DBINARY_STATS
NBODY_OBJ += binary_energy.o binary_properties.o binary_search.o
else ifneq ($(BINARY_STATS),0)
ERROR += "Invalid value for BINARY_STATS : "$(BINARY_STATS)
endif

endif
# ============================================================================


# Kernel used to calculate SPH quantities
# ----------------------------------------------------------------------------
OBJ += kernel.o
ifeq ($(KERNEL),M4)
CFLAGS += -DM4_KERNEL
else ifeq ($(KERNEL),M4TC)
CFLAGS += -DM4_KERNEL -DTC_KERNEL
else ifeq ($(KERNEL),QUINTIC)
CFLAGS += -DQUINTIC_KERNEL
else ifeq ($(KERNEL),QUINTICTC)
CFLAGS += -DQUINTIC_KERNEL -DTC_KERNEL
else ifeq ($(KERNEL),GAUSSIAN_3H)
CFLAGS += -DGAUSSIAN_3H_KERNEL
else ifeq ($(KERNEL),LINEAR)
CFLAGS += -DLINEAR_KERNEL
else
ERROR += "Invalid KERNEL option selected : "$(KERNEL)
endif
ifeq ($(KERNEL_TABLES),1)
CFLAGS += -DKERNEL_TABLES
OBJ += tabulate_kernel_functions.o
else ifneq ($(KERNEL_TABLES),0)
ERROR += "Invalid KERNEL_TABLES option selected : "$(KERNEL_TBALES)
endif


# Gravity flags
# ----------------------------------------------------------------------------
ifeq ($(GRAVITY),1)
CFLAGS += -DGRAVITY
ifeq ($(MEANH_GRAVITY),1)
CFLAGS += -DMEANH_GRAVITY
ifeq ($(INCLUDE_SPH_OBJS),1)
ifeq ($(SPH),GRAD_H_SPH)
OBJ += meanh_zeta.o
endif
endif
endif
endif


# External gravitational force
# ----------------------------------------------------------------------------
ifeq ($(EXTERNAL_FORCE),PLUMMER)
CFLAGS += -DGRAVITY -DEXTERNAL_FORCE -DPLUMMER_POTENTIAL
OBJ += add_external_gravitational_force.o
else ifeq ($(EXTERNAL_FORCE),UDS)
CFLAGS += -DGRAVITY -DEXTERNAL_FORCE -DUDS_POTENTIAL
OBJ += add_external_gravitational_force.o
else ifeq ($(EXTERNAL_FORCE),NFW1996)
CFLAGS += -DGRAVITY -DEXTERNAL_FORCE -DNFW1996_POTENTIAL
OBJ += add_external_gravitational_force.o
else ifneq ($(EXTERNAL_FORCE),0)
ERROR += "Invalid EXTERNAL_FORCE option selected"
endif


# Ewald forces
# ----------------------------------------------------------------------------
ifeq ($(EWALD),1)
CFLAGS += -DEWALD
OBJ += ewald_init.o ewald_force.o
ifeq ($(MPI),1)
ERROR += "EWALD gravity not MPI parallelized"
endif
else ifneq ($(EWALD),0)
ERROR += "Invalid EWALD option selected : "$(EWALD)
endif


# Remove outliers
# ----------------------------------------------------------------------------
ifeq ($(REMOVE_OUTLIERS),1)
CFLAGS += -DREMOVE_OUTLIERS
SPH_OBJ += remove_outlying_particles.o
ifeq ($(MPI),1)
ERROR += "REMOVE_OUTLIERS not MPI parallelized"
endif
else ifneq ($(REMOVE_OUTLIERS),0)
ERROR += "Invalid REMOVE_OUTLIERS option selected : "$(REMOVE_OUTLIERS)
endif


# Sorting algorithm
# ----------------------------------------------------------------------------
ifeq ($(SORT),INSERTION)
CFLAGS += -DINSERTION_SORT
else ifeq ($(SORT),HEAP)
CFLAGS += -DHEAPSORT
else
ERROR += "Invalid SORT option selected : "$(SORT)
endif


# Multiple particle timestep options
# ----------------------------------------------------------------------------
ifeq ($(TIMESTEP),ADAPTIVE)
CFLAGS += -DADAPTIVE_TIMESTEP_LEVELS
else ifeq ($(TIMESTEP),RESTRICTED)
CFLAGS += -DRESTRICTED_TIMESTEP_LEVELS
else ifeq ($(TIMESTEP),FIXED)
CFLAGS += -DFIXED_TIMESTEP_LEVELS
else ifneq ($(TIMESTEP),0)
ERROR += "Invalid TIMESTEP option selected : "$(TIMESTEP)
endif


# Turbulent forcing
# ----------------------------------------------------------------------------
ifeq ($(TURBULENT_FORCING),1)
CFLAGS += -DTURBULENT_FORCING
OBJ += turb_force_module.o turb_force_init.o turb_next_field.o
OBJ += turb_forces.o turb_force_apply.o turb_check_time.o
OBJ += read_turb_fields.o write_turb_fields.o
MPI_OBJ += mpi_share_turb_fields.o
FFTW = 1
else ifneq ($(TURBULENT_FORCING),0)
ERROR += "Invalid TURBULENT_FORCING option selected : "$(TURBULENT_FORCING)
endif


# Timing code
# ----------------------------------------------------------------------------
ifeq ($(TIMING_CODE),1)
CFLAGS += -DTIMING
OBJ += timing.o write_timing_stats.o
else ifneq ($(TIMING_CODE),0)
ERROR += "Invalid TIMING_CODE option selected : "$(TIMING_CODE)
endif


# Test flags and routines
# ----------------------------------------------------------------------------
ifeq ($(TEST),SPIEGEL)
CFLAGS += -DSPIEGEL_TEST -DSTATIC_PARTICLES
else ifeq ($(TEST),FREEFALL)
CFLAGS += -DFREEFALL_TEST
else ifeq ($(TEST),BINARY)
CFLAGS += -DBINARY_TEST
OBJ += write_binary_data.o
else ifeq ($(TEST),PLUMMER)
CFLAGS += -DPLUMMER_TEST -DLAGRANGIAN_RADII
OBJ += write_lagrangian_radii.o
else ifeq ($(TEST),ENTROPY)
CFLAGS += -DENTROPY_CORE_TEST
else ifeq ($(TEST),WIND)
CFLAGS += -DWIND_TEST
OBJ += write_wind_data.o
else ifeq ($(TEST),NBODY_DEFLECTION)
CFLAGS += -DNBODY_DEFLECTION -DSTATIC_PARTICLES
else ifneq ($(TEST),0)
ERROR += "Invalid TEST option selected : "$(TEST)
endif


# FFTW
# ----------------------------------------------------------------------------
ifeq ($(FFTW),1)
ifneq ($(FFTW_DIR),)
INCLUDE_DIR += $(FFTW_DIR)/include
LDFLAGS += -L $(FFTW_DIR)/lib -lfftw3 #-lm
STATIC_LIBS += $(FFTW_DIR)/lib/libfftw3.a
else
LDFLAGS += -lfftw3 #-lm
ifneq ($(wildcard /usr/lib/libfftw3.a),)
STATIC_LIBS += /usr/lib/libfftw3.a
else ifneq ($(wildcard /usr/lib64/libfftw3.a),)
STATIC_LIBS += /usr/lib64/libfftw3.a
else ifneq ($(wildcard /usr/lib/i386-linux-gnu/libfftw3.a),)
STATIC_LIBS += /usr/lib/i386-linux-gnu/libfftw3.a
else ifneq ($(wildcard /usr/lib/x86_64-linux-gnu/libfftw3.a),)
STATIC_LIBS += /usr/lib/x86_64-linux-gnu/libfftw3.a
else ifneq ($(wildcard /lib/i386-linux-gnu/libfftw3.a),)
STATIC_LIBS += /lib/i386-linux-gnu/libfftw3.a
else ifneq ($(wildcard /lib/x86_64-linux-gnu/libfftw3.a),)
STATIC_LIBS += /lib/x86_64-linux-gnu/libfftw3.a
endif
endif
endif


# Debug flags
# ----------------------------------------------------------------------------
ifeq ($(OUTPUT_LEVEL),1)
CFLAGS += -DDEBUG1 -DOUTPUT_LEVEL=1
else ifeq ($(OUTPUT_LEVEL),2)
CFLAGS += -DDEBUG1 -DDEBUG2 -DOUTPUT_LEVEL=2
else ifeq ($(OUTPUT_LEVEL),3)
CFLAGS += -DDEBUG1 -DDEBUG2 -DDEBUG3 -DOUTPUT_LEVEL=3
else ifeq ($(OUTPUT_LEVEL),0)
CFLAGS += -DOUTPUT_LEVEL=0
else
ERROR += "Invalid value for OUTPUT_LEVEL : "$(OUTPUT_LEVEL)
endif


# List of all object files
# ----------------------------------------------------------------------------
OBJ += $(MODULE_OBJ) $(USER_MOD) $(USER_SUB) $(GENERIC_OBJ) $(SPH_OBJ) $(NBODY_OBJ) $(IO_OBJ) $(SETUP_OBJ)

# Add MPI options if required
ifeq ($(MPI),1)
F90 = $(MPIF90)
OBJ += $(MPI_OBJ)
MPI_REQUIRED=mpi_modules.o
endif


# Default goal and phony targets
# ----------------------------------------------------------------------------
.DEFAULT_GOAL = seren
.PHONY : clean check_errors


# Error target
# ----------------------------------------------------------------------------
check_errors :
ifdef ERROR
	@printf "COMPILATION FAILED : Error(s) in Makefile\n"
	@for ERR in $(ERROR); do \
            echo $$ERR ; \
        done
	@exit 2
endif

# Code compilation
# ============================================================================
OPT += $(foreach INCDIR, $(INCLUDE_DIR), -I $(INCDIR))

%.o: %.F90 definitions.o HP_types.o modules.o interface.o $(MPI_REQUIRED)
	$(F90) $(OPT) $(CFLAGS) -c $<

SRC_OBJ = $(OBJ:.o=.F90)


# ----------------------------------------------------------------------------
seren : check_errors $(OBJ) seren.o
ifeq ($(MPI),1)
	$(F90) $(OPT) $(CFLAGS) $(LDFLAGS) -o $(EXEDIR)/seren-mpi $(OBJ) seren.o $(STATIC_LIBS)
else
	$(F90) $(OPT) $(CFLAGS) $(LDFLAGS) -o $(EXEDIR)/seren $(OBJ) seren.o $(STATIC_LIBS)
endif


convert_format :: $(OBJ) convert_format.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/convert_format $(OBJ) convert_format.o

error_norm :: $(OBJ) error_norm.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/error_norm $(OBJ) error_norm.o

ic_BB :: check_errors $(OBJ) $(IC_OBJ) ic_BB.o 
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_BB $(OBJ) $(IC_OBJ) ic_BB.o

ic_binary :: check_errors $(OBJ) ic_binary.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_binary $(OBJ) ic_binary.o

ic_binform :: check_errors $(OBJ) ic_binform.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_binform $(OBJ) ic_binform.o

ic_blob :: check_errors $(OBJ) $(IC_OBJ) ic_blob.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_blob $(OBJ) $(IC_OBJ) ic_blob.o

ic_bondi :: check_errors $(OBJ) $(IC_OBJ) Bondi.o ic_bondi.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_bondi $(OBJ) $(IC_OBJ) Bondi.o ic_bondi.o

ic_bondi_hoyle :: check_errors $(OBJ) $(IC_OBJ) ic_bondi_hoyle.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_bondi_hoyle $(OBJ) $(IC_OBJ) ic_bondi_hoyle.o

ic_core :: check_errors $(OBJ) $(IC_OBJ) ic_core.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_core $(OBJ) $(IC_OBJ) ic_core.o

ic_deflection_cube :: check_errors $(OBJ) $(IC_OBJ) ic_deflection_cube.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_deflection_cube $(OBJ) $(IC_OBJ) ic_deflection_cube.o

ic_entropy_test :: check_errors $(OBJ) $(IC_OBJ) ic_entropy_test.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_entropy_test $(OBJ) $(IC_OBJ) ic_entropy_test.o

ic_hcp :: check_errors $(OBJ) ic_hcp.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_hcp $(OBJ) ic_hcp.o

ic_KH :: check_errors $(OBJ) $(IC_OBJ) ic_KH.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_KH $(OBJ) $(IC_OBJ) ic_KH.o

ic_jeans :: check_errors $(OBJ) ic_jeans.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_jeans $(OBJ) ic_jeans.o

ic_lattice_cube :: check_errors $(OBJ) ic_lattice_cube.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_lattice_cube $(OBJ) ic_lattice_cube.o

ic_mkdynfric :: check_errors $(OBJ) ic_mkdynfric.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_mkdynfric $(OBJ) ic_mkdynfric.o

ic_NTSI :: check_errors $(OBJ) smoothed_velocity.o ic_NTSI.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_NTSI $(OBJ) smoothed_velocity.o ic_NTSI.o

ic_polytrope :: check_errors $(OBJ) $(IC_OBJ) ic_polytrope.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_polytrope $(OBJ) $(IC_OBJ) ic_polytrope.o

ic_plummer :: check_errors $(OBJ) $(IC_OBJ) ic_plummer.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_plummer $(OBJ) $(IC_OBJ) ic_plummer.o

ic_rad_core :: check_errors $(OBJ) ic_radial_core.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_rad_core $(OBJ) ic_radial_core.o

ic_radtest :: check_errors $(OBJ) ic_radtest.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_radtest $(OBJ) ic_radtest.o

ic_random_cube :: check_errors $(OBJ) ic_random_cube.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_random_cube $(OBJ) ic_random_cube.o

ic_replicate_cubes :: check_errors $(OBJ) ic_replicate_cubes.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_replicate_cubes $(OBJ) ic_replicate_cubes.o

ic_RT :: check_errors $(OBJ) $(IC_OBJ) ic_RT.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_RT $(OBJ) $(IC_OBJ) ic_RT.o

ic_sedov :: check_errors $(OBJ) ic_sedov.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_sedov $(OBJ) ic_sedov.o

ic_shear_flow :: check_errors $(OBJ) $(IC_OBJ) ic_shear_flow.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_shear_flow $(OBJ) $(IC_OBJ) ic_shear_flow.o

ic_shocktube :: check_errors $(OBJ) smoothed_velocity.o ic_shocktube.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_shocktube $(OBJ) smoothed_velocity.o ic_shocktube.o

ic_SIS :: check_errors $(OBJ) $(IC_OBJ) ic_SIS.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_SIS $(OBJ) $(IC_OBJ) ic_SIS.o

ic_sphere :: check_errors $(OBJ) ic_sphere.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_sphere $(OBJ) ic_sphere.o

ic_spitzer :: check_errors $(OBJ) ic_spitzer.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_spitzer $(OBJ) ic_spitzer.o

ic_vel_pert :: check_errors $(OBJ) $(IC_OBJ) ic_vel_pert.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_vel_pert $(OBJ) $(IC_OBJ) ic_vel_pert.o

gravtest :: check_errors $(OBJ) direct_sph_gravity.o gravtest.o
	$(F90) $(OPT) $(CFLAGS) -o gravtest $(OBJ) direct_sph_gravity.o gravtest.o

lane_emden :: definitions.o modules.o lane_emden.o
	$(F90) $(OPT) $(CFLAGS) -o lane_emden modules.o lane_emden.o 

nbody_orbits :: definitions.o healpix_types.o modules.o nbody_orbits.o
	$(F90) $(OPT) $(CFLAGS) $(X11LIBS) $(PGPLOTLIBS) -o nbody_orbits definitions.o healpix_types.o modules.o nbody_orbits.o 

radial_average :: check_errors $(OBJ) radial_average.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/radial_average $(OBJ) radial_average.o

ic_subdisk :: check_errors $(OBJ) ic_subdisk.o
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/ic_subdisk $(OBJ) ic_subdisk.o

joinsims ::
	$(F90) $(OPT) $(CFLAGS) -c joinsim/systemargs.F90
	$(F90) $(OPT) $(CFLAGS) -c joinsim/datamodule.F90
	$(F90) $(OPT) $(CFLAGS) -c joinsim/progmodule.F90
	$(F90) $(OPT) $(CFLAGS) -o $(EXEDIR)/joinsims systemargs.o datamodule.o progmodule.o joinsim/join.F90

#joinsims.o : systemargs.o datamodule.o progmodule.o
#progmodule.o : datamodule.o


# User makefile additions
# ----------------------------------------------------------------------------
#include user/user_makefiletail.mk


clean :: 
	\rm -f *.mod
	\rm -f *.o
	\rm -f *__genmod.f90


# Dependencies
# ----------------------------------------------------------------------------
definitions.o : definitions.F90 Makefile macros.h
	$(F90) $(OPT) $(CFLAGS) -c $<

HP_types.o : HP_types.F90
	$(F90) $(OPT) $(CFLAGS) -c $<

modules.o : modules.F90
	$(F90) $(OPT) $(CFLAGS) -c $<

interface.o : interface.F90
	$(F90) $(OPT) $(CFLAGS) -c $<

mpi_modules.o : mpi_modules.F90
	$(F90) $(OPT) $(CFLAGS) -c $<

interface.o : definitions.o HP_types.o modules.o
modules.o : definitions.o HP_types.o
HP_types.o : definitions.o
BH_build_skeleton.o : debug_tests.o
BH_add_particles.o : BH_add_particles.h
turb_force_init.o turb_next_field.o : turb_force_module.o

# MPI dependencies
# ----------------------------------------------------------------------------
ifeq ($(MPI),1)
expanddomainboxes.o nonoverlapbox.o broadcastboundingboxes.o \
mpi_sph_update.o mpi_scatter_ghosts.o mpi_h_gather.o \
transfer_particles.o debug_tests.o get_ghosts.o : domain_comparisons.o

mpi_loadbalance_step.o : debug_tests.o

mpi_modules.o : modules.o

joinsims.o : systemargs.o datamodule.o progmodule.o
progmodule.o : datamodule.o

endif
