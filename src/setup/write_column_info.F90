! WRITE_DEBUG_COLUMN_INFO.F90
! D. A. Hubber - 29/6/2010
! Writes list of column information contained in debug files (if using 
! DEBUG_PLOT_DATA or DEBUG_TRACK_PARTICLE debug flags in Makefile) 
! to file or screen (depending on value of unitno).
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_debug_column_info(unitno)
  use definitions
  implicit none

  integer, intent(in) :: unitno       ! Unit no. to write information to

#if defined(DEBUG_PLOT_DATA) || defined(DEBUG_TRACK_PARTICLE)
  integer :: ncolumns                 ! number of elements in alldata
  integer :: k                        ! dimension counter

! Debug plot data column index 
! ----------------------------------------------------------------------------
  write(unitno,*) "=============================="
  write(unitno,*) "DEBUG PLOT DATA/TRACK PARTICLE"
  write(unitno,*) "=============================="

  ncolumns = 0
  write(unitno,*) "Column  -> Variable name  : Description of variable"

  do k=1,NDIM
     call add_vector_column_info(unitno,ncolumns,k,"r","position")
  end do
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"v","velocity")
  end do
#if defined(SMOOTHED_VELOCITY)
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"v_smooth","smoothed velocity")
  end do
#endif
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"a","acceleration")
  end do
#if defined(DEBUG_FORCES)
#if defined(HYDRO)
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"a_hydro","hydro accel.")
  end do
#endif
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"a_grav","grav. accel.")
  end do
#endif
#if defined(ARTIFICIAL_VISCOSITY)
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"a_visc","viscous accel.")
  end do
#endif
#if defined(STELLAR_WIND)
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"a_wind","wind accel.")
  end do
#endif
#endif
  call add_scalar_column_info(unitno,ncolumns,"h","smoothing length")
  call add_scalar_column_info(unitno,ncolumns,"rho","density")
  call add_scalar_column_info(unitno,ncolumns,"div_v","velocity divergence")
  call add_scalar_column_info(unitno,ncolumns,"dt","ideal timestep")
#if defined(HYDRO)
  call add_scalar_column_info(unitno,ncolumns,"temp","temperature")
  call add_scalar_column_info(unitno,ncolumns,"sound","sound speed")
  call add_scalar_column_info(unitno,ncolumns,"press","thermal pressure")
  call add_scalar_column_info(unitno,ncolumns,"u","specific internal energy")
#if defined(ENERGY_EQN)
  call add_scalar_column_info(unitno,ncolumns,"du_dt","compressional heating rate")
#endif
#if defined(ENTROPIC_FUNCTION)
  call add_scalar_column_info(unitno,ncolumns,"Aent","entropic function")
  call add_scalar_column_info(unitno,ncolumns,"dA_dt","rate of change of A")
#endif
#endif
  call add_scalar_column_info(unitno,ncolumns,"r_rad","radial distance from origin")
  call add_scalar_column_info(unitno,ncolumns,"v_rad","radial velocity")
  call add_scalar_column_info(unitno,ncolumns,"a_rad","radial acceleration")
#if defined(DEBUG_FORCES)
#if defined(HYDRO)
  call add_scalar_column_info(unitno,ncolumns,"ahydro_rad","radial hydro acceleration")
#endif
#if defined(GRAVITY) || defined(EXTERNAL_FORCE)
  call add_scalar_column_info(unitno,ncolumns,"agrav_rad","radial grav. acceleration")
  call add_scalar_column_info(unitno,ncolumns,"gpe","grav. pot. energy")
#endif
#if defined(ARTIFICIAL_VISCOSITY)
  call add_scalar_column_info(unitno,ncolumns,"avisc_rad","radial viscous acceleration")
  call add_scalar_column_info(unitno,ncolumns,"avisc_vel","velocity viscous acceleration")
#endif
#if defined(STELLAR_WIND)
  call add_scalar_column_info(unitno,ncolumns,"awind_rad","radial wind acceleration")
#endif
#endif
#if defined(GRAD_H_SPH)
  call add_scalar_column_info(unitno,ncolumns,"omega","grad-h correction term")
  call add_scalar_column_info(unitno,ncolumns,"hgradh/hp","grad-h correction check")
#if defined(SELF_GRAVITY)
  call add_scalar_column_info(unitno,ncolumns,"zeta","grad-h gravity correction term")
#endif
#endif
  call add_scalar_column_info(unitno,ncolumns,"m","mass")
  call add_scalar_column_info(unitno,ncolumns,"ptype","particle type")
#if defined(RAD_WS) && defined(DEBUG_RAD)
  do k=1,10
     call add_vector_column_info(unitno,ncolumns,k,"rad_info","radiative cooling info")
  end do
#endif
#if defined(DIFFUSION)
  call add_scalar_column_info(unitno,ncolumns,"ueq","radial accel.")
  call add_scalar_column_info(unitno,ncolumns,"dt_therm","tangential accel.")
  call add_scalar_column_info(unitno,ncolumns,"du_dt_diff","radial accel.")
  call add_scalar_column_info(unitno,ncolumns,"k_cond","tangential accel.")
#endif
#if defined(DEBUG_SINK_BOUNDARY_PROPERTIES)
  call add_scalar_column_info(unitno,ncolumns,"a_rad","radial accel.")
  call add_scalar_column_info(unitno,ncolumns,"a_tang","tangential accel.")
#endif
#if defined(VISC_TD)
  call add_scalar_column_info(unitno,ncolumns,"talpha","time-dependent alpha")
#endif
#if defined(VISC_BALSARA)
  call add_scalar_column_info(unitno,ncolumns,"balsara","balsara switch factor")
#endif
#if defined(VISC_PATTERN_REC)
  call add_scalar_column_info(unitno,ncolumns,"pattrec","pattern recognition factor")
#endif
#if defined(DIV_A)
  call add_scalar_column_info(unitno,ncolumns,"div_a","acceleration divergence")
#endif
#if defined(SIGNAL_VELOCITY)
  call add_scalar_column_info(unitno,ncolumns,"vsigmax","maximum signal velocity")
#endif
#if defined(SM2012_SPH)
  call add_scalar_column_info(unitno,ncolumns,"q","internal energy density")
#endif
#if defined(DEBUG_HP_WALK_ALL_RAYS)
  call add_scalar_column_info(unitno,ncolumns,"whichHPlevel","HEALPix level")
#endif
  call add_scalar_column_info(unitno,ncolumns,"time","simulation time")
#endif

  return
END SUBROUTINE write_debug_column_info



! ============================================================================
SUBROUTINE write_sink_column_info(unitno)
  use definitions
  implicit none

  integer, intent(in) :: unitno       ! Unit no. to write information to

#if defined(SINKS)
  integer :: ncolumns                 ! number of elements in alldata
  integer :: k                        ! dimension counter

! Sink data column index 
! ----------------------------------------------------------------------------
  write(unitno,*) "======================="
  write(unitno,*) "SINK OUTPUT COLUMN INFO"
  write(unitno,*) "======================="

  ncolumns = 0
  write(unitno,*) "Column  -> Variable name  : Description of variable"

  call add_scalar_column_info(unitno,ncolumns,"nsteps","total step no. at formation")
  call add_scalar_column_info(unitno,ncolumns,"time","sink formation time")
  do k=1,NDIM
     call add_vector_column_info(unitno,ncolumns,k,"r","position")
  end do
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"v","velocity")
  end do
  call add_scalar_column_info(unitno,ncolumns,"m","mass")
  call add_scalar_column_info(unitno,ncolumns,"h","smoothing length")
  call add_scalar_column_info(unitno,ncolumns,"radius","accretion radius")
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"a","acceleration")
  end do
  do k=1,3
     call add_vector_column_info(unitno,ncolumns,k,"angmom","internal angular momentum")
  end do
  call add_scalar_column_info(unitno,ncolumns,"gpe","grav. potential energy")
  call add_scalar_column_info(unitno,ncolumns,"dmdt","accretion rate")
  call add_scalar_column_info(unitno,ncolumns,"star_radius","protostar radius")
  call add_scalar_column_info(unitno,ncolumns,"luminosity","protostar luminosity")
  call add_scalar_column_info(unitno,ncolumns,"temperature","protostar surface temperature")
  call add_scalar_column_info(unitno,ncolumns,"menc","enclosed SPH mass")
  call add_scalar_column_info(unitno,ncolumns,"taccrete","sink accretion timescale")
#endif

  return
END SUBROUTINE write_sink_column_info



! ============================================================================
SUBROUTINE write_diagnostic_column_info(unitno)
  use definitions
  implicit none

  integer, intent(in) :: unitno       ! Unit no. to write information to

#if defined(DEBUG_DIAGNOSTICS)
  integer :: ncolumns                 ! number of elements in alldata
  integer :: k                        ! dimension counter

! Sink data column index 
! ----------------------------------------------------------------------------
  write(unitno,*) "============================="
  write(unitno,*) "DIAGNOSTIC OUTPUT COLUMN INFO"
  write(unitno,*) "============================="

  ncolumns = 0
  write(unitno,*) "Column  -> Variable name  : Description of variable"

  call add_scalar_column_info(unitno,ncolumns,"nsteps","total no. of steps")
  call add_scalar_column_info(unitno,ncolumns,"time","simulation time")
  call add_scalar_column_info(unitno,ncolumns,"mtot","total mass")
  do k=1,NDIM
     call add_vector_column_info(unitno,ncolumns,k,"rcom","position of COM")
  end do
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"v","velocity of COM")
  end do
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"mom","net momentum")
  end do
  do k=1,3
     call add_vector_column_info(unitno,ncolumns,k,"angmom","net angular momentum")
  end do
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"force","net force")
  end do
#if defined(DEBUG_FORCES)
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"force_grav","net grav.force")
  end do
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"force_hydro","net hydro force")
  end do
  do k=1,VDIM
     call add_vector_column_info(unitno,ncolumns,k,"force_mag","net magnetic force")
  end do
#endif
  call add_scalar_column_info(unitno,ncolumns,"utot","net internal energy")
  call add_scalar_column_info(unitno,ncolumns,"gpetot","net grav. pot. energy")
  call add_scalar_column_info(unitno,ncolumns,"ketot","net kinetic energy")
  call add_scalar_column_info(unitno,ncolumns,"etot","net total energy")
  call add_scalar_column_info(unitno,ncolumns,"Atot","net entropy")
  call add_scalar_column_info(unitno,ncolumns,"error","net energy error")
#endif

  return
END SUBROUTINE write_diagnostic_column_info



! ============================================================================ 
SUBROUTINE add_scalar_column_info(unitno,ncolumns,var_name,message)
  implicit none

  integer, intent(in) :: unitno                ! File unit no.
  integer, intent(inout) :: ncolumns           ! No. of columns
  character(len=*), intent(in) :: var_name     ! Variable name string
  character(len=*), intent(in) :: message      ! Description of variable

  ncolumns = ncolumns + 1
  write(unitno,'(I8,A,A,A,A)') &
       &ncolumns," -> ",trim(var_name),"   : ",trim(message)

  return
END SUBROUTINE add_scalar_column_info



! ============================================================================ 
SUBROUTINE add_vector_column_info(unitno,ncolumns,component,var_name,message)
  implicit none

  integer, intent(in) :: unitno                  ! File unit no.
  integer, intent(inout) :: ncolumns             ! No. of columns
  integer, intent(in) :: component               ! Component of vector
  character(len=*), intent(in) :: var_name       ! Variable name string
  character(len=*), intent(in) :: message        ! Description of variable

  ncolumns = ncolumns + 1
  write(unitno,'(I8,A,A,A,I2,A,A)') &
       &ncolumns," -> ",trim(var_name),"(",component,")   : ",trim(message)

  return
END SUBROUTINE add_vector_column_info
