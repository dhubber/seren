! CALCULATE_SINK_PROPERTIES.F90
! D. A. Hubber - 30/11/2009
! Calculate the stellar, and stellar evolution properties of all sinks.
! For now, compares the sink mass to a look up table (c.f. Crowther et al.) 
! and interpolates the number of UV ionizing photons per second.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE calculate_sink_properties
  use sink_module
  use scaling_module
  use constant_module
  use HP_module
  use stellar_module
  implicit none

  logical :: active                ! Is sink an active feedback source?
  integer :: HPid                  ! HP source id
  integer :: itable                ! ..
  integer :: s                     ! Sink counter
  real(kind=PR) :: intfactor       ! ..
  real(kind=DP) :: intmax_s        ! Maximum ionization integral
  real(kind=DP) :: L_s             ! Total effective luminosity
  real(kind=PR) :: ms              ! Mass of sink s
  real(kind=PR) :: M_loss_s        ! ..
  real(kind=DP) :: N_LyC_s         ! No. of ionizing photons per second
  real(kind=DP) :: radius_s        ! Stellar radius
  real(kind=DP) :: v_wind_s        ! ..

  debug2("Calculate all feedback properties of stars/sinks [calculate_sink_properties.F90]")


! Loop over all sinks
! ============================================================================
  do s=1,stot
     ms = sink(s)%m
     active = .false.
        
     ! Find where sink lies in look-up table for interpolation
     itable = 0
     do
        itable = itable + 1
        if (ms < stellar_table(itable + 1)%mass .or. itable == Ntable - 1) exit
     end do
     if (itable >= Ntable) itable = Ntable - 1
     intfactor = (ms - stellar_table(itable)%mass) / &
          (stellar_table(itable + 1)%mass - stellar_table(itable)%mass)
     if (intfactor > 1.0_PR) intfactor = 1.0_PR
     
     ! Calculated various stellar properties from interpolation of
     ! the look-up table
#if defined(IONIZING_UV_RADIATION)
     N_LyC_s = (1.0_PR - intfactor)*stellar_table(itable)%log_N_LyC &
          & + intfactor*stellar_table(itable + 1)%log_N_LyC
     N_LyC_s = 10.0_DP**(real(N_LyC_s,DP))
     N_LyC_s = N_LyC_s*(tscale*t_SI)
     intmax_s = N_LyC_s*(m_hydrogen/(mscale*m_SI*Xfrac))**2 / &
          &(4.0_PR*PI*a_star)
     if (intmax_s > SMALL_NUMBER) active = .true.
#endif
#if defined(STELLAR_WIND)
     M_loss_s = (1.0_PR - intfactor)*stellar_table(itable)%M_loss &
          & + intfactor*stellar_table(itable + 1)%M_loss
     v_wind_s = (1.0_PR - intfactor)*stellar_table(itable)%v_wind &
          & + intfactor*stellar_table(itable + 1)%v_wind
     if (M_loss_s*v_wind_s > SMALL_NUMBER) active = .true.
#endif
#if defined(STELLAR_LUMINOSITY)
     L_s = (1.0_PR - intfactor)*stellar_table(itable)%log_L &
          & + intfactor*stellar_table(itable + 1)%log_L
     L_s = 10.0_DP**(real(L_s,DP))
     if (L_s > SMALL_NUMBER) active = .true.
#endif
     write(6,*) "Star properties : ",s,itable,intfactor,ms
     !write(6,*) "Luminosity : ",L_s*Lscale,Lunit
     write(6,*) "Ionizing rad    : ",N_LyC_s/tscale/t_SI,intmax_s
     write(6,*) "Wind            : ",M_loss_s*dmdtscale,v_wind_s*vscale
     
     ! If there is some ionizing flux or a strong wind source from the star, 
     ! create a new healpix source from the sink (if required) and 
     ! record all stellar properties for source.
     if (ms > 0.0_PR .and. active) then
        HPid = sink(s)%HPid
        
        ! Create a new source if it doesn't exist
        if (HPid == -1) then
           call create_HP_source(s,sink(s)%r(1:NDIM))
           HPid = HPtot
        end if
 
     end if

#if defined(IONIZING_UV_RADIATION)
     HPsource(HPid)%N_LyC  = N_LyC_s
     HPsource(HPid)%intmax = intmax_s
#endif
#if defined(STELLAR_WINDS)
     HPsource(HPid)%M_loss = Mloss_s
     HPsource(HPid)%v_wind = v_wind_s
#endif
#if defined(STELLAR_LUMINOSITY)
     HPsource(HPid)%L_star = L_s
#endif
     
  end do
! ============================================================================


  return
END SUBROUTINE calculate_sink_properties
