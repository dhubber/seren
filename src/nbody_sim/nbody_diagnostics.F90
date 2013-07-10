! NBODY_DIAGNOSTICS.F90
! D. A. Hubber - 15/9/20120
! ..
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE nbody_diagnostics
  use seren_sim_module
  use particle_module
  use scaling_module
  use type_module
  use hydro_module
  use sink_module
  use Nbody_module
#if defined(DEBUG_DIAGNOSTICS)
  use time_module, only : time
  use filename_module, only : run_dir,run_id
#endif
  implicit none

  integer       :: s                   ! sink counter
  real(kind=DP) :: ap(1:VDIM)          ! acceleration of particle p
  real(kind=DP) :: angmom(1:3)         ! angular momentum
  real(kind=DP) :: etot                ! total energy
  real(kind=DP) :: force(1:VDIM)       ! net force on all particles
  real(kind=DP) :: gpetot              ! total grav. potential energy
  real(kind=DP) :: ketot               ! total kinetic energy
  real(kind=DP) :: mom(1:VDIM)         ! linear momentum
  real(kind=DP) :: mp                  ! mass of particle p
  real(kind=DP) :: mtot_temp           ! total mass
  real(kind=DP) :: rcom_temp(1:NDIM)   ! position of centre of mass
  real(kind=DP) :: rp(1:NDIM)          ! position of particle p
  real(kind=DP) :: vcom_temp(1:VDIM)   ! velocity of centre of mass
  real(kind=DP) :: vp(1:VDIM)          ! velocity of particle p
  real(kind=DP) :: angmomp(1:3)        ! angular momentum of p
  real(kind=DP) :: gpep                ! gpe of particle p
#if defined(DEBUG_DIAGNOSTICS)
  logical :: ex                        ! Does file exist already?
  character(len=256) :: out_file       ! Name of outputted debug file
  real(kind=DP) :: ttemp               ! ..
#endif
#if defined(DEBUG_FORCES)
  real(kind=DP) :: force_grav(1:VDIM)  ! net grav. force
#endif

  debug2("Calculating total energy and momentum [diagnostics.F90]")

! Zero all accumulation variables
! ----------------------------------------------------------------------------
  mtot_temp   = 0.0_DP
  etot   = 0.0_DP
  ketot  = 0.0_DP
  gpetot = 0.0_DP
  force(1:VDIM) = 0.0_DP
  mom(1:VDIM)   = 0.0_DP
  rcom_temp(1:NDIM)  = 0.0_DP
  vcom_temp(1:VDIM)  = 0.0_DP
  angmom(1:3)   = 0.0_DP
#if defined(DEBUG_FORCES)
  force_grav(1:VDIM)  = 0.0_DP
#endif


! Sum up the contributions from all sinks
! ----------------------------------------------------------------------------
  do s=1,stot

     mp = real(star(s)%m,DP)
     rp(1:NDIM) = star(s)%r(1:NDIM)
     vp(1:VDIM) = star(s)%v(1:VDIM)
     ap(1:VDIM) = star(s)%a(1:VDIM)
     angmomp(1:3) = star(s)%angmom(1:3)
     gpep = star(s)%gpe

     ! Sum up the contributions to net force, momentum, energy and mass
     mtot_temp = mtot_temp + mp
     mom(1:VDIM)   = mom(1:VDIM)   + mp*vp(1:VDIM)
     rcom_temp(1:NDIM)  = rcom_temp(1:NDIM)  + mp*rp(1:NDIM)
     vcom_temp(1:VDIM)  = vcom_temp(1:VDIM)  + mp*vp(1:VDIM)
     force(1:VDIM) = force(1:VDIM) + mp*ap(1:VDIM)
     angmom(1:3)   = angmom(1:3)   + angmomp(1:3)
     
     ! Now add orbital angular momentum of sink COM
#if NDIM==3
     angmom(1) = angmom(1) + mp*(rp(2)*vp(3) - rp(3)*vp(2))
     angmom(2) = angmom(2) + mp*(rp(3)*vp(1) - rp(1)*vp(3))
     angmom(3) = angmom(3) + mp*(rp(1)*vp(2) - rp(2)*vp(1))
#elif NDIM==2
     angmom(3) = angmom(3) + mp*(rp(1)*vp(2) - rp(2)*vp(1))
#endif
     ketot  = ketot + mp*dot_product(vp(1:VDIM),vp(1:VDIM))
     gpetot = gpetot + gpep
#if defined(DEBUG_FORCES)
     force_grav(1:VDIM) = force_grav(1:VDIM) + mp*sink(s)%agrav(1:VDIM)
#endif
  end do

! Normalise centre of mass variables
  rcom_temp(1:NDIM) = rcom_temp(1:NDIM) / mtot_temp
  vcom_temp(1:VDIM) = vcom_temp(1:VDIM) / mtot_temp


! Account for constant multipliers, double summation (gravity),
! degrees of freedom (thermal) etc.
! ----------------------------------------------------------------------------
  ketot  = 0.5_DP*ketot
  gpetot = 0.5_DP*gpetot
  etot   = ketot - gpetot


! Output to screen
! ----------------------------------------------------------------------------
5 format(1X,A27,1I10)
10 format(1X,A8,2X,1G18.10)
15 format(1X,A8,3X,1G18.10,3X,1G18.10,3X,1G18.10)
#if NDIM==1
20 format(1X,A8,3X,1G18.10)
#elif NDIM==2
20 format(1X,A8,3X,1G18.10,3X,1G18.10)
#elif NDIM==3
20 format(1X,A8,3X,1G18.10,3X,1G18.10,3X,1G18.10)
#endif

  write(6,5) "Total number of particles : ", ptot
  write(6,5) "Total number of stars     : ", stot
  write(6,10) "mtot_temp   :", mtot_temp*mscale
  write(6,20) "rcom_temp   :", rcom_temp(1:NDIM)*rscale
  write(6,20) "vcom_temp   :", vcom_temp(1:VDIM)*vscale
  write(6,20) "mom    :", mom(1:VDIM)*momscale
  write(6,15) "ang    :", angmom(1:3)*angmomscale
  write(6,20) "force  :", force(1:VDIM)
#if defined(DEBUG_FORCES) && defined(GRAVITY)
  write(6,20) "forceG :", force_grav(1:VDIM)
#endif
  write(6,10) "etot   :", etot*Escale
  write(6,10) "ketot  :", ketot*Escale
  write(6,10) "gpetot :", -gpetot*Escale


! Write to file
! ----------------------------------------------------------------------------
#if defined(DEBUG_DIAGNOSTICS)
  out_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".nbodydiag"

  ! Check if file exists
  inquire(file=out_file,exist=ex)

  ! Read through file to synchronise files with time of simulation.
  if (ex) then
     open(1,file=out_file,status="unknown",&
          &form="formatted",position="rewind")
     do
        read(1,'(1E18.10)',end=50,err=50) ttemp
        if (ttemp > time*tscale) exit
     end do
50   backspace (1,err=100)
  else
     open(1,file=out_file,status="unknown",form="formatted")
  end if
#if defined(DEBUG_DIAGNOSTICS) && defined(DEBUG_FORCES)
100 write(1,'(60E18.10)') time*tscale, mtot_temp*mscale,&
         &rcom_temp(1:NDIM)*rscale,vcom_temp(1:VDIM)*vscale,&
         &mom(1:VDIM)*momscale, angmom(1:3)*angmomscale,force(1:VDIM),&
         &force_grav(1:VDIM),gpetot*Escale,ketot*Escale,etot*Escale
  close(1)
#else
100 write(1,'(60E18.10)') time*tscale,mtot_temp*mscale,&
         &rcom_temp(1:NDIM)*rscale,vcom_temp(1:VDIM)*vscale,&
         &mom(1:VDIM)*momscale, angmom(1:3)*angmomscale,&
         &force(1:VDIM),gpetot*Escale,ketot*Escale,etot*Escale
  close(1)
#endif
#endif

  return
END SUBROUTINE nbody_diagnostics
