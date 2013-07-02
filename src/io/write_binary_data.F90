! WRITE_BINARY_DATA.F90
! D. A. Hubber - ..
! ..
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_binary_data
  use particle_module
  use filename_module
  use scaling_module
  use time_module

  character(len=40) :: out_file     ! Name of outputted debug file
  logical :: ex                     ! Does file exist already?
  integer :: p                      ! ..
  integer :: pfirst                 ! ..
  integer :: plast                  ! ..
  integer :: s                      ! ..
  real(kind=DP) :: dr(1:NDIM)       ! Relative displacement vector
  real(kind=DP) :: drmag            ! Distance
  real(kind=DP) :: drsqd            ! Distance squared
  real(kind=DP) :: dv(1:NDIM)       ! Relative velocity
  real(kind=DP) :: eccent           ! Eccentricity
  real(kind=DP) :: L(1:3)           ! Angular momentum vector
  real(kind=DP) :: Lsqd             ! ..
  real(kind=DP) :: period           ! Period
  real(kind=DP) :: q                ! Mass-ratio
  real(kind=DP) :: rbin(1:NDIM)     ! Position of centre of mass
  real(kind=DP) :: reduced_mass     ! Reduced mass
  real(kind=DP) :: sma              ! Semi-major axis
  real(kind=DP) :: ttemp            ! ..
  real(kind=DP) :: vbin(1:NDIM)     ! Velocity of centre of mass

  type object_node                  ! Stellar object node
     real(kind=DP) :: m             ! Mass
     real(kind=DP) :: rcom(1:NDIM)  ! Position of COM
     real(kind=DP) :: vcom(1:VDIM)  ! Velocity of COM
     real(kind=DP) :: angmom(1:3)   ! Internal angular momentum
     real(kind=DP) :: rmax          ! Maximum radius of particle from COM
  end type object_node
  type(object_node) :: objdata(1:2) ! ..



! ============================================================================
  do s=1,2

     if (s == 1) then
        pfirst = 1
        plast  = ptot/2
     else
        pfirst = ptot/2 + 1
        plast  = ptot
     end if

     ! Zero all variables
     objdata(s)%m            = 0.0_DP
     objdata(s)%rcom(1:NDIM) = 0.0_DP
     objdata(s)%vcom(1:VDIM) = 0.0_DP
     objdata(s)%angmom(1:3)  = 0.0_DP
     objdata(s)%rmax         = 0.0_DP


     ! Find centre of mass of object
     do p=1,ptot
        if (sph(p)%porig < pfirst .or. sph(p)%porig > plast) cycle
        objdata(s)%m = objdata(s)%m + real(sph(p)%m,DP)
        objdata(s)%rcom(1:NDIM) = objdata(s)%rcom(1:NDIM) + &
             & real(sph(p)%m,DP)*real(sph(p)%r(1:NDIM),DP)
        objdata(s)%vcom(1:NDIM) = objdata(s)%vcom(1:NDIM) + &
             & real(sph(p)%m,DP)*real(sph(p)%v(1:NDIM),DP)
     end do
     objdata(s)%rcom(1:NDIM) = objdata(s)%rcom(1:NDIM) / objdata(s)%m
     objdata(s)%vcom(1:VDIM) = objdata(s)%vcom(1:VDIM) / objdata(s)%m


     ! Now find other properties of object
     do p=1,ptot
        if (sph(p)%porig < pfirst .or. sph(p)%porig > plast) cycle
        call distance2_dp(objdata(s)%rcom(1:NDIM),p,dr(1:NDIM),drsqd)
        dv(1:VDIM) = real(sph(p)%v(1:VDIM),DP) - objdata(s)%vcom(1:VDIM)
#if NDIM==3
        objdata(s)%angmom(1) = objdata(s)%angmom(1) + &
             (dr(2)*dv(3) - dr(3)*dv(2))*real(sph(p)%m,DP)
        objdata(s)%angmom(2) = objdata(s)%angmom(2) + &
             (dr(3)*dv(1) - dr(1)*dv(3))*real(sph(p)%m,DP)
        objdata(s)%angmom(3) = objdata(s)%angmom(3) + &
             (dr(1)*dv(2) - dr(2)*dv(1))*real(sph(p)%m,DP)
#elif NDIM == 2
        objdata(s)%angmom(3) = objdata(s)%angmom(3) + &
             (dr(1)*dv(2) - dr(2)*dv(1))*real(sph(p)%m,DP)
#endif
     end do

  end do
! ============================================================================


! Now determine all system properties
  rbin(1:NDIM) = (objdata(1)%m*objdata(1)%rcom(1:NDIM) + &
       &objdata(2)%m*objdata(2)%rcom(1:NDIM))/(objdata(1)%m + objdata(2)%m)
  vbin(1:NDIM) = (objdata(1)%m*objdata(1)%vcom(1:NDIM) + &
       &objdata(2)%m*objdata(2)%vcom(1:NDIM))/(objdata(1)%m + objdata(2)%m)
  dv(1:NDIM) = objdata(2)%vcom(1:NDIM) - objdata(1)%vcom(1:NDIM)
  call distance3_dp(objdata(1)%rcom(1:NDIM),objdata(2)%rcom(1:NDIM),&
       &dr(1:NDIM),drsqd)
  drmag = sqrt(drsqd)

  reduced_mass = objdata(1)%m*objdata(2)%m / (objdata(1)%m + objdata(2)%m)
#if NDIM==3
  L(1) = (dr(2)*dv(3) - dr(3)*dv(2))*reduced_mass
  L(2) = (dr(3)*dv(1) - dr(1)*dv(3))*reduced_mass
  L(3) = (dr(1)*dv(2) - dr(2)*dv(1))*reduced_mass
  Lsqd = L(1)*L(1) + L(2)*L(2) + L(3)*L(3)
#elif NDIM == 2
  L(3) = (dr(1)*dv(2) - dr(2)*dv(1))*reduced_mass
  Lsqd = L(3)*L(3)
#endif

        
! Now calculate all binary parameters
  binen  = 0.5_DP*reduced_mass*dot_product(dv(1:VDIM),dv(1:VDIM)) - &
       &objdata(1)%m*objdata(2)%m/drmag
  sma    = -0.5_DP*objdata(1)%m*objdata(2)%m/binen
  eccent = 1.0_DP - Lsqd/(objdata(1)%m + objdata(2)%m)/sma/(reduced_mass**2)
  eccent = max(0.0_DP,eccent)
  eccent = sqrt(eccent)
!  write(6,*) "ECCENT :",Lsqd,objdata(1)%m,objdata(2)%m,sma,reduced_mass,eccent
!  write(6,*) "L : ",L(1:3),L(1)*L(1),L(2)*L(2),L(3)*L(3)
  period = 2.0_DP*PI*sqrt(sma**3/(objdata(1)%m + objdata(2)%m))
  if (objdata(1)%m > objdata(2)%m) then
     q = objdata(2)%m/objdata(1)%m
  else
     q = objdata(1)%m/objdata(2)%m
  end if

!  write(6,*) "Binary test : ",reduced_mass,objdata(1)%m,objdata(2)%m,&
!       &drmag,L(1:3),Lsqd,dv(1:VDIM),sqrt(dot_product(dv(1:VDIM),dv(1:VDIM))),&
!       &binen,sma,eccent,period,reduced_mass,&
!       &0.5_DP*reduced_mass*dot_product(dv(1:VDIM),dv(1:VDIM)),&
!       &objdata(1)%m*objdata(2)%m/drmag

! ..
! ----------------------------------------------------------------------------
  out_file = trim(adjustl(run_dir))//trim(adjustl(run_id))//".bintest"

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

100 write(1,'(50E18.10)') time*tscale,objdata(1)%m*mscale,&
         &objdata(1)%rcom(1:NDIM)*rscale,objdata(1)%vcom(1:VDIM)*vscale,&
         &objdata(1)%angmom(1:3)*angmomscale,objdata(1)%rmax*rscale,&
         &objdata(2)%m*mscale,objdata(2)%rcom(1:NDIM)*rscale,&
         &objdata(2)%vcom(1:VDIM)*vscale,objdata(2)%angmom(1:3)*angmomscale,&
         &objdata(2)%rmax*rscale,drmag*rscale,sma*rscale,eccent,period*tscale,&
         &L(1:3)*angmomscale,rbin(1:NDIM)*rscale,vbin(1:NDIM)*vscale,&
         &L(1:3) + objdata(1)%angmom(1:3) + objdata(2)%angmom(1:3)

  close(1)


  return
END SUBROUTINE write_binary_data


