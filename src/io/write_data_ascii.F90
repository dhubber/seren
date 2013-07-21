! WRITE_DATA_ASCII.F90
! D. A. Hubber - 09/06/2010
! Writes out a snapshot file in simple column ASCII format.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_data_ascii(out_file)
  use particle_module
  use hydro_module
  use type_module
  use sink_module
  use time_module
  use scaling_module
  implicit none

  character(len=*), intent(in) :: out_file   ! Name of ascii file

  character(len=20) :: auxdata         ! aux. data id variable
  character(len=20) :: data_id(1:100)  ! list of data identifiers 
  integer :: i                         ! aux. loop variable
  integer :: ndata                     ! no. of data columns
  integer :: p                         ! counter to loop over particles
  integer :: auxtype                     ! particle type
  integer :: s                         ! sink id
  real(kind=PR) :: pdata(1:100)        ! aux. particle data array

  debug1("Writing data file : "//trim(out_file)//" [write_data_ascii.F90]")

! Initialise counter
  ndata = 0

! First, read-in ascii column data
  open(1, file="asciicolumns.dat", status="old", form="formatted")
  do
     read(1,*,err=5,end=5) auxdata
     ndata = ndata + 1
     data_id(ndata) = auxdata
  end do
5 close(5)

! Check there are a sufficient no. of data columns, or ptype is not 1st column
  if (ndata <= 1) then
     write(6,*) "Invalid no. of data columns; Cannot write output file"
     return
stop
  end if
  if (data_id(1) /= "ptype") then
     write(6,*) "First column is not particle type id. Cannot write output file"
     return
stop
  end if

  open(1, file=out_file, status="unknown", form="formatted")


! Now read all particle data into main arrays
! ----------------------------------------------------------------------------
  do p=1,ptot

     if (p <= pboundary) then
        auxtype = 6
     else if (p <= pboundary + picm) then
        auxtype = 9
     else if (p <= pboundary + picm + pgas) then
        auxtype = 1
     else if (p <= pboundary + picm + pgas + pcdm) then
        auxtype = 10
     else 
        cycle
     end if
     pdata(1:100) = 0.0_PR

     do i=2,ndata
        if (data_id(i) == "x") then
           pdata(i) = sph(p)%r(1)*real(rscale,PR)
#if NDIM==2 || NDIM==3
        else if (data_id(i) == "y" .and. NDIM > 1) then
           pdata(i) = sph(p)%r(2)*real(rscale,PR)
#endif
#if NDIM==3
        else if (data_id(i) == "z" .and. NDIM == 3) then
           pdata(i) = sph(p)%r(3)*real(rscale,PR)
#endif
        else if (data_id(i) == "m") then
           pdata(i) = sph(p)%m*real(mscale,PR)
        else if (data_id(i) == "h") then
           pdata(i) = sph(p)%h*real(rscale,PR)
        else if (data_id(i) == "vx") then
           pdata(i) = sph(p)%v(1)*real(vscale,PR)
#if VDIM==2 || VDIM==3
        else if (data_id(i) == "vy" .and. VDIM > 1) then
           pdata(i) = sph(p)%v(2)*real(vscale,PR)
#endif
#if VDIM==3
        else if (data_id(i) == "vz" .and. VDIM == 3) then
           pdata(i) = sph(p)%v(3)*real(vscale,PR)
#endif
        else if (data_id(i) == "rho") then
           pdata(i) = sph(p)%rho*real(rhoscale,PR)
#if defined(INTERNAL_ENERGY)
        else if (data_id(i) == "u") then
           pdata(i) = sph(p)%u*real(uscale,PR)
#endif
#if defined(CHEMCOOL)
        else if (data_id(i) == "abundh2") then
           pdata(i) = sph(p)%abundances(1)
        else if (data_id(i) == "abundhp") then
           pdata(i) = sph(p)%abundances(2)
#endif
        end if
     end do
     write(1,'(I3,100E15.7)') auxtype,pdata(2:ndata)

  end do
! ----------------------------------------------------------------------------


! ----------------------------------------------------------------------------
#if defined(SINKS)
  do s=1,stot

     auxtype = -1
     pdata(1:100) = 0.0_PR
     do i=2,ndata
        if (data_id(i) == "x") then
           pdata(i) = sink(s)%r(1)*real(rscale,PR)
#if NDIM==2 || NDIM==3
        else if (data_id(i) == "y" .and. NDIM > 1) then
           pdata(i) = sink(s)%r(2)*real(rscale,PR)
#endif
#if NDIM==3
        else if (data_id(i) == "z" .and. NDIM == 3) then
           pdata(i) = sink(s)%r(3)*real(rscale,PR)
#endif
        else if (data_id(i) == "m") then
           pdata(i) = sink(s)%m*real(mscale,PR)
        else if (data_id(i) == "h") then
           pdata(i) = sink(s)%h*real(rscale,PR)
        else if (data_id(i) == "vx") then
           pdata(i) = sink(s)%v(1)*real(vscale,PR)
#if VDIM==2 || VDIM==3
        else if (data_id(i) == "vy" .and. VDIM > 1) then
           pdata(i) = sink(s)%v(2)*real(vscale,PR)
#endif
#if VDIM==3
        else if (data_id(i) == "vz" .and. VDIM == 3) then
           pdata(i) = sink(s)%v(3)*real(vscale,PR)
#endif
        end if
     end do
     write(1,'(I3,100E15.7)') auxtype,pdata(2:ndata)

  end do
#endif
! ----------------------------------------------------------------------------

  close(1)

  return
END SUBROUTINE write_data_ascii
