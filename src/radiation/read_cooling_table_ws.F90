! READ_COOLING_TABLE_WS.F90
! D. Stamatellos - 3/1/2007
! Reads the table (dens, temp, energy, mu) required to calculate cooling 
! terms in the radiative cooling routine.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_cooling_table_ws
  use Filename_module, only : eos_opa_file
  use constant_module, only : m_sun,kappa_const
  use hydro_module
  use Eos_module
  use scaling_module
#if defined(SPIEGEL_TEST)
  use Tprof_module
#endif
  implicit none
  
  logical :: ex                     ! Does 'eos.dat' file exist?
  integer :: i                      ! Density grid counter 
  integer :: j                      ! Temperature grid counter
  real(kind=DP) :: auxscale         ! Aux. scaling variable
  real(kind=PR) :: rdummy7(1:7)     ! Dummy array for reading arrays
  character(len=256) :: junkstring  ! Junk string variable

  debug1("Reading in eos tables [read_eos.F90]")

! Auxilary scaling variable for specific energy
  auxscale = (Escale*Ecgs) / (mscale*mcgs)

! Open and read opacity tables (if the file exists)
! ----------------------------------------------------------------------------
  inquire(file=eos_opa_file,exist=ex)
  if (ex) then
     open(2,file=eos_opa_file,status='old',form='formatted')
     
     read(2,*) junkstring
     read(2,*) junkstring
     read(2,*) junkstring
     read(2,*) junkstring

     ! Read-in dimensions of EOS table
     read (2,*) dim_dens, dim_temp, fcolumn

     ! Now allocate array sizes according to table dimensions
     allocate(eos_dens(1:dim_dens))
     allocate(eos_temp(1:dim_temp))
     allocate(eos_energy(dim_dens,dim_temp))
     allocate(eos_mu(dim_dens,dim_temp))  
     allocate(kappa(dim_dens,dim_temp))
     allocate(kappar(dim_dens,dim_temp))
     allocate(kappap(dim_dens,dim_temp))

     ! Read-in data from eos.dat
     do i=1,dim_dens
        do j=1,dim_temp
           read (2,*) rdummy7(1:7)
           eos_dens(i)     = real(rdummy7(1),PR)
           eos_temp(j)     = real(rdummy7(2),PR)
           eos_energy(i,j) = real(rdummy7(3),PR)
           eos_mu(i,j)     = real(rdummy7(4),PR)
           kappa(i,j)      = real(rdummy7(5),PR)
           kappar(i,j)     = real(rdummy7(6),PR)
           kappap(i,j)     = real(rdummy7(7),PR)
        end do
     end do

     close(2)
  else
     stop "eos_opa_file not found"
  end if
! ----------------------------------------------------------------------------

#ifdef SPIEGEL_TEST 

#ifndef SPIEGEL_DISPERSION
  write(*,*) "SPIEGEL_TEST: multiplying opacities so that  tau=",ptemp_q,&
       &"(defined as ptemp_q in params.dat)"
!assumes an input sphere of density 1.41ee-19 g/cm3 
! (from Masunaga & Inutsuka test)
    do i=1,dim_dens
        do j=1,dim_temp
        kappa(i,j)=kappa(i,j)*ptemp_q/4.259e-4
        kappap(i,j)=kappap(i,j)*ptemp_q/4.259e-4
        kappar(i,j)=kappar(i,j)*ptemp_q/4.259e-4
        end do
     end do

#else

write(*,*) "SPIEGEL_TEST: multiplying opacities by 2.952e6*10**(-0.2*ptemp_q)(defined as ptemp_q in params.dat)" 

! Assumes an input sphere of density 1.41ee-19 g/cm3 
! (from Masunaga & Inutsuka test)
   do i=1,dim_dens
        do j=1,dim_temp
         kappa(i,j)=kappa(i,j)*2.952e6*10**(-0.2*ptemp_q)
         kappap(i,j)=kappap(i,j)*2.952e6*10**(-0.2*ptemp_q)
         kappar(i,j)=kappar(i,j)*2.952e6*10**(-0.2*ptemp_q)
        end do
   end do

#endif 

#endif
 
! Convert arrays to code units 
  do i=1,dim_dens
     eos_dens(i) = eos_dens(i) / (rhoscale*rhocgs)
     do j=1,dim_temp
        eos_energy(i,j) = eos_energy(i,j) / auxscale
        kappa(i,j) = kappa(i,j) / (kappascale*kappacgs)
        kappar(i,j) = kappar(i,j) / (kappascale*kappacgs)
        kappap(i,j) = kappap(i,j) / (kappascale*kappacgs)
     end do
  end do

! Calculate log factors for fast table referencing
  densmin = eos_dens(1)
  densmax = eos_dens(dim_dens)
  bdens = (real(dim_dens,PR) - 1.0_PR) / (log10(densmax/densmin))

  tempmin = eos_temp(1)
  tempmax = eos_temp(dim_temp)
  btemp = (real(dim_temp,PR) - 1.0_PR) / (log10(tempmax/tempmin))
  
  return
END SUBROUTINE read_cooling_table_ws
