! READ_STELLAR_MODEL_TABLE.F90
! D. A. Hubber - 11/10/2011
! Read-in stellar model table from external file.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_stellar_model_table
  use stellar_module
  use scaling_module
  implicit none

  logical :: ex           ! Does file exist?
  integer :: i            ! Aux. table counter

  debug2("[read_stellar_model_table.F90]")

! Read in the table with stellar properties
  inquire(file='stellar.dat',exist=ex)

  if (ex) then
     open(unit=2,file='stellar.dat',status='old',form='formatted')
     read(2,*) Ntable
     allocate(stellar_table(1:Ntable))
     read(2,*); read(2,*); read(2,*); read(2,*)
     do i=1,Ntable
        read(2,*) stellar_table(i)%mass,stellar_table(i)%log_L,&
             &stellar_table(i)%log_N_LyC,stellar_table(i)%Teff,&
             &stellar_table(i)%M_loss,stellar_table(i)%v_wind
        stellar_table(i)%mass = stellar_table(i)%mass / (mscale)
        stellar_table(i)%log_L = stellar_table(i)%log_L - log10(Lscale)
        stellar_table(i)%Teff = stellar_table(i)%Teff
        stellar_table(i)%M_loss = 1.0E-6*stellar_table(i)%M_loss / (dmdtscale)
        stellar_table(i)%v_wind = stellar_table(i)%v_wind / (vscale)
     end do
     close(2)
  else
     stop 'stellar.dat file not found'
  end if

  return
END SUBROUTINE read_stellar_model_table
