! TRACK_PARTICLE.F90
! C. P. Batty & D. A. Hubber - 7/3/2007
! Write all data for particle p to formatted (i.e. ASCII) file.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE track_particle(p,rorigin)
  use particle_module
  implicit none

  integer, intent(in)      :: p     ! Particle id
  real(kind=PR), intent(in) :: rorigin(1:NDIM) ! origin for output

  integer :: k                      ! Dimension counter
  integer :: nelements              ! Number of elements in alldata
  real(kind=PR) :: alldata(1:100)   ! Array containing all data for output

  open(1, file="track1.dat", position="append", form="formatted")

! Record all data for particle p in array 'alldata'
  call record_particle_data(p,nelements,alldata,rorigin)

! Now write all values to the file
  write(1,'(100E15.7)') (alldata(k),k=1,nelements)

  close(1)

  return
END SUBROUTINE track_particle
