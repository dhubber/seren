! REORDER_ARRAY.F90
! D. A. Hubber - 30/7/2008
! Reorders elements in all particle arrays to the order given in array aorder.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE reorder_particle_arrays(pstart,plast,dummylist)
  use interface_module, only : reorder_array_double,reorder_array_int,&
       &reorder_array_int_2D,reorder_array_logical,reorder_array_long_int,&
       &reorder_array_real,reorder_array_real_2D,reorder_inverse_array_int
  use definitions
  use particle_module
  use hydro_module
  use time_module
  use type_module
  use tree_module
  use Eos_module
  use HP_module
  implicit none

  integer, intent(in) :: pstart             ! id of first particle
  integer, intent(in) :: plast              ! id of last particle
  integer, intent(in) :: dummylist(1:ptot)  ! Array to be reordered

  debug2("Reorder all particle arrays to given order [reorder_particle_arrays]")

! Reorder main sph array here
  call reorder_array_sph(pstart,plast,sph,dummylist)


  return
END SUBROUTINE reorder_particle_arrays



! ============================================================================
SUBROUTINE reorder_array_sph(pstart,ptot,iarray,aorder)
  use definitions
  use particle_module, only : sph_particle
  implicit none

  integer, intent(in)    :: pstart                    ! id of first particle
  integer, intent(in)    :: ptot                      ! No. of particles
  integer, intent(in)    :: aorder(1:ptot)            ! New order of array
  type(sph_particle), intent(inout) :: iarray(1:ptot) ! Array to be reordered

  integer :: p                                        ! Particle counter
  integer :: pold                                     ! Old p
  type(sph_particle), allocatable :: sphtemp(:)       ! Aux storage array

  debug2("[reorder_array_sph.F90]")

  allocate(sphtemp(1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
     sphtemp(p) = iarray(p)
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
    pold = aorder(p)
    iarray(p) = sphtemp(pold)
  end do

  deallocate(sphtemp)

  return
END SUBROUTINE reorder_array_sph




! ============================================================================
SUBROUTINE reorder_array_int_2D(nsize,pstart,ptot,iarray,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: nsize                   ! No. of elements per part
  integer, intent(in)    :: pstart                  ! id of first particle
  integer, intent(in)    :: ptot                    ! No. of particles
  integer, intent(inout) :: iarray(1:nsize,1:ptot)  ! Array to be reordered
  integer, intent(in)    :: aorder(1:ptot)          ! New order of array

  integer :: k                                      ! Element counter
  integer :: p                                      ! Particle counter
  integer :: pold                                   ! Old p
  integer, allocatable :: ptemp(:,:)                ! Aux storage array

  debug2("[reorder_array_int_2D.F90]")

  allocate(ptemp(1:nsize,1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
    do k=1,nsize
      ptemp(k,p) = iarray(k,p)
    end do
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
    pold = aorder(p)
    do k=1,nsize
       iarray(k,p) = ptemp(k,pold)
    end do
  end do

  deallocate(ptemp)

  return
END SUBROUTINE reorder_array_int_2D



! ============================================================================
SUBROUTINE reorder_array_int(pstart,ptot,iarray,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: pstart            ! id of first particle
  integer, intent(in)    :: ptot              ! No. of particles
  integer, intent(inout) :: iarray(1:ptot)    ! Array to be reordered
  integer, intent(in)    :: aorder(1:ptot)    ! New array order

  integer :: p                                ! Particle counter
  integer :: pold                             ! old id
  integer, allocatable :: ptemp(:)            ! Aux. storage array

  debug2("[reorder_array_int.F90]")

  allocate(ptemp(1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
     ptemp(p) = iarray(p)
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
     pold = aorder(p)
     iarray(p) = ptemp(pold)
  end do

  deallocate(ptemp)

  return
END SUBROUTINE reorder_array_int



! ============================================================================
SUBROUTINE reorder_array_long_int(pstart,ptot,iarray,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: pstart                    ! id of first particle
  integer, intent(in)    :: ptot                      ! No. of particles
  integer, intent(in)    :: aorder(1:ptot)            ! New array order
  integer(kind=ILP), intent(inout) :: iarray(1:ptot)  ! Array to be reordered

  integer :: p                                        ! Particle counter
  integer :: pold                                     ! old id
  integer(kind=ILP), allocatable :: ptemp(:)          ! Aux. storage array

  debug2("[reorder_array_long_int.F90]")

  allocate(ptemp(1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
     ptemp(p) = iarray(p)
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
     pold = aorder(p)
     iarray(p) = ptemp(pold)
  end do

  deallocate(ptemp)

  return
END SUBROUTINE reorder_array_long_int



! ============================================================================
SUBROUTINE reorder_array_real_2D(nsize,pstart,ptot,rmain,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: nsize                       ! # elements in array
  integer, intent(in)    :: pstart                      ! id of first particle
  integer, intent(in)    :: ptot                        ! No. of particles
  integer, intent(in)    :: aorder(1:ptot)              ! New array order
  real(kind=PR), intent(inout) :: rmain(1:nsize,1:ptot) ! Array to reorder

  integer :: k                                          ! Dimension counter
  integer :: p                                          ! Particle counter
  integer :: pold                                       ! Old id
  real(kind=PR), allocatable :: rtemp(:,:)              ! Aux. storage array

  debug2("[reorder_array_real_2D.F90]")

  allocate(rtemp(1:nsize,1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
    do k=1,nsize
      rtemp(k,p) = rmain(k,p)
    end do
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
    pold = aorder(p)
    do k=1,nsize
       rmain(k,p) = rtemp(k,pold)
    end do
  end do

  deallocate(rtemp)

  return
END SUBROUTINE reorder_array_real_2D



! ============================================================================
SUBROUTINE reorder_array_real(pstart,ptot,rmain,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: pstart               ! id of first particle
  integer, intent(in)    :: ptot                 ! No. of particles
  integer, intent(in)    :: aorder(1:ptot)       ! New array order
  real(kind=PR), intent(inout) :: rmain(1:ptot)  ! Array to reorder

  integer :: p                                   ! Particle counter
  integer :: pold                                ! Old id
  real(kind=PR), allocatable :: rtemp(:)         ! Aux. storage array

  debug2("[reorder_array_real.F90]")

  allocate(rtemp(1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
     rtemp(p) = rmain(p)
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
     pold = aorder(p)
     rmain(p) = rtemp(pold)
  end do

  deallocate(rtemp)

  return
END SUBROUTINE reorder_array_real



! ============================================================================
SUBROUTINE reorder_array_double(pstart,ptot,rmain,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: pstart               ! id of first particle
  integer, intent(in)    :: ptot                 ! No. of particles
  integer, intent(in)    :: aorder(1:ptot)       ! New array order
  real(kind=DP), intent(inout) :: rmain(1:ptot)  ! Array to reorder

  integer :: p                                   ! Particle counter
  integer :: pold                                ! Old id
  real(kind=DP), allocatable :: rtemp(:)         ! Aux. storage array

  debug2("[reorder_array_double.F90]")

  allocate(rtemp(1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
     rtemp(p) = rmain(p)
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
     pold = aorder(p)
     rmain(p) = rtemp(pold)
  end do

  deallocate(rtemp)

  return
END SUBROUTINE reorder_array_double



! ============================================================================
SUBROUTINE reorder_array_logical(pstart,ptot,larray,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: pstart          ! id of first particle
  integer, intent(in)    :: ptot            ! No. of particles
  logical, intent(inout) :: larray(1:ptot)  ! Array to reorder
  integer, intent(in)    :: aorder(1:ptot)  ! New array order

  integer :: p                              ! Particle counter
  integer :: pold                           ! Old id
  logical, allocatable :: ltemp(:)          ! Aux. storage array

  debug2("[reorder_array_logical.F90]")

  allocate(ltemp(1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
     ltemp(p) = larray(p)
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
     pold = aorder(p)
     larray(p) = ltemp(pold)
  end do

  deallocate(ltemp)

  return
END SUBROUTINE reorder_array_logical



! ============================================================================
SUBROUTINE reorder_inverse_array_int(pstart,ptot,iarray,aorder)
  use definitions
  implicit none

  integer, intent(in)    :: pstart            ! id of first particle
  integer, intent(in)    :: ptot              ! No. of particles
  integer, intent(inout) :: iarray(1:ptot)    ! Array to be reordered
  integer, intent(in)    :: aorder(1:ptot)    ! New array order

  integer :: p                                ! Particle counter
  integer :: pold                             ! old id
  integer, allocatable :: ptemp(:)            ! Aux. storage array

  debug2("[reorder_inverse_array_int.F90]")

  allocate(ptemp(1:ptot))

! Copy main array into temporary array
  do p=pstart,ptot
     ptemp(p) = iarray(p)
  end do

! Now reorder main array according to aorder
  do p=pstart,ptot
     pold = aorder(p)
     iarray(p) = ptemp(pold)
  end do

  deallocate(ptemp)

  return
END SUBROUTINE reorder_inverse_array_int
