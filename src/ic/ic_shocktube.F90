! IC_SHOCKTUBE.F90
! D. A. Hubber - 17/3/2007
! Creates initial conditions file for shocktube tests using initial 
! relaxed particle distributions of uniform periodic box. The parameters 
! for creating the shocktube are read in via the 'sodparams.dat' file.
! outfile                 : Output file name
! out_file_form           : Output file format
! file1                   : Input file 1 name
! file2                   : Input file 2 name
! file1_form              : File 1 format
! file2_form              : File 2 format
! p1, p2                  : No. of particles in file 1, 2
! n1, n2                  : No. of replicas for LHS/RHS
! rho1, rho2              : Density of LHS/RHS layers
! Press1, Press2          : Pressure for LHS/RHS
! x1, x2                  : x
! y1, y2                  : y
! z1, z2                  : z
! v1(1), v2(1)            : vx
! v1(2), v2(2)            : vy
! v1(3), v2(3)            : vz
! B1(1), B2(1)            : Bx
! B1(2), B2(2)            : By
! B1(3), B2(3)            : Bz
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM ic_shocktube
  use particle_module
  use hydro_module
  use filename_module
  use time_module
  use periodic_module
  use mhd_module
  use type_module
  implicit none
 
  character(len=256) :: file1_form       ! File 1 format
  character(len=256) :: file2_form       ! File 2 format
  character(len=256) :: file1            ! File with LHS box
  character(len=256) :: file2            ! File with RHS box
  character(len=256) :: out_file         ! Output file name
  integer :: i                           ! Aux particle counter
  integer :: k                           ! Dimension counter
  integer :: p                           ! Particle counter
  integer :: p1                          ! No. of particles in LHS box
  integer :: p2                          ! No. of particles in RHS box
  integer :: n1                          ! Number of LHS replicas
  integer :: n2                          ! Number of RHS replicas
  real(kind=PR) :: B1(1:3)               ! B of LHS box
  real(kind=PR) :: B2(1:3)               ! B of RHS box
  real(kind=PR) :: offset                ! offset for repeated boxes
  real(kind=PR) :: m1                    ! Mass of LHS particles
  real(kind=PR) :: m2                    ! Mass of RHS particles
  real(kind=PR) :: Press1                ! Pressure of LHS box
  real(kind=PR) :: Press2                ! Pressure of RHS box
  real(kind=PR) :: rcentre(1:NDIM)       ! Position of origin
  real(kind=PR) :: rho1                  ! Density of LHS box
  real(kind=PR) :: rho2                  ! Density of RHS box
  real(kind=PR) :: T1                    ! Temp. of LHS box
  real(kind=PR) :: T2                    ! Temp. of RHS box
  real(kind=PR) :: volume                ! Volume of box
  real(kind=PR) :: v1(1:3)               ! x-velocity of LHS box
  real(kind=PR) :: v2(1:3)               ! x-velocity of RHS box
  real(kind=PR) :: x1                    ! x-length of LHS box
  real(kind=PR) :: x2                    ! x-length of RHS box
  real(kind=PR) :: y1                    ! y-length of LHS box
  real(kind=PR) :: y2                    ! y-length of RHS box
  real(kind=PR) :: z1                    ! z-length of LHS box
  real(kind=PR) :: z2                    ! z-length of RHS box
  real(kind=PR), allocatable :: r1(:,:)  ! Position of LHS particles
  real(kind=PR), allocatable :: r2(:,:)  ! Position of RHS particles
  real(kind=PR), allocatable :: vs(:,:)  ! Smoothed velocity

! Reading parameter file
  call default_parameters

! Read in sod parameters file 
  write(6,*) "Opening sodparams file"
  open(unit=1,file="sodparams.dat",status='old')
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*); read(1,*); read(1,*); read(1,*)
  read(1,*) out_file
  read(1,*) out_file_form
  read(1,*) file1
  read(1,*) file2
  read(1,*) file1_form
  read(1,*) file2_form
  read(1,*); read(1,*); read(1,*)
  read(1,*) p1, p2
  read(1,*) n1, n2
  read(1,*) gamma
  read(1,*) rho1, rho2
  read(1,*) Press1, Press2
  read(1,*) x1, x2
  read(1,*) y1, y2
  read(1,*) z1, z2
  read(1,*) v1(1), v2(1)
  read(1,*) v1(2), v2(2)
  read(1,*) v1(3), v2(3)
  read(1,*) B1(1), B2(1)
  read(1,*) B1(2), B2(2)
  read(1,*) B1(3), B2(3)
  close(1)

#if defined(IDEAL_MHD)
  if (out_file_form=='dragon_form' .OR. out_file_form=='dragon_unform') then
     stop 'Dragon format incompatible with MHD'
  end if
#endif
  file1         = trim(adjustl(file1))
  file1_form    = trim(adjustl(file1_form))
  file2         = trim(adjustl(file2))
  file2_form    = trim(adjustl(file2_form))
  out_file      = trim(adjustl(out_file))
  out_file_form = trim(adjustl(out_file_form))

! Initialise some variables using parameters
  call initialize_seren_variables_1
  call set_default_particle_types

! Checking compiler flags and parameter values
  call sanitycheck

! Setting up scaling units for simulation
  call units

! Initialize periodic variables
  periodic_min(1) = -real(n1,PR)*x1
  periodic_min(2:3) = 0.0_PR
  periodic_max(1) = real(n2,PR)*x2
  periodic_max(2) = y1
  periodic_max(3) = z1
  periodic_size(1:3) = periodic_max(1:3) - periodic_min(1:3)
  periodic_half(1:3) = 0.5_PR*periodic_size(1:3)

! Calculate temperature of gas
  T1 = Press1 / rho1
  T2 = Press2 / rho2


! 1-D shock-tube
! ============================================================================
#if NDIM==1
  ptot = n1*p1 + n2*p2
  call allocate_memory(.FALSE.)
  allocate(vs(1:VDIM,1:ptot))

! Prepare LHS first
  write(6,*) "Preparing LHS"
  m1 = rho1*x1 / real(p1,PR)
  do p=1,n1*p1
     sph(p)%r(1)      = x1*(real(p,PR) - 0.5_PR)/real(p1,PR) + periodic_min(1)
     sph(p)%m         = m1
     sph(p)%v(1:VDIM) = v1(1:VDIM)
     sph(p)%temp      = T1
     sph(p)%porig     = p
#if defined(INTERNAL_ENERGY)
     sph(p)%u = Pconst*sph(p)%temp/(gamma - 1.0_PR)
#endif
#if defined(IDEAL_MHD)
     B(1:BDIM,p)      = B1(1:BDIM)
#endif
  end do

! Now prepare RHS
  write(6,*) "Preparing RHS"
  m2 = rho2*x2 / real(p2,PR)
  do p=1,n2*p2
     i = n1*p1 + p
     sph(i)%r(1)      = x2*(real(p,PR) - 0.5_PR) / real(p2,PR) + &
          & 0.5_PR*(periodic_min(1) + periodic_max(1))
     sph(i)%m         = m2
     sph(i)%v(1:VDIM) = v2(1:VDIM)
     sph(i)%temp      = T2
     sph(i)%porig     = i
#if defined(INTERNAL_ENERGY)
     sph(i)%u = Pconst*sph(p)%temp/(gamma - 1.0_PR)
#endif
#ifdef IDEAL_MHD
     B(1:BDIM,i)      = B2(1:BDIM)
#endif
  end do


! 2-D or 3-D shock-tubes
! ============================================================================
#else


! Read in file for LHS
! ---------------------------------------------------------------------------- 
  call read_data(file1,file1_form)
  allocate(r1(1:NDIM,1:n1*ptot))

  ! Record info into temp file 
  do i=0,n1-1
     do p=1,ptot
        do k=1,NDIM
           if (k==1) then 
              offset = x1*real(i,PR) + periodic_min(1)
           else 
              offset = 0
           end if
           r1(k,i*ptot + p) = sph(p)%r(k) + offset
        end do
     end do
  end do

  p1 = ptot
#if NDIM==1
  volume = x1
#elif NDIM==2
  volume = x1*y1
#elif NDIM==3
  volume = x1*y1*z1
#endif
  m1 = rho1*volume/real(p1,PR)
  call clean_up

! Now read in data for RHS
! ----------------------------------------------------------------------------
  call read_data(file2,file2_form)
  allocate(r2(1:NDIM,1:n2*ptot))

  ! Record info into temp file 
  do i=0,n2-1
     do p=1,ptot
        do k=1,NDIM
           if (k==1) then 
              offset = x2*real(i,PR) + &
                   &0.5_PR*(periodic_min(1) + periodic_max(1))
           else 
              offset = 0
           end if
           r2(k,i*ptot + p) = sph(p)%r(k) + offset
        end do
     end do
  end do

  p2 = ptot
#if NDIM==1
  volume = x2
#elif NDIM==2
  volume = x2*y2
#elif NDIM==3
  volume = x2*y2*z2
#endif
  m2 = rho2*volume/real(p2,PR)
  call clean_up

! Combine into main arrays and add velocity information 
! ----------------------------------------------------------------------------
  ptot = n1*p1 + n2*p2
  call allocate_memory(.FALSE.)
  allocate(vs(1:VDIM,1:ptot))


! First write LHS of shock tube info 
! ----------------------------------------------------------------------------
  do p=1,n1*p1
     sph(p)%r(1:NDIM) = r1(1:NDIM,p)
     sph(p)%m         = m1
     sph(p)%h         = 1.0_PR
     sph(p)%v(1:VDIM) = v1(1:VDIM)
     sph(p)%temp      = T1
     sph(p)%porig     = p
#if defined(IDEAL_MHD)
     B(1:BDIM,p)      = B1(1:BDIM)
#endif
  end do


! Now write RHS of shock tube info 
! ----------------------------------------------------------------------------
  do p=1,n2*p2
     i = p1*n1 + p
     sph(i)%r(1:NDIM) = r2(1:NDIM,p)
     sph(i)%m         = m2
     sph(i)%h         = 1.0_PR
     sph(i)%v(1:VDIM) = v2(1:VDIM)
     sph(i)%temp      = T2
     sph(i)%porig     = i
#if defined(IDEAL_MHD)
     B(1:BDIM,i)      = B2(1:BDIM)
#endif
  end do

#endif
! ============================================================================

! Setting up for different particle types
  pboundary = 0
  picm      = 0
  pgas      = ptot
  pcdm      = 0
  pdust     = 0
  pion      = 0
  call types
  write(6,*) "No. of particles : ",ptot

#if defined(KERNEL_TABLES)
  call tabulate_kernel_functions
#endif
  call initialize_seren_variables_2
  call initialize_sph_variables_1

! Build and stock trees for first time
  call tree_update(nbuild,nstock)

! Make initial guesses of h either using tree or a global estimate
#if defined(BH_TREE) && !defined(CONSTANT_H)
  if (.not. restart) call BHhydro_hguess
#elif !defined(CONSTANT_H)
  if (.not. restart) call h_guess
#endif

! Calculating initial SPH quantities (h, neibs, rho, etc.)
  call sph_update
  call initialize_thermal_properties

! Calculate smoothed velocities of all particles
  do p=1,ptot
     call smoothed_velocity(p,vs(1:VDIM,p))
  end do
  sph(1:ptot)%v(1) = vs(1,1:ptot)

! Calculate hydro forces on all SPH particles
  call sph_hydro_forces

! Write everything to file 
  call write_data(out_file,out_file_form)
#if defined(DEBUG_PLOT_DATA)
  rcentre(1:NDIM) = 0.0_PR
  call write_data_debug("ICSHOCK.debug.dat",rcentre(1:NDIM))
#endif

#if NDIM==2 || NDIM==3
  deallocate(r1)
  deallocate(r2)
#endif

  deallocate(vs)
  call clean_up


  stop
END PROGRAM ic_shocktube
