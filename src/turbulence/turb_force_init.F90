! TURB_FORCE_INIT.F90
! A.McLeod - 24/11/2012
! Initialise turbulent forcing routine, and generate first two turbulent fields
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE turb_force_init
   use definitions
   use time_module
   use filename_module
   use turbulence_module
   use turb_force_module
#if defined(USE_MPI)
   use mpi_communication_module
#endif
   implicit none

   integer       :: i, j, k                ! Loop variables
   integer       :: k_vec(1:NDIM)          ! Wavevector
   
   integer              :: n_seed          ! Length of random seed
   integer              :: clock           ! Integer clock time
   integer, allocatable :: seed(:)         ! Random seed
   
   real(kind=PR)        :: power_norm      ! Normalization from power spectrum
   real(kind=PR)        :: proj_norm       ! Normalization from projection
   real(kind=PR)        :: OU_norm         ! Normalization for OU process
   
   debug2("Initializing turbulence [turb_force_init.F90]")

   ! Tasks to always be done (including MPI non-root tasks)
   ! ---------------------------------------------------------------------------
   
   ! Allocate grids
   allocate(afield_last(1:NDIM,0:TURB_GS-1,0:TURB_GS-1,0:TURB_GS-1))
   allocate(afield_next(1:NDIM,0:TURB_GS-1,0:TURB_GS-1,0:TURB_GS-1))
   allocate(afield_now(1:NDIM,0:TURB_GS-1,0:TURB_GS-1,0:TURB_GS-1))
   
   ! Set grid spacing
   turb_gs_real = real(TURB_GS, PR)
   turb_space = (turb_max - turb_min) / turb_gs_real

   ! Set turbulence update time from autocorrelation time and number of substeps
   turb_dt = turb_T / real(turb_Ndt,DP)
   
#if defined(USE_MPI)
   if (rank/=0) then
      call mpi_share_turb_fields
      return
   end if
#endif

   ! Tasks that do not need to be done by MPI non-root tasks
   ! ---------------------------------------------------------------------------

   ! Initialize turb_changed
   turb_changed = .TRUE.

   ! Set up random seed (from gfortran docs)
   call random_seed(size=n_seed)
   allocate(seed(n_seed))

   call system_clock(count=clock)

   seed = clock + 37 * (/(i-1,i=1,n_seed)/)
   call random_seed(put=seed)

   deallocate(seed)

   ! Allocate grids
   allocate(turb_last(1:NDIM,0:TURB_GS-1,0:TURB_GS-1,0:TURB_GS-1))
   allocate(turb_next(1:NDIM,0:TURB_GS-1,0:TURB_GS-1,0:TURB_GS-1))
   allocate(power_spec(0:TURB_GS-1,0:TURB_GS-1,0:TURB_GS-1))
   
   ! Set decay fraction per timestep dt
   turb_decay_frac = turb_dt / turb_T ! == 1 / turbNdt

   ! Set solenoidal fraction from compressive fraction
   sol_frac = 1.0_PR - comp_frac

   ! Set turbulent field time
   turb_next_time = time ! will be updated later

   ! Set initial power distribution (should be parameterised in some fashion)

   do k=0,TURB_GS-1
      if (k > TURB_GS / 2) then
         k_vec(3) = k - TURB_GS
      else
         k_vec(3) = k
      end if
      do j=0,TURB_GS-1
         if (j > TURB_GS / 2) then
            k_vec(2) = j - TURB_GS
         else
            k_vec(2) = j
         end if
         do i=0,TURB_GS-1
            if (i > TURB_GS / 2) then
               k_vec(1) = i - TURB_GS
            else
               k_vec(1) = i
            end if
            power_spec(i,j,k) = power_spectrum(k_vec)
         end do
      end do
   end do
   
   ! Calculate turbulent normalization
   ! Power normalization comes from FFT of initial power spectrum
   power_norm = power_rms_norm(power_spec)
   ! Projection normalization was empirically estimated and fitted
   proj_norm = proj_rms_norm(sol_frac)
   ! OU norm comes from standard deviation of OU process
   OU_norm = sqrt(turb_T / 2.0_PR)
   ! Combination of all factored (reciprocal for easy multiplication)
   turb_norm = 1.0_PR / (power_norm * proj_norm * OU_norm)

   if (restart) then
      ! Restart - load turbulent fields from files and perform FFTs
      call read_turb_fields
      call FFT_3D(turb_last, afield_last)
      call FFT_3D(turb_next, afield_next)
#if defined(USE_MPI)
      call mpi_share_turb_fields
#endif
   else
      ! Not a restart - set up initial field and perform FFT
      turb_next = cmplx(0.0_PR, 0.0_PR)
      call add_turbulence(turb_next, turb_dt)
   
      ! Fourier transform
      call FFT_3D(turb_next, afield_next)
      afield_next = afield_next * turb_norm * turb_rms

      ! Call turb_next_field to create second field
      call turb_next_field
   end if

   return
  
END SUBROUTINE turb_force_init