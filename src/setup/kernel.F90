! KERNEL.F90
! D. A. Hubber - 14/02/2011
! All SPH kernel functions.
! i)  M4 cubic spline kernel (ref??), 
! ii) Quintic spline kernel (Morris 1996) 
! ============================================================================

#include "macros.h"

! ============================================================================
! W0
! Main kernel function
! ============================================================================
FUNCTION w0(s)
  use definitions
  use kernel_module
  implicit none

  real(kind=PR) :: w0                    ! kernel value
  real(kind=PR), intent(in) :: s         ! r / h
#if defined(KERNEL_TABLES)
  integer :: ikern                       ! kernel table element
#endif

  ! Use tabulated value
  ! --------------------------------------------------------------------------
#if defined(KERNEL_TABLES)
  ikern = int(HALFKERNTOT*s)
  ikern = min(ikern,KERNTOT)
  ikern = max(ikern,0)
  w0    = w0table(ikern)

  ! M4 kernel
  ! --------------------------------------------------------------------------
#elif defined(M4_KERNEL)
  if (s < 1.0_PR) then
     w0 = KERNNORM*(1.0_PR - 1.5_PR*s*s + 0.75_PR*s**3)
  else if (s < 2.0_PR) then
     w0 = 0.25_PR*KERNNORM*(2.0_PR - s)**3
  else
     w0 = 0.0_PR
  end if

  ! Quintic kernel
  ! --------------------------------------------------------------------------
#elif defined(QUINTIC_KERNEL)
  if (s < 1.0_PR) then
     w0 = KERNNORM*(66.0_PR - 60.0_PR*s*s + 30.0_PR*s**4 - 10.0_PR*s**5)
  else if (s < 2.0_PR) then
     w0 = KERNNORM*(51.0_PR+ 75.0_PR*s - 210.0_PR*s*s + 150.0_PR*s**3 - &
          & 45.0_PR*s**4 + 5.0_PR*s**5)
  else if (s < 3.0_PR) then
     w0 = KERNNORM*(243.0_PR - 405.0_PR*s + 270.0_PR*s*s - &
          & 90.0_PR*s**3 + 15.0_PR*s**4 - s**5)
  else
     w0 = 0.0_PR
  end if

  ! Top-hat kernel
  ! --------------------------------------------------------------------------
#elif defined(LINEAR_KERNEL)
  if (s < KERNRANGE) then
     w0 = KERNNORM*(KERNRANGE - s)
  else
     w0 = 0.0_PR
  end if

#endif
  ! --------------------------------------------------------------------------

END FUNCTION w0



! ============================================================================
! W1
! 1st spatial derivative (dW/dr) of main kernel function.  Contains flat 
! inner-region of Thomas & Couchman (199?) to prevent clumping instability.
! ============================================================================
FUNCTION w1(s)
  use definitions
  use kernel_module
  implicit none

  real(kind=PR) :: w1                    ! kernel value
  real(kind=PR), intent(in) :: s         ! r / h
#if defined(KERNEL_TABLES)
  integer :: ikern                       ! kernel table element
#endif

  ! Use tabulated value
  ! --------------------------------------------------------------------------
#if defined(KERNEL_TABLES)
  ikern = int(HALFKERNTOT*s)
  ikern = min(ikern,KERNTOT)
  ikern = max(ikern,0)
  w1    = w1table(ikern)

  ! M4 kernel
  ! --------------------------------------------------------------------------
#elif defined(M4_KERNEL) && defined(TC_KERNEL)
  if (s < TWOTHIRDS) then
     w1 = -KERNNORM
  else if (s < 1.0_PR) then
     w1 = KERNNORM*(-3.0_PR*s + 2.25_PR*s*s)
  else if (s < 2.0_PR) then
     w1 = -0.75_PR*KERNNORM*(2.0_PR - s)**2
  else
     w1 = 0.0_PR
  end if

  ! M4 kernel
  ! --------------------------------------------------------------------------
#elif defined(M4_KERNEL)
  if (s < 1.0_PR) then
     w1 = KERNNORM*(-3.0_PR*s + 2.25_PR*s*s)
  else if (s < 2.0_PR) then
     w1 = -0.75_PR*KERNNORM*(2.0_PR - s)**2
  else
     w1 = 0.0_PR
  end if

  ! Quintic kernel
  ! --------------------------------------------------------------------------
#elif defined(QUINTIC_KERNEL) && defined(TC_KERNEL)
  if (s < 0.75929848074_PR) then
     w1 = -55.2040_PR*KERNNORM
  else if (s < 1.0_PR) then
     w1 = KERNNORM*(-120.0_PR*s + 120.0_PR*s**3 - 50.0_PR*s**4)
  else if (s < 2.0_PR) then
     w1 = KERNNORM*(75.0_PR - 420.0_PR*s + 450.0_PR*s*s - &
             & 180.0_PR*s**3 + 25.0_PR*s**4)
  else if (s < 3.0_PR) then
     w1 = KERNNORM*(-405.0_PR + 540.0_PR*s - 270.0_PR*s*s + &
             & 60.0_PR*s**3 - 5.0_PR*s**4)
  else
     w1 = 0.0_PR
  end if

  ! Quintic kernel
  ! --------------------------------------------------------------------------
#elif defined(QUINTIC_KERNEL)
  if (s < 1.0_PR) then
     w1 = KERNNORM*(-120.0_PR*s + 120.0_PR*s**3 - 50.0_PR*s**4)
  else if (s < 2.0_PR) then
     w1 = KERNNORM*(75.0_PR - 420.0_PR*s + 450.0_PR*s*s - &
             & 180.0_PR*s**3 + 25.0_PR*s**4)
  else if (s < 3.0_PR) then
     w1 = KERNNORM*(-405.0_PR + 540.0_PR*s - 270.0_PR*s*s + &
             & 60.0_PR*s**3 - 5.0_PR*s**4)
  else
     w1 = 0.0_PR
  end if

  ! Top-hat kernel
  ! --------------------------------------------------------------------------
#elif defined(LINEAR_KERNEL)
  if (s < KERNRANGE) then
     w1 = -KERNNORM
  else
     w1 = 0.0_PR
  end if

#endif
  ! --------------------------------------------------------------------------

END FUNCTION w1



! ============================================================================
! W2
! Unmodified 1st spatial derivative (dW/dr) of main kernel function.
! ============================================================================
FUNCTION w2(s)
  use definitions
  use kernel_module
  implicit none

  real(kind=PR) :: w2                    ! kernel value
  real(kind=PR), intent(in) :: s         ! r / h
#if defined(KERNEL_TABLES)
  integer :: ikern                       ! kernel table element
#endif

  ! Use tabulated value
  ! --------------------------------------------------------------------------
#if defined(KERNEL_TABLES)
  ikern = int(HALFKERNTOT*s)
  ikern = min(ikern,KERNTOT)
  ikern = max(ikern,0)
  w2    = w2table(ikern)

  ! M4 kernel
  ! --------------------------------------------------------------------------
#elif defined(M4_KERNEL)
  if (s < 1.0_PR) then
     w2 = KERNNORM*(-3.0_PR*s + 2.25_PR*s*s)
  else if (s < 2.0_PR) then
     w2 = -0.75_PR*KERNNORM*(KERNRANGE - s)**2
  else
     w2 = 0.0_PR
  end if

  ! Quintic kernel
  ! --------------------------------------------------------------------------
#elif defined(QUINTIC_KERNEL)
  if (s < 1.0_PR) then
     w2 = KERNNORM*(-120.0_PR*s + 120.0_PR*s**3 - 50.0_PR*s**4)
  else if (s < 2.0_PR) then
     w2 = KERNNORM*(75.0_PR - 420.0_PR*s + 450.0_PR*s*s - &
             & 180.0_PR*s**3 + 25.0_PR*s**4)
  else if (s < 3.0_PR) then
     w2 = KERNNORM*(-405.0_PR + 540.0_PR*s - 270.0_PR*s*s + &
             & 60.0_PR*s**3 - 5.0_PR*s**4)
  else
     w2 = 0.0_PR
  end if

  ! Top-hat kernel
  ! --------------------------------------------------------------------------
#elif defined(LINEAR_KERNEL)
  if (s < KERNRANGE) then
     w2 = -KERNNORM
  else
     w2 = 0.0_PR
  end if

#endif
  ! --------------------------------------------------------------------------

END FUNCTION w2



! ============================================================================
! WOMEGA
! Derivative of kernel w.r.t. h for computing grad-h omega correction factor.
! ============================================================================
FUNCTION womega(s)
  use definitions
  use kernel_module
  implicit none

  real(kind=PR) :: womega                ! kernel function
  real(kind=PR), intent(in) :: s         ! r / h
#if defined(KERNEL_TABLES)
  integer :: ikern                       ! kernel table element
#endif

  ! Use tabulated value
  ! --------------------------------------------------------------------------
#if defined(KERNEL_TABLES)
  ikern  = int(HALFKERNTOT*s)
  ikern  = min(ikern,KERNTOT)
  ikern  = max(ikern,0)
  womega = womegatable(ikern)

  ! M4 kernel
  ! --------------------------------------------------------------------------
#elif defined(M4_KERNEL)
  if (s < 1.0_PR) then
     womega = KERNNORM*(-NDIMPR + 1.5_PR*(NDIMPR + 2.0_PR)*s*s - &
          &0.75_PR*(NDIMPR + 3.0_PR)*s**3)
  else if (s < 2.0_PR) then
     womega = KERNNORM*(-2.0_PR*NDIMPR + 3.0_PR*(NDIMPR + 1.0_PR)*s - &
          & 1.50_PR*(NDIMPR + 2.0_PR)*s*s + 0.25_PR*(NDIMPR + 3.0_PR)*s**3)
  else
     womega = 0.0_PR
  end if

  ! Quintic kernel
  ! --------------------------------------------------------------------------
#elif defined(QUINTIC_KERNEL)
  if (s < 1.0_PR) then
     womega = KERNNORM*(-66.0_PR*NDIMPR + 60.0_PR*(NDIMPR + 2.0_PR)*s*s - &
          & 30.0_PR*(NDIMPR + 4.0_PR)*s**4 + 10.0_PR*(NDIMPR + 5)*s**5)
  else if (s < 2.0_PR) then
     womega = KERNNORM*(-51.0_PR*NDIMPR - 75.0_PR*(NDIMPR + 1.0_PR)*s + &
          & 210.0_PR*(NDIMPR + 2.0_PR)*s**2 - &
          & 150.0_PR*(NDIMPR + 3.0_PR)*s**3 + &
          & 45.0_PR*(NDIMPR + 4.0_PR)*s**4 - 5.0_PR*(NDIMPR + 5.0_PR)*s**5)
  else if (s < 3.0_PR) then
     womega = KERNNORM*(-243.0_PR*NDIMPR + 405.0_PR*(NDIMPR + 1.0_PR)*s - &
          & 270.0_PR*(NDIMPR + 2.0_PR)*s**2 + &
          & 90.0_PR*(NDIMPR + 3.0_PR)*s**3 - &
          & 15.0_PR*(NDIMPR + 4.0_PR)*s**4 + (NDIMPR + 5.0_PR)*s**5)
  else
     womega = 0.0_PR
  end if

  ! Top-hat kernel
  ! --------------------------------------------------------------------------
#elif defined(LINEAR_KERNEL)
  if (s < KERNRANGE) then
     womega = -KERNNORM*(NDIMPR*KERNRANGE - (NDIMPR + 1.0_PR)*s)
  else
     womega = 0.0_PR
  end if

#endif
  ! --------------------------------------------------------------------------

END FUNCTION womega



! ============================================================================
! WZETA
! Derivative of potential kernel function w.r.t. h for computing grad-h 
! zeta correction factor for gravitational forces.
! ============================================================================
FUNCTION wzeta(s)
  use definitions
  use kernel_module
  implicit none

  real(kind=PR) :: wzeta                 ! kernel value
  real(kind=PR), intent(in) :: s         ! r / h
#if defined(KERNEL_TABLES)
  integer :: ikern                       ! kernel table element
#endif

  ! Use tabulated value
  ! --------------------------------------------------------------------------
#if defined(KERNEL_TABLES)
  ikern = int(HALFKERNTOT*s)
  ikern = min(ikern,KERNTOT)
  ikern = max(ikern,0)
  wzeta = wzetatable(ikern)

  ! M4 kernel
  ! --------------------------------------------------------------------------
#elif defined(M4_KERNEL)
  if (s < 1.0_PR) then
     wzeta = 1.4_PR - 2.0_PR*s*s + 1.5_PR*s**4 - 0.6_PR*s**5
  else if (s < 2.0_PR) then
     wzeta = 1.6_PR - 4.0_PR*s*s + 4.0_PR*s**3 - 1.5_PR*s**4 + 0.2_PR*s**5
  else
     wzeta = 0.0_PR
  end if

  ! Quintic kernel
  ! --------------------------------------------------------------------------
#elif defined(QUINTIC_KERNEL)
  if (s < 1.0_PR) then
     wzeta = 33.0_PR*s*s - 15.0_PR*s**4 + 5.0_PR*s**6 &
          & - 1.42857142857_PR*s**7 - 34.14285714_PR
  else if (s < 2.0_PR) then
     wzeta = 25.5_PR*s*s + 25.0_PR*s**3 - 52.5_PR*s**4 + 30.0_PR*s**5 &
          & - 7.5_PR*s**6 + 0.7142857143_PR*s**7 - 33.785714286_PR
  else if (s < 3.0_PR) then
     wzeta = 121.5_PR*s*s - 135.0_PR*s**3 + 67.5_PR*s**4 - 18.0_PR*s**5 &
          & + 2.5_PR*s**6 - 0.142857143_PR*s**7 - 52.07142857_PR
  else
     wzeta = 0.0_PR
  end if

  ! Top-hat kernel
  ! --------------------------------------------------------------------------
#elif defined(LINEAR_KERNEL)
  if (s < KERNRANGE) then
     wzeta = -(6.0_PR*s*s/(KERNRANGE**3) - &
          &4.0_PR*s**3/(KERNRANGE**4) - 2.0_PR/KERNRANGE)
  else
     wzeta = 0.0_PR
  end if

#endif
  ! --------------------------------------------------------------------------

END FUNCTION wzeta



! ============================================================================
! WGRAV
! Gravitational force kernel (Price & Monaghan 2007)
! ============================================================================
FUNCTION wgrav(s)
  use definitions
  use kernel_module
  implicit none

  real(kind=PR) :: wgrav                 ! kernel value
  real(kind=PR), intent(in) :: s         ! r / h
#if defined(KERNEL_TABLES)
  integer :: ikern                       ! kernel table element
#endif

  ! Use tabulated value
  ! --------------------------------------------------------------------------
#if defined(KERNEL_TABLES)
  ikern = int(HALFKERNTOT*s)
  ikern = min(ikern,KERNTOT)
  ikern = max(ikern,0)
  wgrav = wgravtable(ikern)

  ! M4 kernel
  ! --------------------------------------------------------------------------
#elif defined(M4_KERNEL)
  if (s < 1.0_PR) then
     wgrav = 1.33333333333_PR*s - 1.2_PR*s**3 + 0.5_PR*s**4
  else if (s < 2.0_PR) then
     wgrav = 2.666666666667_PR*s - 3.0_PR*s**2 + 1.2_PR*s**3 &
          & - 0.166666666667_PR*s**4 - 0.066666666667_PR/s**2
  else
     wgrav = 1.0_PR/(s*s)
  end if

  ! Quintic kernel
  ! --------------------------------------------------------------------------
#elif defined(QUINTIC_KERNEL)
  if (s < 1.0_PR) then
     wgrav = (12.0_PR/359.0_PR)*(22.0_PR*s - 12.0_PR*s**3 + &
          & (30.0_PR/7.0_PR)*s**5 - (5.0_PR/4.0_PR)*s**6)
  else if (s < 2.0_PR) then
     wgrav = (12.0_PR/359.0_PR)*(17.0_PR*s + (75.0_PR/4.0_PR)*s**2 - &
          & 42.0_PR*s**3 + 25.0_PR*s**4 - (45.0_PR/7.0_PR)*s**5 + &
          & (5.0_PR/8.0_PR)*s**6 + (5.0_PR/56.0_PR)/s**2)
  else if (s < 3.0_PR) then
     wgrav = (12.0_PR/359.0_PR)*(81.0_PR*s - (405.0_PR/4.0_PR)*s**2 + &
          & 54.0_PR*s**3 - 15.0_PR*s**4 + (15.0_PR/7.0_PR)*s**5 - &
          & (1.0_PR/8.0_PR)*s**6 - (507.0_PR/56.0_PR)/s**2)
  else
     wgrav = 1.0_PR/(s*s)
  end if

  ! Top-hat kernel
  ! --------------------------------------------------------------------------
#elif defined(LINEAR_KERNEL)
  if (s < KERNRANGE) then
     wgrav = 4.0_PR*s/(KERNRANGE**3) - 3.0_PR*s*s/(KERNRANGE**4)
  else
     wgrav = 1.0_PR/(s**2)
  end if

#endif
  ! --------------------------------------------------------------------------

END FUNCTION wgrav



! ============================================================================
! WPOT
! Gravitational potential kernel (Price & Monaghan 2007)
! ============================================================================
FUNCTION wpot(s)
  use definitions
  use kernel_module
  implicit none

  real(kind=PR) :: wpot                  ! kernel value
  real(kind=PR), intent(in) :: s         ! r / h
#if defined(KERNEL_TABLES)
  integer :: ikern                       ! kernel table element
#endif

  ! Use tabulated value
  ! --------------------------------------------------------------------------
#if defined(KERNEL_TABLES)
  ikern = int(HALFKERNTOT*s)
  ikern = min(ikern,KERNTOT)
  ikern = max(ikern,0)
  wpot  = wpottable(ikern)

  ! M4 kernel
  ! --------------------------------------------------------------------------
#elif defined(M4_KERNEL)
  if (s < 1.0_PR) then
     wpot = 1.4_PR - 0.66666666666666_PR*s**2 + 0.3_PR*s**4 - 0.1_PR*s**5  
  else if (s < 2.0_PR) then
     wpot = -1.0_PR/(15.0_PR*s) + 1.6_PR - 1.3333333333333_PR*s**2 + &
          & s**3 - 0.3_PR*s**4 + (1.0_PR/30.0_PR)*s**5
  else
     wpot = 1.0_PR/s
  end if

  ! Quintic kernel
  ! --------------------------------------------------------------------------
#elif defined(QUINTIC_KERNEL)
  if (s < 1.0_PR) then
     wpot = (12.0_PR/359.0_PR)*(-11.0_PR*s*s + 3.0_PR*s**4 - &
          & (5.0_PR/7.0_PR)*s**6 + (5.0_PR/28.0_PR)*s**7 + (478.0_PR/14.0_PR))
  else if (s < 2.0_PR) then
     wpot = (12.0_PR/359.0_PR)*(-(17.0_PR/2.0_PR)*s*s - &
          & (25.0_PR/4.0_PR)*s**3 + (21.0_PR/2.0_PR)*s**4 - 5.0_PR*s**5 + &
          & (15.0_PR/14.0_PR)*s**6 - (5.0_PR/56.0_PR)*s**7 + &
          & (473.0_PR/14.0_PR) + (5.0_PR/56.0_PR)/s)
  else if (s < 3.0_PR) then
     wpot = (12.0_PR/359.0_PR)*(-(81.0_PR/2.0_PR)*s*s + &
          & (135.0_PR/4.0_PR)*s**3 - (27.0_PR/2.0_PR)*s**4 + 3.0_PR*s**5 - &
          & (5.0_PR/14.0_PR)*s**6 + (1.0_PR/56.0_PR)*s**7 + &
          & (729.0_PR/14.0_PR) - (507.0_PR/56.0_PR)/s)
  else
     wpot = 1.0_PR/s
  end if

  ! Top-hat kernel
  ! --------------------------------------------------------------------------
#elif defined(LINEAR_KERNEL)
  if (s < KERNRANGE) then
     wpot = 2.0_PR/KERNRANGE - 2.0_PR*s*s/(KERNRANGE**3) + &
          &1.0_PR*s**3/(KERNRANGE**4) 
  else
     wpot = 1.0_PR/s
  end if

#endif
  ! --------------------------------------------------------------------------

END FUNCTION wpot
