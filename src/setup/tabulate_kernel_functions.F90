! TABULATE_KERNEL_FUNCTIONS.F90
! C. P. Batty & D. A. Hubber - 8/12/2006
! Tabulates kernel tables for 
! i)   M4 cubic spline kernel (ref??), 
! ii)  Thomas-Couchman modification of M4 kernel (Thomas & Couchman 199?), or
! iii) Quintic spline kernel (Morris 1996) 
! Tables contain KERNTOT + 1 entries going from i=0 (s=0) to i=KERNTOT (s=2) 
! with an equal spacing of 2/KERNTOT.  
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE tabulate_kernel_functions
  use kernel_module
  implicit none

  integer :: i                  ! counter of kernel table elements
  integer :: kern1              ! piecewise boundary in table
  real(kind=PR) :: invkerntot   ! ( 1 / KERNTOT )
  real(kind=PR) :: s            ! ( r / h ) argument to calculate W
  real(kind=PR) :: saux         ! ( 2 - s )
#if defined(TC_KERNEL)
  integer :: kerntc             ! piecewise boundary for Thomas-Couchman
#endif
#if defined(QUINTIC_KERNEL)
  integer :: kern2              ! 2nd piecewise boundary in table
#endif

  debug1("Calculating kernel tables [kernel.F90]")

! Limits of different kernel regions
  invkerntot = 1.0_PR / real(KERNTOT,PR)
#if defined(M4_KERNEL)
  kern1 = int(HALFKERNTOT)
#if defined(TC_KERNEL)
  kerntc = int(real(KERNTOT,PR) / 3.0_PR)
#endif
#elif defined(QUINTIC_KERNEL)
  kern1 = int(ONETHIRD*real(KERNTOT,PR))
  kern2 = int(2.*ONETHIRD*real(KERNTOT,PR))
#if defined(TC_KERNEL)
  kerntc = int(0.75929848074_PR*real(KERNTOT,PR) / 3.0_PR)
#endif
#endif

! Loop over and set values to all kernel table elements
! ============================================================================
  do i=0,KERNTOT
     s = (KERNRANGE*real(i,PR) + 1.0_PR)*invkerntot  
     saux = KERNRANGE - s

#if defined(M4_KERNEL)
     ! -----------------------------------------------------------------------
     if (i < kern1) then
#if defined(TC_KERNEL)
        ! Thomas-Couchman (1992) kernel gradient
        if (i < kerntc) then
           w1table(i) = -1.0_PR
        else
           w1table(i) = -3.0_PR*s + 2.25_PR*s*s
        end if
#else
        ! Standard M4 kernel gradient
        w1table(i) = -3.0_PR*s + 2.25_PR*s*s
#endif
        w0table(i) = 1.0_PR - 1.5_PR*s*s + 0.75_PR*s**3
        w2table(i) = -3.0_PR*s + 2.25_PR*s*s
        wgravtable(i) = (4.0_PR/3.0_PR)*s - 1.2_PR*s**3 + 0.5_PR*s**4
        wpottable(i) = 7.0_PR/5.0_PR - (2.0_PR/3.0_PR)*s**2 + &
             & 0.3_PR*s**4 - 0.1_PR*s**5  
#if defined(GRAD_H_SPH)
        ! *** alternatively/equivalently: womegatable(i) = -NDIM*w0table(i) - s*w2table(i)
        !                                 wzetatable(i) = -1/h^2 * (
        womegatable(i) = -NDIMPR + 1.5*(NDIMPR + 2.0_PR)*s*s - &
             & 0.75_PR*(NDIMPR + 3.0_PR)*s**3
        wzetatable(i) = 1.4_PR - 2.0_PR*s*s + 1.5_PR*s**4 - 0.6_PR*s**5
#endif
     ! -----------------------------------------------------------------------
     else
        w0table(i) = 0.25_PR*saux**3
        w1table(i) = -0.75_PR*saux*saux
        w2table(i) = -0.75_PR*saux*saux
        wgravtable(i) = (8.0_PR/3.0_PR)*s - 3.0_PR*s**2 + 1.2_PR*s**3 - &
             & (1.0_PR/6.0_PR)*s**4 - 1.0_PR/(15.0_PR*s**2)
        wpottable(i) = -1.0_PR/(15.0_PR*s) + 8.0_PR/5.0_PR - (4.0_PR/3.0_PR)*s**2 + &
             & s**3 - 0.3_PR*s**4 + (1.0_PR/30.0_PR)*s**5
#if defined(GRAD_H_SPH)
        ! *** alternatively/equivalently: womegatable(i) = -NDIM*w0table(i) - s*w2table(i)
        womegatable(i) = -2.0_PR*NDIMPR + 3.0_PR*(NDIMPR + 1.0_PR)*s - &
             & 1.50_PR*(NDIMPR + 2.0_PR)*s*s + 0.25_PR*(NDIMPR + 3.0_PR)*s**3
        wzetatable(i) = 1.6_PR - 4.0_PR*s*s + 4.0_PR*s**3 - 1.5_PR*s**4 + 0.2_PR*s**5
#endif
     end if
     ! -----------------------------------------------------------------------

     ! Normalise kernel tables depending on dimensionality
#if NDIM == 1
     w0table(i) = w0table(i)*TWOTHIRDS
     w1table(i) = w1table(i)*TWOTHIRDS
     w2table(i) = w2table(i)*TWOTHIRDS
#if defined(GRAD_H_SPH)
     womegatable(i) = womegatable(i)*TWOTHIRDS
#endif
#elif NDIM == 2
     w0table(i) = w0table(i)*INVPI*10.0_PR/7.0_PR
     w1table(i) = w1table(i)*INVPI*10.0_PR/7.0_PR
     w2table(i) = w2table(i)*INVPI*10.0_PR/7.0_PR
#if defined(GRAD_H_SPH)
     womegatable(i) = womegatable(i)*INVPI*10.0_PR/7.0_PR
#endif
#else
     w0table(i) = w0table(i)*INVPI
     w1table(i) = w1table(i)*INVPI
     w2table(i) = w2table(i)*INVPI
#if defined(GRAD_H_SPH)
     womegatable(i) = womegatable(i)*INVPI
#endif
#endif
#endif


#if defined(QUINTIC_KERNEL)
     ! -----------------------------------------------------------------------
     if (i < kern1) then
#if defined(TC_KERNEL)
        ! Thomas-Couchman (1992) kernel gradient
        if (i < kerntc) then
           w1table(i) = -55.2040_PR
        else
           w1table(i) = -120.0_PR*s + 120.0_PR*s**3 - 50.0_PR*s**4
        end if
#else
        ! Standard Quintic kernel gradient
        w1table(i) = -120.0_PR*s + 120.0_PR*s**3 - 50.0_PR*s**4
#endif
        w0table(i) = 66.0_PR - 60.0_PR*s*s + 30.0_PR*s**4 - 10.0_PR*s**5
        w2table(i) = -120.0_PR*s + 120.0_PR*s**3 - 50.0_PR*s**4
        wgravtable(i) = 22.0_PR*s - 12.0_PR*s**3 + (30.0_PR/7.0_PR)*s**5 - &
             & (5.0_PR/4.0_PR)*s**6
        wpottable(i) = -11.0_PR*s*s + 3.0_PR*s**4 - (5.0_PR/7.0_PR)*s**6 + &
             & (5.0_PR/28.0_PR)*s**7 + (478.0_PR/14.0_PR)
#if defined(GRAD_H_SPH)
        womegatable(i) = -66.0_PR*NDIMPR + 60.0_PR*(NDIMPR + 2.0_PR)*s*s - &
             & 30.0_PR*(NDIMPR + 4.0_PR)*s**4 + 10.0_PR*(NDIMPR + 5)*s**5
        wzetatable(i) = 33.0_PR*s*s - 15.0_PR*s**4 + 5.0_PR*s**6 - &
             (10.0_PR/7.0_PR)*s**7 - (478.0_PR/14.0_PR)
#endif
     ! -----------------------------------------------------------------------
     else if (i < kern2) then
        w0table(i) = 51.0_PR+ 75.0_PR*s - 210.0_PR*s*s + 150.0_PR*s**3 - &
             & 45.0_PR*s**4 + 5.0_PR*s**5
        w1table(i) = 75.0_PR - 420.0_PR*s + 450.0_PR*s*s - &
             & 180.0_PR*s**3 + 25.0_PR*s**4
        w2table(i) = 75.0_PR - 420.0_PR*s + 450.0_PR*s*s - &
             & 180.0_PR*s**3 + 25.0_PR*s**4
        wgravtable(i) = 17.0_PR*s + (75.0_PR/4.0_PR)*s**2 - 42.0_PR*s**3 + &
             & 25.0_PR*s**4 - (45.0_PR/7.0_PR)*s**5 &
             &+ (5.0_PR/8.0_PR)*s**6 + (5.0_PR/56.0_PR)/s**2
        wpottable(i) = -(17.0_PR/2.0_PR)*s*s - (25.0_PR/4.0_PR)*s**3 + &
             & (21.0_PR/2.0_PR)*s**4 - 5.0_PR*s**5 + &
             & (15.0_PR/14.0_PR)*s**6 - (5.0_PR/56.0_PR)*s**7 + &
             & (473.0_PR/14.0_PR) + (5.0_PR/56.0_PR)/s
#if defined(GRAD_H_SPH)
        womegatable(i) = -51.0_PR*NDIMPR - 75.0_PR*(NDIMPR + 1.0_PR)*s + &
             & 210.0_PR*(NDIMPR + 2.0_PR)*s**2 - &
             & 150.0_PR*(NDIMPR + 3.0_PR)*s**3 + &
             & 45.0_PR*(NDIMPR + 4.0_PR)*s**4 - &
             & 5.0_PR*(NDIMPR + 5.0_PR)*s**5
        wzetatable(i) = (51.0_PR/2.0_PR)*s*s + 25.0_PR*s**3 - (105.0_PR/2.0_PR)*s**4 &
             & + 30.0_PR*s**5 - (15.0_PR/2.0_PR)*s**6 + (5.0_PR/7.0_PR)*s**7 &
             & - (473.0_PR/14.0_PR)
#endif
    ! ------------------------------------------------------------------------
    else
        w0table(i) = 243.0_PR - 405.0_PR*s + 270.0_PR*s*s - &
             & 90.0_PR*s**3 + 15.0_PR*s**4 - s**5
        w1table(i) = -405.0_PR + 540.0_PR*s - 270.0_PR*s*s + &
             & 60.0_PR*s**3 - 5.0_PR*s**4
        w2table(i) = -405.0_PR + 540.0_PR*s - 270.0_PR*s*s + &
             & 60.0_PR*s**3 - 5.0_PR*s**4
        wgravtable(i) = 81.0_PR*s - (405.0_PR/4.0_PR)*s**2 + 54.0_PR*s**3 - &
             & 15.0_PR*s**4 + (15.0_PR/7.0_PR)*s**5 - &
             & (1.0_PR/8.0_PR)*s**6 - (507.0_PR/56.0_PR)/s**2
        wpottable(i) = -(81.0_PR/2.0_PR)*s*s + (135.0_PR/4.0_PR)*s**3 - &
             & (27.0_PR/2.0_PR)*s**4 + 3.0_PR*s**5 - &
             & (5.0_PR/14.0_PR)*s**6 + (1.0_PR/56.0_PR)*s**7 + &
             & (729.0_PR/14.0_PR) - (507.0_PR/56.0_PR)/s
#if defined(GRAD_H_SPH)
        womegatable(i) = -243.0_PR*NDIMPR + 405.0_PR*(NDIMPR + 1.0_PR)*s - &
             & 270.0_PR*(NDIMPR + 2.0_PR)*s**2 + &
             & 90.0_PR*(NDIMPR + 3.0_PR)*s**3 - &
             & 15.0_PR*(NDIMPR + 4.0_PR)*s**4 + (NDIMPR + 5.0_PR)*s**5
        wzetatable(i) = (243.0_PR/2.0_PR)*s*s - 135.0_PR*s**3 + &
             & (135.0_PR/2.0_PR)*s**4 &
             & - 18.0_PR*s**5 + (5.0_PR/2.0_PR)*s**6 - (1.0_PR/7.0_PR)*s**7 &
             & - (729.0_PR/14.0_PR)
#endif
     end if
     ! -----------------------------------------------------------------------


     ! Normalise kernel tables depending on dimensionality
#if NDIM == 1
     w0table(i) = w0table(i)*120.0_PR
     w1table(i) = w1table(i)*120.0_PR
     w2table(i) = w2table(i)*120.0_PR
#if defined(GRAD_H_SPH)
     womegatable(i) = womegatable(i)*120.0_PR
#endif
#elif NDIM == 2
     w0table(i) = w0table(i)*INVPI*7.0_PR/478.0_PR
     w1table(i) = w1table(i)*INVPI*7.0_PR/478.0_PR
     w2table(i) = w2table(i)*INVPI*7.0_PR/478.0_PR
#if defined(GRAD_H_SPH)
     womegatable(i) = womegatable(i)*INVPI*7.0_PR/478.0_PR
#endif
#else
     w0table(i) = w0table(i)*INVPI*3.0_PR/359.0_PR
     w1table(i) = w1table(i)*INVPI*3.0_PR/359.0_PR
     w2table(i) = w2table(i)*INVPI*3.0_PR/359.0_PR
     wgravtable(i) = wgravtable(i)*12.0_PR/359.0_PR
     wpottable(i) = wpottable(i)*12.0_PR/359.0_PR
#if defined(GRAD_H_SPH)
     womegatable(i) = womegatable(i)*INVPI*3.0_PR/359.0_PR
#endif
#if defined(QUINTIC_KERNEL) && defined(GRAD_H_SPH)
     wzetatable(i) = -wzetatable(i)*12.0_PR/359.0_PR
#endif
#endif
#endif

  end do
! ============================================================================

! Output kernel tables to file "kernel.dat"
#if defined(DEBUG_KERNEL)
  write(6,*) "Writing kernel tables to file"
  open(1, file="kernel.dat", status="unknown", form="formatted")
  do i=0,KERNTOT
     s = (KERNRANGE * real(i,PR) + 1.) * invkerntot
#if defined(GRAD_H_SPH)
     write(1,'(8E15.7)') s, w0table(i), w1table(i), w2table(i), &
          &wgravtable(i), wpottable(i), womegatable(i), wzetatable(i)
#else
     write(1,'(6E15.7)') s, w0table(i), w1table(i), w2table(i), &
          &wgravtable(i), wpottable(i)
#endif
  end do
  close(1)
#endif

  return
END SUBROUTINE tabulate_kernel_functions
