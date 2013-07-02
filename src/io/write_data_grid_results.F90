! WRITE_DATA_GRID_RESULTS.F90
! D. A. Hubber - 12/12/2006
! ..
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE write_data_grid_results(filename)
  use interface_module, only : distance2,w0,w2
  use neighbour_module
  use hydro_module
  use kernel_module
  use particle_module
  use periodic_module
  use tree_module
  implicit none

  character(len=*), intent(in) :: filename  ! ..

#if defined(DEBUG_GRID_RESULTS)
  integer :: igrid                       ! ..
  integer :: jgrid                       ! ..
  integer, parameter :: ngridx = 1280    ! ..
  integer, parameter :: ngridy = 128     ! ..
  integer :: i                           ! auxilary neighbour counter 
  integer :: pp                          ! neighbouring particles (p')
  integer :: pp_pot                      ! number of neighbours 
  integer, allocatable :: pp_templist(:) ! temp. list of neighbours
  real(kind=PR) :: div_v_p               ! Local copy of velocity divergence
  real(kind=PR) :: dr(1:NDIM)            ! relative position vector
  real(kind=PR) :: dv(1:VDIM)            ! Relative velocity vector
  real(kind=PR) :: dvdr                  ! Scalar product of dv and dr
  real(kind=PR) :: drmag                 ! magnitude of separation
  real(kind=PR) :: drsqd                 ! separation squared
  real(kind=PR) :: hfactor               ! invhp ^ NDIM
  real(kind=PR) :: hpp                   ! Smoothing length of particle p
  real(kind=PR) :: invdrmag              ! ( 1 / drmag )
  real(kind=PR) :: invhpp                ! ( 1 / hp )
  real(kind=PR) :: mpp                   ! mass of neighbour pp
  real(kind=PR) :: rhotemp               ! local value of density
  real(kind=PR) :: rp(1:NDIM)            ! position of particle p
  real(kind=PR) :: vp(1:VDIM)            ! Local copy of velocity of p
  real(kind=PR) :: normsum               ! Unity summation for normalisation
#if defined(SINKS) && defined(SELF_GRAVITY)
  real(kind=PR) :: gpotmin               ! Min. potential of neighbours 
#endif
#if defined(VISC_BALSARA)
  real(kind=PR) :: curl_v_p(1:3)         ! Curl of vel. vector (always 3D)
  real(kind=PR) :: dvXdr(1:3)            ! Cross product of dv and dr
#endif

  !debug3("Calculating SPH quantities [all_sph.F90] for particle ",p)

  allocate(pp_templist(1:ptot))
  open(1, file=filename, status="unknown", form="formatted")
  write(1,'(3I8)') NDIM,ngridx,ngridy
  write(1,'(A)') "export"
  write(1,'(A)') "X Y RHO"
#if defined(DEBUG1)
  write(6,*) "Grid output file :",trim(filename),"   (grid file)"
#endif


! Loop over x-dimension
! ============================================================================
  do jgrid=1,ngridy

     ! Loop over y-dimension
     ! =======================================================================
     do igrid=1,ngridx

        rp(1) = periodic_min(1) + &
             &(real(igrid,PR) - 0.5_PR)*periodic_size(1)/real(ngridx,PR)
        rp(2) = periodic_min(2) + &
             &(real(jgrid,PR) - 0.5_PR)*periodic_size(2)/real(ngridy,PR)
#if NDIM==3
        rp(3) = 0.0_PR
#endif
        
        rhotemp = 0.0_PR
        div_v_p = 0.0_PR
        normsum = 0.0_PR
#if defined(SINKS)
        gpotmin = 0.0_PR
#endif
#if defined(VISC_BALSARA)
        curl_v_p(1:3) = 0.0_PR
        dvXdr(1:3) = 0.0_PR
#endif

        ! Determine neighbour list
#if defined(BH_TREE)
        pp_pot = 0
        call BHhydro_walk(rp(1:NDIM),0.0_PR,pp_pot,ptot,&
             &pp_templist(1:ptot),ctot_hydro,BHhydro(0:ctot_hydro))
#endif

        ! Now loop over all neighbours to find contributions
        ! --------------------------------------------------------------------
#if defined(BH_TREE)
        do i=1,pp_pot
           pp = pp_templist(i)
#else
        do pp=1,ptot
#endif
           mpp = sph(pp)%m
           hpp = sph(pp)%h
           call distance2(rp(1:NDIM),pp,dr(1:NDIM),drsqd)
           dv(1:NDIM) = sph(pp)%v(1:VDIM)
           dvdr  = dot_product(dv(1:NDIM),dr(1:NDIM))
           drmag = sqrt(drsqd) + SMALL_NUMBER
           if (drmag >= KERNRANGE*hpp) cycle
           invdrmag = 1.0_PR / drmag
           invhpp = 1.0_PR / hpp
           hfactor = invhpp**(NDIMPR)
           rhotemp = rhotemp + mpp*hfactor*w0(drmag*invhpp)
           div_v_p = div_v_p - mpp*dvdr*w2(drmag*invhpp)*hfactor*invhpp*invdrmag
           normsum = normsum + mpp*hfactor*w0(drmag*invhpp)/sph(pp)%rho
#if defined(VISC_BALSARA)
#if NDIM==3
           dvXdr(1) = dv(2)*dr(3) - dv(3)*dr(2)
           dvXdr(2) = dv(3)*dr(1) - dv(1)*dr(3)
#endif
#if NDIM==2 || NDIM==3
           dvXdr(3) = dv(1)*dr(2) - dv(2)*dr(1)
           curl_v_p(1:3) = curl_v_p(1:3) - mpp*dvXdr(1:3)*invdrmag* &
2                & w2(drmag*invhpp)*hfactor*invhpp
#endif
#endif
        end do
        ! --------------------------------------------------------------------


        ! Normalise and store SPH quantities in main arrays
        div_v_p   = div_v_p / rhotemp

        ! Store magnitude of curl in balsara array for now, and calculate 
        ! balsara factor after thermal once we know sound speed.
#if defined(VISC_BALSARA)
        curl_v_p(1:3) = curl_v_p(1:3) / rhotemp
#endif  
        if (normsum > 0.0_PR) then
           write(1,'(3E18.9)') rp(1),rp(2),rhotemp/normsum
        else
           write(1,'(3E18.9)') rp(1),rp(2),0.0_PR
        end if

     end do
     ! =======================================================================

  end do
! ============================================================================

  close(1)
  deallocate(pp_templist)

#endif


  return
END SUBROUTINE write_data_grid_results
