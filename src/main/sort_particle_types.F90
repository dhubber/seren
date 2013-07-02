! SORT_PARTICLE_TYPES.F90
! D. A. Hubber - 16/03/2011
! Sort all particles into the correct 'order' in the main arrays according 
! to their types given by the ptype array.
! ============================================================================
SUBROUTINE sort_particle_types(ptype)
  use particle_module, only : ptot
  use hydro_module
  use type_module
  use type_module
  implicit none

  integer, intent(in) :: ptype(1:ptot)         ! Particle types

  integer :: boundaryslot                      ! Aux boundary ptcl counter
  integer :: cdmslot                           ! ..
  integer :: dustslot                          ! ..
  integer :: gasslot                           ! Aux gas ptcl counter
  integer :: icmslot                           ! Aux icm ptcl counter
  integer :: ionslot                           ! ..
  integer :: p                                 ! Particle counter
  integer, allocatable :: porder(:)            ! New particle order

  write(6,*) "Sorting particles types [sort_particle_types.F90]"

! Calculate average value of smoothing length and gas cloud radius  
! to know where to select boundary particles.
  pboundary = 0
  picm      = 0
  pgas      = 0
  pcdm      = 0
  pdust     = 0
  pion      = 0
  allocate(porder(1:ptot))
  

! Find how many particle there are of each type
! ----------------------------------------------------------------------------
  do p=1,ptot
     if (ptype(p) == boundaryid) then
        pboundary = pboundary + 1
     else if (ptype(p) == icmid) then
        picm = picm + 1
     else if (ptype(p) == gasid) then
        pgas = pgas + 1
     else if (ptype(p) == cdmid) then
        pcdm = pcdm + 1
     else if (ptype(p) == dustid) then
        pdust = pdust + 1
     else if (ptype(p) == ionid) then
        pion = pion + 1
     else
        write(6,*) "Unidentified particle type : ",p,ptype(p)
        stop
     end if
  end do

  write(6,*) "pboundary        :",pboundary
  write(6,*) "picm             :",picm
  write(6,*) "pgas             :",pgas
  write(6,*) "pcdm             :",pcdm
  write(6,*) "pdust            :",pdust
  write(6,*) "pion             :",pion
  write(6,*) "ptot             :",ptot


! Now make new ordered list of particles
! ----------------------------------------------------------------------------
  boundaryslot = 0
  icmslot      = pboundary
  gasslot      = pboundary + picm
  cdmslot      = pboundary + picm + pgas
  dustslot     = pboundary + picm + pgas + pcdm
  ionslot      = pboundary + picm + pgas + pcdm + pdust

  do p=1,ptot
     if (ptype(p) == boundaryid) then
        boundaryslot = boundaryslot + 1
        porder(boundaryslot) = p
     else if (ptype(p) == icmid) then
        icmslot = icmslot + 1
        porder(icmslot) = p
     else if (ptype(p) == gasid) then
        gasslot = gasslot + 1
        porder(gasslot) = p
     else if (ptype(p) == cdmid) then
        cdmslot = cdmslot + 1
        porder(cdmslot) = p
     else if (ptype(p) == dustid) then
        dustslot = dustslot + 1
        porder(dustslot) = p
     else if (ptype(p) == ionid) then
        ionslot = ionslot + 1
        porder(ionslot) = p
     end if
  end do

! Now reorder all arrays
  call reorder_particle_arrays(1,ptot,porder(1:ptot))
  
  deallocate(porder)

  return
END SUBROUTINE sort_particle_types
