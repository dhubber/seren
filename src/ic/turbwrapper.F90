! turbwrapper.F90
! Andrew McLeod 03/01/2008
! Almost certainly doesn't do what it is supposed to
! which is generate a turbulent velocity field from a
! customizable power spectrum
! ----------------------------------------------------------------------------

program turbwrapper
  use systemargs
  implicit none

  integer, parameter   :: PR=selected_real_kind(p=15) ! double precision
  real, parameter      :: pi = 3.14159265359         ! Pi
  character (len=40)   :: filename                   ! Output filename
  character (LEN=3)    :: yn                         ! Yes/No input
  integer              :: ierr                       ! error variable
  integer              :: i,j,k
  integer              :: n1,n2,n3,k1,k2,k3          ! Real space and fourier space loop counters
  integer              :: Kmax                       ! Fourier space extent
  integer              :: Nmax                       ! Output grid size
  integer              :: d                          ! Dimension counter
  integer, parameter   :: NDIM=3 ! don't change this ! Number of dimensions
  real(kind=PR)        :: alpha                      ! Power spectrum index (=0 for flat)
  real(kind=PR), allocatable :: xreal(:,:,:,:)       ! Output velocity field
  real(kind=PR)        :: rmswanted             ! r.m.s. velocity

  logical              :: commandline        ! Use the command line arguments
  logical              :: commandyn          ! We have the overwrite flag from the command line
  integer              :: nargs              ! Number of command line arguments

  alpha = -2.0 ! Slope of turbulent spectrum
  Kmax = 128     ! Extent of grid in Fourier space
  Nmax = 128     ! Resolution of velocity grid
! ----------------------------------------------------------------------------------------------

  nargs = num_of_args()
  commandline = .FALSE.
  commandyn = .FALSE.

  if (nargs == 2 .OR. nargs == 3) then
    commandline = .TRUE.
    call get_cmd_arg(1,rmswanted)
    call get_cmd_arg(2,filename)
    filename = trim(adjustl(filename))
    if (nargs == 3) then
      commandyn = .TRUE.
      call get_cmd_arg(3,yn)
      if (.NOT.(yn=="yes".OR.yn=="y".OR.yn=="no".OR.yn=="n")) then
        write (6,*) "Invalid command line flag! (third argument should be 'y' or 'n')"
        stop
      end if
    end if
  else if (nargs /= 0) then
    write (6,*) "Invalid number (",nargs,") of command line arguments!"
    stop
  else
    write (6,*) "Enter desired r.m.s. of quantity:"
    read (5,*) rmswanted

    write (6,*) "Enter output filename:"
    read (5,*) filename
  end if

  ! Open output file
  open(unit=10,file=filename,status="NEW",form="FORMATTED",iostat=ierr)
  if (ierr /= 0) then
    if (.NOT. commandyn) then
      write (6,*) "File already exists. Overwrite? (y/n)"
      read (5,*) yn
    end if
    if (yn(1:1) == "y") then
      open(unit=10,file=filename,status="REPLACE",form="FORMATTED",iostat=ierr)
      if (ierr /= 0) then
        write (6,*) "File error, iostat=", ierr
        stop
      end if
    else
      close (10)
      stop
    end if
  end if

! ----------------------------------------------------------------------------------------------

  allocate(xreal(1:NDIM,1:Nmax,1:Nmax,1:Nmax))

  call turbsub(2,rmswanted,alpha,Kmax,Nmax,xreal)

  ! Write velocity grid to file

  do n1=1,Nmax
    do n2=1,Nmax
      do n3=1,Nmax
        write (10,"(I0,A,I0,A,I0,A,G22.15,1X,G22.15,1X,G22.15)") n1," ",n2," ",n3," ",xreal(1:3,n1,n2,n3)
      end do
    end do
  end do

  close(10)

  stop
end program turbwrapper
