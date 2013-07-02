! NBODY_ORBITS.F90
! D. A. Hubber - 5/9/2008
! Reads data files produced by sink/star particles and produces plots or 
! animations using pgplot
! ============================================================================

#include "macros.h"

! ============================================================================
PROGRAM nbody_orbits
  use definitions
  use sink_module, only : stot
  use filename_module, only : run_id
  use Nbody_module, only : star
  implicit none

  logical :: file_exists                  ! Flag if file exists or not
  character(len=256) :: file_ext          ! ..
  character(len=256) :: out_file          ! ..
  character(len=256) :: filenames(1:SMAX) ! ..
  character(len=256) :: label             ! ..
  integer :: filestatus(1:SMAX)           ! 1 = sink, 2 = star, -1 = done
  integer :: idum                         ! ..
  integer :: linestyle(1:SMAX)            ! PGPLOT line style for stars
  integer :: nline                        ! Integer time for next line
  integer :: nmark                        ! Integer time for next mark
  integer :: nmarklast                    ! Integer time of last mark
  integer :: s                            ! sink/star counter
  integer :: unit_no                      ! File unit no.
  real(kind=DP) :: dt                     ! ..
  real(kind=DP) :: line_interval          ! Interval between plotting line
  real(kind=DP) :: line_first             ! ..
  real(kind=DP) :: mark_interval          ! Interval between plotting marks
  real(kind=DP) :: rdum(1:21)             ! Aux. array for reading file
  real(kind=SP) :: rtemp(1:NDIM)          ! ..
  real(kind=DP) :: tline                  ! Next line point
  real(kind=DP) :: tmark                  ! Next mark point
  real(kind=DP) :: tmin                   ! Minimum time for plotting
  real(kind=DP) :: tmax                   ! Maximum time for plotting
  real(kind=SP) :: xlabel                 ! ..
  real(kind=PR) :: xmin                   ! ..
  real(kind=PR) :: xmax                   ! ..
  real(kind=SP) :: ylabel                 ! ..
  real(kind=PR) :: ymin                   ! ..
  real(kind=PR) :: ymax                   ! ..
  real(kind=DP) :: rold(1:NDIM,1:SMAX)    ! ..
  real(kind=DP) :: rnew(1:NDIM,1:SMAX)    ! ..
  real(kind=DP) :: told(1:SMAX)           ! ..
  real(kind=DP) :: tnew(1:SMAX)           ! ..
  real(kind=DP) :: vold(1:NDIM,1:SMAX)    ! ..
  real(kind=DP) :: aold(1:NDIM,1:SMAX)    ! ..
  real(kind=DP) :: vnew(1:NDIM,1:SMAX)    ! ..
  real(kind=DP) :: anew(1:NDIM,1:SMAX)    ! ..
  real(kind=DP) :: adot(1:NDIM)           ! ..
  real(kind=DP) :: adot2(1:NDIM)          ! ..
  integer, allocatable :: nlines(:)       ! No. of line points
  real(kind=PR), allocatable :: x(:,:)    ! x-position of line points
  real(kind=PR), allocatable :: y(:,:)    ! y-position of line points

  write(6,*) "[nbody_orbits.F90]"

! Read in input parameters
  open(unit=1,file='nbodyorbits.dat')
  read(1,*) run_id
  read(1,*) tmin
  read(1,*) tmax
  read(1,*) mark_interval
  read(1,*) line_interval
  read(1,*) xmin
  read(1,*) xmax
  read(1,*) ymin
  read(1,*) ymax
  read(1,*) label
  close(1)

  linestyle(1:SMAX) = 0
  filestatus(1:SMAX) = 0
  s = 0
  rold = 0.0_PR; rnew = 0.0_PR
  vold = 0.0_PR; vnew = 0.0_PR
  aold = 0.0_PR; anew = 0.0_PR
  told = 0.0_PR; tnew = 0.0_PR
  allocate(star(1:SMAX))
 

! Search for any star files with matching runid
! ----------------------------------------------------------------------------
  do
     s = s + 1
     if (s>=100) then
        write(file_ext,"(I3)") s
     else if (s>=10) then
        write(file_ext,"(I2)") s
     else
        write(file_ext,"(I1)") s
     end if

     ! First search for the 'run_id.sinkX' files.  If they don't exist, 
     ! next search for the 'run_id.starX' files
     out_file = trim(adjustl(run_id))//".sink"//trim(adjustl(file_ext))
     write(6,*) "Looking for ",out_file
     inquire(file=out_file,exist=file_exists)
     if (file_exists) then
        filestatus(s) = 1
     else
        out_file = trim(adjustl(run_id))//".star"//trim(adjustl(file_ext))
        write(6,*) "Looking for ",out_file
        inquire(file=out_file,exist=file_exists)
        if (file_exists) filestatus(s) = 2
     end if

     ! If either the star or sink file is found, open and read the first 
     ! record for initializing the plot
     if (filestatus(s) > 0) then
        stot = s
        unit_no = 10 + s
        open(unit=unit_no,file=out_file)
        read(unit_no,'(I8,21E15.7)') idum,rdum(1:21) 
        star(s)%tcreate = rdum(1)
        told(s) = rdum(1)
        tnew(s) = rdum(1)
        rold(1:NDIM,s) = rdum(2:4)
        vold(1:NDIM,s) = rdum(5:7)
        aold(1:NDIM,s) = rdum(11:13)
        rnew(1:NDIM,s) = rdum(2:4)
        vnew(1:NDIM,s) = rdum(5:7)
        anew(1:NDIM,s) = rdum(11:13)
     else
        exit
     end if
  end do
! ----------------------------------------------------------------------------

! If no files exist, exit program
  if (stot == 0) then 
     write(6,*) "No files found"
     deallocate(star)
     stop
  end if

  xlabel = real(xmin + 0.1*(xmax - xmin),SP)
  ylabel = real(ymin + 0.9*(ymax - ymin),SP)

! Set-up pgplot routines
  call pgbegin(0,'?',1,1)
!  call pgslw(1)
!  call pgsch(1)
!  call pgsci(1)
  call pgenv(real(xmin,SP),real(xmax,SP),real(ymin,SP),real(ymax,SP),1,-2)
!  call pgslw(1.0)
  call pgaxis('',real(xmin,SP),0.,real(xmax,SP),0., &
       &real(xmin,SP), real(xmax,SP), 1., 1,  0.25, 0.25, 1., 0., 90.)
  call pgaxis('',0.,real(ymin,SP),0.,real(ymax,SP), &
       &real(ymin,SP), real(ymax,SP), 1., 1,  0.25, 0.25, 1., 0., 90.)
!  call pgslw(4)
!  call pgsch(3)
!  call pgsci(1)
!              call pgsch(20)
!              call pgslw(20)
!              call pgsci(1)
!              call pgpt1(0.0,0.0,17)
  call pgmtxt('T',-1.,0.05,0.,label)
!  call pglabel('x','y','')

! Set up real and integer times for next line/mark plots
  nline = int(tmin/line_interval)
  tline = nline*line_interval
  tmark = minval(told(1:stot))
  nmark = int(tline/mark_interval) + 1
  tmark = nmark*mark_interval
  nmarklast = -1

!  write(6,*) "tmin :",tmin,"    tmax :",tmax
!  write(6,*) "tline : ",tline,"   tmark :",tmark

  allocate(x(1:100000,1:stot))
  allocate(y(1:100000,1:stot))
  allocate(nlines(1:stot))
  nlines(1:stot) = 0


! Main loop
! ============================================================================
  do 

     ! Loop over all sinks/stars
     ! -----------------------------------------------------------------------
     do s=1,stot
        unit_no = 10 + s

        ! Continuously loop reading the file until the current line can
        ! be drawn.  
        ! --------------------------------------------------------------------
        do while (told(s) < tline .AND. tnew(s) < tline .AND. filestatus(s)>0)
           rold(1:NDIM,s) = rnew(1:NDIM,s)
           vold(1:NDIM,s) = vnew(1:NDIM,s)
           aold(1:NDIM,s) = anew(1:NDIM,s)
           told(s) = tnew(s)

           if (filestatus(s) > 0) then
              read(unit_no,'(I8,21E15.7)',end=120) idum,rdum(1:21) 
              tnew(s) = rdum(1)
              star(s)%r(1:NDIM) = rdum(2:4)
              star(s)%v(1:NDIM) = rdum(5:7)
              star(s)%m = rdum(8)
              star(s)%h = rdum(9)
              star(s)%radius = rdum(10)
              star(s)%a(1:NDIM) = rdum(11:13)
              star(s)%angmom(1:3) = rdum(14:16)   
              star(s)%gpe = rdum(17)
              star(s)%dmdt = rdum(18)
              rnew(1:NDIM,s) = star(s)%r(1:NDIM)
              vnew(1:NDIM,s) = star(s)%v(1:NDIM)
              anew(1:NDIM,s) = star(s)%a(1:NDIM)
           end if
           cycle

120        if (filestatus(s)==2) then
              close(unit_no)
              filestatus(s) = -1
           end if

           if (filestatus(s)==1) then
              close(unit_no)
              if (s>=100) then
                 write(file_ext,"(I3)") s
              else if (s>=10) then
                 write(file_ext,"(I2)") s
              else
                 write(file_ext,"(I1)") s
              end if
              
              out_file = trim(adjustl(run_id))//".star"//trim(adjustl(file_ext))
              write(6,*) "Looking for ",out_file
              inquire(file=out_file,exist=file_exists)
              if (file_exists) then
                 filestatus(s) = 2
                 open(unit=unit_no,file=out_file)
              end if
           end if

        end do
        ! --------------------------------------------------------------------


        ! Interpolate the end position of the current line and store in 
        ! memory to be drawn later.
        if (tline > told(s) .AND. tline <= tnew(s)) then
           dt = tnew(s) - told(s)
           adot(1:NDIM) = 2.*(-3.*(vold(1:NDIM,s) - vnew(1:NDIM,s) - &
                &2.*(anew(1:NDIM,s) + aold(1:NDIM,s))*dt))/dt**2
           adot2(1:NDIM) = 6.*(2.*(vold(1:NDIM,s) - vnew(1:NDIM,s) + &
                &(anew(1:NDIM,s) +aold(1:NDIM,s))*dt))/dt**3
           dt = tline - told(s)
           rtemp(1:NDIM) = real(rold(1:NDIM,s) + vold(1:NDIM,s)*dt &
                &+ 0.5*aold(1:NDIM,s)*dt*dt,SP) !+ &
               ! &(1./6.)*adot(1:NDIM)*dt**3 + &
               ! &(1./24.)*adot2(1:NDIM)*dt**4,SP)
           nlines(s) = nlines(s) + 1
           if (nlines(s) > 100000) stop 'oops, nlines too big'
           x(nlines(s),s) = rtemp(1)
           y(nlines(s),s) = rtemp(2)

           ! If current time range covers mark, plot it here
           if (tmark > told(s) .AND. tmark <= tnew(s)) then
              dt = tmark - told(s)
              rtemp(1:NDIM) = real(rold(1:NDIM,s) + vold(1:NDIM,s)*dt &
                   &+ 0.5*aold(1:NDIM,s)*dt*dt,SP)  !+ &
                 !  &(1./6.)*adot(1:NDIM)*dt**3 + &
                 !  &(1./24.)*adot2(1:NDIM)*dt**4,SP)
!              call pgsch(20)
!              call pgslw(2)
              call pgsci(1)
              call pgpt1(rtemp(1),rtemp(2),17)
              nmarklast = nmark
           end if
        end if

     end do
     ! -----------------------------------------------------------------------

     ! Move to next mark if done
     if (nmarklast == nmark) then
        nmark = nmark + 1
        tmark = nmark*mark_interval
     else if (tmark < tline) then
        nmark = int(tline/mark_interval) + 1
        tmark = nmark*mark_interval
     end if

     nline = nline + 1
     tline = nline*line_interval

     if (tline > tmax) exit

  end do
! ============================================================================

  linestyle(1) = 4
  linestyle(2) = 2
  linestyle(3) = 1

  ! Plot tracks of all sinks/stars
  do s=1,stot
     if (nlines(s) > 1) then
        call pgslw(1)
        call pgsch(1)
        call pgsls(linestyle(s))
        call pgline(nlines(s),x(1:nlines(s),s),y(1:nlines(s),s))
     end if
  end do


! Finish plotting, close all open files and free up memory
  call pgend
  do s=1,stot
     unit_no = 10 + s
     close(unit_no)
  end do
  deallocate(star)

  stop
END PROGRAM nbody_orbits

