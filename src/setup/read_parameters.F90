! READ_PARAMETERS.F90
! D. A. Hubber & K. Rawiraswattana - 20/08/2010
! Read parameters file from a given input file.  The params file is read 
! in by a text parser which looks for lines of the form : 
! 'Comments  : variable_name = value' or simply 'variable_name = value'
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_parameters(filename)
  use definitions
  use filename_module
  use sink_module
  implicit none

  character(len=*), intent(in) :: filename   ! params file name

  character(len=500) :: line                 ! line
  character(len=100) :: var_name             ! variable name
  character(len=100) :: var_value            ! variable value
  logical :: alldone                         ! all parameters read?
  integer :: colon_pos                       ! position of ':' in string
  integer :: equal_pos                       ! position of '=' in string
  integer :: i                               ! parameter counter

  debug1("Reading file : "//trim(filename)//" [read_parameters.F90]")

  open(1, file=filename, status="old", form="formatted")

  params(1:nparams)%done = .false.
  alldone = .true.


! Scan through file and parse text to read parameter values
! ============================================================================
  do

     ! Read entire line into
     read(1,'(A500)',err=10,end=10) line

     ! Find positions of ':' and '=' in line
     colon_pos = index(line,':')
     equal_pos = colon_pos + index(line(colon_pos+1:),'=')

     ! If the ':' is after the '=', or neither exist, skip line
     if (colon_pos >= equal_pos) cycle

     ! Find string inbetween ':' and '=', i.e. the variable name
     var_name  = trim(adjustl(line(colon_pos+1:equal_pos-1)))
     var_value = trim(adjustl(line(equal_pos+1:)))


     ! Now search through all possible parameters and find match
     ! -----------------------------------------------------------------------
     do i=1,nparams
        if (trim(adjustl(var_name)) == trim(params(i)%var_name)) then

           if (params(i)%done) then
              write(6,*) "Fatal error; parameter already read : ",&
                   trim(adjustl(var_name))
              stop
           else
              params(i)%done = .true.
           end if

           select case(params(i)%var_type)
              case('i'); read(var_value(:),*,err=5) params(i)%var_i
              case('j'); read(var_value(:),*,err=5) params(i)%var_j
              case('r'); read(var_value(:),*,err=5) params(i)%var_r
              case('d'); read(var_value(:),*,err=5) params(i)%var_d
              case('c'); read(var_value(:),*,err=5) params(i)%var_c
              case('u'); read(var_value(:),*,err=5) params(i)%var_u
              case('l'); read(var_value(:),*,err=5) params(i)%var_l
              case default; stop
           end select

           cycle
           
5          write(6,*) "Problem reading variable : ",var_name,var_value
           stop
           
        end if
     end do
     ! -----------------------------------------------------------------------

  end do
! ============================================================================

10 close (1)

! Now check through parameters to make sure all have been read
  if (.not. incomplete_params) then
     do i=1,nparams
        if (.not. params(i)%done) then
           alldone = .false.
           write(6,*) "Fatal error!  Parameter not read : ",&
                &trim(adjustl(params(i)%var_name))
        end if
     end do
     if (.not. alldone) stop
  end if

  return
END SUBROUTINE read_parameters
