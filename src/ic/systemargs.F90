! For command line arguments, following the method of Daniel Price in Splash
! Andrew McLeod

module systemargs
#ifdef NAG
  use f90_unix
#endif
  implicit none

  INTERFACE get_cmd_arg
    module procedure get_cmd_arg_char, get_cmd_arg_int, get_cmd_arg_real, get_cmd_arg_double
  END INTERFACE get_cmd_arg

  contains

#ifdef FORTRAN2003
  function num_of_args()
    integer :: num_of_args

    num_of_args = COMMAND_ARGUMENT_COUNT()

  end function num_of_args

  subroutine get_cmd_arg_char(i,cmd_arg)
    integer, intent(in)             :: i
    character (LEN=*), intent(out)  :: cmd_arg

    call GET_COMMAND_ARGUMENT(i,cmd_arg)

  end subroutine get_cmd_arg_char

  subroutine get_cmd_arg_int(i,int_arg)
    integer, intent(in)             :: i
    integer, intent(out)            :: int_arg
    character (LEN=10)              :: cmd_arg

    call GET_COMMAND_ARGUMENT(i,cmd_arg)

    read (cmd_arg,*) int_arg

  end subroutine get_cmd_arg_int

  subroutine get_cmd_arg_real(i,real_arg)
    integer, intent(in)             :: i
    real, intent(out)               :: real_arg
    character (LEN=10)              :: cmd_arg

    call GET_COMMAND_ARGUMENT(i,cmd_arg)

    read (cmd_arg,*) real_arg

  end subroutine get_cmd_arg_real

  subroutine get_cmd_arg_double(i,double_arg)
    integer, intent(in)             :: i
    double precision, intent(out)   :: double_arg
    character (LEN=10)              :: cmd_arg

    call GET_COMMAND_ARGUMENT(i,cmd_arg)

    read (cmd_arg,*) double_arg

  end subroutine get_cmd_arg_double
#else

  function num_of_args()
    integer :: num_of_args
#ifndef NAG
    integer :: iargc
#endif

    num_of_args = iargc()

  end function num_of_args

  subroutine get_cmd_arg_char(i,cmd_arg)
    integer, intent(in)             :: i
    character (LEN=*), intent(out)  :: cmd_arg

    call getarg(i,cmd_arg)

  end subroutine get_cmd_arg_char

  subroutine get_cmd_arg_int(i,int_arg)
    integer, intent(in)             :: i
    integer, intent(out)            :: int_arg
    character (LEN=10)               :: cmd_arg

    call getarg(i,cmd_arg)

    read (cmd_arg,*) int_arg

  end subroutine get_cmd_arg_int

  subroutine get_cmd_arg_real(i,real_arg)
    integer, intent(in)              :: i
    real, intent(out)                :: real_arg
    character (LEN=10)               :: cmd_arg

    call getarg(i,cmd_arg)

    read (cmd_arg,*) real_arg

  end subroutine get_cmd_arg_real

  subroutine get_cmd_arg_double(i,double_arg)
    integer, intent(in)              :: i
    double precision, intent(out)    :: double_arg
    character (LEN=10)               :: cmd_arg

    call getarg(i,cmd_arg)

    read (cmd_arg,*) double_arg

  end subroutine get_cmd_arg_double
#endif

end module systemargs
