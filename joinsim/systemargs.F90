! For command line arguments, following the method of Daniel Price in Splash
! Andrew McLeod

module systemargs
#ifdef NAG
  use f90_unix
#endif
  implicit none

  INTERFACE get_cmd_arg
    module procedure get_cmd_arg_char, get_cmd_arg_int
  END INTERFACE get_cmd_arg

  contains

#ifdef F2003ARGS
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
    character (LEN=100)             :: cmd_arg

    call GET_COMMAND_ARGUMENT(i,cmd_arg)

    read (cmd_arg,*) int_arg

  end subroutine get_cmd_arg_int
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
    character (LEN=100)             :: cmd_arg

    call getarg(i,cmd_arg)

    read (cmd_arg,*) int_arg

  end subroutine get_cmd_arg_int
#endif

end module systemargs

subroutine new_pause()
  ! Because pause is deleted from the standard
  character (LEN=100) :: dummytext

  write (6,*) "Press enter to continue..."
  read (5,*) dummytext

  return
end subroutine new_pause
