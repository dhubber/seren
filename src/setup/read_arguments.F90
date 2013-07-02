! READ_ARGUMENTS.F90
! D. A. Hubber - 21/02/2010
! Process command line arguments for the parameters file.  
! If no argument is given, then use the default 'params.dat' file.
! ============================================================================

#include "macros.h"

! ============================================================================
SUBROUTINE read_arguments
  use interface_module, only : write_debug_column_info,write_makefile_options
  use filename_module
  implicit none

  integer :: narguments            ! No. of command line arguments
  character(len=256) :: auxstring  ! Aux. string variable

  narguments = command_argument_count()
MPI_ROOT
  write(6,*) "No. of command line arguments : ",narguments
MPI_END

! Terminate program if there are too many command line arguments
! ----------------------------------------------------------------------------
  if (narguments > 1) then
MPI_ROOT
     stop 'Too many command line arguments!!'
MPI_END

! Use default "params.dat" file if there is no argument
! ----------------------------------------------------------------------------
  else if (narguments == 0) then
     param_file = "params.dat"

! Process individual argument
! ----------------------------------------------------------------------------
  else if (narguments == 1) then
     call get_command_argument(1,auxstring)
     select case(auxstring)

     ! -----------------------------------------------------------------------
     case ("-H","-h","-help")
MPI_ROOT
        write(6,*)
        write(6,*) "Available command line options in Seren"
        write(6,*) "---------------------------------------"
        write(6,*) "-d, -D, --debug          : Data columns in debug files"
        write(6,*) "--diag                   : Data columns in diagnostic file"
        write(6,*) "-i, -I, --incomplete     : Accept incomplete params file"
        write(6,*) "-h, -H, --help           : help"
        write(6,*) "-m, -M, --makefile       : &
             &Makefile options used to compile Seren"
        write(6,*) "-s, -S, --sinks, --stars : Data columns in sink files"
        write(6,*) "-v, -V, --version        : version number"
        write(6,*) "'paramsfile'             : &
             &Reads named parameter file instead of default"
#ifdef USE_MPI
        ! Gracefully end MPI tasks (do not print walltimes)
        call mpi_finish(.true.)
#endif
        stop
MPI_END

     ! -----------------------------------------------------------------------
     case ("-I","-i","--incomplete")
        incomplete_params = .true.

     ! -----------------------------------------------------------------------
     case ("-V","-v","--version")
MPI_ROOT
        write(6,*) "Seren version " // SEREN_VERSION
MPI_END
#ifdef USE_MPI
        ! Gracefully end MPI tasks (do not print walltimes)
        call mpi_finish(.true.)
#endif
        stop

     ! -----------------------------------------------------------------------
     case ("-M","-m","--makefile")
MPI_ROOT
        call write_makefile_options(6)
MPI_END
#ifdef USE_MPI
        ! Gracefully end MPI tasks (do not print walltimes)
        call mpi_finish(.true.)
#endif
        stop

     ! -----------------------------------------------------------------------
     case ("-D","-d","--debug")
MPI_ROOT
        call write_debug_column_info(6)
MPI_END
#ifdef USE_MPI
        ! Gracefully end MPI tasks (do not print walltimes)
        call mpi_finish(.true.)
#endif
        stop

     ! -----------------------------------------------------------------------
     case ("-S","-s","--sinks","--stars")
MPI_ROOT
        call write_sink_column_info(6)
MPI_END
#ifdef USE_MPI
        ! Gracefully end MPI tasks (do not print walltimes)
        call mpi_finish(.true.)
#endif
        stop

     ! -----------------------------------------------------------------------
     case ("--diag")
MPI_ROOT
        call write_diagnostic_column_info(6)
MPI_END
#ifdef USE_MPI
        ! Gracefully end MPI tasks (do not print walltimes)
        call mpi_finish(.true.)
#endif
        stop

     ! -----------------------------------------------------------------------
     case default
        param_file = trim(adjustl(auxstring))

     end select
     ! -----------------------------------------------------------------------

  end if
! ----------------------------------------------------------------------------

  write(6,*) "param_file : ",trim(adjustl(param_file))

  return
END SUBROUTINE read_arguments
