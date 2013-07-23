! SEREN_SPLASH_SCREEN.F90
! D. A. Hubber - 5/03/2012 
! Write splash screen to standard output
! ============================================================================
SUBROUTINE seren_splash_screen
  use definitions
  implicit none

  write(6,*) "**********************************************************************"
  write(6,*) "*            ****     ******    *****     ******   *     *           *"
  write(6,*) "*           *    *    *         *    *    *        **    *           *"
  write(6,*) "*           *         *         *    *    *        * *   *           *"
  write(6,*) "*            ****     *****     *****     ******   *  *  *           *"
  write(6,*) "*                *    *         *    *    *        *   * *           *"
  write(6,*) "*           *    *    *         *    *    *        *    **           *"
  write(6,*) "*            ****     ******    *    *    ******   *     *           *"
  write(6,*) "*                                                                    *"
  write(6,*) "*                            Version 1.6.0                           *"
  write(6,*) "*                              17/07/2013                            *"
  write(6,*) "*                                                                    *"
  write(6,*) "*   Coders           : David Hubber, Chris Batty & Andrew McLeod     *"
  write(6,*) "*   Contributions by : Thomas Bisbas, Simon Goodwin,                 *"
  write(6,*) "*                      Krisada Rawiraswattana, Dimitri Stamatellos,  *"
  write(6,*) "*                      Stefanie Walch, Anthony Whitworth             *"
  write(6,*) "*                                                                    *"
  write(6,*) "*                   https://github.com/dhubber/seren                 *"
  write(6,*) "**********************************************************************"
#endif

  return
END SUBROUTINE seren_splash_screen
