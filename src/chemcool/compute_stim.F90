subroutine compute_stim(a10, e10, rad_temp, b10)

use chemistry_constants
implicit none
double precision  a10, e10, rad_temp, b10
double precision x

x = e10 / (KBOLTZ * rad_temp)
if (x .lt. 5d0) then
   b10 = a10 / (dexp(x) - 1d0)
else
   b10 = 0d0
endif
end subroutine
