!>Module randdp: double precision random number generation
!!
!!   Nick Maclaren's double precision random number generator.
!!   This version, which is compatible with Lahey's ELF90 compiler,
!!   is by Alan Miller ( alan @ vic.cmis.csiro.au, www.vic.cmis.csiro.au/~alan )
!!
!!   Latest revision - 18 December 1997
!!
!!   Copyright (C) 1992  N.M. Maclaren
!!   Copyright (C) 1992  The University of Cambridge
!!
!!   This software may be reproduced and used freely, provided that all
!!   users of it agree that the copyright holders are not liable for any
!!   damage or injury caused by use of this software and that this
!!   condition is passed onto all subsequent recipients of the software,
!!   whether modified or not.

module mod_randdp
  !Replaces COMMON /randpp/ and defines dp
  implicit none

  integer, parameter, public :: dp=SELECTED_REAL_KIND(15,60)
  REAL(dp), save, public :: poly(101), other, offset
  integer, save, public :: index

contains

  !*******************************************************************************
  !> Subroutine sdprnd - seeds random number generator in this processor using integer input
  !Input: iseed
  !Output: poly(1:101), other, offset, index
  !*******************************************************************************
  subroutine sdprnd (iseed)

    integer, intent(in) :: iseed
    !Local variables
    real(dp) :: x
    real(dp), parameter :: xmod = 1000009711.0_dp, ymod = 33554432.0_dp
    integer :: ix, iy, iz, i
    logical, save :: inital = .true.

    !ISEED should be set to an integer between 0 and 9999 inclusive;
    !a value of 0 will initialise the generator only if it has not already been done.
    if(inital .or. iseed /=0) then
      inital = .false.
    else
      return  !<The following statements are no longer executed
    end if

    !INDEX must be initialised to an integer between 1 and 101 inclusive,
    !POLY(1...N) to integers between 0 and 1000009710 inclusive (not all 0),
    !and OTHER to a non-negative proper fraction with denominator 33554432.
    !It uses the Wichmann-Hill generator to do this.

    ix = mod(abs(iseed), 10000) + 1
    iy = 2*ix + 1
    iz = 3*ix + 1
    do i = -10,101
      if(i >= 1) then
        poly(i) = aint(xmod*x)
      end if
      ix = mod(171*ix, 30269)
      iy = mod(172*iy, 30307)
      iz = mod(170*iz, 30323)
      x = mod(dble(ix)/30269.0_dp + dble(iy)/30307.0_dp + dble(iz)/30323.0_dp, 1.0_dp)
    end do
    other = aint(ymod*x)/ymod
    offset = 1.0_dp/ymod
    index = 1

    return
  end subroutine sdprnd

  !*******************************************************************************
  !>Function dprand - generates a double precision random number between 0 and 1
  !*******************************************************************************
  function dprand() result(fn_val)

    real(dp) :: fn_val  !Output
    !Local variables
    real(dp) :: x, y
    !N.B. ymod has been removed from the previous DATA statement; it caused a fatal error as it is not used.
    real(dp), parameter :: xmod = 1000009711.0_dp, xmod2 = 2000019422.0_dp, xmod4 = 4000038844.0_dp
    real(dp), parameter :: tiny = 1.0E-17_dp,  zero = 0.0_dp, one = 1.0_dp
    integer :: n
    logical, save :: inital = .true.

    !This returns a uniform (0,1) random number, with extremely good uniformity properties.
    !It assumes that real(dp) provides at least 33 bits of accuracy, and uses a power of two base.
    if(inital) then
      call sdprnd(0)
      inital = .false.
    end if

    !See [Knuth] for why this implements the algorithm described in the paper.
    !Note that this code is tuned for machines with fast real(dp), but slow multiply and divide; many, many other options are possible.
    n = index - 64  !n=-63?
    if(n <= 0) then
      n = n + 101 !n=38?
    end if
    x = poly(index) + poly(index)
    x = xmod4 - poly(n) - poly(n) - x - x - poly(index)
    if(x < zero) then
      if(x < -xmod) then
        x = x + xmod2
      end if
      if(x < zero) then
        x = x + xmod
      end if
    else
      if(x >= xmod2) then
        x = x - xmod2
        if(x >= xmod) then
          x = x - xmod
        end if
      end if
      if(x >= xmod) then
        x = x - xmod
      end if
    end if
    poly(index) = x
    index = index + 1
    if(index > 101) then
      index = index - 101
    end if

    !Add in the second generator modulo 1, and force to be non-zero.
    !The restricted ranges largely cancel themselves out.
    do
      y = 37.0_dp*other + offset
      other = y - aint(y)
      if(other /= zero) then
        exit
      end if
    end do

    x = x/xmod + other
    if(x >= one) then
      x = x - one
    end if
    fn_val = x + tiny

    return
  end function dprand

end module mod_randdp

