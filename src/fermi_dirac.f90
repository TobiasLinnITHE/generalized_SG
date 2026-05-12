! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

module fermi_dirac_m
  !! Contains several functions and subroutines related to the Fermi-Dirac integral:
  !!   - Fermi-Dirac integral of half integer order (-9/2, -7/2, -5/2, -3/2, -1/2, 1/2)
  !!   - Inverse function of Fermi-Dirac integral with order 1/2
  !!   - Incomplete Fermi-Dirac integral for half integer order
  !!
  !! References:
  !!
  !!   Fukushima, T. (2015). Precise and fast computation of Fermi–Dirac integral of integer and half integer order by
  !!   piecewise minimax rational approximation. Applied Mathematics and Computation, 259, 708-729.
  !!   coefficients extracted from https://www.researchgate.net/publication/267330765_xfdhtxt_Fortran_program_package_to_compute_Fermi-Dirac_integrals_of_half_integer_orders_k-92-72-52212_and_of_integer_orders_k01210
  !!
  !!   Fukushima, T. (2015). Precise and fast computation of inverse Fermi–Dirac integral of order 1/2 by minimax
  !!   rational function approximation. Applied Mathematics and Computation, 259, 698-707.
  !!   coefficients extracted from https://www.researchgate.net/publication/275275569_xfdi1htxt_Fortran_program_to_compute_the_inverse_Fermi-Dirac_integral_of_order_k12
  !!

  use, intrinsic :: iso_fortran_env, only: real64

  use fermi_dirac_parameters_m

  implicit none

  private

  public :: fd_m9h, fd_m7h, fd_m5h, fd_m3h, fd_m1h, fd_1h, fd_kh
  public :: inv_fd_1h

contains

  elemental function fd_m9h(x) result(f)
    real(real64), intent(in) :: x
    real(real64)             :: f

    f = fd_kh(x, -9)
  end function

  elemental function fd_m7h(x) result(f)
    real(real64), intent(in) :: x
    real(real64)             :: f

    f = fd_kh(x, -7)
  end function

  elemental function fd_m5h(x) result(f)
    real(real64), intent(in) :: x
    real(real64)             :: f

    f = fd_kh(x, -5)
  end function

  elemental function fd_m3h(x) result(f)
    real(real64), intent(in) :: x
    real(real64)             :: f

    f = fd_kh(x, -3)
  end function

  elemental function fd_m1h(x) result(f)
    real(real64), intent(in) :: x
    real(real64)             :: f

    f = fd_kh(x, -1)
  end function

  elemental function fd_1h(x) result(f)
    real(real64), intent(in) :: x
    real(real64)             :: f

    f = fd_kh(x, 1)
  end function

  pure function fd_kh(x, k) result(f)
    !! piecewise approximation of normalized Fermi-Dirac integral of half-integer order k / 2 by rational functions
    real(real64), intent(in) :: x
      !! argument
    integer,      intent(in) :: k
      !! must be one of the values {-9, -7, -5, -3, -1, 1}
    real(real64)             :: f
      !! return value of Fermi-Dirac integral

    integer      :: l, i, izero
    real(real64) :: s, t, w

    ! get ascending index with l = -9 -> 1, -7 -> 2, -5 -> 3, ...
    l = (k + 9) / 2 + 1

    ! get interval such that x < X_INTERVAL(i) (or i == 8 when x out of bounds)
    i = find_interval(x, X_INTERVAL)

    ! if any zero is in this interval, get its index
    izero = findloc(ZEROS_INTERVAL(IZEROS(l):IZEROS(l+1)-1), i, dim = 1)
    if (izero > 0) izero = izero - 1 + IZEROS(l)

    ! select relevant numerator and denominator polynomial coefficients based on interval i and function index l
    associate (num => COEFF(FD_NUM(1,i,l):FD_NUM(2,i,l)), &
        &      den => COEFF(FD_DEN(1,i,l):FD_DEN(2,i,l)))
      ! evaluate approximation with special treatment of first and last intervals
      if (i == 1) then ! x < -2
        w = exp(x)
        t = w * E2
        f = w * (GAMMA_FACTOR(l) + w * horner(t, num) / horner1(t, den))
      else if (i > size(X_INTERVAL)) then ! x > 40
        w = 1 / (x * x)
        t = X_INTERVAL(7)**2 * w
        if (k /= -1) then
          t = 1 - t
          f = x**(0.5 * (k + 2)) / (0.5 * (k + 2)) * (1 + w * horner(t, num) / horner1(t, den))
        else
          ! for k = -1 use a simple polynomial instead of a rational function
          f = x**(0.5 * (k + 2)) / (0.5 * (k + 2)) * (1 + w * horner(t, num))
        end if
      else ! -2 < x < 40
        ! normalize x relative to interval
        t = (x - X_INTERVAL(i-1)) / (X_INTERVAL(i) - X_INTERVAL(i-1))

        if (i == 2) then
          ! for interval -2 .. 0 flip variable around for denominator polynomial
          s = 1 - t

          ! for k == -9, -7, -5 use flipped variable in numerator polynomial as well
          if (k <= -5) t = s
        else
          ! for other intervals, use unflipped variable for numerator and denominator
          s = t
        end if

        ! evaluate rational function
        f = horner(t, num) / horner1(s, den)

        ! apply zero factor if it exists
        if (izero > 0) f = f * (x - ZEROS(izero))
      end if
    end associate

    ! scale by inverse Gamma factor
    f = SCALING_FACTOR(l) * f
  end function


  elemental function inv_fd_1h(F) result(x)
    !! piecewise approximation of inverse of normalized Fermi-Dirac integral of order 1/2 by rational functions
    real(real64), intent(in) :: F
      !! value of normalized Fermi-Dirac integral of order 1/2
    real(real64)             :: x
      !! return argument such that fd1h(x) == F

    integer      :: i
    real(real64) :: s, t, u, w, z

    ! denormalize function value
    u = F * GAMMA_FACTOR(6)

    ! find interval
    i = find_interval(u, F1H_INTERVAL)

    ! select relevant numerator and denominator polynomial coefficients based on interval i
    associate (num => COEFF(IFD1H_NUM(1,i):IFD1H_NUM(2,i)), &
      &        den => COEFF(IFD1H_DEN(1,i):IFD1H_DEN(2,i)))
      ! evaluate approximation with special treatment of first and last intervals
      if (i == 1) then
        t = u / F1H_INTERVAL(1)
        x = log(t * horner(t, num) / horner1(t, den))
      else if (i > size(F1H_INTERVAL)) then
        s = (F1H_INTERVAL(5) / u)**(4.0_real64 / 3.0_real64)
        t = 1 - s
        x = sqrt(horner1(t, num) / (s * horner(t, den)))
      else
        t = (u - F1H_INTERVAL(i-1)) / (F1H_INTERVAL(i) - F1H_INTERVAL(i-1))
        x = horner(t, num) / horner1(t, den)
      end if
    end associate
  end function


  pure function find_interval(x, xi) result(i)
    !! helper function to find index of x interval
    real(real64), intent(in) :: x
      !! argument to be searched
    real(real64), intent(in) :: xi(:)
      !! x interval boundaries in ascending order, should have small size
    integer                  :: i
      !! return index such that x < xi(i)

    ! linear search efficient due to small size of xi
    do i = 1, size(xi)
      if (x < xi(i)) exit
    end do
  end function

  pure function horner(x, a) result(p)
    !! helper function to evaluate polynomial given in horner form
    real(real64), intent(in) :: x
      !! argument
    real(real64), intent(in) :: a(:)
      !! polynomial coefficients in descending order
    real(real64)             :: p
      !! return value of polynomial

    integer :: i, n

    n = size(a)

    p = a(1)
    do i = 2, n
      p = p * x + a(i)
    end do
  end function

  pure function horner1(x, a) result(p)
    !! helper function to evaluate polynomial given in horner form
    !! assumes that the highest coefficient is 1
    real(real64), intent(in) :: x
      !! argument
    real(real64), intent(in) :: a(:)
      !! polynomial coefficients in descending order, starting at the second highest (a(0) assumed to be 1)
    real(real64)             :: p
      !! return value of polynomial

    integer :: i, n

    n = size(a)

    p = x + a(1) ! = a(0) * x + a(1)
    do i = 2, n
      p = p * x + a(i)
    end do
  end function

end module
