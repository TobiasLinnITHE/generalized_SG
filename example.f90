! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

program example

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, IEEE_NEGATIVE_INF, IEEE_POSITIVE_INF
  use, intrinsic :: iso_c_binding,   only: c_double, c_int
  use, intrinsic :: iso_fortran_env, only: real64, real128

  use generalized_SG_m,     only: get_current_tab, get_current
  use distribution_table_m, only: distribution_table
  use util_m,               only: expm1

  implicit none

  logical                  :: status
  real(real64)             :: F(2), dpot, j, j_SG, djdF(2), djddpot
  type(distribution_table) :: tab_FD, tab_GF

  F    = [1.0, 42.0]
  dpot = 7.0

  j_SG = scharfetter_gummel(F, dpot)
  call get_current(maxwell_boltzmann, inv_maxwell_boltzmann, F, dpot, j, djdF, djddpot)
  print "(A)", "Maxwell-Boltzmann"
  print "(A,2ES25.16E3)", "F       = ", F
  print "(A,ES25.16E3)",  "dpot    = ", dpot
  print "(A,ES25.16E3)",  "j_SG    = ", j_SG
  print "(A,2ES25.16E3)", "j, err  = ", j, abs(j - j_SG) / abs(j_SG)
  print "(A,2ES25.16E3)", "djdF    = ", djdF
  print "(A,ES25.16E3)",  "djddpot = ", djddpot
  print *

  call tab_FD%load("tab_FD.bin", status)
  if (.not. status) then
    call tab_FD%init(sqrt_dos, 0.0, ieee_value(1.0, IEEE_POSITIVE_INF), fermi_dirac_density, .false., -50.0, 50.0, 3)
    call tab_FD%save("tab_FD.bin")
  end if
  call get_current_tab(tab_FD, F, dpot, j, djdF, djddpot)
  print "(A)", "Fermi-Dirac"
  print "(A,2ES25.16E3)", "F       = ", F
  print "(A,ES25.16E3)",  "dpot    = ", dpot
  print "(A,ES25.16E3)",  "j       = ", j
  print "(A,2ES25.16E3)", "djdF    = ", djdF
  print "(A,ES25.16E3)",  "djddpot = ", djddpot
  print *

  F    = [0.2, 0.8]
  dpot = 10.0

  call tab_GF%load("tab_GF.bin", status)
  if (.not. status) then
    call tab_GF%init(gaussian_dos, ieee_value(1.0, IEEE_NEGATIVE_INF), ieee_value(1.0, IEEE_POSITIVE_INF), fermi_dirac_density, .true., -40.0, 40.0, 3)
    call tab_GF%save("tab_GF.bin")
  end if
  call get_current_tab(tab_GF, F, dpot, j, djdF, djddpot)
  print "(A)", "Gauss-Fermi"
  print "(A,2ES25.16E3)", "F       = ", F
  print "(A,ES25.16E3)",  "dpot    = ", dpot
  print "(A,ES25.16E3)",  "j       = ", j
  print "(A,2ES25.16E3)", "djdF    = ", djdF
  print "(A,ES25.16E3)",  "djddpot = ", djddpot

contains

  subroutine maxwell_boltzmann(eta, k, F, dFdeta)
    !! k-th derivative of cumulative distribution function
    real(c_double), value, intent(in)  :: eta
      !! normalized chemical potential
    integer(c_int), value, intent(in)  :: k
      !! derivative (can be 0)
    real(c_double),        intent(out) :: F
      !! output k-th derivative of distribution function
    real(c_double),        intent(out) :: dFdeta
      !! output derivative of F wrt eta (used for Newton iteration, might be slightly different from F with k + 1)

    F      = exp(eta)
    dFdeta = F
  end subroutine

  subroutine inv_maxwell_boltzmann(F, eta, detadF)
    !! inverse of cumulative distribution function
    real(c_double), value, intent(in)  :: F
      !! value of cumulative distribution function
    real(c_double),        intent(out) :: eta
      !! output normalized chemical potential corresponding to F
    real(c_double),        intent(out) :: detadF
      !! output derivative of eta wrt F

    eta    = log(F)
    detadF = 1 / F
  end subroutine

  function scharfetter_gummel(n, dpot) result(j)
    real(real64), intent(in) :: n(2)
    real(real64), intent(in) :: dpot
    real(real64)             :: j

    j = ber(-dpot) * n(1) - ber(dpot) * n(2)
  end function

  elemental function ber(x) result(b)
    !! Bernoulli function ber(x) = x / (exp(x) - 1)
    real(real64), intent(in) :: x
    real(real64)             :: b

    if (x < -37) then
      b = - x
    elseif (abs(x) < 1e-9) then
      b = 1.0 - 0.5 * x
    elseif (abs(x) < 1e-4) then
      b = 1.0 + x * (-0.5 + x / 12.0)
    elseif (x < 712.5) then
      b = 0.5 * x * exp(-0.5 * x) / sinh(0.5 * x)
    else
      b = 0.0
    end if
  end function

  function sqrt_dos(t) result(Z) bind(c)
    !! square-root density of states in quad-precision
    real(real128), value, intent(in) :: t
    real(real128)                    :: Z

    real(real128), parameter :: PI = 4 * atan(1.0_16)

    if (t < 0) then
      Z = 0
    else
      Z = 2 * sqrt(t / PI)
    end if
  end function

  function gaussian_dos(t) result(Z) bind(c)
    !! gaussian density of states in quad-precision
    real(real128), value, intent(in) :: t
    real(real128)                    :: Z

    real(real128), parameter :: PI = 4 * atan(1.0_16), SIGMA = 4.0_16

    Z = 1.0_16 / (sqrt(2 * PI) * SIGMA) * exp(-t**2 / (2 * SIGMA**2))
  end function

  function fermi_dirac_density(u, k) result(f) bind(c)
    !! Fermi-Dirac distribution (density) in quad-precision
    real(real128), value, intent(in) :: u
    integer,       value, intent(in) :: k
    real(real128)                    :: f

    real(real128) :: e, f0

    e  = exp(u)
    f0 = 1 / (1 + e)
    f  = 0

    select case (k)
    case (0)
      f = f0
    case (1)
      if (ieee_is_finite(e)) f = (- e * f0) * f0
    case (2)
      if (ieee_is_finite(e)) f = ((e * f0) * (expm1(u) * f0)) * f0
    case (3)
      if (ieee_is_finite(e)) f = (e * f0) * (- 1 + 6 * (e * f0) * (1 - e * f0)) * f0
    case (4)
      if (ieee_is_finite(e)) f = (- (e * f0) + 14 * (e * f0)**2 - 36 * (e * f0)**3 + 24 * (e * f0)**4) * f0
    end select
  end function

end program
