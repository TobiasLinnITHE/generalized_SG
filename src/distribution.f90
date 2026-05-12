! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

module distribution_m

  use, intrinsic :: iso_c_binding,   only: c_double, c_float128, c_int
  use, intrinsic :: iso_fortran_env, only: real64, real128
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_negative_inf, ieee_positive_inf, ieee_is_finite

  use cheb_table_m,     only: cheb_table
  use fermi_dirac_m,    only: fd_1h, fd_m1h, fd_m3h, fd_m5h, fd_m7h, fd_m9h, inv_fd_1h
  use integral_table_m, only: integral_table
  use quad_m,           only: quad
  use util_m,           only: expm1, linspace, PI, PI_128, LOGGING_INFO, int2str

  implicit none

  private
  public :: fermi_dirac_density
  public :: maxwell_boltzmann, inv_maxwell_boltzmann, int_maxwell_boltzmann
  public :: init_fermi_dirac, fermi_dirac_1h, inv_fermi_dirac_1h, int_fermi_dirac_1h
  public :: init_band_tailing, band_tailing, inv_band_tailing, int_band_tailing
  public :: init_gauss_fermi, gauss_fermi, inv_gauss_fermi, int_gauss_fermi

  type(integral_table) :: itab_fd_1h(-5:5)

  type band_tailing_params
    !! private helper class to store tables and parameters for the band-tailing Fermi-Dirac integral evaluation

    real(real64) :: t0
      !! band-tailing cutoff energy
    real(real64) :: uc
      !! critical value with 1.0 + exp(-uc) = 1.0 and exp(uc) + 1.0 = exp(uc) in double precision

    type(cheb_table)     :: ctab(0:5)
      !! Chebyshev lookup table for central region (eta_small < eta < eta_large)
    type(integral_table) :: itab(-5:5)
      !! evaluate integrals int_{eta1}^{eta2} F(eta)^k d{eta}

    real(real64) :: A, B, C, D
      !! parameters to evaluate integral for eta <= t0-uc and eta >= t0+uc

    real(real64) :: eta_small
      !! t0 - uc
    real(real64) :: eta_large
      !! t0 + uc
    real(real64) :: F_small
      !! F(eta_small)
    real(real64) :: F_large
      !! F(eta_large)
  end type
  type(band_tailing_params) :: bt

  type gauss_fermi_params
    !! private helper class to store tables and parameters for the Gauss-Fermi integral evaluation

    real(real64) :: sigma
      !! width of gaussian distribution

    type(cheb_table)     :: ctab(0:5)
      !! Chebyshev lookup table
    type(integral_table) :: itab(-5:5)
      !! evaluate integrals int_{eta1}^{eta2} F(eta)^k d{eta}

    real(real64) :: eta_min
      !! minimum supported eta
    real(real64) :: eta_max
      !! maximum supported eta
    real(real64) :: F_min
      !! F(eta_min)
    real(real64) :: F_max
      !! F(eta_max)
  end type
  type(gauss_fermi_params) :: gf

contains

  function fermi_dirac_density(u, k) result(f) bind(c, name = "ITHE_fermi_dirac_density")
    !! k-th derivative of Fermi-Dirac distribution density in quad precision
    real(c_float128), value, intent(in) :: u
      !! argument u = t - eta
    integer(c_int), value, intent(in) :: k
      !! derivative index
    real(c_float128)             :: f
      !! output f^{(k)}

    real(real128), parameter :: ZERO_3 = acosh(2.0_real128)
      !! zeros of third derivative (take +/-  of ZERO3)
    real(real128), parameter :: ZERO_4 = acosh(5.0_real128)
      !! non-trivial zeros of fourth derivative (take +/- of ZERO4)
    real(real128), parameter :: ZEROS_5(2) = [acosh((13 + sqrt(105.0_real128))/2), acosh((13 - sqrt(105.0_real128))/2)]
      !! zeros of fifth derivative (take +/- of each one)
    real(real128), parameter :: BIG = log(huge(1.0_real128))
      !! largest value for which exponential function does not produce an overflow

    real(real128) :: c, s1, s2, s3, s4

    select case (k)
    case (0)
      if (u <= BIG) then
        f = 1 / (1 + exp(u))
      else
        f = 0
      end if
    case (1)
      if (abs(u) <= BIG) then
        c = cosh(0.5 * u)
        f = - 0.25 / c**2
      else
        f = 0
      end if
    case (2)
      if (abs(u) < 1e-16) then
        f = u*(0.125 - u**2 / 24)
      else
        if (abs(u) <= BIG) then
          ! trivial zero is included implicitly: tanh(0.5 * u) = sinh(0.5 * (u - 0)) / cosh(0.5 * u)
          c = cosh(0.5 * u)
          f = 0.25 * tanh(0.5 * u) / c**2
        else
          f = 0
        end if
      end if
    case (3)
      if (abs(u) <= BIG) then
        s1 = sinh(0.5 * (u + ZERO_3))
        s2 = sinh(0.5 * (u - ZERO_3))
        c  = cosh(0.5 * u)
        f  = -0.25 * (s1 / c) * (s2 / c) / c**2
      else
        f = 0
      end if
    case (4)
      if (abs(u) <= BIG) then
        ! trivial zero is included implicitly: tanh(0.5 * u) = sinh(0.5 * (u - 0)) / cosh(0.5 * u)
        s1 = sinh(0.5 * (u + ZERO_4))
        s2 = sinh(0.5 * (u - ZERO_4))
        c  = cosh(0.5 * u)
        f  = 0.25 * (s1 / c) * (s2 / c) * tanh(0.5 * u) / c**2
      else
        f = 0
      end if
    case (5)
      if (abs(u) <= BIG) then
        s1 = sinh(0.5 * (u + ZEROS_5(1)))
        s2 = sinh(0.5 * (u - ZEROS_5(1)))
        s3 = sinh(0.5 * (u + ZEROS_5(2)))
        s4 = sinh(0.5 * (u - ZEROS_5(2)))
        c  = cosh(0.5 * u)
        f  = - 0.25 * (s1 / c) * (s2 / c) * (s3 / c) * (s4 / c) / c**2
      else
        f = 0
      end if

    case default
      error stop "fermi-dirac not implemented for k = " // int2str(k)
    end select
  end function

  function maxwell_boltzmann(eta, k) result(F) bind(c, name = "ITHE_maxwell_boltzmann")
    !! k-th derivative of Maxwell-Boltzmann distribution function
    real(c_double), value, intent(in) :: eta
      !! normalized chemical potential
    integer(c_int), value, intent(in) :: k
      !! derivative (can be 0)
    real(c_double)                    :: F
      !! return k-th derivative of distribution function

    F = exp(eta)
  end function

  subroutine inv_maxwell_boltzmann(F, eta, deta_dF) bind(c, name = "ITHE_inv_maxwell_boltzmann")
    !! inverse of Maxwell-Boltzmann distribution function
    real(c_double), value, intent(in)  :: F
      !! value of Maxwell-Boltzmann distribution function
    real(c_double),        intent(out) :: eta
      !! output normalized chemical potential corresponding to F
    real(c_double),        intent(out) :: deta_dF
      !! output derivative of eta wrt F

    eta     = log(F)
    deta_dF = 1 / F
  end subroutine

  function int_maxwell_boltzmann(eta1, eta2, k) result(I) bind(c, name = "ITHE_int_maxwell_boltzmann")
    !! integrate Maxwell-Boltzmann distribution function exponentiated by k
    real(c_double), value, intent(in) :: eta1
      !! lower bound
    real(c_double), value, intent(in) :: eta2
      !! upper bound
    integer(c_int), value, intent(in) :: k
      !! integrand exponent
    real(c_double)                    :: I
      !! return value of integral

    if (k == 0) then
      I = eta2 - eta1
    else
      I = exp(k * eta1) * expm1(k * (eta2 - eta1)) / k
    end if
  end function



  subroutine init_fermi_dirac(eta_min, eta_max, logging) bind(c, name = "ITHE_init_fermi_dirac")
    real(c_double), value, intent(in) :: eta_min
      !! minimum supported eta
    real(c_double), value, intent(in) :: eta_max
      !! maximum supported eta
    integer(c_int), value, intent(in) :: logging
      !! logging level

    integer      :: k
    real(real64) :: par(1)

    if (logging >= LOGGING_INFO) print "(A)", "init_fermi_dirac"

    do k = lbound(itab_fd_1h, 1), ubound(itab_fd_1h, 1)
      if (k == 0) cycle

      par(1) = k
      call itab_fd_1h(k)%init(eta_min, eta_max, fermi_dirac_1h_integrand, par, 5e-15, 0.0, 1e-10, logging)
    end do
  end subroutine

  function fermi_dirac_1h(eta, k) result(F) bind(c, name = "ITHE_fermi_dirac_1h")
    !! k-th derivative of Fermi-Dirac integral of order 1/2
    real(c_double), value, intent(in) :: eta
      !! normalized chemical potential
    integer(c_int), value, intent(in) :: k
      !! derivative index (can be 0)
    real(c_double)                    :: F
      !! return k-th derivative of distribution function

    select case (k)
    case (0)
      F = fd_1h(eta)
    case (1)
      F = fd_m1h(eta)
    case (2)
      F = fd_m3h(eta)
    case (3)
      F = fd_m5h(eta)
    case (4)
      F = fd_m7h(eta)
    case (5)
      F = fd_m9h(eta)
    end select
  end function

  subroutine inv_fermi_dirac_1h(F, eta, deta_dF) bind(c, name = "ITHE_inv_fermi_dirac_1h")
    !! inverse of Fermi-Dirac integral of order 1/2
    real(c_double), value, intent(in)  :: F
      !! value of cumulative distribution function
    real(c_double),        intent(out) :: eta
      !! output normalized chemical potential corresponding to F
    real(c_double),        intent(out) :: deta_dF
      !! output derivative of eta wrt F

    eta     = inv_fd_1h(F)
    deta_dF = 1 / fd_m1h(eta)
  end subroutine

  function int_fermi_dirac_1h(eta1, eta2, k) result(I) bind(c, name = "ITHE_int_fermi_dirac_1h")
    !! integrate Fermi-Dirac integral of order 1/2 exponentiated by k
    real(c_double), value, intent(in) :: eta1
      !! lower bound
    real(c_double), value, intent(in) :: eta2
      !! upper bound
    integer(c_int), value, intent(in) :: k
      !! integrand exponent
    real(c_double)                    :: I
      !! return value of integral

    if (k == 0) then
      I = eta2 - eta1
    else
      I = itab_fd_1h(k)%get(fermi_dirac_1h_integrand, eta1, eta2, 5e-15)
    end if
  end function

  function fermi_dirac_1h_integrand(eta, num_par, par) result(F) bind(c)
    !! Fermi-Dirac integral of order 1/2 as integrand for itab_fd1h
    real(c_double), value, intent(in) :: eta
      !! inner argument
    integer(c_int), value, intent(in) :: num_par
      !! number of optional parameters (can be 0)
    real(c_double),        intent(in) :: par(num_par)
      !! optional parameters
    real(c_double)                    :: F
      !! return integrand

    integer :: k

    ! get F(eta)
    F = fermi_dirac_1h(eta, 0)

    ! apply exponent (stored in parameter array)
    k = nint(par(1))
    if (k /= 1) F = F**k
  end function



  subroutine init_band_tailing(t0, eta_min, eta_max, logging) bind(c, name = "ITHE_init_band_tailing")
    !! initialize band tailing distribution function
    real(c_double), value, intent(in) :: t0
      !! band-tailing energy
    real(c_double), value, intent(in) :: eta_min
      !! minimum supported eta
    real(c_double), value, intent(in) :: eta_max
      !! maximum supported eta
    integer(c_int), value, intent(in) :: logging
      !! logging level

    real(real64),  parameter :: RTOL      = 5e-15_real64
    real(real64),  parameter :: XTOL_CTAB = 1e-6_real64
    real(real64),  parameter :: XTOL_ITAB = 1e-10_real64

    integer      :: k
    real(real64) :: par(2)
    real(real128) :: intgrl

    if (logging >= LOGGING_INFO) print "(A)", "init_band_tailing"

    bt%t0 = t0
    bt%uc = - log(epsilon(1.0_real64))
    bt%eta_small = bt%t0 - bt%uc
    bt%eta_large = bt%t0 + bt%uc

    par(1) = t0
    do k = 0, 5
      par(2) = k
      call bt%ctab(k)%init(band_tailing_center, par, bt%eta_small, bt%eta_large, XTOL_CTAB, 2.0, RTOL, 16, logging)
    end do

    call quad(band_tailing_exp_integrand, real(-bt%uc, kind = real128), real(bt%uc, kind = real128), real([t0], kind = real128), intgrl)
    bt%A = 2 * sqrt(t0 / (PI * exp(1.0_real64))) * ( &
      2 * t0 * exp(-bt%uc/(2*t0)) * (1 - exp(-bt%uc)/(1 + 2 * t0)) + real(intgrl) - 2*t0/(1-2*t0) * exp(bt%uc*((1 - 2*t0)/(2*t0))))
    bt%B = 2 / (1 - 2 * t0) * sqrt(t0 / PI) * exp(-t0) + erfc(sqrt(t0))
    bt%C = 8 * sqrt(t0)**3 / (3 * sqrt(PI))
    bt%D = 2 * exp(t0) / (2 * t0 + 1) * sqrt(t0/PI) - erfi(sqrt(t0))

    bt%F_small = band_tailing(bt%eta_small, 0)
    bt%F_large = band_tailing(bt%eta_large, 0)

    do k = lbound(itab_fd_1h, 1), ubound(itab_fd_1h, 1)
      if (k == 0) cycle
      par(2) = k
      call bt%itab(k)%init(eta_min, eta_max, band_tailing_integrand, par, RTOL, 0.0, XTOL_ITAB, logging)
    end do
  end subroutine

  function band_tailing(eta, k) result(F) bind(c, name = "ITHE_band_tailing")
    !! k-th derivative of Fermi-Dirac integral with band-tailing
    real(c_double), value, intent(in) :: eta
      !! normalized chemical potential
    integer(c_int), value, intent(in) :: k
      !! derivative index (can be 0)
    real(c_double)                    :: F
      !! return k-th derivative of distribution function

    if (eta <= bt%eta_small) then
      F = bt%A * exp(eta / (2 * bt%t0)) / (2 * bt%t0)**k + bt%B * exp(eta)
    elseif (eta < bt%eta_large) then
      F = bt%ctab(k)%eval(eta)
    else ! eta >= eta_large
      F = fermi_dirac_1h(eta, k)
      if (mod(k, 2) == 0) then
        F = F + bt%D * exp(-eta)
      else
        F = F - bt%D * exp(-eta)
      end if
      if (k == 0) F = F + bt%C
    end if
  end function

  subroutine inv_band_tailing(F, eta, deta_dF) bind(c, name = "ITHE_inv_band_tailing")
    !! inverse of Fermi-Dirac integral with band-tailing
    real(c_double), value, intent(in)  :: F
      !! value of cumulative distribution function
    real(c_double),        intent(out) :: eta
      !! output normalized chemical potential corresponding to F
    real(c_double),        intent(out) :: deta_dF
      !! output derivative of eta wrt F

    integer,      parameter :: MIN_ITERATIONS = 1
    integer,      parameter :: MAX_ITERATIONS = 10
    real(real64), parameter :: RTOL = 5e-15

    integer      :: iteration
    real(real64) :: deta, err, G, dG_deta

    ! initial guess
    if (F <= bt%F_small) then
      ! eta <= eta_small
      if (bt%t0 < 0.5) then
        eta = log(F / bt%B)
      else
        eta = 2 * bt%t0 * log(F / bt%A)
      end if
    elseif (F < bt%F_large) then
      ! eta_small < eta < eta_large
      eta = bt%eta_small + (bt%eta_large - bt%eta_small) * log(F / bt%F_small) / log(bt%F_large / bt%F_small)
    else
      ! eta >= eta_large
      eta = inv_fd_1h(F - bt%C)
    end if

    ! Newton iteration
    iteration = 0
    err = huge(1.0_real64)
    do while ((iteration < MIN_ITERATIONS) .or. (err > RTOL * abs(eta)))
      iteration = iteration + 1
      if (iteration > MAX_ITERATIONS) then
        error stop "no convergence after "//int2str(MAX_ITERATIONS)//" iterations"
      end if

      G       = band_tailing(eta, 0)
      dG_deta = band_tailing(eta, 1)

      deta = (F - G) / dG_deta
      eta  = eta + deta
      err  = abs(deta)
    end do

    ! implicit differentiation
    deta_dF = 1 / dG_deta
  end subroutine

  function int_band_tailing(eta1, eta2, k) result(I) bind(c, name = "ITHE_int_band_tailing")
    !! integrate Fermi-Dirac integral with band-tailing exponentiated by k
    real(c_double), value, intent(in) :: eta1
      !! lower bound
    real(c_double), value, intent(in) :: eta2
      !! upper bound
    integer(c_int), value, intent(in) :: k
      !! integrand exponent
    real(c_double)                    :: I
      !! return value of integral

    if (k == 0) then
      I = eta2 - eta1
    else
      I = bt%itab(k)%get(band_tailing_integrand, eta1, eta2, 5e-15)
    end if
  end function

  function band_tailing_integrand(eta, num_par, par) result(F) bind(c)
    !! Fermi-Dirac integral with band-tailing as integrand for bt%itab
    real(c_double), value, intent(in) :: eta
      !! inner argument
    integer(c_int), value, intent(in) :: num_par
      !! number of optional parameters (can be 0)
    real(c_double),        intent(in) :: par(num_par)
      !! optional parameters
    real(c_double)                    :: F
      !! return integrand

    integer :: k

    ! get F(eta)
    F = band_tailing(eta, 0)

    ! apply exponent (stored in parameter array)
    k = nint(par(2))
    if (k /= 1) F = F**k
  end function

  function erfi(x) result(e)
    !! imaginary error function
    real(real64), intent(in) :: x
    real(real64)             :: e

    real(real128) :: x128, e128, par(0)

    x128 = x
    call quad(erfi_integrand, 0.0_real128, x128, par, e128)
    e = real(e128, real64)
  end function

  function erfi_integrand(x, n_par, par) result(f) bind(c)
    !! integrand for the imaginary error function
    real(c_float128), value, intent(in) :: x
      !! argument
    integer(c_int),   value, intent(in) :: n_par
      !! number of additional parameters: 0
    real(c_float128),        intent(in) :: par(n_par)
      !! additional parameters: none
    real(c_float128)                    :: f
      !! output integrand

    f = 2 * exp(x * x) / sqrt(PI_128)
  end function

  function band_tailing_exp_integrand(u, num_par, par) result(f) bind(c)
    real(c_float128), value, intent(in) :: u
      !! inner argument
    integer(c_int),   value, intent(in) :: num_par
      !! number of optional parameters (can be 0)
    real(c_float128),        intent(in) :: par(num_par)
      !! optional parameters
    real(c_float128) :: f
      !! return integrand

    real(c_float128) :: t0

    t0 = par(1)

    f = exp(u / (2 * t0)) / (exp(u) + 1)
  end function

  function band_tailing_center(eta, n_par, par) result(F)
    !! slow evaluation of F(eta) for t0 - uc < eta < t0 + uc using numerical quadrature
    real(c_double), value, intent(in) :: eta
      !! eta
    integer(c_int), value, intent(in) :: n_par
      !! number of additional parameters
    real(c_double),        intent(in) :: par(n_par)
      !! additional parameters
    real(c_double)                    :: F
      !! return F(eta)

    integer(c_int) :: k
    real(real128) :: t0, minf, pinf, qpar(3), Fe, F0

    t0 = par(1)
    k  = nint(par(2))

    ! perform numerical quadrature in quad precision
    minf = ieee_value(minf, ieee_negative_inf)
    pinf = ieee_value(pinf, ieee_positive_inf)
    qpar(1) = eta
    qpar(2) = t0
    qpar(3) = k
    call quad(band_tailing_Zf_exp,  minf,   t0, qpar, Fe)
    call quad(band_tailing_Zf_sqrt,   t0, pinf, qpar, F0)

    ! add both parts and convert down to double precision
    F = real(Fe + F0, c_double)
  end function

  function band_tailing_Zf_exp(t, n_par, par) result(Zf) bind(c)
    !! density of states times distribution density for exponential density of states
    real(c_float128), value, intent(in) :: t
      !! integration variable
    integer(c_int),   value, intent(in) :: n_par
      !! number of additional parameters: 3
    real(c_float128),        intent(in) :: par(n_par)
      !! additional parameters: [eta, t0, k]
    real(c_float128)                    :: Zf
      !! return Z(t) * f(t - eta)

    integer(c_int) :: k
    real(real128) :: eta, t0, Z, f

    eta = par(1)
    t0  = par(2)
    k   = nint(par(3))

    ! density of states
    Z = 2 * sqrt(t0 / PI_128) * exp((t - t0) / (2 * t0))

    ! k-th derivative of distribution density wrt u = t - eta
    f = fermi_dirac_density(t - eta, k)

    ! d/deta = (-1)**k * d/du with u
    if (mod(k, 2) /= 0) f = - f

    ! combine
    Zf = Z * f
  end function

  function band_tailing_Zf_sqrt(t, n_par, par) result(Zf) bind(c)
    !! density of states times distribution density for square root density of states
    real(c_float128), value, intent(in) :: t
      !! integration variable
    integer(c_int),   value, intent(in) :: n_par
      !! number of additional parameters: 3
    real(c_float128),        intent(in) :: par(n_par)
      !! additional parameters: [eta, t0, k]
    real(c_float128)                    :: Zf
      !! return Z(t) * f(t - eta)

    integer       :: k
    real(real128) :: eta, t0, Z, f

    eta = par(1)
    t0  = par(2)
    k   = nint(par(3))

    ! density of states
    Z = 2 * sqrt(t / PI_128)

    ! k-th derivative of distribution density wrt u = t - eta
    f = fermi_dirac_density(t - eta, k)

    ! d/deta = (-1)**k * d/du with u
    if (mod(k, 2) /= 0) f = - f

    ! combine
    Zf = Z * f
  end function


  subroutine init_gauss_fermi(sigma, eta_min, eta_max, logging) bind(c, name = "ITHE_init_gauss_fermi")
    !! initialize Gauss-Fermi distribution function
    real(c_double), value, intent(in) :: sigma
      !! width of gaussian distribution
    real(c_double), value, intent(in) :: eta_min
      !! minimum supported eta
    real(c_double), value, intent(in) :: eta_max
      !! maximum supported eta
    integer(c_int), value, intent(in) :: logging
      !! logging level

    real(real64), parameter :: RTOL = 5e-15_real64
    real(real64), parameter :: XTOL = 1e-6_real64

    integer      :: k
    real(real64) :: par(2)

    if (logging >= LOGGING_INFO) print "(A)", "init_gauss_fermi"

    gf%sigma   = sigma
    gf%eta_min = eta_min
    gf%eta_max = eta_max

    par(1) = sigma
    do k = 0, 5
      par(2) = k
      call gf%ctab(k)%init(band_tailing_center, par, eta_min, eta_max, XTOL, 2.0, RTOL, 16, logging)
    end do

    do k = lbound(itab_fd_1h, 1), ubound(itab_fd_1h, 1)
      if (k == 0) cycle
      par(2) = k
      call gf%itab(k)%init(eta_min, eta_max, band_tailing_integrand, par, RTOL, 0.0, 1e-10, logging)
    end do
  end subroutine

  function gauss_fermi(eta, k) result(F) bind(c, name = "ITHE_gauss_fermi")
    !! k-th derivative of Gauss-Fermi integral
    real(c_double), value, intent(in) :: eta
      !! normalized chemical potential
    integer(c_int), value, intent(in) :: k
      !! derivative index (can be 0)
    real(c_double)                    :: F
      !! return k-th derivative of distribution function

    F = gf%ctab(k)%eval(eta)
  end function

  subroutine inv_gauss_fermi(F, eta, deta_dF) bind(c, name = "ITHE_inv_gauss_fermi")
    !! inverse of Gauss-Fermi integral
    real(c_double), value, intent(in)  :: F
      !! value of cumulative distribution function
    real(c_double),        intent(out) :: eta
      !! output normalized chemical potential corresponding to F
    real(c_double),        intent(out) :: deta_dF
      !! output derivative of eta wrt F

    integer,      parameter :: MIN_ITERATIONS = 1
    integer,      parameter :: MAX_ITERATIONS = 10
    real(real64), parameter :: RTOL = 5e-15

    integer      :: iteration
    real(real64) :: deta, err, G, dG_deta

    ! initial guess
    if (F <= 0.5_real64) then
      eta = log(F)
    else
      eta = 10
    end if

    ! Newton iteration
    iteration = 0
    err = huge(1.0_real64)
    do while ((iteration < MIN_ITERATIONS) .or. (err > RTOL * abs(eta)))
      iteration = iteration + 1
      if (iteration > MAX_ITERATIONS) then
        error stop "no convergence after "//int2str(MAX_ITERATIONS)//" iterations"
      end if

      G       = gauss_fermi(eta, 0)
      dG_deta = gauss_fermi(eta, 1)

      deta = (F - G) / dG_deta
      eta  = eta + deta
      err  = abs(deta)
    end do

    ! implicit differentiation
    deta_dF = 1 / dG_deta
  end subroutine

  function int_gauss_fermi(eta1, eta2, k) result(I) bind(c, name = "ITHE_int_gauss_fermi")
    !! integrate Fermi-Dirac integral with band-tailing exponentiated by k
    real(c_double), value, intent(in) :: eta1
      !! lower bound
    real(c_double), value, intent(in) :: eta2
      !! upper bound
    integer(c_int), value, intent(in) :: k
      !! integrand exponent
    real(c_double)                    :: I
      !! return value of integral

    if (k == 0) then
      I = eta2 - eta1
    else
      I = gf%itab(k)%get(gauss_fermi_integrand, eta1, eta2, 5e-15)
    end if
  end function

  function gauss_fermi_integrand(eta, num_par, par) result(F) bind(c)
    !! Gauss-Fermi integral with band-tailing as integrand for gf%itab
    real(c_double), value, intent(in) :: eta
      !! inner argument
    integer(c_int), value, intent(in) :: num_par
      !! number of optional parameters (can be 0)
    real(c_double),        intent(in) :: par(num_par)
      !! optional parameters
    real(c_double)                    :: F
      !! return integrand

    integer :: k

    ! get F(eta)
    F = gauss_fermi(eta, 0)

    ! apply exponent (stored in parameter array)
    k = nint(par(2))
    if (k /= 1) F = F**k
  end function

end module
