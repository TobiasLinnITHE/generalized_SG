! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

module quad_m
  !! general purpose integration module (double exponential integration)

  use, intrinsic :: iso_c_binding,   only: c_double, c_float128, c_int
  use, intrinsic :: iso_fortran_env, only: real64, real128
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  #:set precisions = [ 64, 128 ]
  #:set kinds      = { 64: "real64",   128: "real128"    }
  #:set c_kinds    = { 64: "c_double", 128: "c_float128" }

  private
  public quad
  #:for prec in precisions
  public integrand_${prec}$
  public integrand_derivatives_${prec}$
  #:endfor

  interface
    #:for prec in precisions
    #:set c_kind = c_kinds[prec]
    function integrand_${prec}$(x, n_par, par) result(f) bind(c)
      import ${c_kinds[prec]}$, c_int

      real(${c_kind}$), value, intent(in)  :: x
        !! argument
      integer(c_int),   value, intent(in)  :: n_par
        !! number of additional parameters
      real(${c_kind}$),        intent(in)  :: par(n_par)
        !! additional parameters
      real(${c_kind}$)                     :: f
        !! return integrand
    end function

    subroutine integrand_derivatives_${prec}$(x, n_par, par, f, df_dx, df_dpar) bind(c)
      import ${c_kinds[prec]}$, c_int

      real(${c_kind}$), value, intent(in)  :: x
        !! argument
      integer(c_int),   value, intent(in)  :: n_par
        !! number of additional parameters
      real(${c_kind}$),        intent(in)  :: par(n_par)
        !! additional parameters
      real(${c_kind}$),        intent(out) :: f
        !! output integrand
      real(${c_kind}$),        intent(out) :: df_dx
        !! output derivative of f wrt x
      real(${c_kind}$),        intent(out) :: df_dpar(n_par)
        !! output derivatives of f wrt p
    end subroutine
    #:endfor
  end interface

  interface quad
    #:for prec in precisions
    module procedure :: quad_${prec}$
    module procedure :: quad_derivatives_${prec}$
    #:endfor
  end interface

  integer, parameter :: MODE_TANH_SINH = 0
    !! tanh-sinh integration mode: x = (a + b)/2 + (b - a)/2 * tanh(sinh(t))
  integer, parameter :: MODE_EXP_SINH  = 1
    !! exp-sinh integration mode: x = a + exp(sinh(t)) or x = b - exp(sinh(t))
  integer, parameter :: MODE_SINH_SINH = 2
    !! sinh-sinh integration mode: x = sinh(sinh(t))

  integer,       parameter :: DEFAULT_MIN_LEVELS = 10
    !! default value for minimum number of refinement values
  integer,       parameter :: DEFAULT_MAX_LEVELS = 18
    !! default value for maximum number of refinement values
  real(real64),  parameter :: DEFAULT_RTOL_64 = 1e-12
    !! default value for relative error tolerance for double precision
  real(real128), parameter :: DEFAULT_RTOL_128 = real(epsilon(1.0_real64), real128)
    !! default value for relative error tolerance for quadruple precision (exact in double precision)

  integer, parameter :: N_INIT = 5
    !! 2 * N_INIT + 1 is the minimum number of quadrature points used for level 0

contains

  #:for prec in precisions
  #:set kind = kinds[prec]
  subroutine quad_${prec}$(intg, a, b, par, I, rtol, err, min_levels, max_levels, n_calls)
    !! general purpose integration routine with derivatives (double exponential integration)
    procedure(integrand_${prec}$)         :: intg
    real(${kind}$),           intent(in)  :: a
      !! lower bound (can be infinite)
    real(${kind}$),           intent(in)  :: b
      !! upper bound (can be infinite)
    real(${kind}$),           intent(in)  :: par(:)
      !! additional parameters passed to the integrand
    real(${kind}$),           intent(out) :: I
      !! output value of the integral
    real(${kind}$), optional, intent(in)  :: rtol
      !! relative error tolerance, overestimate (default: DEFAULT_RTOL)
    real(${kind}$), optional, intent(out) :: err
      !! output error estimate
    integer,        optional, intent(in)  :: min_levels
      !! minimum number of refinement levels (default: DEFAULT_MIN_LEVELS)
    integer,        optional, intent(in)  :: max_levels
      !! maximum number of refinement levels (default: DEFAULT_MAX_LEVELS)
    integer,        optional, intent(out) :: n_calls
      !! output number of calls to integrand

    integer        :: min_levels_, max_levels_, n_calls_
    integer        :: sgn, mode, level, dir, n0, n, n_step, n_max, const_counter
    real(${kind}$) :: bnd(2), rtol_, scale, h, g, g_max
    real(${kind}$) :: sum, sum_err, x, x_old, x_ref, x_sgn
    real(${kind}$) :: f, f_tmp
    real(${kind}$) :: part, part_old, eh, enh, enh_step, w, r, s

    ! optional arguments
    min_levels_ = DEFAULT_MIN_LEVELS
    max_levels_ = DEFAULT_MAX_LEVELS
    rtol_       = DEFAULT_RTOL_${prec}$
    if (present(min_levels)) min_levels_ = min_levels
    if (present(max_levels)) max_levels_ = max_levels
    if (present(rtol)) rtol_ = rtol

    ! input sanity checks
    if (rtol_ <= 0) error stop "rtol <= 0"
    if (min_levels_ <= 0) error stop "min_levels <= 0"
    if (max_levels_ <= min_levels_) error stop "max_levels <= min_levels"

    ! sort integration bounds
    if (a < b) then
      bnd = [a, b]
      sgn = 1
    else
      bnd = [b, a]
      sgn = -1
    end if

    ! determine integration mode and first node
    if (ieee_is_finite(bnd(1)) .and. ieee_is_finite(bnd(2))) then
      mode  = MODE_TANH_SINH
      scale = 0.5 * (bnd(2) - bnd(1))
      x     = 0.5 * (bnd(1) + bnd(2)) ! x = (a + b) / 2 + (b - a) / 2 * tanh(sinh(0)) = (a + b) / 2
    elseif (ieee_is_finite(bnd(1))) then
      mode  = MODE_EXP_SINH
      scale = 1
      x_sgn = 1
      x_ref = bnd(1)
      x     = x_ref + x_sgn ! x = a + exp(sinh(0)) = a + 1
    elseif (ieee_is_finite(bnd(2))) then
      mode  = MODE_EXP_SINH
      scale = 1
      x_sgn = -1
      x_ref = bnd(2)
      x     = x_ref + x_sgn ! x = b - exp(sinh(0)) = b - 1
    else
      mode  = MODE_SINH_SINH
      scale = 1
      x     = 0 ! x = sinh(sinh(0)) = 0
    end if

    ! start sum with first point, w(t=0) = 1 regardless of integration mode
    f        = intg(x, size(par), par)
    n_calls_ = 1
    sum      = f
    sum_err  = huge(sum_err)

    ! loop over levels
    level  = 0
    h      = 2      ! twice initial step size
    n_step = 1      ! n stepsize (1 for level == 0, 2 for level > 0)
    g_max  = abs(f) ! initial maximum
    n_max  = 0      ! initial maximum position
    do while ((level < min_levels_) .or. ((sum_err > rtol_ * abs(sum)) .and. (level <= max_levels_)))
      h  = 0.5 * h
      eh = exp(h)

      ! start n at maximum of previous level
      n0 = 2 * n_max

      ! reset partial sum value
      part = 0

      ! loop over two directions
      do dir = -1, 1, 2
        n   = n0 + dir ! n +/- 1
        enh = eh**n    ! exp(n * h)

        ! step size (positive or negative depending on direction)
        if (level == 0) then
          n_step = dir
        else
          ! jump over every second point for higher levels (these are already included in the previous level)
          n_step = 2 * dir
        end if

        ! exp(n*h) update factor
        enh_step = eh**n_step

        ! reset temporary values
        x = 0
        f = 0

        ! evaluate partial sum
        const_counter = 0
        do while (((level == 0) .and. (abs(n - n0) < N_INIT)) .or. (const_counter < 2))
          x_old = x

          ! get nodes and weights depending on integration mode
          select case (mode)
          case (MODE_TANH_SINH)
            r = exp(1 / enh - enh)            ! r = exp(-2 * sinh(n * h))
            s = 2 * r / (1 + r)               ! s = 1 - tanh(sinh(n * h))
            w = (enh + 1 / enh) * s / (1 + r) ! w = cosh(n * h) / cosh(sinh(n * h))**2
            x = bnd(2) - scale * s            ! x = (a+b)/2 + (b-a)/2 * tanh(sinh(n * h))
          case (MODE_EXP_SINH)
            r = exp(0.5 * (enh - 1 / enh))    ! r = exp(sinh(n * h))
            w = 0.5 * (enh + 1 / enh) * r     ! w = cosh(n * h) * exp(sinh(n * h))
            x = x_ref + x_sgn * r             ! x = a[b] +[-] exp(sinh(n * h))
          case (MODE_SINH_SINH)
            r = exp(0.5 * (enh - 1 / enh))    ! r = exp(sinh(n * h))
            x = 0.5 * (r - 1 / r)             ! x = sinh(sinh(n * h))
            r = 0.5 * (r + 1 / r)             ! r = cosh(sinh(n * h))
            w = 0.5 * (enh + 1 / enh) * r     ! w = cosh(n * h) * cosh(sinh(n * h))
          end select

          ! evaluate function (reuse value if too close to prev. x)
          if (x /= x_old) then
            f_tmp = intg(x, size(par), par)
            n_calls_ = n_calls_ + 1
            if (ieee_is_finite(f_tmp)) f = f_tmp
          end if

          ! actual integrand including weight from transformation
          g = w * f

          ! update maximum integrand
          if (abs(g) > g_max) then
            g_max = abs(g)
            n_max = n
          end if

          ! update partial sum
          part_old  = part
          part      = part + g
          if (part == part_old) const_counter = const_counter + 1

          ! update n and exp(n*h)
          n   = n + n_step
          enh = enh * enh_step
        end do
      end do

      ! update sum and estimate error (error is based on difference to prev. level => overestimate)
      sum_err = abs(sum - part) ! corresponds to |I_level - I_level-1| >= rtol_ |I_level|
      sum     = sum + part

      ! go to next level
      level = level + 1
    end do

    ! calculate integral from raw sum
    I = sgn * scale * sum * h

    ! optional output
    if (present(err)) err = sum_err / (abs(sum) + 1e-16)
    if (present(n_calls)) n_calls = n_calls_
  end subroutine
  #:endfor

  #:for prec in precisions
  #:set kind = kinds[prec]
  subroutine quad_derivatives_${prec}$(intg, a, b, par, I, dI_da, dI_db, dI_dpar, rtol, err, min_levels, max_levels, n_calls)
    !! general purpose integration routine with derivatives (double exponential integration)
    procedure(integrand_derivatives_${prec}$) :: intg
    real(${kind}$),               intent(in)  :: a
      !! lower bound (can be infinite)
    real(${kind}$),               intent(in)  :: b
      !! upper bound (can be infinite)
    real(${kind}$),               intent(in)  :: par(:)
      !! additional parameters passed to the integrand
    real(${kind}$),               intent(out) :: I
      !! output value of the integral
    real(${kind}$),               intent(out) :: dI_da
      !! output derivative of I wrt a
    real(${kind}$),               intent(out) :: dI_db
      !! output derivative of I wrt b
    real(${kind}$),               intent(out) :: dI_dpar(:)
      !! output derivatives of I wrt par
    real(${kind}$), optional,     intent(in)  :: rtol
      !! relative error tolerance, overestimate (default: DEFAULT_RTOL)
    real(${kind}$), optional,     intent(out) :: err
      !! output error estimate
    integer,        optional,     intent(in)  :: min_levels
      !! minimum number of refinement levels (default: DEFAULT_MIN_LEVELS)
    integer,        optional,     intent(in)  :: max_levels
      !! maximum number of refinement levels (default: DEFAULT_MAX_LEVELS)
    integer,        optional,     intent(out) :: n_calls
      !! output number of calls to integrand

    integer        :: min_levels_, max_levels_, n_calls_
    integer        :: sgn, mode, level, dir, n0, n, n_step, n_max, const_counter
    real(${kind}$) :: bnd(2), rtol_, dI_dbnd(2), scale, h, dpart_dbnd(2), dpart_dpar(size(par)), g, g_max
    real(${kind}$) :: dsum_dbnd(2), dsum_dpar(size(par)), sum, sum_err, x, x_old, dx_dbnd(2), x_ref, x_sgn
    real(${kind}$) :: f, f_tmp, df_dx, df_dbnd(2), df_dpar(size(par)), df_dpar_tmp(size(par))
    real(${kind}$) :: part, part_old, eh, enh, enh_step, w, r, s

    ! optional arguments
    min_levels_ = DEFAULT_MIN_LEVELS
    max_levels_ = DEFAULT_MAX_LEVELS
    rtol_       = DEFAULT_RTOL_${prec}$
    if (present(min_levels)) min_levels_ = min_levels
    if (present(max_levels)) max_levels_ = max_levels
    if (present(rtol)) rtol_ = rtol

    ! input sanity checks
    if (rtol_ <= 0) error stop "rtol <= 0"
    if (min_levels_ <= 0) error stop "min_levels <= 0"
    if (max_levels_ <= min_levels_) error stop "max_levels <= min_levels"

    ! sort integration bounds
    if (a < b) then
      bnd = [a, b]
      sgn = 1
    else
      bnd = [b, a]
      sgn = -1
    end if

    ! determine integration mode and first node
    if (ieee_is_finite(bnd(1)) .and. ieee_is_finite(bnd(2))) then
      mode    = MODE_TANH_SINH
      scale   = 0.5 * (bnd(2) - bnd(1))
      x       = 0.5 * (bnd(1) + bnd(2)) ! x = (a + b) / 2 + (b - a) / 2 * tanh(sinh(0)) = (a + b) / 2
      dx_dbnd = [0.5, 0.5]
    elseif (ieee_is_finite(bnd(1))) then
      mode    = MODE_EXP_SINH
      scale   = 1
      x_sgn   = 1
      x_ref   = bnd(1)
      x       = x_ref + x_sgn ! x = a + exp(sinh(0)) = a + 1
      dx_dbnd = [1.0, 0.0]
    elseif (ieee_is_finite(bnd(2))) then
      mode    = MODE_EXP_SINH
      scale   = 1
      x_sgn   = -1
      x_ref   = bnd(2)
      x       = x_ref + x_sgn ! x = b - exp(sinh(0)) = b - 1
      dx_dbnd = [0.0, 1.0]
    else
      mode    = MODE_SINH_SINH
      scale   = 1
      x       = 0 ! x = sinh(sinh(0)) = 0
      dx_dbnd = 0
    end if

    ! start sum with first point, w(t=0) = 1 regardless of integration mode
    call intg(x, size(par), par, f, df_dx, df_dpar)
    n_calls_  = 1
    sum       = f
    dsum_dbnd = df_dx * dx_dbnd
    dsum_dpar = df_dpar
    sum_err   = huge(sum_err)

    ! loop over levels
    level  = 0
    h      = 2      ! twice initial step size
    n_step = 1      ! n stepsize (1 for level == 0, 2 for level > 0)
    g_max  = abs(f) ! initial maximum
    n_max  = 0      ! initial maximum position
    do while ((level < min_levels_) .or. ((sum_err > rtol_ * abs(sum)) .and. (level <= max_levels_)))
      h  = 0.5 * h
      eh = exp(h)

      ! start n at maximum of previous level
      n0 = 2 * n_max

      ! reset temporary values
      part       = 0
      dpart_dbnd = 0
      dpart_dpar = 0

      ! loop over two directions
      do dir = -1, 1, 2
        n   = n0 + dir ! n +/- 1
        enh = eh**n    ! exp(n * h)

        ! step size (positive or negative depending on direction)
        if (level == 0) then
          n_step = dir
        else
          ! jump over every second point for higher levels (these are already included in the previous level)
          n_step = 2 * dir
        end if

        ! exp(n*h) update factor
        enh_step = eh**n_step

        ! reset temporary values
        x       = 0
        f       = 0
        df_dbnd = 0
        df_dpar = 0

        ! evaluate partial sum
        const_counter = 0
        do while (((level == 0) .and. (abs(n - n0) < N_INIT)) .or. (const_counter < 2))
          x_old = x

          ! get nodes and weights depending on integration mode
          select case (mode)
          case (MODE_TANH_SINH)
            r       = exp(1 / enh - enh)            ! r = exp(-2 * sinh(n * h))
            s       = 2 * r / (1 + r)               ! s = 1 - tanh(sinh(n * h))
            w       = (enh + 1 / enh) * s / (1 + r) ! w = cosh(n * h) / cosh(sinh(n * h))**2
            x       = bnd(2) - scale * s            ! x = (a+b)/2 + (b-a)/2 * tanh(sinh(n * h))
            dx_dbnd = [0.0, 1.0] + [0.5, -0.5] * s
          case (MODE_EXP_SINH)
            r = exp(0.5 * (enh - 1 / enh))          ! r = exp(sinh(n * h))
            w = 0.5 * (enh + 1 / enh) * r           ! w = cosh(n * h) * exp(sinh(n * h))
            x = x_ref + x_sgn * r                   ! x = a[b] +[-] exp(sinh(n * h))
          case (MODE_SINH_SINH)
            r = exp(0.5 * (enh - 1 / enh))          ! r = exp(sinh(n * h))
            x = 0.5 * (r - 1 / r)                   ! x = sinh(sinh(n * h))
            r = 0.5 * (r + 1 / r)                   ! r = cosh(sinh(n * h))
            w = 0.5 * (enh + 1 / enh) * r           ! w = cosh(n * h) * cosh(sinh(n * h))
          end select

          ! evaluate function (reuse value if too close to prev. x)
          if (x /= x_old) then
            call intg(x, size(par), par, f_tmp, df_dx, df_dpar_tmp)
            n_calls_ = n_calls_ + 1
            if (ieee_is_finite(f_tmp)) then
              f       = f_tmp
              df_dbnd = df_dx * dx_dbnd
              df_dpar = df_dpar_tmp
            end if
          end if

          ! actual integrand including weight from transformation
          g = w * f

          ! update maximum integrand
          if (abs(g) > g_max) then
            g_max = abs(g)
            n_max = n
          end if

          ! update partial sum
          part_old  = part
          part      = part + g
          if (part == part_old) const_counter = const_counter + 1
          dpart_dbnd = dpart_dbnd + w * df_dbnd
          dpart_dpar = dpart_dpar + w * df_dpar

          ! update n and exp(n*h)
          n   = n + n_step
          enh = enh * enh_step
        end do
      end do

      ! update sum and estimate error (error is based on difference to prev. level => overestimate)
      sum_err   = abs(sum - part) ! corresponds to |I_level - I_level-1| >= rtol_ |I_level|
      sum       = sum   + part
      dsum_dbnd = dsum_dbnd + dpart_dbnd
      dsum_dpar = dsum_dpar + dpart_dpar

      ! go to next level
      level = level + 1
    end do

    ! calculate integral from raw sum
    I       = sgn * scale * sum * h
    dI_dbnd = sgn * scale * dsum_dbnd * h
    if (mode == MODE_TANH_SINH) then
      dI_dbnd = dI_dbnd + sgn * [-0.5, 0.5] * sum * h
    end if
    dI_dpar = sgn * scale * dsum_dpar * h
    if (sgn == 1) then
      dI_da = dI_dbnd(1)
      dI_db = dI_dbnd(2)
    else
      dI_da = dI_dbnd(2)
      dI_db = dI_dbnd(1)
    end if

    ! optional output
    if (present(err)) err = sum_err / (abs(sum) + 1e-16)
    if (present(n_calls)) n_calls = n_calls_
  end subroutine
  #:endfor

end module
