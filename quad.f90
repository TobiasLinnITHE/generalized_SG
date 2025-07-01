! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

module quad_m
  !! general purpose integration module (tanh-sinh quadrature)

  use ieee_arithmetic, only: ieee_is_finite

  implicit none

  private
  public quad

  interface
    subroutine integrand(x, p, f, dfdx, dfdp)
      real, intent(in)  :: x
        !! argument
      real, intent(in)  :: p(:)
        !! parameters
      real, intent(out) :: f
        !! output integrand
      real, intent(out) :: dfdx
        !! output derivative of f wrt x
      real, intent(out) :: dfdp(:)
        !! output derivatives of f wrt p
    end subroutine

    subroutine integrand_16(x, p, f, dfdx, dfdp)
      real(kind=16), intent(in)  :: x
        !! argument
      real(kind=16), intent(in)  :: p(:)
        !! parameters
      real(kind=16), intent(out) :: f
        !! output integrand
      real(kind=16), intent(out) :: dfdx
        !! output derivative of f wrt x
      real(kind=16), intent(out) :: dfdp(:)
        !! output derivatives of f wrt p
    end subroutine
  end interface

  interface quad
    module procedure :: quad_8, quad_16
  end interface

  integer, parameter :: MODE_TANH_SINH = 0
  integer, parameter :: MODE_EXP_SINH  = 1
  integer, parameter :: MODE_SINH_SINH = 2

contains

  subroutine quad_8(func, a, b, p, I, dIda, dIdb, dIdp, rtol, err, min_levels, max_levels, ncalls)
    !! general purpose integration routine using tanh-sinh (finite interval), exp-sinh (one integration bound is infinite)
    !! or sinh-sinh (both integration bounds are infinite) quadrature using quad precision
    procedure(integrand) :: func
      !! function to integrate
    real,              intent(in)  :: a
      !! lower integration bound
    real,              intent(in)  :: b
      !! upper integration bound
    real,              intent(in)  :: p(:)
      !! additional parameters passed to the integrand
    real,              intent(out) :: I
      !! output value of integral
    real,              intent(out) :: dIda
      !! output derivative of I wrt lower bound
    real,              intent(out) :: dIdb
      !! output derivative of I wrt upper bound
    real,              intent(out) :: dIdp(:)
      !! output derivatives of I wrt parameters
    real,    optional, intent(in)  :: rtol
      !! error tolerance, overestimate (default: 1e-13)
    real,    optional, intent(out) :: err
      !! output error estimate
    integer, optional, intent(in)  :: min_levels
      !! minimum number of levels (default: 2)
    integer, optional, intent(in)  :: max_levels
      !! maximal number of levels (default: 16)
    integer, optional, intent(out) :: ncalls
      !! output number of function calls

    integer, parameter :: N_INIT = 5
      !! 2 * N_INIT + 1 is the minimum number of points for level 0

    integer :: min_levels_, max_levels_, ncalls_, sgn, mode, level, dir, n0, n, nstep, nmax, const_counter
    real    :: bnd(2), rtol_, dIdbnd(2), scale, h, dpartdbnd(2), dpartdp(size(p)), g, gmax
    real    :: dsumdbnd(2), dsumdp(size(p)), sum, sum_err, x, x_old, dxdbnd(2), xref
    real    :: f, f_tmp, dfdx, dfdbnd(2), dfdp(size(p)), dfdp_tmp(size(p))
    real    :: part, part_old, eh, enh, enh_step, w, r, s

    ! optional arguments
    min_levels_ = 2
    max_levels_ = 16
    rtol_       = 1e-13
    if (present(min_levels)) min_levels_ = min_levels
    if (present(max_levels)) max_levels_ = max_levels
    if (present(rtol)) rtol_ = rtol

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
      mode   = MODE_TANH_SINH
      scale  = 0.5 * (bnd(2) - bnd(1))
      x      = 0.5 * (bnd(1) + bnd(2)) ! x = (a + b) / 2 + (b - a) / 2 * tanh(sinh(0)) = (a + b) / 2
      dxdbnd = [0.5, 0.5]
    elseif (ieee_is_finite(bnd(1))) then
      mode   = MODE_EXP_SINH
      scale  = 1
      xref   = bnd(1)
      x      = xref + scale ! x = a + exp(sinh(0)) = a + 1
      dxdbnd = [1.0, 0.0]
    elseif (ieee_is_finite(bnd(2))) then
      mode   = MODE_EXP_SINH
      scale  = -1
      xref   = bnd(2)
      x      = xref + scale ! x = b - exp(sinh(0)) = b - 1
      dxdbnd = [0.0, 1.0]
    else
      mode   = MODE_SINH_SINH
      scale  = 1
      x      = 0 ! x = sinh(sinh(0)) = 0
      dxdbnd = 0
    end if

    ! start sum with first point, w(t=0) = 1 regardless of integration mode
    call func(x, p, f, dfdx, dfdp)
    ncalls_  = 1
    sum      = f
    dsumdbnd = dfdx * dxdbnd
    dsumdp   = dfdp
    sum_err  = huge(1.0)

    ! loop over levels
    level   = 0
    h       = 2      ! twice initial step size
    nstep   = 1      ! n stepsize (1 for level == 0, 2 for level > 0)
    gmax    = abs(f) ! initial maximum
    nmax    = 0      ! initial maximum position
    do while ((level < min_levels_) .or. ((sum_err > rtol_ * abs(sum)) .and. (level <= max_levels_)))
      h  = 0.5 * h
      eh = exp(h)

      ! start n at maximum of previous level
      n0 = 2 * nmax

      ! reset temporary values
      part      = 0
      dpartdbnd = 0
      dpartdp   = 0

      ! loop over two directions
      do dir = -1, 1, 2
        n   = n0 + dir ! n +/- 1
        enh = eh**n    ! exp(n * h)

        ! step size (positive or negative depending on direction)
        if (level == 0) then
          nstep = dir
        else
          ! jump over every second point for higher levels (these are already included in the previous level)
          nstep = 2 * dir
        end if

        ! exp(n*h) update factor
        enh_step = eh**nstep

        ! reset temporary values
        x      = 0
        f      = 0
        dfdbnd = 0
        dfdp   = 0

        ! evaluate partial sum
        const_counter = 0
        do while (((level == 0) .and. (abs(n - n0) < N_INIT)) .or. (const_counter < 2))
          x_old = x

          ! get nodes and weights depending on integration mode
          select case (mode)
          case (MODE_TANH_SINH)
            r      = exp(1 / enh - enh)            ! r = exp(-2 * sinh(n * h))
            s      = 2 * r / (1 + r)               ! s = 1 - tanh(sinh(n * h))
            w      = (enh + 1 / enh) * s / (1 + r) ! w = cosh(n * h) / cosh(sinh(n * h))**2
            x      = bnd(2) - scale * s            ! x = (a+b)/2 + (b-a)/2 * tanh(sinh(n * h))
            dxdbnd = [0.0, 1.0] + [0.5, -0.5] * s
          case (MODE_EXP_SINH)
            r = exp(0.5 * (enh - 1 / enh))         ! r = exp(sinh(n * h))
            w = 0.5 * (enh + 1 / enh) * r          ! w = cosh(n * h) * exp(sinh(n * h))
            x = xref + scale * r                   ! x = a[b] +[-] exp(sinh(n * h))
          case (MODE_SINH_SINH)
            r = exp(0.5 * (enh - 1 / enh))         ! r = exp(sinh(n * h))
            x = 0.5 * (r - 1 / r)                  ! x = sinh(sinh(n * h))
            r = 0.5 * (r + 1 / r)                  ! r = cosh(sinh(n * h))
            w = 0.5 * (enh + 1 / enh) * r          ! w = cosh(n * h) * cosh(sinh(n * h))
          end select

          ! evaluate function (reuse value if too close to prev. x)
          if (x /= x_old) then
            call func(x, p, f_tmp, dfdx, dfdp_tmp)
            ncalls_ = ncalls_ + 1
            if (ieee_is_finite(f_tmp)) then
              f      = f_tmp
              dfdbnd = dfdx * dxdbnd
              dfdp   = dfdp_tmp
            end if
          end if

          ! actual integrand including weight from transformation
          g = w * f

          ! update maximum integrand
          if (abs(g) > gmax) then
            gmax = abs(g)
            nmax = n
          end if

          ! update partial sum
          part_old  = part
          part      = part + g
          if (part == part_old) const_counter = const_counter + 1
          dpartdbnd = dpartdbnd + w * dfdbnd
          dpartdp   = dpartdp   + w * dfdp

          ! update n and exp(n*h)
          n   = n + nstep
          enh = enh * enh_step
        end do
      end do

      ! update sum and estimate error (error is based on difference to prev. level => overestimate)
      sum_err  = abs(sum - part) ! corresponds to |I_{level} - I_{level-1}| >= rtol_ |I_{level}|
      sum      = sum   + part
      dsumdbnd = dsumdbnd + dpartdbnd
      dsumdp   = dsumdp   + dpartdp

      ! go to next level
      level = level + 1
    end do

    ! calculate integral from raw sum
    I      = sgn * scale * sum * h
    dIdbnd = sgn * scale * dsumdbnd * h
    if (mode == MODE_TANH_SINH) then
      dIdbnd = dIdbnd + sgn * [-0.5, 0.5] * sum * h
    end if
    dIdp   = sgn * scale * dsumdp * h
    if (sgn == 1) then
      dIda = dIdbnd(1)
      dIdb = dIdbnd(2)
    else
      dIda = dIdbnd(2)
      dIdb = dIdbnd(1)
    end if

    ! optional output
    if (present(err)) err = sum_err / (abs(sum) + 1e-16)
    if (present(ncalls)) ncalls = ncalls_
  end subroutine

  subroutine quad_16(func, a, b, p, I, dIda, dIdb, dIdp, rtol, err, min_levels, max_levels, ncalls)
    !! general purpose integration routine using tanh-sinh (finite interval), exp-sinh (one integration bound is infinite)
    !! or sinh-sinh (both integration bounds are infinite) quadrature using quad precision
    procedure(integrand_16)              :: func
      !! function to integrate
    real(kind=16),           intent(in)  :: a
      !! lower integration bound
    real(kind=16),           intent(in)  :: b
      !! upper integration bound
    real(kind=16),           intent(in)  :: p(:)
      !! additional parameters passed to the integrand
    real(kind=16),           intent(out) :: I
      !! output value of integral
    real(kind=16),           intent(out) :: dIda
      !! output derivative of I wrt lower bound
    real(kind=16),           intent(out) :: dIdb
      !! output derivative of I wrt upper bound
    real(kind=16),           intent(out) :: dIdp(:)
      !! output derivatives of I wrt parameters
    real(kind=16), optional, intent(in)  :: rtol
      !! error tolerance, overestimate (default: 1e-15)
    real(kind=16), optional, intent(out) :: err
      !! output error estimate
    integer,       optional, intent(in)  :: min_levels
      !! minimum number of levels (default: 2)
    integer,       optional, intent(in)  :: max_levels
      !! maximal number of levels (default: 16)
    integer,       optional, intent(out) :: ncalls
      !! output number of function calls

    integer, parameter :: N_INIT = 5
      !! 2 * N_INIT + 1 is the minimum number of points for level 0

    integer       :: min_levels_, max_levels_, ncalls_, sgn, mode, level, dir, n0, n, nstep, nmax, const_counter
    real(kind=16) :: bnd(2), rtol_, dIdbnd(2), scale, h, dpartdbnd(2), dpartdp(size(p)), g, gmax
    real(kind=16) :: dsumdbnd(2), dsumdp(size(p)), sum, sum_err, x, x_old, dxdbnd(2), xref
    real(kind=16) :: f, f_tmp, dfdx, dfdbnd(2), dfdp(size(p)), dfdp_tmp(size(p))
    real(kind=16) :: part, part_old, eh, enh, enh_step, w, r, s

    ! optional arguments
    min_levels_ = 2
    max_levels_ = 16
    rtol_       = 1e-15
    if (present(min_levels)) min_levels_ = min_levels
    if (present(max_levels)) max_levels_ = max_levels
    if (present(rtol)) rtol_ = rtol

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
      mode   = MODE_TANH_SINH
      scale  = 0.5 * (bnd(2) - bnd(1))
      x      = 0.5 * (bnd(1) + bnd(2)) ! x = (a + b) / 2 + (b - a) / 2 * tanh(sinh(0)) = (a + b) / 2
      dxdbnd = [0.5, 0.5]
    elseif (ieee_is_finite(bnd(1))) then
      mode   = MODE_EXP_SINH
      scale  = 1
      xref   = bnd(1)
      x      = xref + scale ! x = a + exp(sinh(0)) = a + 1
      dxdbnd = [1.0, 0.0]
    elseif (ieee_is_finite(bnd(2))) then
      mode   = MODE_EXP_SINH
      scale  = -1
      xref   = bnd(2)
      x      = xref + scale ! x = b - exp(sinh(0)) = b - 1
      dxdbnd = [0.0, 1.0]
    else
      mode   = MODE_SINH_SINH
      scale  = 1
      x      = 0 ! x = sinh(sinh(0)) = 0
      dxdbnd = 0
    end if

    ! start sum with first point, w(t=0) = 1 regardless of integration mode
    call func(x, p, f, dfdx, dfdp)
    ncalls_  = 1
    sum      = f
    dsumdbnd = dfdx * dxdbnd
    dsumdp   = dfdp
    sum_err  = huge(1.0_16)

    ! loop over levels
    level   = 0
    h       = 2      ! twice initial step size
    nstep   = 1      ! n stepsize (1 for level == 0, 2 for level > 0)
    gmax    = abs(f) ! initial maximum
    nmax    = 0      ! initial maximum position
    do while ((level < min_levels_) .or. ((sum_err > rtol_ * abs(sum)) .and. (level <= max_levels_)))
      h  = 0.5 * h
      eh = exp(h)

      ! start n at maximum of previous level
      n0 = 2 * nmax

      ! reset temporary values
      part      = 0
      dpartdbnd = 0
      dpartdp   = 0

      ! loop over two directions
      do dir = -1, 1, 2
        n   = n0 + dir ! n +/- 1
        enh = eh**n    ! exp(n * h)

        ! step size (positive or negative depending on direction)
        if (level == 0) then
          nstep = dir
        else
          ! jump over every second point for higher levels (these are already included in the previous level)
          nstep = 2 * dir
        end if

        ! exp(n*h) update factor
        enh_step = eh**nstep

        ! reset temporary values
        x      = 0
        f      = 0
        dfdbnd = 0
        dfdp   = 0

        ! evaluate partial sum
        const_counter = 0
        do while (((level == 0) .and. (abs(n - n0) < N_INIT)) .or. (const_counter < 2))
          x_old = x

          ! get nodes and weights depending on integration mode
          select case (mode)
          case (MODE_TANH_SINH)
            r      = exp(1 / enh - enh)            ! r = exp(-2 * sinh(n * h))
            s      = 2 * r / (1 + r)               ! s = 1 - tanh(sinh(n * h))
            w      = (enh + 1 / enh) * s / (1 + r) ! w = cosh(n * h) / cosh(sinh(n * h))**2
            x      = bnd(2) - scale * s            ! x = (a+b)/2 + (b-a)/2 * tanh(sinh(n * h))
            dxdbnd = [0.0, 1.0] + [0.5, -0.5] * s
          case (MODE_EXP_SINH)
            r = exp(0.5 * (enh - 1 / enh))         ! r = exp(sinh(n * h))
            w = 0.5 * (enh + 1 / enh) * r          ! w = cosh(n * h) * exp(sinh(n * h))
            x = xref + scale * r                   ! x = a[b] +[-] exp(sinh(n * h))
          case (MODE_SINH_SINH)
            r = exp(0.5 * (enh - 1 / enh))         ! r = exp(sinh(n * h))
            x = 0.5 * (r - 1 / r)                  ! x = sinh(sinh(n * h))
            r = 0.5 * (r + 1 / r)                  ! r = cosh(sinh(n * h))
            w = 0.5 * (enh + 1 / enh) * r          ! w = cosh(n * h) * cosh(sinh(n * h))
          end select

          ! evaluate function (reuse value if too close to prev. x)
          if (x /= x_old) then
            call func(x, p, f_tmp, dfdx, dfdp_tmp)
            ncalls_ = ncalls_ + 1
            if (ieee_is_finite(f_tmp)) then
              f      = f_tmp
              dfdbnd = dfdx * dxdbnd
              dfdp   = dfdp_tmp
            end if
          end if

          ! actual integrand including weight from transformation
          g = w * f

          ! update maximum integrand
          if (abs(g) > gmax) then
            gmax = abs(g)
            nmax = n
          end if

          ! update partial sum
          part_old  = part
          part      = part + g
          if (part == part_old) const_counter = const_counter + 1
          dpartdbnd = dpartdbnd + w * dfdbnd
          dpartdp   = dpartdp   + w * dfdp

          ! update n and exp(n*h)
          n   = n + nstep
          enh = enh * enh_step
        end do
      end do

      ! update sum and estimate error (error is based on difference to prev. level => overestimate)
      sum_err  = abs(sum - part) ! corresponds to |I_{level} - I_{level-1}| >= rtol_ |I_{level}|
      sum      = sum   + part
      dsumdbnd = dsumdbnd + dpartdbnd
      dsumdp   = dsumdp   + dpartdp

      ! go to next level
      level = level + 1
    end do

    ! calculate integral from raw sum
    I      = sgn * scale * sum * h
    dIdbnd = sgn * scale * dsumdbnd * h
    if (mode == MODE_TANH_SINH) then
      dIdbnd = dIdbnd + sgn * [-0.5, 0.5] * sum * h
    end if
    dIdp   = sgn * scale * dsumdp * h
    if (sgn == 1) then
      dIda = dIdbnd(1)
      dIdb = dIdbnd(2)
    else
      dIda = dIdbnd(2)
      dIdb = dIdbnd(1)
    end if

    ! optional output
    if (present(err)) err = sum_err / (abs(sum) + 1e-32)
    if (present(ncalls)) ncalls = ncalls_
  end subroutine

end module
