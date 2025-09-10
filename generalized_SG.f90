! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

module generalized_SG_m
  !! get edge-current by solving integral equation

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_next_after
  use, intrinsic :: iso_c_binding,   only: c_double, c_int
  use, intrinsic :: iso_fortran_env, only: real64, real128

  use distribution_table_m, only: distribution_table
  use quad_m,               only: quad
  use util_m,               only: int2str, expm1

  implicit none

  private
  public :: distribution, inv_distribution, DEBUG_GENERALIZED_SG, get_current_tab, get_current

  interface
    subroutine distribution(eta, k, F, dFdeta)  bind(c)
      import :: c_double, c_int

      !! k-th derivative of cumulative distribution function
      real(c_double), value, intent(in)  :: eta
        !! normalized chemical potential
      integer(c_int), value, intent(in)  :: k
        !! derivative (can be 0)
      real(c_double),        intent(out) :: F
        !! output k-th derivative of distribution function
      real(c_double),        intent(out) :: dFdeta
        !! output derivative of F wrt eta (used for Newton iteration, might be slightly different from F with k + 1)
    end subroutine

    subroutine inv_distribution(F, eta, detadF)  bind(c)
      import :: c_double

      !! inverse of cumulative distribution function
      real(c_double), value, intent(in)  :: F
        !! value of cumulative distribution function
      real(c_double),        intent(out) :: eta
        !! output normalized chemical potential corresponding to F
      real(c_double),        intent(out) :: detadF
        !! output derivative of eta wrt F
    end subroutine
  end interface

  logical :: DEBUG_GENERALIZED_SG = .false.

  integer,      parameter :: CASE0A  = 1
  integer,      parameter :: CASE0B  = 2
  integer,      parameter :: CASE0C  = 3
  integer,      parameter :: CASE1A  = 4
  integer,      parameter :: CASE1B  = 5
  integer,      parameter :: CASE2A  = 6
  integer,      parameter :: CASE2B  = 7
  character(2), parameter :: CASENAME(7) = ["0A", "0B", "0C", "1a", "1b", "2a", "2b"]
  real(real64), parameter :: CASETOL = 1e-3

  integer,      parameter :: MIN_IT = 6
  integer,      parameter :: MAX_IT = 50
  real(real64), parameter :: RTOL   = 2e-14
  real(real64), parameter :: ATOL   = 1e-16

contains

  subroutine get_current_tab(tab, n, dpot, j, djdn, djddpot)
    !! get edge current, distribution function given by table
    type(distribution_table), intent(in)  :: tab
      !! lookup table for cummulative distribution function
    real(real64),             intent(in)  :: n(2)
      !! normalized densities at left/right edge end-point (same as F)
    real(real64),             intent(in)  :: dpot
      !! normalized potential drop
    real(real64),             intent(out) :: j
      !! output normalized current density on edge
    real(real64),             intent(out) :: djdn(2)
      !! output derivative of j wrt n
    real(real64),             intent(out) :: djddpot
      !! output derivative of j wrt dpot

    call get_current(dist, idist, n, dpot, j, djdn, djddpot)

  contains

    subroutine dist(eta, k, F, dFdeta)
      real(c_double), value, intent(in)  :: eta
      integer(c_int), value, intent(in)  :: k
      real(c_double),        intent(out) :: F
      real(c_double),        intent(out) :: dFdeta

      call tab%get(eta, k, F, dFdeta)
    end subroutine

    subroutine idist(F, eta, detadF)
      real(c_double), value, intent(in)  :: F
      real(c_double),        intent(out) :: eta
      real(c_double),        intent(out) :: detadF

      call tab%inv(F, eta, detadF)
    end subroutine

  end subroutine

  recursive subroutine get_current(dist, idist, n, dpot, j, djdn, djddpot)
    !! get edge current
    procedure(distribution)     :: dist
      !! cumulative distribution function
    procedure(inv_distribution) :: idist
      !! inverse of cumulative distribution function
    real(real64),   intent(in)  :: n(2)
      !! normalized densities at left/right edge end-point (same as F)
    real(real64),   intent(in)  :: dpot
      !! normalized potential drop
    real(real64),   intent(out) :: j
      !! output normalized current density on edge
    real(real64),   intent(out) :: djdn(2)
      !! output derivative of j wrt n
    real(real64),   intent(out) :: djddpot
      !! output derivative of j wrt dpot

    integer      :: cs, it
    logical      :: lfinite, rfinite
    real(real64) :: eta(2), detadn(2), deta, djdeta(2), tmp, dtmp(2)
    real(real64) :: jj, jjmin, jjmax, jj_old, djj, err, dpot2
    real(real64) :: res, res_old, dresdjj, dresdeta(2), dresddpot

    ! flip edge direction if potential drop is negative
    if (dpot < 0) then
      call get_current(dist, idist, [n(2), n(1)], - dpot, j, djdn, djddpot)

      ! flip current (dj/ddpot = d(-j)/d(-dpot) unchanged)
      j      = - j
      djdn = - djdn(2:1:-1)
      return
    end if

    ! get eta
    call idist(n(1), eta(1), detadn(1))
    call idist(n(2), eta(2), detadn(2))
    deta = eta(2) - eta(1)

    ! get case
    if (abs(dpot) > CASETOL) then
      if (deta - dpot > dpot * CASETOL) then
        cs = CASE0A
      elseif (deta - dpot >= - dpot * CASETOL) then
        cs = CASE1A
      elseif (deta > CASETOL) then
        cs = CASE0B
      elseif (deta >= - CASETOL) then
        cs = CASE1B
      else
        cs = CASE0C
      end if
    else
      if (abs(deta) > CASETOL) then
        cs = CASE2A
      else
        cs = CASE2B
      end if
    end if

    ! debug info
    if (DEBUG_GENERALIZED_SG) then
      print *
      print "(A,ES25.16E3)", "n1   = ", n(1)
      print "(A,ES25.16E3)", "n2   = ", n(2)
      print "(A,ES25.16E3)", "eta1 = ", eta(1)
      print "(A,ES25.16E3)", "eta2 = ", eta(2)
      print "(A,ES25.16E3)", "dpot = ", dpot
      print "(A,ES25.16E3)", "deta = ", deta
      print "(2A)",          "case ", CASENAME(cs)
    end if

    ! calculate current based on case
    select case (cs)
    case (CASE0A, CASE0B, CASE0C) ! general case
      ! range for jj (jj = j scaled by 1/dpot) by mean value theorem
      jjmin = abs(1 - deta / dpot) * min(n(1), n(2))
      jjmax = abs(1 - deta / dpot) * max(n(1), n(2))
      if (deta > dpot) then
        ! flip sign and order
        jj    = jjmin
        jjmin = - jjmax
        jjmax = - jj
      end if

      ! further reduce range if possible
      if (cs == CASE0A) then
        call integrate_dist(dist, [eta(2) - dpot, eta(2)], -1, tmp, dtmp)
        jjmin = max(jjmin, (dpot - deta) / tmp)
      elseif (cs == CASE0B) then
        jjmax = min(jjmax, n(1))
      elseif (cs == CASE0C) then
        jjmin = max(jjmin, n(1))
      end if

      ! initial value: enhanced diffusion approximation
      if (abs(deta) < 1e-6) then
        dpot2 = dpot / (n(1) * detadn(1))
      else
        dpot2 = dpot * log(n(2) / n(1)) / deta
      end if
      jj = - n(1) / expm1(-dpot2) - n(2) / expm1(dpot2)
      if ((jj < jjmin) .or. (jj > jjmax)) jj = 0.5 * (jjmin + jjmax)

      if (DEBUG_GENERALIZED_SG) then
        print "(A,ES25.16E3)", "jjmin = ", jjmin
        print "(A,ES25.16E3)", "jjmax = ", jjmax
        print "(A,ES25.16E3)", "jj    = ", jj
      end if

      ! Newton iteration
      err = huge(1.0)
      it  = 0
      do while ((it < MIN_IT) .or. (err > RTOL * abs(jj) + ATOL / abs(dpot)))
        it = it + 1
        if (it > MAX_IT) call error("No convergence after " // int2str(MAX_IT) // " iterations")

        ! evaluate residual and get Newton update
        call residual(jj, res, dresdjj, dresdeta, dresddpot)

        ! treat singularity of res at interval endpoints
        if (.not. ieee_is_finite(res)) then
          jj_old  = jj
          res_old = res

          jj = jjmin
          call residual(jj, res, dresdjj, dresdeta, dresddpot)
          lfinite = ieee_is_finite(res)

          jj = jjmax
          call residual(jj, res, dresdjj, dresdeta, dresddpot)
          rfinite = ieee_is_finite(res)

          if (lfinite .and. rfinite) then
            ! reduce range (singularity inside of interval)
            if (abs(jj - jjmin) < abs(jj - jjmax)) then
              jjmin   = jj
              lfinite = .false.
            else
              jjmax   = jj
              rfinite = .false.
            end if
          end if
          if (.not. lfinite .and. .not. rfinite) call error("residual not finite at both endpoints")

          jj = jj_old
          res = res_old
          do while (.not. ieee_is_finite(res))
            if (rfinite) then
              jj = ieee_next_after(jj, huge(1.0))
              jjmin = jj
            else
              jj = ieee_next_after(jj, -huge(1.0))
              jjmax = jj
            end if
            call residual(jj, res, dresdjj, dresdeta, dresddpot)
          end do
        end if

        ! Newton update
        djj = - res / dresdjj
        err = abs(djj)

        ! update bounds (assume monotonic behaviour)
        if (djj > 0) then
          jjmin = jj
        else
          jjmax = jj
        end if

        ! update solution
        jj_old = jj
        jj     = jj + djj

        ! bisection
        if ((jj < jjmin) .or. (jj > jjmax) .or. ((jj_old == jjmin) .and. (jj == jjmax))) then
          jj = 0.5 * (jjmin + jjmax)
          err = min(err, jjmax - jjmin)
          if (DEBUG_GENERALIZED_SG) then
            print "(A)", "bisection"
            print "(A,ES25.16E3)", "  jj_old = ", jj_old + djj
            print "(A,ES25.16E3)", "  jj_new = ", jj
          end if
        end if

        if (DEBUG_GENERALIZED_SG) then
          print "(I6,A,ES25.16E3,A,ES25.16E3,A,ES25.16E3)", it, ": jj = ", jj, "  +/-", err, "  ,", err / jj
        end if
      end do

      ! get current and derivatives with implicit differentiation
      call residual(jj, res, dresdjj, dresdeta, dresddpot)
      j       = dpot * jj
      djdeta  =    - dpot * dresdeta  / dresdjj
      djddpot = jj - dpot * dresddpot / dresdjj
    case (CASE1A)
      call case_1a(j, djdeta, djddpot)
    case (CASE1B)
      call case_1b(j, djdeta, djddpot)
    case (CASE2A)
      call case_2a(j, djdeta, djddpot)
    case (CASE2B)
      call case_2b(j, djdeta, djddpot)
    end select

    ! derivatives wrt densities
    djdn = djdeta * detadn

  contains

    subroutine case_1a(j, djdeta, djddpot)
      !! dpot >> 0 and deta ~ dpot
      real,       intent(out) :: j
        !! output normalized current density on edge
      real,       intent(out) :: djdeta(2)
        !! output derivative of j wrt eta
      real,       intent(out) :: djddpot
        !! output derivative of j wrt dpot

      integer       :: k
      real(real64)  :: I(-5:0), dI(2,-5:0), jc(1:5), djc(2,1:5)
      real(real128) :: I_16(-5:0), jc_16(1:5), djc_16(2,1:5), tmp_16, dtmp_16(2), j_16

      ! get I_k = integral_eta1^eta2 F(eta)^k deta for k = -5 to 0
      do k = -5, -1
        call integrate_dist(dist, eta, k, I(k), dI(:,k))
      end do
      I(0)    = deta
      dI(:,0) = [-1.0, 1.0]
      I_16 = I

      ! first Taylor coefficient
      jc_16(1)    = I_16(0) / I_16(-1)
      djc_16(:,1) = (dI(:,0) - dI(:,-1) * I_16(0) / I_16(-1)) / I_16(-1)

      ! second Taylor coefficient
      tmp_16   = I_16(-1)**2 - I_16(-2)*I_16(0)
      dtmp_16  = 2*I_16(-1)*dI(:,-1) - dI(:,-2)*I_16(0) - I_16(-2)*dI(:,0)
      jc_16(2) = tmp_16 / I_16(-1)**3
      djc_16(:,2) = (dtmp_16 - 3 * tmp_16 * dI(:,-1) / I_16(-1)) / I_16(-1)**3

      ! third Taylor coefficient
      tmp_16   = 2*I_16(-2)**2*I_16(0) - I_16(-1)**2*I_16(-2) - I_16(-1)*I_16(-3)*I_16(0)
      dtmp_16  = 4*I_16(-2)*dI(:,-2)*I_16(0) + 2*I_16(-2)**2*dI(:,0) - 2*I_16(-1)*dI(:,-1)*I_16(-2) - I_16(-1)**2*dI(:,-2) &
        &      - dI(:,-1)*I_16(-3)*I_16(0) - I_16(-1)*dI(:,-3)*I_16(0) - I_16(-1)*I_16(-3)*dI(:,0)
      jc_16(3) = tmp_16 / I_16(-1)**5
      djc_16(:,3) = (dtmp_16 - 5 * tmp_16 * dI(:,-1) / I_16(-1)) / I_16(-1)**5

      ! fourth Taylor coefficient
      tmp_16   = 2*I_16(-1)**2*I_16(-2)**2 - I_16(-1)**3*I_16(-3) - 5*I_16(-2)**3*I_16(0) &
        &      - I_16(-1)**2*I_16(-4)*I_16(0) + 5*I_16(-1)*I_16(-2)*I_16(-3)*I_16(0)
      dtmp_16  = 4*I_16(-1)*dI(:,-1)*I_16(-2)**2 + 4*I_16(-1)**2*I_16(-2)*dI(:,-2) - 3*I_16(-1)**2*dI(:,-1)*I_16(-3) &
        &      - I_16(-1)**3*dI(:,-3) - 15*I_16(-2)**2*dI(:,-2)*I_16(0) - 5*I_16(-2)**3*dI(:,0) &
        &      - 2*I_16(-1)*dI(:,-1)*I_16(-4)*I_16(0) - I_16(-1)**2*dI(:,-4)*I_16(0) - I_16(-1)**2*I_16(-4)*dI(:,0) &
        &      + 5*dI(:,-1)*I_16(-2)*I_16(-3)*I_16(0) + 5*I_16(-1)*dI(:,-2)*I_16(-3)*I_16(0) &
        &      + 5*I_16(-1)*I_16(-2)*dI(:,-3)*I_16(0) + 5*I_16(-1)*I_16(-2)*I_16(-3)*dI(:,0)
      jc_16(4) = tmp_16 / I_16(-1)**7
      djc_16(:,4) = (dtmp_16 - 7 * tmp_16 * dI(:,-1) / I_16(-1)) / I_16(-1)**7

      ! fifth Taylor coefficient
      tmp_16   = 14*I_16(0)*I_16(-2)**4 - I_16(-1)**4*I_16(-4) - 5*I_16(-1)**2*I_16(-2)**3 &
        &      - I_16(0)*I_16(-1)**3*I_16(-5) + 5*I_16(-1)**3*I_16(-2)*I_16(-3) + 3*I_16(0)*I_16(-1)**2*I_16(-3)**2 &
        &      - 21*I_16(0)*I_16(-1)*I_16(-2)**2*I_16(-3) + 6*I_16(0)*I_16(-1)**2*I_16(-2)*I_16(-4)
      dtmp_16  = 14*dI(:,0)*I_16(-2)**4 + 56*I_16(0)*I_16(-2)**3*dI(:,-2) - 4*I_16(-1)**3*dI(:,-1)*I_16(-4) &
        &      - I_16(-1)**4*dI(:,-4) - 10*I_16(-1)*dI(:,-1)*I_16(-2)**3 - 15*I_16(-1)**2*I_16(-2)**2*dI(:,-2) &
        &      - dI(:,0)*I_16(-1)**3*I_16(-5) - 3*I_16(0)*I_16(-1)**2*dI(:,-1)*I_16(-5) - I_16(0)*I_16(-1)**3*dI(:,-5) &
        &      + 15*I_16(-1)**2*dI(:,-1)*I_16(-2)*I_16(-3) + 5*I_16(-1)**3*dI(:,-2)*I_16(-3) &
        &      + 5*I_16(-1)**3*I_16(-2)*dI(:,-3) + 3*dI(:,0)*I_16(-1)**2*I_16(-3)**2 &
        &      + 6*I_16(0)*I_16(-1)*dI(:,-1)*I_16(-3)**2 + 6*I_16(0)*I_16(-1)**2*I_16(-3)*dI(:,-3) &
        &      - 21*dI(:,0)*I_16(-1)*I_16(-2)**2*I_16(-3) - 21*I_16(0)*dI(:,-1)*I_16(-2)**2*I_16(-3) &
        &      - 42*I_16(0)*I_16(-1)*I_16(-2)*dI(:,-2)*I_16(-3) - 21*I_16(0)*I_16(-1)*I_16(-2)**2*dI(:,-3) &
        &      + 6*dI(:,0)*I_16(-1)**2*I_16(-2)*I_16(-4) + 12*I_16(0)*I_16(-1)*dI(:,-1)*I_16(-2)*I_16(-4) &
        &      + 6*I_16(0)*I_16(-1)**2*dI(:,-2)*I_16(-4) + 6*I_16(0)*I_16(-1)**2*I_16(-2)*dI(:,-4)
      jc_16(5) = tmp_16 / I_16(-1)**9
      djc_16(:,5) = (dtmp_16 - 9 * tmp_16 * dI(:,-1) / I_16(-1)) / I_16(-1)**9

      ! evaluate Taylor series
      j_16     = 0
      djdeta  = 0
      djddpot = 0
      jc      = real(jc_16, kind = real64)
      djc     = real(djc_16, kind = real64)
      do k = 1, 5
        if (.not. ieee_is_finite(jc_16(k))) cycle
        j_16     = j_16 + jc_16(k) * (dpot - I_16(0))**k
        djdeta  = djdeta + (djc(:,k) * (dpot - I(0)) - k * jc(k) * dI(:,0)) * (dpot - I(0))**(k-1)
        djddpot = djddpot + jc(k) * k * (dpot - I(0))**(k-1)
      end do
      j = real(j_16, kind = real64)
    end subroutine

    subroutine case_1b(j, djdeta, djddpot)
      !! |dpot| >> 0 and deta ~ 0
      real(real64), intent(out) :: j
        !! output normalized current density on edge
      real(real64), intent(out) :: djdeta(2)
        !! output derivative of j wrt eta
      real(real64), intent(out) :: djddpot
        !! output derivative of j wrt dpot

      integer      :: k
      real(real64) :: etam, F(0:3), dF(0:3), G(2), dG(2), e, dedetam, deddpot, u, dudetam, duddpot, v, dvdetam, dvddpot
      real(real64) :: jc(0:3), djcdetam(0:3), djcddpot(0:3)

      ! evaluate distribution and derivatives at mean eta
      etam = 0.5 * (eta(1) + eta(2))
      do k = 0, 3
        call dist(etam, k, F(k), dF(k))
      end do

      ! F(1) and F(2) scaled by 1/F(0)
      do k = 1, 2
        G(k) = F(k) / F(0)
        dG(k) = (dF(k) - G(k) * dF(0)) / F(0)
      end do

      ! abbreviation
      e = 1 / tanh(0.5 * G(1) * dpot)
      u = 1 / sinh(0.5 * G(1) * dpot)**2
      dedetam = - 0.5 * dpot * dG(1) * u
      deddpot = - 0.5 * G(1) * u

      ! zero-th Taylor coefficient
      jc(0)       = dpot * F(0)
      djcdetam(0) = dpot * dF(0)
      djcddpot(0) = F(0)

      ! first Taylor coefficient
      jc(1)       = - 0.5 * F(1) * dpot * e
      djcdetam(1) = - 0.5 * dpot * (dF(1) * e + F(1) * dedetam)
      djcddpot(1) = - 0.5 * F(1) * (e + dpot * deddpot)

      ! second Taylor coefficient = 0.25 * dpot * (0.5 * F(2) + (F(2) - F(1)**2/F(0)) * (1 - e**2) * (0.5*dpot*F(1)/F(0)*e - 1))
      u           = (F(2) - F(1) * G(1)) * (1 - e**2)
      dudetam     = (dF(2) - dF(1) * G(1) - F(1) * dG(1)) * (1 - e**2) - 2 * (F(2) - F(1) * G(1)) * e * dedetam
      duddpot     = - 2 * (F(2) - F(1) * G(1)) * e * deddpot
      v           = 0.5 * dpot * G(1) * e - 1
      dvdetam     = 0.5 * dpot * (dG(1) * e + G(1) * dedetam)
      dvddpot     = 0.5 * G(1) * (e + dpot * deddpot)
      jc(2)       = 0.25 * dpot * (0.5 * F(2) + u * v)
      djcdetam(2) = 0.25 * dpot * (0.5 * dF(2) + dudetam * v + u * dvdetam)
      djcddpot(2) = 0.25 * (0.5 * F(2) + u * v + dpot * (duddpot * v + u * dvddpot))

      ! third Taylor coefficient
      u           = e * (1 - e**2)
      dudetam     = (1 - 3*e**2) * dedetam
      duddpot     = (1 - 3*e**2) * deddpot
      v           = G(1)**2 * F(1) + 0.25 * (F(2)/F(1))*F(2) - 1.5 * G(1) * F(2)
      dvdetam     = G(1) * (2*F(1)*dG(1) + G(1)*dF(1)) + 0.25 * (F(2)/F(1)) * (2*dF(2) - (F(2)/F(1)) * dF(1)) - 1.5 * (dG(1) * F(2) + G(1) * dF(2))
      jc(3)       = 0.25 * dpot * u * v
      djcdetam(3) = 0.25 * dpot * (dudetam * v + u * dvdetam)
      djcddpot(3) = 0.25 * (u + dpot * duddpot) * v
      dudetam     = dudetam * (1 - 2 * e**2) - 4 * u * e * dedetam
      duddpot     = duddpot * (1 - 2 * e**2) - 4 * u * e * deddpot
      u           = u * (1 - 2 * e**2)
      v           = G(1) * (G(1)**2*F(2) - 0.5*G(2)*F(2) - 0.5*G(1)**3*F(1))
      dvdetam     = G(1)**2*(3*F(2)*dG(1) + G(1)*dF(2)) - 0.5*(dG(1)*G(2)*F(2) + G(1)*dG(2)*F(2) + G(1)*G(2)*dF(2)) - 0.5*G(1)**3*(4*dG(1)*F(1) + G(1)*dF(1))
      jc(3)       = jc(3) + dpot**3/16 * u * v
      djcdetam(3) = djcdetam(3) + dpot**3/16 * (dudetam * v + u * dvdetam)
      djcddpot(3) = djcddpot(3) + dpot**2/16 * (3 * u + dpot * duddpot) * v
      u           = (1 - e**2) * (1 - 5*e**2)
      dudetam     = - (12 - 20*e**2) * e * dedetam
      duddpot     = - (12 - 20*e**2) * e * deddpot
      v           = G(2) * F(2)
      dvdetam     = dG(2) * F(2) + G(2) * dF(2)
      jc(3)       = jc(3) + dpot**2/32 * u * v
      djcdetam(3) = djcdetam(3) + dpot**2/32 * (dudetam * v + u * dvdetam)
      djcddpot(3) = djcddpot(3) + dpot/32 * (2 * u + dpot * duddpot) * v
      u           = (1 - e**2) * (1 - 4*e**2)
      dudetam     = - (10 - 16*e**2) * e * dedetam
      duddpot     = - (10 - 16*e**2) * e * deddpot
      v           = G(1)**3 * F(1)
      dvdetam     = G(1)**2 * (3 * F(1) * dG(1) + G(1) * dF(1))
      jc(3)       = jc(3) + dpot**2/16 * u * v
      djcdetam(3) = djcdetam(3) + dpot**2/16 * (dudetam * v + u * dvdetam)
      djcddpot(3) = djcddpot(3) + dpot/16 * (2 * u + dpot * duddpot) * v
      u           = e * (2 - 3*e**2)
      dudetam     = (2 - 9*e**2) * dedetam
      duddpot     = (2 - 9*e**2) * deddpot
      jc(3)       = jc(3) + dpot/48 * u * F(3)
      djcdetam(3) = djcdetam(3) + dpot/48 * (dudetam * F(3) + u * dF(3))
      djcddpot(3) = djcddpot(3) + (u + dpot * duddpot) * F(3) / 48
      u           = (1 - e**2) * (1.5 - 7*e**2)
      dudetam     = - (17 - 28*e**2) * e * dedetam
      duddpot     = - (17 - 28*e**2) * e * deddpot
      v           = G(1)**2 * F(2)
      dvdetam     = G(1) * (2 * F(2) * dG(1) + G(1) * dF(2))
      jc(3)       = jc(3) - dpot**2/16 * u * v
      djcdetam(3) = djcdetam(3) - dpot**2/16 * (dudetam * v + u * dvdetam)
      djcddpot(3) = djcddpot(3) - dpot/16 * (2 * u + dpot * duddpot) * v
      u           = e**2 * (1 - e**2)
      dudetam     = 2 * e * (1 - 2*e**2) * dedetam
      duddpot     = 2 * e * (1 - 2*e**2) * deddpot
      v           = G(1) * F(3)
      dvdetam     = dG(1) * F(3) + G(1) * dF(3)
      jc(3)       = jc(3) - dpot**2/32 * u * v
      djcdetam(3) = djcdetam(3) - dpot**2/32 * (dudetam * v + u * dvdetam)
      djcddpot(3) = djcddpot(3) - dpot/32 * (2 * u + dpot * duddpot) * v

      ! evaluate Taylor series
      j = jc(0) + jc(1) * deta + jc(2) * deta**2 + jc(3) * deta**3
      djdeta(1) = 0.5 * (djcdetam(0) + djcdetam(1) * deta + djcdetam(2) * deta**2 + djcdetam(3) * deta**3) - jc(1) - 2 * jc(2) * deta - 3 * jc(3) * deta**2
      djdeta(2) = 0.5 * (djcdetam(0) + djcdetam(1) * deta + djcdetam(2) * deta**2 + djcdetam(3) * deta**3) + jc(1) + 2 * jc(2) * deta + 3 * jc(3) * deta**2
      djddpot   = djcddpot(0) + djcddpot(1) * deta + djcddpot(2) * deta**2 + djcddpot(3) * deta**3
    end subroutine

    subroutine case_2a(j, djdeta, djddpot)
      !! dpot ~ 0 and |deta| >> 0
      real(real64), intent(out) :: j
        !! output normalized current density on edge
      real(real64), intent(out) :: djdeta(2)
        !! output derivative of j wrt eta
      real(real64), intent(out) :: djddpot
        !! output derivative of j wrt dpot

      integer      :: k
      real(real64) :: I(5), dI(2,5), jc(0:4), djc(2,0:4), tmp, dtmp(2)

      ! get I_k = integral_eta1^eta2 F(eta)^k deta for k = 0 to 5
      do k = 1, 5
        call integrate_dist(dist, eta, k, I(k), dI(:,k))
      end do

      ! normalize I_{2..5} wrt I_1 to improve numerical stability (avoid underflows etc.)
      do k = 2, 5
        I(   k) = I(k) / I(1)
        dI(:,k) = (dI(:,k) - I(k) * dI(:,1)) / I(1)
      end do

      ! zero-th Taylor coefficient
      jc(0)    = - I(1)
      djc(:,0) = - dI(:,1)

      ! first Taylor coefficient
      jc(1)    = I(2)
      djc(:,1) = dI(:,2)

      ! second Taylor coefficient
      tmp      = I(2)**2 - I(3)
      dtmp     = 2*I(2)*dI(:,2) - dI(:,3)
      jc(2)    = tmp / I(1)
      djc(:,2) = (dtmp - jc(2) * dI(:,1)) / I(1)

      ! third Taylor coefficient
      tmp      = I(4) + 2*I(2)**3 - 3*I(2)*I(3)
      dtmp     = dI(:,4) + 6*I(2)**2*dI(:,2) - 3*dI(:,2)*I(3) - 3*I(2)*dI(:,3)
      jc(3)    = tmp / I(1)**2
      djc(:,3) = (dtmp - 2 * tmp * dI(:,1) / I(1)) / I(1)**2

      ! fourth Taylor coefficient
      tmp      = 5*I(2)**4 - I(5) + 2*I(3)**2 - 10*I(2)**2*I(3) + 4*I(2)*I(4)
      dtmp     = 20*I(2)**3*dI(:,2) - dI(:,5) + 4*I(3)*dI(:,3) - 20*I(2)*dI(:,2)*I(3) - 10*I(2)**2*dI(:,3) + 4*dI(:,2)*I(4) + 4*I(2)*dI(:,4)
      jc(4)    = tmp / I(1)**3
      djc(:,4) = (dtmp - 3 * tmp * dI(:,1) / I(1)) / I(1)**3

      ! evaluate Taylor series
      j       = jc(0)
      djdeta  = djc(:,0)
      djddpot = 0
      do k = 1, 4
        j       = j + jc(k) * dpot**k
        djdeta  = djdeta + djc(:,k) * dpot**k
        djddpot = djddpot + jc(k) * k * dpot**(k-1)
      end do
    end subroutine

    subroutine case_2b(j, djdeta, djddpot)
      !! dpot ~ 0 and deta ~ 0
      real(real64), intent(out) :: j
        !! output normalized current density on edge
      real(real64), intent(out) :: djdeta(2)
        !! output derivative of j wrt eta
      real(real64), intent(out) :: djddpot
        !! output derivative of j wrt dpot

      integer      :: k
      real(real64) :: etam, F(0:2), dF(0:2), tmp, dtmpdetam, dtmpddeta, dtmpddpot, djdetam, djddeta

      ! evaluate distribution and derivatives at mean eta
      etam = 0.5 * (eta(1) + eta(2))
      do k = 0, 2
        call dist(etam, k, F(k), dF(k))
      end do

      ! factor
      tmp       = F(0) - (F(1)/F(0)) * F(1) / 12 * deta * dpot + F(2) / 24 * deta**2
      dtmpdetam = dF(0) - (2 * dF(1) - F(1) * dF(0) / F(0)) * F(1) / (12 * F(0)) * deta * dpot
      dtmpddeta = - F(1)**2 / (12 * F(0)) * dpot + F(2) / 12 * deta
      dtmpddpot = - F(1)**2 / (12 * F(0)) * deta

      ! current
      j = (dpot - deta) * tmp
      djdetam = (dpot - deta) * dtmpdetam
      djddeta = - tmp + (dpot - deta) * dtmpddeta
      djdeta  = djdetam * [0.5, 0.5] + djddeta * [-1.0, 1.0]
      djddpot = tmp + (dpot - deta) * dtmpddpot
    end subroutine

    subroutine residual(jj, res, dresdjj, dresdeta, dresddpot)
      real(real64), intent(in)  :: jj
      real(real64), intent(out) :: res
      real(real64), intent(out) :: dresdjj
      real(real64), intent(out) :: dresdeta(2)
      real(real64), intent(out) :: dresddpot

      real(real64) :: dresdjj1(1)

      ! integrate using tanh-sinh
      call quad(integrand, eta(1), eta(2), [jj], res, dresdeta(1), dresdeta(2), dresdjj1, max_levels = 8)
      dresdjj = dresdjj1(1)

      ! subtract delta potential
      res       = res - dpot
      dresddpot = -1
    end subroutine

    subroutine integrand(eta, jj, u, dudeta, dudjj)
      real(real64), intent(in)  :: eta
      real(real64), intent(in)  :: jj(:)
      real(real64), intent(out) :: u
      real(real64), intent(out) :: dudeta
      real(real64), intent(out) :: dudjj(:)

      real(real64) :: t

      call dist(eta, 0, u, dudeta)
      t        = u - jj(1)
      u        = u / t
      dudeta   = - (jj(1) / t) * (dudeta / t)
      dudjj(1) = u / t
    end subroutine

    subroutine error(msg)
      character(*), intent(in) :: msg

      print "(A,ES25.16E3)", "n1        = ", n(1)
      print "(A,ES25.16E3)", "n2        = ", n(2)
      print "(A,ES25.16E3)", "eta1      = ", eta(1)
      print "(A,ES25.16E3)", "eta2      = ", eta(2)
      print "(A,ES25.16E3)", "dpot      = ", dpot
      print "(A,ES25.16E3)", "jj        = ", jj
      print "(2A)",          "case      = ", CASENAME(cs)

      error stop msg
    end subroutine

  end subroutine

  subroutine integrate_dist(dist, eta, k, I, dIdeta)
    !! integral_eta1^eta2 dist(eta)^k deta
    procedure(distribution)   :: dist
      !! cumulative distribution function
    real(real64), intent(in)  :: eta(2)
      !! integration bounds
    integer,      intent(in)  :: k
      !! exponent
    real(real64), intent(out) :: I
      !! output integral over dist^k
    real(real64), intent(out) :: dIdeta(2)
      !! output derivatives of I wrt eta

    real(real64) :: dum(0), dum2(0)

    if (k == 0) then
      ! I = integral_eta1^eta2 deta = eta2 - eta1
      I = eta(2) - eta(1)
      dIdeta = [-1.0, 1.0]
    else
      ! integrate using tanh-sinh
      call quad(dist_k, eta(1), eta(2), dum, I, dIdeta(1), dIdeta(2), dum2, max_levels = 8)
    end if

  contains

    subroutine dist_k(eta, p, F, dFdeta, dFdp)
      real(real64), intent(in)  :: eta
      real(real64), intent(in)  :: p(:)
      real(real64), intent(out) :: F
      real(real64), intent(out) :: dFdeta
      real(real64), intent(out) :: dFdp(:)

      call dist(eta, 0, F, dFdeta)

      if (k /= 1) then
        dFdeta = k * F**(k-1) * dFdeta
        F      = F**k
      end if
    end subroutine

  end subroutine


end module
