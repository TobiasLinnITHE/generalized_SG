! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

module generalized_SG_m

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_next_after
  use, intrinsic :: iso_c_binding,   only: c_double, c_int
  use, intrinsic :: iso_fortran_env, only: real64, real128

  use dual_m
  use quad_m, only: quad
  use util_m, only: ber, expm1, int2str, log1p, LOGGING_SILENT, LOGGING_ERROR, LOGGING_INFO

  implicit none

  private
  public :: distribution, inv_distribution, int_distribution
  public :: STATUS_SUCCESS, STATUS_INVALID_INPUT, STATUS_ITERATION_FAILED, STATUS_UNKNOWN_ERROR
  public :: get_edge_current

  interface

    function distribution(eta, k) result(F) bind(c)
      !! k-th derivative of cumulative distribution function

      import :: c_double, c_int

      real(c_double), value, intent(in) :: eta
        !! normalized chemical potential
      integer(c_int), value, intent(in) :: k
        !! derivative index with k in {0, 1, 2, 3, 4}
      real(c_double)                    :: F
        !! return k-th derivative of distribution function
    end function

    function inv_distribution(F) result(eta) bind(c)
      !! inverse of cumulative distribution function

      import :: c_double

      real(c_double), value, intent(in)  :: F
        !! value of cumulative distribution function
      real(c_double)                     :: eta
        !! output normalized chemical potential corresponding to F
    end function

    function int_distribution(eta1, eta2, k) result(I) bind(c)
      !! integrate cumulative distribution function

      import :: c_double, c_int

      real(c_double), value, intent(in) :: eta1
        !! lower bound
      real(c_double), value, intent(in) :: eta2
        !! upper bound
      integer(c_int), value, intent(in) :: k
        !! integrand exponent
      real(c_double)                    :: I
        !! return value of integral
    end function

  end interface

  enum, bind(c)
    enumerator :: STATUS_SUCCESS          = 0
      !! successful calculation
    enumerator :: STATUS_INVALID_INPUT    = 1
      !! invalid input parameters (e.g. non-positive density)
    enumerator :: STATUS_ITERATION_FAILED = 2
      !! iteration failed to converge
    enumerator :: STATUS_UNKNOWN_ERROR    = 3
      !! unknown error occured, result is not finite
  end enum

  enum, bind(c)
    !! cases encode relative ordering and sizes of potential and chemical potential drops

    enumerator :: CASE0A  = 1
      !! dpot > 0, dpot < deta
    enumerator :: CASE0B  = 2
      !! dpot > 0, 0 < deta < dpot
    enumerator :: CASE0C  = 3
      !! dpot > 0, deta < 0
    enumerator :: CASE1A  = 4
      !! dpot > 0, deta ~ dpot
    enumerator :: CASE1B  = 5
      !! dpot > 0, deta ~ 0
    enumerator :: CASE2A  = 6
      !! dpot ~ 0, deta > 0
    enumerator :: CASE2B  = 7
      !! dpot ~ 0, deta ~ 0
  end enum

  character(2), parameter :: CASENAME(7) = ["0a", "0b", "0c", "1a", "1b", "2a", "2b"]
    !! case names for logging purposes
  real(real64), parameter :: CASETOL = 1e-3_real64
    !! tolerance for approximates in case selection, chosen in conjunction with taylor series accuracy

  real(real64), parameter :: RTOL = 2e-14_real64
    !! relative error tolerance for Newton iteration in case 0
  real(real64), parameter :: ATOL = 1e-24_real64
    !! absolute error tolerance for Newton iteration in case 0
  integer,      parameter :: MIN_ITERATIONS = 2
    !! minimum number of Newton iterations in case 0
  integer,      parameter :: MAX_ITERATIONS = 20
    !! maximum number of Newton iterations in case 0
  integer,      parameter :: MIN_QUAD_LEVELS = 4
    !! maximum number of levels in tanh-sinh quadrature (8 gives good performance and is accurate enough)
  integer,      parameter :: MAX_QUAD_LEVELS = 8
    !! maximum number of levels in tanh-sinh quadrature (8 gives good performance and is accurate enough)
  real(real64), parameter :: DETA_TOL = 1e-6_real64
    !! tolerance for chemical potential drop for diffusion enhancement factor calculation in case 0
  real(real64), parameter :: T_LIMIT = log(1e300_real64)
      !! limit for t to avoid overflow of exp(t), maximal double precision value is ~ 1.79e308
  real(real64), parameter :: EPS  = epsilon(1.0_real64)
    !! machine precision used to determine cutoff for analytical integration
  real(real64), parameter :: CUTOFF_SAFETY = 0.5
    !! cutoff safety factor for analytical integration (error estimate using second derivative is not exact)

  type edge_state
    !! internally used to store commonly used edge quantities

    real(real64) :: F(2)
      !! normalized densities on edge endpoints (F(eta_L), F(eta_R))
    real(real64) :: dpot
      !! potential drop along edge
    real(real64) :: eta(2), deta_dF(2)
      !! chemical potential on edge endpoints
    real(real64) :: deta
      !! chemical potential drop along edge
    real(real64) :: deta_sgn
      !! sign of deta
    integer      :: case
      !! relative ordering of dpot, deta
  contains
    procedure :: init       => edge_state_init
    procedure :: print_info => edge_state_print_info
  end type

  type analytical_integration
    !! internally used to store parameters for analytical part of the residual calculation

    real(real64) :: cutoff, dcutoff_detaL
      !! cut-off where distribution function can be approximated by linear function around etaL
    real(real64) :: scaling, dscaling_detaL, dscaling_ddpot
      !! scaling factor: 1 / (dpot**2 * F'(etaL))
    real(real64) :: shift, dshift_detaL
      !! shift inside logarithm: cutoff * F'(etaL) / F(etaL)
    real(real64) :: const, dconst_detaL, dconst_ddpot
      !! constant part of integral: sign(deta) * cutoff / dpot
  end type

contains

  subroutine get_edge_current(dist, inv_dist, int_dist, F, dpot, j, dj_dF, dj_ddpot, logging, status)
    !! get normalized edge current by solving integral equation
    procedure(distribution)     :: dist
      !! cumulative distribution function
    procedure(inv_distribution) :: inv_dist
      !! inverse of cumulative distribution function
    procedure(int_distribution) :: int_dist
      !! integral over cumulative distribution function
    real(real64),   intent(in)  :: F(2)
      !! normalized densities at left/right edge end-point (= F(etaL), F(etaR))
    real(real64),   intent(in)  :: dpot
      !! normalized potential drop
    real(real64),   intent(out) :: j
      !! output normalized current density on edge
    real(real64),   intent(out) :: dj_dF(2)
      !! output derivative of j wrt F
    real(real64),   intent(out) :: dj_ddpot
      !! output derivative of j wrt dpot
    integer,        intent(in)  :: logging
      !! logging level
    integer,        intent(out) :: status
      !! output status code

    logical          :: flip_edge
    type(edge_state) :: es

    if (any(.not. ieee_is_finite(F)) .or. any(F <= 0) .or. .not. ieee_is_finite(dpot)) then
      if (logging >= LOGGING_ERROR) then
        print "(A,ES25.16E3)", "nL   = ", F(1)
        print "(A,ES25.16E3)", "nR   = ", F(2)
        print "(A,ES25.16E3)", "dpot = ", dpot
      end if
      status = STATUS_INVALID_INPUT
      return
    end if

    ! flip edge direction if potential drop is negative
    flip_edge = (dpot < 0)

    ! initialize edge state
    if (flip_edge) then
      ! keep potential drop positive by flipping edge direction
      call es%init(dist, inv_dist, [F(2), F(1)], - dpot)
    else
      call es%init(dist, inv_dist, F, dpot)
    end if

    ! debug info
    if (logging >= LOGGING_INFO) call es%print_info()

    ! perform actual calculation of current
    call calc_current(es, dist, int_dist, j, dj_dF, dj_ddpot, logging, status)
    if (status /= STATUS_SUCCESS) return

    ! use symmetry: j(nL, nR, dpot) = - j (nR, nL, -dpot)
    if (flip_edge) then
      ! j is negated, dj/dn negated + reversed, dj/ddpot unchanged (negated twice)
      j     = - j
      dj_dF = - dj_dF(2:1:-1)
    end if

    ! perform sanity check of result
    if (.not. ieee_is_finite(j) .or. any(.not. ieee_is_finite(dj_dF)) .or. .not. ieee_is_finite(dj_ddpot)) then
      if (logging >= LOGGING_ERROR) then
        print "(A,ES25.16E3)", "j        = ", j
        print "(A,ES25.16E3)", "dj_dF    = ", dj_dF
        print "(A,ES25.16E3)", "dj_ddpot = ", dj_ddpot
      end if
      status = STATUS_UNKNOWN_ERROR
    end if
  end subroutine


  subroutine edge_state_init(this, dist, inv_dist, F, dpot)
    !! initialize edge state
    class(edge_state), intent(out) :: this
    procedure(distribution)        :: dist
    procedure(inv_distribution)    :: inv_dist
      !! inverse distribution function
    real(real64),      intent(in)  :: F(2)
      !! normalized densities at left/right edge end-point (same as F)
    real(real64),      intent(in)  :: dpot
      !! normalized potential drop

    this%F    = F
    this%dpot = dpot

    ! get chemical potential at edge endpoints
    this%eta(1) = inv_dist(this%F(1))
    this%eta(2) = inv_dist(this%F(2))
    this%deta_dF(1) = 1 / dist(this%eta(1), 1)
    this%deta_dF(2) = 1 / dist(this%eta(2), 1)

    ! chemical potential drop
    this%deta = this%eta(2) - this%eta(1)
    this%deta_sgn = sign(1.0_real64, this%deta)

    ! get the case depending on whether dpot, deta or both are small or close to each other or not
    if (abs(this%dpot) > CASETOL) then
      if (this%deta - this%dpot > this%dpot * CASETOL) then
        this%case = CASE0A
      elseif (this%deta - this%dpot >= - this%dpot * CASETOL) then
        this%case = CASE1A
      elseif (this%deta > CASETOL) then
        this%case = CASE0B
      elseif (this%deta >= - CASETOL) then
        this%case = CASE1B
      else
        this%case = CASE0C
      end if
    else
      if (abs(this%deta) > CASETOL) then
        this%case = CASE2A
      else
        this%case = CASE2B
      end if
    end if
  end subroutine

  subroutine edge_state_print_info(this)
    class(edge_state), intent(in) :: this

    print *
    print "(A,ES24.16E3)", "FL   = ", this%F(1)
    print "(A,ES24.16E3)", "FR   = ", this%F(2)
    print "(A,ES24.16E3)", "etaL = ", this%eta(1)
    print "(A,ES24.16E3)", "etaR = ", this%eta(2)
    print "(A,ES24.16E3)", "dpot = ", this%dpot
    print "(A,ES24.16E3)", "deta = ", this%deta
    print "(2A)",          "case = ", CASENAME(this%case)
  end subroutine


  subroutine calc_current(es, dist, int_dist, j, dj_dF, dj_ddpot, logging, status)
    !! calculate current
    type(edge_state), intent(in)  :: es
    procedure(distribution)       :: dist
      !! cumulative distribution function
    procedure(int_distribution)   :: int_dist
      !! integral over cumulative distribution function
    real(real64),     intent(out) :: j
      !! output edge current
    real(real64),     intent(out) :: dj_dF(2)
      !! output derivatives of j wrt n
    real(real64),     intent(out) :: dj_ddpot
      !! output derivative of j wrt dpot
    integer,          intent(in)  :: logging
      !! logging level
    integer,          intent(out) :: status
      !! output status code

    real(real64) :: dj_deta(2)

    status = STATUS_SUCCESS

    ! calculate current based on case
    select case (es%case)
    case (CASE0A, CASE0B, CASE0C)
      ! general case with newton iteration
      call case_0(es, dist, int_dist, j, dj_deta, dj_ddpot, logging, status)

    case (CASE1A)
      ! dpot > 0, deta ~ dpot
      call case_1a(es, int_dist, j, dj_deta, dj_ddpot)

    case (CASE1B)
      ! dpot > 0, deta ~ 0
      call case_1b(es, dist, j, dj_deta, dj_ddpot)

    case (CASE2A)
      ! dpot ~ 0, deta > 0
      call case_2a(es, int_dist, j, dj_deta, dj_ddpot)

    case (CASE2B)
      ! dpot ~ 0, deta ~ 0
      call case_2b(es, dist, j, dj_deta, dj_ddpot)

    end select

    ! derivatives wrt densities
    dj_dF = dj_deta * es%deta_dF
  end subroutine

  subroutine case_0(es, dist, int_dist, j, dj_deta, dj_ddpot, logging, status)
    type(edge_state), intent(in)  :: es
    procedure(distribution)       :: dist
      !! cumulative distribution function
    procedure(int_distribution)   :: int_dist
      !! integral over cumulative distribution function
    real(real64),     intent(out) :: j
      !! output edge current
    real(real64),     intent(out) :: dj_deta(2)
      !! output derivatives of j wrt eta
    real(real64),     intent(out) :: dj_ddpot
      !! output derivative of j wrt dpot
    integer,          intent(in)  :: logging
      !! logging level
    integer,          intent(out) :: status
      !! output status code

    integer                      :: iteration
    real(real64)                 :: j_min, j_max, dj_dt, t_min, t_max, t, t_old, t_update, err
    real(real64)                 :: res, dres_dt, dres_deta(2), dres_ddpot
    type(analytical_integration) :: an

    ! get initial range for j
    call get_j_range(es, int_dist, j_min, j_max)

    ! initial guess for j: enhanced diffusion approximation
    call get_j_approx(es, j)

    ! transform initial range and approximation with: j = dpot * n(1) * (1 - deta_sgn * exp(-t))
    call transform_j_to_t(es, j_min, t_min)
    call transform_j_to_t(es, j_max, t_max)
    if (t_min > t_max) then
      t = t_min
      t_min = t_max
      t_max = t
    end if
    call transform_j_to_t(es, j, t)

    ! adjust initial guess if not in range
    if ((t < t_min) .or. (t > t_max)) t = 0.5 * (t_min + t_max)

    ! initialize params for analytical part of integration
    call init_analytical_integration(es, dist, an)

    if (logging >= LOGGING_INFO) then
      print "(A,ES24.16E3)", "t_min = ", t_min
      print "(A,ES24.16E3)", "t_max = ", t_max
      print "(A,ES24.16E3)", "t     = ", t
    end if

    ! Newton iteration
    err = huge(1.0)
    iteration = 0
    if (logging >= LOGGING_INFO) then
      print "(A)", "Newton iteration:"
    end if
    do while (((iteration < MIN_ITERATIONS) .or. (err > RTOL * abs(j) + ATOL)) .and. (t_min < t_max))
      iteration = iteration + 1
      if (iteration > MAX_ITERATIONS) then
        if (logging >= LOGGING_ERROR) then
          call es%print_info()
          print "(A)", "No convergence after " // int2str(MAX_ITERATIONS) // " iterations"
        end if
        status = STATUS_ITERATION_FAILED
        return
      end if

      ! evaluate residual and get Newton update
      call residual(es, an, dist, t, res, dres_dt, dres_deta, dres_ddpot, j, dj_dt, dj_deta, dj_ddpot)

      ! Newton update
      t_update = - res / dres_dt
      err      = abs(t_update)

      ! update bounds (assume monotonic behaviour)
      if (t_update > 0) then
        t_min = min(t, t_max)
      else
        t_max = max(t, t_min)
      end if

      ! update solution
      t_old = t
      t    = t + t_update

      ! bisection
      if ((t < t_min) .or. (t > t_max) .or. ((t_old == t_min) .and. (t == t_max))) then
        if (logging >= LOGGING_INFO) then
          print "(A)", "  performing bisection"
          print "(A,ES24.16E3)", "    t_min = ", t_min
          print "(A,ES24.16E3)", "    t_max = ", t_max
          print "(A,ES24.16E3)", "    t_old = ", t_old
          print "(A,ES24.16E3)", "    t     = ", t
          print "(A,ES24.16E3)", "    t_new = ", 0.5 * (t_min + t_max)
        end if
        t   = 0.5 * (t_min + t_max)
        err = min(err, t_max - t_min)
      end if

      ! update current error
      err = abs(dj_dt) * err

      if (logging >= LOGGING_INFO) then
        print "(I6,A,ES24.16E3,A,ES24.16E3,A,ES24.16E3)", iteration, ": t = ", t, ", j = ", j, ", err/j = ", err / abs(j)
      end if
    end do

    ! get current + derivatives with implicit differentiation
    call residual(es, an, dist, t, res, dres_dt, dres_deta, dres_ddpot, j, dj_dt, dj_deta, dj_ddpot)
    dj_deta  = dj_deta  - dres_deta  * dj_dt / dres_dt
    dj_ddpot = dj_ddpot - dres_ddpot * dj_dt / dres_dt
  end subroutine

  subroutine get_j_range(es, int_dist, j_min, j_max)
    !! get the valid range for j for the general case (CASE0A, CASE0B, CASE0C)
    type(edge_state), intent(in)  :: es
    procedure(int_distribution)   :: int_dist
      !! integral over cumulative distribution function
    real(real64),     intent(out) :: j_min
      !! output lower bound of j
    real(real64),     intent(out) :: j_max
      !! output upper bound of j

    real(real64) :: dum(2)

    select case (es%case)
    case (CASE0A)
      j_min = int_dist(es%eta(2) - es%dpot, es%eta(2), -1)
      j_min = es%dpot * (es%dpot - es%deta) / j_min
      j_max = es%F(1) * (es%dpot - es%deta)

    case (CASE0B)
      j_min = es%F(1) * (es%dpot - es%deta)
      j_max = min(es%F(1) * es%dpot, es%F(2) * (es%dpot - es%deta))

    case (CASE0C)
      j_min = max(es%F(1) * es%dpot, es%F(2) * (es%dpot - es%deta))
      j_max = es%F(1) * (es%dpot - es%deta)

    end select
  end subroutine

  subroutine get_j_approx(es, j)
    !! get enhanced diffusion approximation
    type(edge_state), intent(in)  :: es
    real(real64),     intent(out) :: j
      !! output approximation for current

    real(real64) :: g

    ! get diffusion enhancement factor
    if (abs(es%deta) < DETA_TOL) then
      ! avoid NaN from 0 / 0 => lim_{deta->0} deta / log(F_R / F_L) = F_L * deta_dF(1)
      g = es%F(1) * es%deta_dF(1)
    else
      g = es%deta / log(es%F(2) / es%F(1))
    end if

    ! enhanced diffusion scheme
    j = g * (ber(-es%dpot / g) * es%F(1) - ber(es%dpot / g) * es%F(2))
  end subroutine

  subroutine transform_j_to_t(es, j, t)
    !! transform edge current to solution variable t
    type(edge_state), intent(in)  :: es
    real(real64),     intent(in)  :: j
      !! edge current
    real(real64),     intent(out) :: t
      !! output transformed variable

    if (es%deta_sgn > 0) then
      if (j < es%dpot * es%F(1)) then
        t = - log1p(- j / (es%dpot * es%F(1)))
      else
        t = T_LIMIT
      end if
    else
      if (j > es%dpot * es%F(1)) then
        t = - log(j / (es%dpot * es%F(1)) - 1)
      else
        t = T_LIMIT
      end if
    end if
  end subroutine

  subroutine transform_t_to_j(es, t, j, dj_deta, dj_ddpot, dj_dt)
    !! transform solution variable t to edge current j
    type(edge_state), intent(in)  :: es
      !! common parameters
    real(real64),     intent(in)  :: t
      !! transformed solution variable
    real(real64),     intent(out) :: j
      !! output edge current
    real(real64),     intent(out) :: dj_deta(2)
      !! output derivatives of j wrt eta
    real(real64),     intent(out) :: dj_ddpot
      !! output derivative of j wrt dpot
    real(real64),     intent(out) :: dj_dt
      !! output derivative of j wrt t

    real(real64) :: e, e1

    ! transform: j = dpot * FL * (1 - deta_sgn * exp(-t))
    e = exp(-t)
    if (es%deta_sgn > 0) then
      e1 = - expm1(-t)
    else
      e1 = e + 1
    end if
    j          = es%dpot * es%F(1) * e1
    dj_deta(1) = es%dpot * e1 / es%deta_dF(1)
    dj_deta(2) = 0
    dj_ddpot   = es%F(1) * e1
    dj_dt      = es%dpot * es%F(1) * es%deta_sgn * e
  end subroutine

  subroutine residual(es, an, dist, t, res, dres_dt, dres_deta, dres_ddpot, j, dj_dt, dj_deta, dj_ddpot)
    !! calculate residual for newton iteration of case 0
    type(edge_state),             intent(in)  :: es
    type(analytical_integration), intent(in)  :: an
      !! parameters and coefficients for analytical part of integral
    procedure(distribution)                   :: dist
    real(real64),                 intent(in)  :: t
      !! transformed solution variables
    real(real64),                 intent(out) :: res
      !! output residual
    real(real64),                 intent(out) :: dres_dt
      !! output derivative of res wrt t
    real(real64),                 intent(out) :: dres_deta(2)
      !! output derivatives of res wrt eta
    real(real64),                 intent(out) :: dres_ddpot
      !! output derivatives of res wrt dpot
    real(real64),                 intent(out) :: j
      !! output current density corresponding to t
    real(real64),                 intent(out) :: dj_dt
      !! output derivative of j wrt t
    real(real64),                 intent(out) :: dj_deta(2)
      !! output derivatives of j wrt eta
    real(real64),                 intent(out) :: dj_ddpot
      !! output derivative of j wrt dpot

    real(real64) :: e, par(2), eta_cutoff, dres_dpar(2)

    ! transform: j = dpot * FL * (1 - deta_sgn * exp(-t))
    call transform_t_to_j(es, t, j, dj_deta, dj_ddpot, dj_dt)

    ! tanh-sinh quadrature, leave out part close to eta(1) (might have singularity)
    par = [j, es%dpot]
    eta_cutoff = es%eta(1) + es%deta_sgn * an%cutoff
    call quad(integrand, eta_cutoff, es%eta(2), par, res, dres_deta(1), dres_deta(2), dres_dpar, min_levels = MIN_QUAD_LEVELS, max_levels = MAX_QUAD_LEVELS)
    dres_dt      = dres_dpar(1) * dj_dt
    dres_deta(1) = dres_deta(1) * (1 + es%deta_sgn * an%dcutoff_detaL) + dres_dpar(1) * dj_deta(1)
    dres_ddpot   = dres_dpar(2) + dres_dpar(1) * dj_ddpot

    ! analytical integration for missing part (distribution function approximated by first order Taylor series)
    call residual_analytical(an, t, j, dj_dt, dj_deta, dj_ddpot, res, dres_dt, dres_deta, dres_ddpot)

    ! res = integral - 1 == 0
    res = res - 1

  contains

    subroutine integrand(eta, n_par, par, u, du_deta, du_dpar) bind(c)
      real(c_double), value, intent(in)  :: eta
        !! chemical potential
      integer(c_int), value, intent(in)  :: n_par
        !! number of additional parameters
      real(c_double),        intent(in)  :: par(n_par)
        !! par = [j, dpot]
      real(c_double),        intent(out) :: u
        !! output integrand u(eta)
      real(c_double),        intent(out) :: du_deta
        !! output derivative of u wrt eta
      real(c_double),        intent(out) :: du_dpar(n_par)
        !! du_dpar = [du_dj, du_ddpot]

      real(real64) :: F, dF_deta, q

      associate (j => par(1), dpot => par(2), du_dj => du_dpar(1), du_ddpot => du_dpar(2))
        ! u = F / (dpot * F - j)
        F        = dist(eta, 0)
        dF_deta  = dist(eta, 1)
        q        = dpot * F - j
        u        = F / q
        du_deta  = (1 - u * dpot) / q * dF_deta
        du_dj    = u / q
        du_ddpot = - u**2
      end associate
    end subroutine

  end subroutine



  subroutine init_analytical_integration(es, dist, an)
    !! initialize parameters for the analytical part of the residual calculation
    type(edge_state),             intent(in)  :: es
    procedure(distribution)                   :: dist
      !! cumulative distribution function
    type(analytical_integration), intent(out) :: an

    real(real64) :: FL2, FL3

    ! get cutoff such that F can be approximated by a linear function around etaL +/- cutoff
    ! use second derivative to estimate error of first order Taylor
    FL2 = dist(es%eta(1), 2)
    FL3 = dist(es%eta(1), 3)
    an%cutoff        = CUTOFF_SAFETY * sqrt(2 * EPS * es%F(1) / abs(FL2))
    an%dcutoff_detaL = (1.0 / es%deta_dF(1) - es%F(1) * FL3 / FL2) * EPS / (4 * an%cutoff * abs(FL2))

    ! scaling factor for logarithm part
    an%scaling        = es%deta_dF(1) / es%dpot**2 ! = 1 / (dpot**2 * dndeta(1))
    an%dscaling_detaL = - FL2 * (es%deta_dF(1) / es%dpot)**2
    an%dscaling_ddpot = - 2 * es%deta_dF(1) / es%dpot**3

    ! argument shift for logarithm part
    an%shift          = an%cutoff / (es%deta_dF(1) * es%F(1))
    an%dshift_detaL   = (an%dcutoff_detaL / es%deta_dF(1) + an%cutoff * (FL2 - 1 / (es%deta_dF(1)**2 * es%F(1)))) / es%F(1)

    ! constant part of integral
    an%const          = es%deta_sgn * an%cutoff / es%dpot
    an%dconst_detaL   = es%deta_sgn * an%dcutoff_detaL / es%dpot
    an%dconst_ddpot   = - es%deta_sgn * an%cutoff / es%dpot**2
  end subroutine

  subroutine residual_analytical(an, t, j, dj_dt, dj_deta, dj_ddpot, res, dres_dt, dres_deta, dres_ddpot)
    !! evaluate analytical part of residual calculation
    type(analytical_integration), intent(in)    :: an
    real(real64),                 intent(in)    :: t
      !! transformed solution variable
    real(real64),                 intent(in)    :: j
      !! normalized current density
    real(real64),                 intent(in)    :: dj_dt
      !! derivative of j wrt t
    real(real64),                 intent(in)    :: dj_deta(2)
      !! derivatives of j wrt eta
    real(real64),                 intent(in)    :: dj_ddpot
      !! derivative of j wrt dpot
    real(real64),                 intent(inout) :: res
      !! update residual
    real(real64),                 intent(inout) :: dres_dt
      !! update derivative of res wrt t
    real(real64),                 intent(inout) :: dres_deta(2)
      !! update derivatives of res wrt eta
    real(real64),                 intent(inout) :: dres_ddpot
      !! update derivative of res wrt dpot

    real(real64) :: e, logarithm, dlogarithm_dshift, dlogarithm_dt

    e = exp(-t)
    if (t <= - log(an%shift)) then
      logarithm = log1p(an%shift / e)
    else
      logarithm = t + log(e + an%shift)
    end if
    dlogarithm_dshift  = 1 / (an%shift + e)
    dlogarithm_dt      = an%shift * dlogarithm_dshift

    res          = res          + an%scaling * j * logarithm + an%const
    dres_dt      = dres_dt      + an%scaling * (dj_dt * logarithm + j * dlogarithm_dt)
    dres_deta(1) = dres_deta(1) + (an%dscaling_detaL * j + an%scaling * dj_deta(1)) * logarithm &
      &                         + an%scaling * j * dlogarithm_dshift * an%dshift_detaL + an%dconst_detaL
    dres_ddpot   = dres_ddpot   + (an%dscaling_ddpot * j + an%scaling * dj_ddpot) * logarithm + an%dconst_ddpot
  end subroutine

  subroutine case_1a(es, int_dist, j, dj_deta, dj_ddpot)
    !! dpot >> 0 and deta ~ dpot
    !! numerically difficult to evaluate, need to use quad precision
    type(edge_state), intent(in)  :: es
    procedure(int_distribution)   :: int_dist
      !! integral over cumulative distribution function
    real(real64),     intent(out) :: j
      !! output normalized current density on edge
    real(real64),     intent(out) :: dj_deta(2)
      !! output derivative of j wrt eta
    real(real64),     intent(out) :: dj_ddpot
      !! output derivative of j wrt dpot

    integer               :: k
    real(real64)          :: Itmp, dItmp(2)
    type(dual_3(real128)) :: I(-5:0), jc(5), dpot, jd, tmp, tmp2

    ! get I_k = integral_etaL^etaR F(eta)^k deta for k = -5 to 0
    do k = -5, -1
      Itmp = int_dist(es%eta(1), es%eta(2), k)
      call init_dual(I(k), real(Itmp, kind = real128), dx = real([- es%F(1)**k, es%F(2)**k, 0.0], kind = real128))
    end do
    call init_dual(I(0), real(es%deta, kind = real128), dx = [-1.0_real128, 1.0_real128, 0.0_real128])

    ! first series coefficient
    jc(1) = I(0)

    ! second series coefficient
    jc(2) = I(-1)**2 - I(-2)*I(0)

    ! third series coefficient
    jc(3) = 2*I(-2)**2*I(0) - I(-1)**2*I(-2) - I(-1)*I(-3)*I(0)

    ! fourth series coefficient
    jc(4) = 2*I(-1)**2*I(-2)**2 - I(-1)**3*I(-3) - 5*I(-2)**3*I(0) &
      &     - I(-1)**2*I(-4)*I(0) + 5*I(-1)*I(-2)*I(-3)*I(0)

    ! fifth series coefficient
    jc(5) = 14*I(0)*I(-2)**4 - I(-1)**4*I(-4) - 5*I(-1)**2*I(-2)**3 &
      &     - I(0)*I(-1)**3*I(-5) + 5*I(-1)**3*I(-2)*I(-3) + 3*I(0)*I(-1)**2*I(-3)**2 &
      &     - 21*I(0)*I(-1)*I(-2)**2*I(-3) + 6*I(0)*I(-1)**2*I(-2)*I(-4)

    ! scale coefficients by 1 / I(-1)**(2 * k - 1)
    tmp  = 1 / I(-1)
    tmp2 = tmp * tmp
    jc(1) = jc(1) * tmp
    do k = 2, 5
      tmp   = tmp * tmp2
      jc(k) = jc(k) * tmp
    end do

    ! evaluate series
    call init_dual(dpot, real(es%dpot, kind = real128), dx = [0.0_real128, 0.0_real128, 1.0_real128])
    tmp  = dpot - I(0)
    tmp2 = tmp
    jd = jc(1) * tmp
    do k = 2, 5
      ! skip coefficient if not finite (might lead to precision loss)
      if (.not. ieee_is_finite(jc(k)%x)) cycle

      ! add next term, set tmp = (dpot - I(0))**k
      tmp = tmp * tmp2
      jd  = jd + jc(k) * tmp
    end do

    ! convert to double precision and extract derivatives
    j        = real(jd%x,       kind = real64)
    dj_deta  = real(jd%dx(1:2), kind = real64)
    dj_ddpot = real(jd%dx(3),   kind = real64)
  end subroutine

  subroutine case_1b(es, dist, j, dj_deta, dj_ddpot)
    !! |dpot| >> 0 and deta ~ 0
    type(edge_state), intent(in)  :: es
    procedure(distribution)       :: dist
      !! cumulative distribution function
    real(real64),     intent(out) :: j
      !! output normalized current density on edge
    real(real64),     intent(out) :: dj_deta(2)
      !! output derivative of j wrt eta
    real(real64),     intent(out) :: dj_ddpot
      !! output derivative of j wrt dpot

    integer              :: k
    real(real64)         :: Ftmp(0:4)
    type(dual_3(real64)) :: eta3, dpot, deta, e, F(0:3), G(2), jc(0:3), jc3(3), jd
    type(dual_3(real64)) :: e2, u0, u1, u2, u3, u4, u5, u6, u7, G12, G13

    ! express mean chemical potential and potential drop as dual numbers
    call init_dual(eta3, 0.5 * (es%eta(1) + es%eta(2)), dx = [0.5, 0.5, 0.0])
    call init_dual(dpot, es%dpot, dx = [ 0.0, 0.0, 1.0])
    call init_dual(deta, es%deta, dx = [-1.0, 1.0, 0.0])

    ! evaluate distribution and derivatives at mean chemical potential
    do k = 0, 4
      Ftmp(k) = dist(eta3%x, k)
    end do
    do k = 0, 3
      call init_dual(F(k), Ftmp(k), dx = Ftmp(k+1) * [0.5, 0.5, 0.0])
    end do

    ! scale F1 and F2 by F0
    G = F(1:2) / F(0)

    ! zero-th coefficient
    jc(0) = dpot * F(0)

    ! first coefficient
    e     = 1 / tanh(0.5 * G(1) * dpot)
    jc(1) = - 0.5 * F(1) * dpot * e

    ! precompute common subexpressions for jc(2) and jc(3)
    e2  = e * e
    u0  = 1 - e2
    u1  = e * u0
    u2  =  e * (2 - 3 * e2)
    u3  = u0 * (1 - 5 * e2)
    u4  = u0 * (1 - 4 * e2)
    u5  = u0 * (1.5_real64 - 7 * e2)
    u6  = e2 * u0
    u7  = u1 * (1 - 2 * e2)
    G12 = G(1) * G(1)
    G13 = G12 * G(1)

    ! second coefficient
    jc(2) = 0.25 * dpot * (0.5 * F(2) + (F(2) - F(1) * G(1)) * u0 * (0.5 * dpot * G(1) * e - 1))

    ! third coefficient evaluated using horner scheme
    jc3(1) = 4 * u1 * (G12 * F(1) + 0.25 * (F(2) / F(1)) * F(2) - 1.5 * G(1) * F(2)) + u2 * F(3) / 3
    jc3(2) = 0.5 * u3 * G(2) * F(2) + u4 * G13 * F(1) - u5 * G12 * F(2) - 0.5 * u6 * G(1) * F(3)
    jc3(3) = u7 * G(1) * (G12 * F(2) - 0.5 * G(2) * F(2) - 0.5 * G13 * F(1))
    jc(3)  = dpot * (jc3(1) + dpot * (jc3(2) + dpot * jc3(3))) / 16

    ! evaluate series
    jd = jc(0)
    do k = 1, 3
      jd = jd + jc(k) * deta**k
    end do

    ! extract value and derivatives
    j        = jd%x
    dj_deta  = jd%dx(1:2)
    dj_ddpot = jd%dx(3)
  end subroutine

  subroutine case_2a(es, int_dist, j, dj_deta, dj_ddpot)
    !! dpot ~ 0 and |deta| >> 0
    type(edge_state), intent(in)  :: es
    procedure(int_distribution)   :: int_dist
      !! integral over cumulative distribution function
    real(real64),     intent(out) :: j
      !! output normalized current density on edge
    real(real64),     intent(out) :: dj_deta(2)
      !! output derivative of j wrt eta
    real(real64),     intent(out) :: dj_ddpot
      !! output derivative of j wrt dpot

    integer              :: k
    real(real64)         :: Itmp
    type(dual_3(real64)) :: dpot, I(1:5), jc(0:4), jd

    call init_dual(dpot, es%dpot, dx = [0.0, 0.0, 1.0])

    ! get I_k = integral_etaL^etaR F(eta)^k deta for k = 0 to 5
    do k = 1, 5
      Itmp = int_dist(es%eta(1), es%eta(2), k)
      call init_dual(I(k), Itmp, dx = [-es%F(1)**k, es%F(2)**k, 0.0_real64])
    end do

    ! normalize I_{2..5} wrt I_1 to improve numerical stability (avoid underflows etc.)
    I(2:5) = I(2:5) / I(1)

    ! coefficients
    jc(0) = - I(1)
    jc(1) = I(2)
    jc(2) = (I(2)**2 - I(3)) / I(1)
    jc(3) = (I(4) + 2*I(2)**3 - 3*I(2)*I(3)) / I(1)**2
    jc(4) = (5*I(2)**4 - I(5) + 2*I(3)**2 - 10*I(2)**2*I(3) + 4*I(2)*I(4)) / I(1)**3

    ! evaluate series
    jd = jc(0)
    do k = 1, 4
      jd = jd + jc(k) * dpot**k
    end do

    ! extract value and derivatives
    j        = jd%x
    dj_deta  = jd%dx(1:2)
    dj_ddpot = jd%dx(3)
  end subroutine

  subroutine case_2b(es, dist, j, dj_deta, dj_ddpot)
    !! dpot ~ 0 and deta ~ 0
    type(edge_state), intent(in)  :: es
    procedure(distribution)        :: dist
      !! cumulative distribution function
    real(real64),     intent(out) :: j
      !! output normalized current density on edge
    real(real64),     intent(out) :: dj_deta(2)
      !! output derivative of j wrt eta
    real(real64),     intent(out) :: dj_ddpot
      !! output derivative of j wrt dpot

    integer              :: k
    real(real64)         :: Ftmp(0:3)
    type(dual_3(real64)) :: eta3, dpot, deta, F(0:2), jd

    ! express mean chemical potential and potential drop as dual numbers
    call init_dual(eta3, 0.5 * (es%eta(1) + es%eta(2)), dx = [0.5, 0.5, 0.0])
    call init_dual(dpot, es%dpot, dx = [ 0.0, 0.0, 1.0])
    call init_dual(deta, es%deta, dx = [-1.0, 1.0, 0.0])

    ! evaluate distribution and derivatives at mean chemical potential
    do k = 0, 3
      Ftmp(k) = dist(eta3%x, k)
    end do
    do k = 0, 2
      call init_dual(F(k), Ftmp(k), dx = Ftmp(k+1) * [0.5, 0.5, 0.0])
    end do

    ! evaluate current
    jd = (dpot - deta) * (F(0) + deta * (0.5 * F(2) * deta - (F(1) / F(0)) * F(1) * dpot) / 12)

    ! extract value and derivatives
    j        = jd%x
    dj_deta  = jd%dx(1:2)
    dj_ddpot = jd%dx(3)
  end subroutine

end module
