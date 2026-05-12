module dual_m

  use, intrinsic :: iso_fortran_env, only: real64, real128

  use util_m, only: log1p, expm1

  implicit none

#:set min_dim = 3
#:set max_dim = 3
#:set real_kinds = ["real64", "real128"]

#:for DIM in range(min_dim, max_dim + 1)
  public dual_${DIM}$
#:endfor

  public operator(+), operator(-), operator(*), operator(/), operator(**)
  public abs, sqrt, exp, expm1, log, log1p, sin, sinh, cos, cosh, tan, tanh
  public sum, dot_product

#:for DIM in range(min_dim, max_dim + 1)
  type :: dual_${DIM}$(kind)
    !! dual number type for automatic differentiation in forward mode
    !! with ${DIM}$ infinitesimal components

    integer, kind :: kind
      !! working precision

    real(kind) :: x
      !! principal part
    real(kind) :: dx(${DIM}$)
      !! infinitesimal part
  end type

#:endfor

  interface init_dual
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: init_dual_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface operator(+)
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_add_dual_${DIM}$_${REAL_KIND}$
    module procedure :: dual_add_real_${DIM}$_${REAL_KIND}$
    module procedure :: dual_add_int_${DIM}$_${REAL_KIND}$
    module procedure :: real_add_dual_${DIM}$_${REAL_KIND}$
    module procedure :: int_add_dual_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface operator(-)
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_neg_${DIM}$_${REAL_KIND}$
    module procedure :: dual_sub_dual_${DIM}$_${REAL_KIND}$
    module procedure :: dual_sub_real_${DIM}$_${REAL_KIND}$
    module procedure :: dual_sub_int_${DIM}$_${REAL_KIND}$
    module procedure :: real_sub_dual_${DIM}$_${REAL_KIND}$
    module procedure :: int_sub_dual_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface operator(*)
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_mul_dual_${DIM}$_${REAL_KIND}$
    module procedure :: dual_mul_real_${DIM}$_${REAL_KIND}$
    module procedure :: dual_mul_int_${DIM}$_${REAL_KIND}$
    module procedure :: real_mul_dual_${DIM}$_${REAL_KIND}$
    module procedure :: int_mul_dual_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface operator(/)
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_div_dual_${DIM}$_${REAL_KIND}$
    module procedure :: dual_div_real_${DIM}$_${REAL_KIND}$
    module procedure :: dual_div_int_${DIM}$_${REAL_KIND}$
    module procedure :: real_div_dual_${DIM}$_${REAL_KIND}$
    module procedure :: int_div_dual_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface operator(**)
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_pow_dual_${DIM}$_${REAL_KIND}$
    module procedure :: dual_pow_real_${DIM}$_${REAL_KIND}$
    module procedure :: dual_pow_int_${DIM}$_${REAL_KIND}$
    module procedure :: real_pow_dual_${DIM}$_${REAL_KIND}$
    module procedure :: int_pow_dual_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface abs
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_abs_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface sqrt
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_sqrt_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface exp
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_exp_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface expm1
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_expm1_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface log
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_log_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface log1p
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_log1p_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface sin
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_sin_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface sinh
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_sinh_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface cos
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_cos_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface cosh
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_cosh_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface tan
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_tan_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface tanh
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_tanh_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface sum
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_sum_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

  interface dot_product
#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
    module procedure :: dual_dot_dual_${DIM}$_${REAL_KIND}$
    module procedure :: dual_dot_real_${DIM}$_${REAL_KIND}$
    module procedure :: real_dot_dual_${DIM}$_${REAL_KIND}$
#:endfor
#:endfor
  end interface

contains

#:for DIM in range(min_dim, max_dim + 1)
#:for REAL_KIND in real_kinds
  subroutine init_dual_${DIM}$_${REAL_KIND}$(this, x, dx)
    !! initialize dual number with ${DIM}$ components (${REAL_KIND}$ version)
    class(dual_${DIM}$(${REAL_KIND}$)), intent(out) :: this
    real(kind=${REAL_KIND}$),           intent(in)  :: x
      !! principal value
    real(kind=${REAL_KIND}$), optional, intent(in)  :: dx(${DIM}$)
      !! mark this as i-th dependency

    this%x = x
    if (present(dx)) then
      this%dx = dx
    else
      this%dx = 0
    end if
  end subroutine

  elemental function dual_add_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! add two dual numbers with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x  + y%x
    r%dx = x%dx + y%dx
  end function

  elemental function dual_add_real_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! add real to dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    real(kind=${REAL_KIND}$),          intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x  + y
    r%dx = x%dx
  end function

  elemental function dual_add_int_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! add integer to dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x  + y
    r%dx = x%dx
  end function

  elemental function real_add_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! add dual to real number with ${DIM}$ components (${REAL_KIND}$)
    real(kind=${REAL_KIND}$),          intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x + y%x
    r%dx =     y%dx
  end function

  elemental function int_add_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! add dual to integer with ${DIM}$ components (${REAL_KIND}$)
    integer,                           intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x + y%x
    r%dx =     y%dx
  end function

  elemental function dual_neg_${DIM}$_${REAL_KIND}$(x) result(r)
    !! negate dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = - x%x
    r%dx = - x%dx
  end function

  elemental function dual_sub_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! subtract two dual numbers with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x  - y%x
    r%dx = x%dx - y%dx
  end function

  elemental function dual_sub_real_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! subtract real from dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    real(kind=${REAL_KIND}$),          intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x  - y
    r%dx = x%dx
  end function

  elemental function dual_sub_int_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! subtract integer from dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x  - y
    r%dx = x%dx
  end function

  elemental function real_sub_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! subtract dual from real number with ${DIM}$ components (${REAL_KIND}$)
    real(kind=${REAL_KIND}$),          intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x - y%x
    r%dx =   - y%dx
  end function

  elemental function int_sub_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! subtract dual from integer with ${DIM}$ components (${REAL_KIND}$)
    integer,                           intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x - y%x
    r%dx =   - y%dx
  end function

  elemental function dual_mul_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! multiply two dual numbers with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x * y%x
    r%dx = x%dx * y%x + x%x * y%dx
  end function

  elemental function dual_mul_real_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! multiply dual with real number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    real(kind=${REAL_KIND}$),          intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x  * y
    r%dx = x%dx * y
  end function

  elemental function dual_mul_int_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! multiply dual with integer with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x  * y
    r%dx = x%dx * y
  end function

  elemental function real_mul_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! multiply real with dual number with ${DIM}$ components (${REAL_KIND}$)
    real(kind=${REAL_KIND}$),          intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x * y%x
    r%dx = x * y%dx
  end function

  elemental function int_mul_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! multiply integer with dual number with ${DIM}$ components (${REAL_KIND}$)
    integer,                           intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x * y%x
    r%dx = x * y%dx
  end function

  elemental function dual_div_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! divide two dual numbers with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x / y%x
    r%dx = (x%dx - (x%x / y%x) * y%dx) / y%x
  end function

  elemental function dual_div_real_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! divide dual by real number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    real(kind=${REAL_KIND}$),          intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x  / y
    r%dx = x%dx / y
  end function

  elemental function dual_div_int_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! divide dual by integer with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x  / y
    r%dx = x%dx / y
  end function

  elemental function real_div_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! divide real by dual number with ${DIM}$ components (${REAL_KIND}$)
    real(kind=${REAL_KIND}$),          intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x / y%x
    r%dx = - (r%x * y%dx) / y%x
  end function

  elemental function int_div_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! divide integer by dual number with ${DIM}$ components (${REAL_KIND}$)
    integer,                           intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x / y%x
    r%dx = - (r%x * y%dx) / y%x
  end function

  elemental function dual_pow_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! raise a dual to the power of a dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x**y%x
    r%dx = y%x * x%x ** (y%x-1) * x%dx

    if ((x%x /= 0) .or. (y%x <= 0)) then
      r%dx = r%dx + r%x * log(x%x) * y%dx
    end if
  end function

  elemental function dual_pow_real_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! raise a dual to the power of a real number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    real(kind=${REAL_KIND}$),          intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x%x**y
    r%dx = y * x%x**(y-1) * x%dx
  end function

  elemental function dual_pow_int_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! raise a dual number to the power of an integer with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    integer                            :: p
    type(dual_${DIM}$(${REAL_KIND}$)) :: r2

    ! for negative integers: compute inverse first
    if (y < 0) then
      r = 1.0_${REAL_KIND}$ / x
      p = - y
    elseif (y > 0) then
      r = x
      p = y
    else
      r%x  = 1.0_${REAL_KIND}$
      r%dx = 0.0_${REAL_KIND}$
      return
    end if

    ! compute low integer powers by repeated multiplication; high powers by real power function
    select case (p)
      case (2)
        r = r * r
      case (3)
        r = (r * r) * r
      case (4)
        r2 = r * r
        r  = r2 * r2
      case default
        r = r ** real(p, kind=${REAL_KIND}$)
    end select
  end function

  elemental function real_pow_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! raise real to the power of a dual number with ${DIM}$ components (${REAL_KIND}$)
    real(kind=${REAL_KIND}$),          intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = x**y%x
    if ((x == 0) .and. (y%x > 0)) then
      r%dx = 0 ! avoid log(0)
    else
      r%dx = r%x * log(x) * y%dx
    end if
  end function

  elemental function int_pow_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! raise integer to the power of a dual number with ${DIM}$ components (${REAL_KIND}$)
    integer,                           intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r = real(x, kind=${REAL_KIND}$) ** y
  end function

  elemental function dual_abs_${DIM}$_${REAL_KIND}$(x) result(r)
    !! absolute value of dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    if (x%x >= 0) then
      r =  x
    else
      r = -x
    end if
  end function

  elemental function dual_sqrt_${DIM}$_${REAL_KIND}$(x) result(r)
    !! square root of dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = sqrt(x%x)
    r%dx = x%dx / (2 * r%x)
  end function

  elemental function dual_exp_${DIM}$_${REAL_KIND}$(x) result(r)
    !! exponential of dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = exp(x%x)
    r%dx = r%x * x%dx
  end function

  elemental function dual_expm1_${DIM}$_${REAL_KIND}$(x) result(r)
    !! exp(x) - 1 for dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = expm1(x%x)
    r%dx = exp(x%x) * x%dx
  end function

  elemental function dual_log_${DIM}$_${REAL_KIND}$(x) result(r)
    !! natural logarithm of dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = log(x%x)
    r%dx = x%dx / x%x
  end function

  elemental function dual_log1p_${DIM}$_${REAL_KIND}$(x) result(r)
    !! log(1 + x) for dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = log1p(x%x)
    r%dx = x%dx / (1.0_${REAL_KIND}$ + x%x)
  end function

  elemental function dual_sin_${DIM}$_${REAL_KIND}$(x) result(r)
    !! sine of dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = sin(x%x)
    r%dx = cos(x%x) * x%dx
  end function

  elemental function dual_sinh_${DIM}$_${REAL_KIND}$(x) result(r)
    !! sine of dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = sinh(x%x)
    r%dx = cosh(x%x) * x%dx
  end function

  elemental function dual_cos_${DIM}$_${REAL_KIND}$(x) result(r)
    !! cosine of dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = cos(x%x)
    r%dx = - sin(x%x) * x%dx
  end function

  elemental function dual_cosh_${DIM}$_${REAL_KIND}$(x) result(r)
    !! cosine of dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = cosh(x%x)
    r%dx = sinh(x%x) * x%dx
  end function

  elemental function dual_tan_${DIM}$_${REAL_KIND}$(x) result(r)
    !! tangent of dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = tan(x%x)
    r%dx = x%dx / cos(x%x)**2
  end function

  elemental function dual_tanh_${DIM}$_${REAL_KIND}$(x) result(r)
    !! tangent of dual number with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r%x  = tanh(x%x)
    r%dx = x%dx / cosh(x%x)**2
  end function

  function dual_sum_${DIM}$_${REAL_KIND}$(x) result(r)
    !! sum of dual number array with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x(:)
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    integer :: i

    call init_dual(r, 0.0_${REAL_KIND}$)
    do i = 1, size(x)
      r = r + x(i)
    end do
  end function

  function dual_dot_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! dot product of two dual vectors with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x(:)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y(:)
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    integer :: i

    call init_dual(r, 0.0_${REAL_KIND}$)
    do i = 1, size(x)
      r = r + x(i) * y(i)
    end do
  end function

  function real_dot_dual_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! dot product of real and dual vectors with ${DIM}$ components (${REAL_KIND}$)
    real(kind=${REAL_KIND}$),          intent(in) :: x(:)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: y(:)
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    integer :: i

    call init_dual(r, 0.0_${REAL_KIND}$)
    do i = 1, size(x)
      r = r + x(i) * y(i)
    end do
  end function

  function dual_dot_real_${DIM}$_${REAL_KIND}$(x, y) result(r)
    !! dot product of dual and real vectors with ${DIM}$ components (${REAL_KIND}$)
    type(dual_${DIM}$(${REAL_KIND}$)), intent(in) :: x(:)
    real(kind=${REAL_KIND}$),          intent(in) :: y(:)
    type(dual_${DIM}$(${REAL_KIND}$))             :: r

    r = dot_product(y, x)
  end function

#:endfor
#:endfor

end module
