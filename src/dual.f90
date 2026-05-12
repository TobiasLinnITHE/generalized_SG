module dual_m

  use, intrinsic :: iso_fortran_env, only: real64, real128

  use util_m, only: log1p, expm1

  implicit none


  public dual_3

  public operator(+), operator(-), operator(*), operator(/), operator(**)
  public abs, sqrt, exp, expm1, log, log1p, sin, sinh, cos, cosh, tan, tanh
  public sum, dot_product

  type :: dual_3(kind)
    !! dual number type for automatic differentiation in forward mode
    !! with 3 infinitesimal components

    integer, kind :: kind
      !! working precision

    real(kind) :: x
      !! principal part
    real(kind) :: dx(3)
      !! infinitesimal part
  end type


  interface init_dual
    module procedure :: init_dual_3_real64
    module procedure :: init_dual_3_real128
  end interface

  interface operator(+)
    module procedure :: dual_add_dual_3_real64
    module procedure :: dual_add_real_3_real64
    module procedure :: dual_add_int_3_real64
    module procedure :: real_add_dual_3_real64
    module procedure :: int_add_dual_3_real64
    module procedure :: dual_add_dual_3_real128
    module procedure :: dual_add_real_3_real128
    module procedure :: dual_add_int_3_real128
    module procedure :: real_add_dual_3_real128
    module procedure :: int_add_dual_3_real128
  end interface

  interface operator(-)
    module procedure :: dual_neg_3_real64
    module procedure :: dual_sub_dual_3_real64
    module procedure :: dual_sub_real_3_real64
    module procedure :: dual_sub_int_3_real64
    module procedure :: real_sub_dual_3_real64
    module procedure :: int_sub_dual_3_real64
    module procedure :: dual_neg_3_real128
    module procedure :: dual_sub_dual_3_real128
    module procedure :: dual_sub_real_3_real128
    module procedure :: dual_sub_int_3_real128
    module procedure :: real_sub_dual_3_real128
    module procedure :: int_sub_dual_3_real128
  end interface

  interface operator(*)
    module procedure :: dual_mul_dual_3_real64
    module procedure :: dual_mul_real_3_real64
    module procedure :: dual_mul_int_3_real64
    module procedure :: real_mul_dual_3_real64
    module procedure :: int_mul_dual_3_real64
    module procedure :: dual_mul_dual_3_real128
    module procedure :: dual_mul_real_3_real128
    module procedure :: dual_mul_int_3_real128
    module procedure :: real_mul_dual_3_real128
    module procedure :: int_mul_dual_3_real128
  end interface

  interface operator(/)
    module procedure :: dual_div_dual_3_real64
    module procedure :: dual_div_real_3_real64
    module procedure :: dual_div_int_3_real64
    module procedure :: real_div_dual_3_real64
    module procedure :: int_div_dual_3_real64
    module procedure :: dual_div_dual_3_real128
    module procedure :: dual_div_real_3_real128
    module procedure :: dual_div_int_3_real128
    module procedure :: real_div_dual_3_real128
    module procedure :: int_div_dual_3_real128
  end interface

  interface operator(**)
    module procedure :: dual_pow_dual_3_real64
    module procedure :: dual_pow_real_3_real64
    module procedure :: dual_pow_int_3_real64
    module procedure :: real_pow_dual_3_real64
    module procedure :: int_pow_dual_3_real64
    module procedure :: dual_pow_dual_3_real128
    module procedure :: dual_pow_real_3_real128
    module procedure :: dual_pow_int_3_real128
    module procedure :: real_pow_dual_3_real128
    module procedure :: int_pow_dual_3_real128
  end interface

  interface abs
    module procedure :: dual_abs_3_real64
    module procedure :: dual_abs_3_real128
  end interface

  interface sqrt
    module procedure :: dual_sqrt_3_real64
    module procedure :: dual_sqrt_3_real128
  end interface

  interface exp
    module procedure :: dual_exp_3_real64
    module procedure :: dual_exp_3_real128
  end interface

  interface expm1
    module procedure :: dual_expm1_3_real64
    module procedure :: dual_expm1_3_real128
  end interface

  interface log
    module procedure :: dual_log_3_real64
    module procedure :: dual_log_3_real128
  end interface

  interface log1p
    module procedure :: dual_log1p_3_real64
    module procedure :: dual_log1p_3_real128
  end interface

  interface sin
    module procedure :: dual_sin_3_real64
    module procedure :: dual_sin_3_real128
  end interface

  interface sinh
    module procedure :: dual_sinh_3_real64
    module procedure :: dual_sinh_3_real128
  end interface

  interface cos
    module procedure :: dual_cos_3_real64
    module procedure :: dual_cos_3_real128
  end interface

  interface cosh
    module procedure :: dual_cosh_3_real64
    module procedure :: dual_cosh_3_real128
  end interface

  interface tan
    module procedure :: dual_tan_3_real64
    module procedure :: dual_tan_3_real128
  end interface

  interface tanh
    module procedure :: dual_tanh_3_real64
    module procedure :: dual_tanh_3_real128
  end interface

  interface sum
    module procedure :: dual_sum_3_real64
    module procedure :: dual_sum_3_real128
  end interface

  interface dot_product
    module procedure :: dual_dot_dual_3_real64
    module procedure :: dual_dot_real_3_real64
    module procedure :: real_dot_dual_3_real64
    module procedure :: dual_dot_dual_3_real128
    module procedure :: dual_dot_real_3_real128
    module procedure :: real_dot_dual_3_real128
  end interface

contains

  subroutine init_dual_3_real64(this, x, dx)
    !! initialize dual number with 3 components (real64 version)
    class(dual_3(real64)), intent(out) :: this
    real(kind=real64),           intent(in)  :: x
      !! principal value
    real(kind=real64), optional, intent(in)  :: dx(3)
      !! mark this as i-th dependency

    this%x = x
    if (present(dx)) then
      this%dx = dx
    else
      this%dx = 0
    end if
  end subroutine

  elemental function dual_add_dual_3_real64(x, y) result(r)
    !! add two dual numbers with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x  + y%x
    r%dx = x%dx + y%dx
  end function

  elemental function dual_add_real_3_real64(x, y) result(r)
    !! add real to dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    real(kind=real64),          intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x  + y
    r%dx = x%dx
  end function

  elemental function dual_add_int_3_real64(x, y) result(r)
    !! add integer to dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x  + y
    r%dx = x%dx
  end function

  elemental function real_add_dual_3_real64(x, y) result(r)
    !! add dual to real number with 3 components (real64)
    real(kind=real64),          intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x + y%x
    r%dx =     y%dx
  end function

  elemental function int_add_dual_3_real64(x, y) result(r)
    !! add dual to integer with 3 components (real64)
    integer,                           intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x + y%x
    r%dx =     y%dx
  end function

  elemental function dual_neg_3_real64(x) result(r)
    !! negate dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64))             :: r

    r%x  = - x%x
    r%dx = - x%dx
  end function

  elemental function dual_sub_dual_3_real64(x, y) result(r)
    !! subtract two dual numbers with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x  - y%x
    r%dx = x%dx - y%dx
  end function

  elemental function dual_sub_real_3_real64(x, y) result(r)
    !! subtract real from dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    real(kind=real64),          intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x  - y
    r%dx = x%dx
  end function

  elemental function dual_sub_int_3_real64(x, y) result(r)
    !! subtract integer from dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x  - y
    r%dx = x%dx
  end function

  elemental function real_sub_dual_3_real64(x, y) result(r)
    !! subtract dual from real number with 3 components (real64)
    real(kind=real64),          intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x - y%x
    r%dx =   - y%dx
  end function

  elemental function int_sub_dual_3_real64(x, y) result(r)
    !! subtract dual from integer with 3 components (real64)
    integer,                           intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x - y%x
    r%dx =   - y%dx
  end function

  elemental function dual_mul_dual_3_real64(x, y) result(r)
    !! multiply two dual numbers with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x * y%x
    r%dx = x%dx * y%x + x%x * y%dx
  end function

  elemental function dual_mul_real_3_real64(x, y) result(r)
    !! multiply dual with real number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    real(kind=real64),          intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x  * y
    r%dx = x%dx * y
  end function

  elemental function dual_mul_int_3_real64(x, y) result(r)
    !! multiply dual with integer with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x  * y
    r%dx = x%dx * y
  end function

  elemental function real_mul_dual_3_real64(x, y) result(r)
    !! multiply real with dual number with 3 components (real64)
    real(kind=real64),          intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x * y%x
    r%dx = x * y%dx
  end function

  elemental function int_mul_dual_3_real64(x, y) result(r)
    !! multiply integer with dual number with 3 components (real64)
    integer,                           intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x * y%x
    r%dx = x * y%dx
  end function

  elemental function dual_div_dual_3_real64(x, y) result(r)
    !! divide two dual numbers with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x / y%x
    r%dx = (x%dx - (x%x / y%x) * y%dx) / y%x
  end function

  elemental function dual_div_real_3_real64(x, y) result(r)
    !! divide dual by real number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    real(kind=real64),          intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x  / y
    r%dx = x%dx / y
  end function

  elemental function dual_div_int_3_real64(x, y) result(r)
    !! divide dual by integer with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x  / y
    r%dx = x%dx / y
  end function

  elemental function real_div_dual_3_real64(x, y) result(r)
    !! divide real by dual number with 3 components (real64)
    real(kind=real64),          intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x / y%x
    r%dx = - (r%x * y%dx) / y%x
  end function

  elemental function int_div_dual_3_real64(x, y) result(r)
    !! divide integer by dual number with 3 components (real64)
    integer,                           intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x / y%x
    r%dx = - (r%x * y%dx) / y%x
  end function

  elemental function dual_pow_dual_3_real64(x, y) result(r)
    !! raise a dual to the power of a dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x**y%x
    r%dx = y%x * x%x ** (y%x-1) * x%dx

    if ((x%x /= 0) .or. (y%x <= 0)) then
      r%dx = r%dx + r%x * log(x%x) * y%dx
    end if
  end function

  elemental function dual_pow_real_3_real64(x, y) result(r)
    !! raise a dual to the power of a real number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    real(kind=real64),          intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x%x**y
    r%dx = y * x%x**(y-1) * x%dx
  end function

  elemental function dual_pow_int_3_real64(x, y) result(r)
    !! raise a dual number to the power of an integer with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_3(real64))             :: r

    integer                            :: p
    type(dual_3(real64)) :: r2

    ! for negative integers: compute inverse first
    if (y < 0) then
      r = 1.0_real64 / x
      p = - y
    elseif (y > 0) then
      r = x
      p = y
    else
      r%x  = 1.0_real64
      r%dx = 0.0_real64
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
        r = r ** real(p, kind=real64)
    end select
  end function

  elemental function real_pow_dual_3_real64(x, y) result(r)
    !! raise real to the power of a dual number with 3 components (real64)
    real(kind=real64),          intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r%x  = x**y%x
    if ((x == 0) .and. (y%x > 0)) then
      r%dx = 0 ! avoid log(0)
    else
      r%dx = r%x * log(x) * y%dx
    end if
  end function

  elemental function int_pow_dual_3_real64(x, y) result(r)
    !! raise integer to the power of a dual number with 3 components (real64)
    integer,                           intent(in) :: x
    type(dual_3(real64)), intent(in) :: y
    type(dual_3(real64))             :: r

    r = real(x, kind=real64) ** y
  end function

  elemental function dual_abs_3_real64(x) result(r)
    !! absolute value of dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64))             :: r

    if (x%x >= 0) then
      r =  x
    else
      r = -x
    end if
  end function

  elemental function dual_sqrt_3_real64(x) result(r)
    !! square root of dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64))             :: r

    r%x  = sqrt(x%x)
    r%dx = x%dx / (2 * r%x)
  end function

  elemental function dual_exp_3_real64(x) result(r)
    !! exponential of dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64))             :: r

    r%x  = exp(x%x)
    r%dx = r%x * x%dx
  end function

  elemental function dual_expm1_3_real64(x) result(r)
    !! exp(x) - 1 for dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64))             :: r

    r%x  = expm1(x%x)
    r%dx = exp(x%x) * x%dx
  end function

  elemental function dual_log_3_real64(x) result(r)
    !! natural logarithm of dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64))             :: r

    r%x  = log(x%x)
    r%dx = x%dx / x%x
  end function

  elemental function dual_log1p_3_real64(x) result(r)
    !! log(1 + x) for dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64))             :: r

    r%x  = log1p(x%x)
    r%dx = x%dx / (1.0_real64 + x%x)
  end function

  elemental function dual_sin_3_real64(x) result(r)
    !! sine of dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64))             :: r

    r%x  = sin(x%x)
    r%dx = cos(x%x) * x%dx
  end function

  elemental function dual_sinh_3_real64(x) result(r)
    !! sine of dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64))             :: r

    r%x  = sinh(x%x)
    r%dx = cosh(x%x) * x%dx
  end function

  elemental function dual_cos_3_real64(x) result(r)
    !! cosine of dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64))             :: r

    r%x  = cos(x%x)
    r%dx = - sin(x%x) * x%dx
  end function

  elemental function dual_cosh_3_real64(x) result(r)
    !! cosine of dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64))             :: r

    r%x  = cosh(x%x)
    r%dx = sinh(x%x) * x%dx
  end function

  elemental function dual_tan_3_real64(x) result(r)
    !! tangent of dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64))             :: r

    r%x  = tan(x%x)
    r%dx = x%dx / cos(x%x)**2
  end function

  elemental function dual_tanh_3_real64(x) result(r)
    !! tangent of dual number with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x
    type(dual_3(real64))             :: r

    r%x  = tanh(x%x)
    r%dx = x%dx / cosh(x%x)**2
  end function

  function dual_sum_3_real64(x) result(r)
    !! sum of dual number array with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x(:)
    type(dual_3(real64))             :: r

    integer :: i

    call init_dual(r, 0.0_real64)
    do i = 1, size(x)
      r = r + x(i)
    end do
  end function

  function dual_dot_dual_3_real64(x, y) result(r)
    !! dot product of two dual vectors with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x(:)
    type(dual_3(real64)), intent(in) :: y(:)
    type(dual_3(real64))             :: r

    integer :: i

    call init_dual(r, 0.0_real64)
    do i = 1, size(x)
      r = r + x(i) * y(i)
    end do
  end function

  function real_dot_dual_3_real64(x, y) result(r)
    !! dot product of real and dual vectors with 3 components (real64)
    real(kind=real64),          intent(in) :: x(:)
    type(dual_3(real64)), intent(in) :: y(:)
    type(dual_3(real64))             :: r

    integer :: i

    call init_dual(r, 0.0_real64)
    do i = 1, size(x)
      r = r + x(i) * y(i)
    end do
  end function

  function dual_dot_real_3_real64(x, y) result(r)
    !! dot product of dual and real vectors with 3 components (real64)
    type(dual_3(real64)), intent(in) :: x(:)
    real(kind=real64),          intent(in) :: y(:)
    type(dual_3(real64))             :: r

    r = dot_product(y, x)
  end function

  subroutine init_dual_3_real128(this, x, dx)
    !! initialize dual number with 3 components (real128 version)
    class(dual_3(real128)), intent(out) :: this
    real(kind=real128),           intent(in)  :: x
      !! principal value
    real(kind=real128), optional, intent(in)  :: dx(3)
      !! mark this as i-th dependency

    this%x = x
    if (present(dx)) then
      this%dx = dx
    else
      this%dx = 0
    end if
  end subroutine

  elemental function dual_add_dual_3_real128(x, y) result(r)
    !! add two dual numbers with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x  + y%x
    r%dx = x%dx + y%dx
  end function

  elemental function dual_add_real_3_real128(x, y) result(r)
    !! add real to dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    real(kind=real128),          intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x  + y
    r%dx = x%dx
  end function

  elemental function dual_add_int_3_real128(x, y) result(r)
    !! add integer to dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x  + y
    r%dx = x%dx
  end function

  elemental function real_add_dual_3_real128(x, y) result(r)
    !! add dual to real number with 3 components (real128)
    real(kind=real128),          intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x + y%x
    r%dx =     y%dx
  end function

  elemental function int_add_dual_3_real128(x, y) result(r)
    !! add dual to integer with 3 components (real128)
    integer,                           intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x + y%x
    r%dx =     y%dx
  end function

  elemental function dual_neg_3_real128(x) result(r)
    !! negate dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128))             :: r

    r%x  = - x%x
    r%dx = - x%dx
  end function

  elemental function dual_sub_dual_3_real128(x, y) result(r)
    !! subtract two dual numbers with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x  - y%x
    r%dx = x%dx - y%dx
  end function

  elemental function dual_sub_real_3_real128(x, y) result(r)
    !! subtract real from dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    real(kind=real128),          intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x  - y
    r%dx = x%dx
  end function

  elemental function dual_sub_int_3_real128(x, y) result(r)
    !! subtract integer from dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x  - y
    r%dx = x%dx
  end function

  elemental function real_sub_dual_3_real128(x, y) result(r)
    !! subtract dual from real number with 3 components (real128)
    real(kind=real128),          intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x - y%x
    r%dx =   - y%dx
  end function

  elemental function int_sub_dual_3_real128(x, y) result(r)
    !! subtract dual from integer with 3 components (real128)
    integer,                           intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x - y%x
    r%dx =   - y%dx
  end function

  elemental function dual_mul_dual_3_real128(x, y) result(r)
    !! multiply two dual numbers with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x * y%x
    r%dx = x%dx * y%x + x%x * y%dx
  end function

  elemental function dual_mul_real_3_real128(x, y) result(r)
    !! multiply dual with real number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    real(kind=real128),          intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x  * y
    r%dx = x%dx * y
  end function

  elemental function dual_mul_int_3_real128(x, y) result(r)
    !! multiply dual with integer with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x  * y
    r%dx = x%dx * y
  end function

  elemental function real_mul_dual_3_real128(x, y) result(r)
    !! multiply real with dual number with 3 components (real128)
    real(kind=real128),          intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x * y%x
    r%dx = x * y%dx
  end function

  elemental function int_mul_dual_3_real128(x, y) result(r)
    !! multiply integer with dual number with 3 components (real128)
    integer,                           intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x * y%x
    r%dx = x * y%dx
  end function

  elemental function dual_div_dual_3_real128(x, y) result(r)
    !! divide two dual numbers with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x / y%x
    r%dx = (x%dx - (x%x / y%x) * y%dx) / y%x
  end function

  elemental function dual_div_real_3_real128(x, y) result(r)
    !! divide dual by real number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    real(kind=real128),          intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x  / y
    r%dx = x%dx / y
  end function

  elemental function dual_div_int_3_real128(x, y) result(r)
    !! divide dual by integer with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x  / y
    r%dx = x%dx / y
  end function

  elemental function real_div_dual_3_real128(x, y) result(r)
    !! divide real by dual number with 3 components (real128)
    real(kind=real128),          intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x / y%x
    r%dx = - (r%x * y%dx) / y%x
  end function

  elemental function int_div_dual_3_real128(x, y) result(r)
    !! divide integer by dual number with 3 components (real128)
    integer,                           intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x / y%x
    r%dx = - (r%x * y%dx) / y%x
  end function

  elemental function dual_pow_dual_3_real128(x, y) result(r)
    !! raise a dual to the power of a dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x**y%x
    r%dx = y%x * x%x ** (y%x-1) * x%dx

    if ((x%x /= 0) .or. (y%x <= 0)) then
      r%dx = r%dx + r%x * log(x%x) * y%dx
    end if
  end function

  elemental function dual_pow_real_3_real128(x, y) result(r)
    !! raise a dual to the power of a real number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    real(kind=real128),          intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x%x**y
    r%dx = y * x%x**(y-1) * x%dx
  end function

  elemental function dual_pow_int_3_real128(x, y) result(r)
    !! raise a dual number to the power of an integer with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    integer,                           intent(in) :: y
    type(dual_3(real128))             :: r

    integer                            :: p
    type(dual_3(real128)) :: r2

    ! for negative integers: compute inverse first
    if (y < 0) then
      r = 1.0_real128 / x
      p = - y
    elseif (y > 0) then
      r = x
      p = y
    else
      r%x  = 1.0_real128
      r%dx = 0.0_real128
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
        r = r ** real(p, kind=real128)
    end select
  end function

  elemental function real_pow_dual_3_real128(x, y) result(r)
    !! raise real to the power of a dual number with 3 components (real128)
    real(kind=real128),          intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r%x  = x**y%x
    if ((x == 0) .and. (y%x > 0)) then
      r%dx = 0 ! avoid log(0)
    else
      r%dx = r%x * log(x) * y%dx
    end if
  end function

  elemental function int_pow_dual_3_real128(x, y) result(r)
    !! raise integer to the power of a dual number with 3 components (real128)
    integer,                           intent(in) :: x
    type(dual_3(real128)), intent(in) :: y
    type(dual_3(real128))             :: r

    r = real(x, kind=real128) ** y
  end function

  elemental function dual_abs_3_real128(x) result(r)
    !! absolute value of dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128))             :: r

    if (x%x >= 0) then
      r =  x
    else
      r = -x
    end if
  end function

  elemental function dual_sqrt_3_real128(x) result(r)
    !! square root of dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128))             :: r

    r%x  = sqrt(x%x)
    r%dx = x%dx / (2 * r%x)
  end function

  elemental function dual_exp_3_real128(x) result(r)
    !! exponential of dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128))             :: r

    r%x  = exp(x%x)
    r%dx = r%x * x%dx
  end function

  elemental function dual_expm1_3_real128(x) result(r)
    !! exp(x) - 1 for dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128))             :: r

    r%x  = expm1(x%x)
    r%dx = exp(x%x) * x%dx
  end function

  elemental function dual_log_3_real128(x) result(r)
    !! natural logarithm of dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128))             :: r

    r%x  = log(x%x)
    r%dx = x%dx / x%x
  end function

  elemental function dual_log1p_3_real128(x) result(r)
    !! log(1 + x) for dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128))             :: r

    r%x  = log1p(x%x)
    r%dx = x%dx / (1.0_real128 + x%x)
  end function

  elemental function dual_sin_3_real128(x) result(r)
    !! sine of dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128))             :: r

    r%x  = sin(x%x)
    r%dx = cos(x%x) * x%dx
  end function

  elemental function dual_sinh_3_real128(x) result(r)
    !! sine of dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128))             :: r

    r%x  = sinh(x%x)
    r%dx = cosh(x%x) * x%dx
  end function

  elemental function dual_cos_3_real128(x) result(r)
    !! cosine of dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128))             :: r

    r%x  = cos(x%x)
    r%dx = - sin(x%x) * x%dx
  end function

  elemental function dual_cosh_3_real128(x) result(r)
    !! cosine of dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128))             :: r

    r%x  = cosh(x%x)
    r%dx = sinh(x%x) * x%dx
  end function

  elemental function dual_tan_3_real128(x) result(r)
    !! tangent of dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128))             :: r

    r%x  = tan(x%x)
    r%dx = x%dx / cos(x%x)**2
  end function

  elemental function dual_tanh_3_real128(x) result(r)
    !! tangent of dual number with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x
    type(dual_3(real128))             :: r

    r%x  = tanh(x%x)
    r%dx = x%dx / cosh(x%x)**2
  end function

  function dual_sum_3_real128(x) result(r)
    !! sum of dual number array with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x(:)
    type(dual_3(real128))             :: r

    integer :: i

    call init_dual(r, 0.0_real128)
    do i = 1, size(x)
      r = r + x(i)
    end do
  end function

  function dual_dot_dual_3_real128(x, y) result(r)
    !! dot product of two dual vectors with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x(:)
    type(dual_3(real128)), intent(in) :: y(:)
    type(dual_3(real128))             :: r

    integer :: i

    call init_dual(r, 0.0_real128)
    do i = 1, size(x)
      r = r + x(i) * y(i)
    end do
  end function

  function real_dot_dual_3_real128(x, y) result(r)
    !! dot product of real and dual vectors with 3 components (real128)
    real(kind=real128),          intent(in) :: x(:)
    type(dual_3(real128)), intent(in) :: y(:)
    type(dual_3(real128))             :: r

    integer :: i

    call init_dual(r, 0.0_real128)
    do i = 1, size(x)
      r = r + x(i) * y(i)
    end do
  end function

  function dual_dot_real_3_real128(x, y) result(r)
    !! dot product of dual and real vectors with 3 components (real128)
    type(dual_3(real128)), intent(in) :: x(:)
    real(kind=real128),          intent(in) :: y(:)
    type(dual_3(real128))             :: r

    r = dot_product(y, x)
  end function


end module
