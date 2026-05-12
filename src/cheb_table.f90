! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

include "mkl_dfti.f90"

module cheb_table_m

  use, intrinsic :: iso_c_binding,   only: c_double, c_int
  use, intrinsic :: iso_fortran_env, only: real64

  use mkl_dfti
  use util_m, only: bin_search, int2str, linspace, LOGGING_ERROR, LOGGING_INFO, PI, qsort

  implicit none

  private
  public funct, cheb_table

  interface

    function funct(x, n_par, par) result(f) bind(c)
      !! function to interpolate

      import :: c_double, c_int

      real(c_double), value, intent(in) :: x
        !! argument
      integer(c_int), value, intent(in) :: n_par
        !! number of additional parameters
      real(c_double),        intent(in) :: par(n_par)
        !! additional parameters
      real(c_double)                    :: f
        !! value of function
    end function

  end interface

  type cheb_table
    !! Chebyshev interpolation table.
    !!
    !! This class is used to approximate an expensive function using piecewise defined Chebyshev polynomials up to a
    !! maximum degree. The intervals are adaptively chosen to minimize interpolation error. Only the necessary
    !! coefficients are stored, so the actual degree can be smaller depending on the interval.
    !!
    !! Some functions behave exponentially within certain intervals, so wherever it is advantageous instead of the
    !! function itself the logarithm of the absolute value of the function is interpolated. This can drastically reduce
    !! the number of intervals and coefficients.

    real(real64), allocatable :: x(:)
      !! interval boundaries
    real(real64), allocatable :: f_scale(:)
      !! function scaling factors
    logical,      allocatable :: use_log(:)
      !! use logarithmic instead of direct interpolation
    real(real64), allocatable :: coeff(:)
      !! flat coefficient array for all intervals
    integer,      allocatable :: icoeff(:)
      !! the coefficients for interval i are given by coeff(icoeff(i):icoeff(i+1))
  contains
    procedure :: init => cheb_table_init
    procedure :: eval => cheb_table_eval

    procedure, private :: pack => cheb_table_pack
  end type

  integer, parameter :: MAX_DEGREE_LIMIT = 32
    !! upper limit for  maximum interpolation degree

  type interval
    !! private helper class to keep track of interval bounds, function values and chebyshev coefficients

    real(real64) :: x0, x1
      !! lower and upper interval boundaries
    real(real64) :: f(0:MAX_DEGREE_LIMIT)
      !! function values evaluated at chebyshev nodes of 2nd kind
    real(real64) :: f_scale
      !! scaling factor, f_scale = max(|f|)
    real(real64) :: cheb_coeff(0:MAX_DEGREE_LIMIT)
      !! coefficients
    integer      :: degree
      !! degree for direct interpolation

    logical      :: use_log
      !! use logarithmic interpolation instead of direct
    integer      :: f_sgn
      !! sign of f if it is the same for all f, otherwise 0
    real(real64) :: log_f(0:MAX_DEGREE_LIMIT)
      !! log(|f| / f_scale), only used if all f have the same sign
    real(real64) :: log_cheb_coeff(0:MAX_DEGREE_LIMIT)
      !! coefficients for logarithmic interpolation
    integer      :: log_degree
      !! degree for logarithmic interpolation
  contains
    procedure :: set_boundaries => interval_set_boundaries
  end type

  type interval_buffer
    !! private helper class to store intervals to be processed (input) and new intervals for the next iteration (refine)

    integer :: input_index
      !! active buffer index (0 or 1)
    integer :: refine_index
      !! refinement buffer index (1 or 0, refine_index = 1 - input_index)

    integer :: n_intvl
      !! number of intervals in active buffer
    integer :: n_refine
      !! number of intervals in next buffer (refinement)

    type(interval), allocatable :: buf(:,:)
      !! internal interval buffer (capacity x 0:1)

    type(interval), pointer :: input(:)  => null()
      !! pointer to input buffer slice
    type(interval), pointer :: refine(:) => null()
      !! pointer to refinement buffer slice
  contains
    procedure :: init    => interval_buffer_init
    procedure :: reserve => interval_buffer_reserve
    procedure :: switch  => interval_buffer_switch
  end type

  type finished_buffer
    !! private helper class to store finished interval data

    integer :: n_intvl
      !! number of intervals added so far
    integer :: n_coeff
      !! number of coefficients added so far

    real(real64), allocatable :: x0(:)
    real(real64), allocatable :: f_scale(:)
    logical,      allocatable :: use_log(:)
    integer,      allocatable :: coeff_begin(:)
    integer,      allocatable :: coeff_end(:)
    real(real64), allocatable :: coeff(:)
  contains
    procedure :: init => finished_buffer_init
    procedure :: push => finished_buffer_push
  end type

contains

  subroutine cheb_table_init(this, func, par, xmin, xmax, dx_min, dx_max, rtol, max_degree, logging)
    !! initialize Chebyshev lookup table
    class(cheb_table), intent(out) :: this
    procedure(funct)               :: func
      !! function to interpolate
    real(real64),      intent(in)  :: par(:)
      !! additional parameters passed to function
    real(real64),      intent(in)  :: xmin
      !! minimal supported x
    real(real64),      intent(in)  :: xmax
      !! maximal supported x
    real(real64),      intent(in)  :: dx_min
      !! minimal interval size
    real(real64),      intent(in)  :: dx_max
      !! maximal interval size
    real(real64),      intent(in)  :: rtol
      !! relative error tolerance
    integer,           intent(in)  :: max_degree
      !! maximum polynomial degree (e.g. 16)
    integer,           intent(in)  :: logging
      !! logging level

    integer                   :: i, level, max_level, n0
    logical                   :: imprecise
    real(real64), allocatable :: xcheb(:), xinit(:), finit(:), fft(:)
    type(interval_buffer)     :: intervals
    type(finished_buffer)     :: finished

    ! input sanity checks
    if (xmin >= xmax) error stop "xmin >= xmax"
    if (dx_min <= 0) error stop "dx_min <= 0"
    if (dx_max <= 0) error stop "dx_max <= 0"
    if (dx_max < dx_min) error stop "dx_max < dx_min"
    if (dx_max > (xmax - xmin)) error stop "dx_max > (xmax - xmin)"
    if (rtol <= 0) error stop "rtol <= 0"
    if (rtol < 2 * epsilon(1.0)) error stop "rtol < 2 * epsilon"
    if (max_degree < 3) error stop "max_degree < 3"
    if (max_degree > MAX_DEGREE_LIMIT) error stop "max_degree > "//int2str(MAX_DEGREE_LIMIT)

    if (logging >= LOGGING_INFO) print "(A)", "Generating Chebyshev lookup table..."

    ! chebyshev nodes
    allocate (xcheb(0:max_degree))
    xcheb = cos((PI / max_degree) * [(i, i = 0, max_degree)])

    ! number of initial intervals
    n0 = ceiling((xmax - xmin) / dx_max)

    ! determine maximum refinement level based on minimal interval size
    max_level = ceiling(log(dx_max / dx_min) / log(2.0))

    ! evaluate function on initial interval boundaries in parallel
    allocate (xinit(n0 + 1), finit(n0 + 1))
    xinit = linspace(xmin, xmax, n0 + 1)
    !$omp parallel do schedule(dynamic) default(none) &
    !$omp private(i) &
    !$omp shared(par,xinit,finit,n0)
    do i = 1, n0 + 1
      finit(i) = func(xinit(i), size(par), par)
    end do
    !$omp end parallel do

    ! allocate interval and finished intervals buffer
    call intervals%init(n0)
    call finished%init(8 * n0, max_degree)

    ! set boundary values for initial intervals
    do i = 1, n0
      call intervals%input(i)%set_boundaries(max_degree, xinit(i), finit(i), xinit(i+1), finit(i+1))
    end do

    ! main loop
    do level = 0, max_level
      if (intervals%n_intvl == 0) exit
      if (logging >= LOGGING_INFO) then
        print "(A,I2,A,I2)",  "  Level: ", level, " / ", max_level
        print "(A,I0)",       "    #intervals: ", intervals%n_intvl
        print "(A,ES12.5E3)", "    dx:         ", (xmax - xmin) / (n0 * 2**level)
      end if

      ! evaluate function on chebyshev nodes inside of intervals
      call sample_intervals(func, par, xcheb, intervals%input)

      ! perform FFT on all intervals at once in batch mode
      call compute_fft(intervals%input, max_degree, fft)

      ! prepare interval buffer for refinement
      call intervals%reserve()

      ! process coefficients, select degree, log/linear, finished or refine
      call process_coeffs(max_degree, fft, rtol, intervals, finished, level < max_level, imprecise)

      ! switch to other interval buffer for next iteration
      call intervals%switch()
    end do

    if (imprecise .and. (logging >= LOGGING_ERROR)) print "(A)", "Warning: Could not reach required precision"

    ! pack and store coefficients in table
    call this%pack(finished, xmax)

    if (logging >= LOGGING_INFO) then
      print "(A,I0)", "  Total #intervals:    ", size(this%x) - 1
      print "(A,I0)", "  Total #use_log:      ", count(this%use_log)
      print "(A,I0)", "  Total #coefficients: ", size(this%coeff)
      print "(A,I0)", "  Avg.  #coefficients: ", nint(real(size(this%coeff),real64) / (size(this%x) - 1))
    end if
  end subroutine

  subroutine cheb_table_pack(this, finished, xmax)
    !! sort, pack and store coefficients in table
    class(cheb_table),     intent(inout) :: this
    type(finished_buffer), intent(in)    :: finished
      !! finished interval data
    real(real64),          intent(in)    :: xmax
      !! upper bound of last interval

    integer              :: i, i0, i1, j0, j1
    integer, allocatable :: perm(:)

    ! allocate memory
    allocate (this%x(finished%n_intvl + 1), this%f_scale(finished%n_intvl), this%use_log(finished%n_intvl), &
      &       this%icoeff(finished%n_intvl + 1), this%coeff(finished%n_coeff))

    ! sort meta-data by interval boundaries
    allocate (perm(finished%n_intvl))
    this%x(1:finished%n_intvl) = finished%x0(1:finished%n_intvl)
    call qsort(this%x(1:finished%n_intvl), perm)
    this%x(finished%n_intvl+1) = xmax
    this%f_scale = finished%f_scale(perm)
    this%use_log = finished%use_log(perm)

    ! reorganize coefficients
    j1 = 0
    do i = 1, finished%n_intvl
      ! lower and upper bound in original coefficient array
      i0 = finished%coeff_begin(perm(i))
      i1 = finished%coeff_end(perm(i))

      ! lower and upper bound in reorganized coefficient array
      j0 = j1 + 1
      j1 = j0 + i1 - i0

      ! store start index for this interval
      this%icoeff(i) = j0

      ! store coefficients in reversed order (leads to simpler evaluation)
      this%coeff(j0:j1) = finished%coeff(i1:i0:-1)
    end do

    ! last index
    this%icoeff(finished%n_intvl + 1) = j1 + 1
  end subroutine

  function cheb_table_eval(this, x) result(f)
    !! evaluate function at x using interpolation table
    class(cheb_table), intent(in) :: this
    real(real64),      intent(in) :: x
      !! argument, must be between xmin and xmax
    real(real64)                  :: f
      !! return interpolated function value

    integer      :: i, i0, i1, k
    real(real64) :: t, P, P_prev, P_prev2

    ! input sanity check
    if ((x < this%x(1)) .or. (x > this%x(size(this%x)))) then
      print "(A,ES25.16E3)", "xmin = ", this%x(1)
      print "(A,ES25.16E3)", "xmax = ", this%x(size(this%x))
      print "(A,ES25.16E3)", "x    = ", x
      error stop "x out of range"
    end if

    ! search for interval
    i = bin_search(this%x, x)
    if (i == size(this%x)) i = i - 1

    ! get coefficient indices
    i0 = this%icoeff(i)
    i1 = this%icoeff(i+1)-1

    ! scale argument to interval [-1, +1]
    t = (2 * x - (this%x(i) + this%x(i+1))) / (this%x(i+1) - this%x(i))

    ! evaluate polynomial by clenshaw algorithm
    P_prev = 0
    P      = 0
    do k = i0, i1
      P_prev2 = P_prev
      P_prev  = P
      P       = 2 * t * P_prev - P_prev2 + this%coeff(k)
    end do
    P = P - t * P_prev

    ! retrieve function value, for logarithmic interpolation we need to exponentiate
    if (this%use_log(i)) then
      f = this%f_scale(i) * exp(P)
    else
      f = this%f_scale(i) * P
    end if
  end function



  subroutine interval_set_boundaries(this, max_degree, x0, f0, x1, f1)
    !! set boundary values for interval
    class(interval), intent(inout) :: this
    integer,         intent(in)    :: max_degree
    real(real64),    intent(in)    :: x0
    real(real64),    intent(in)    :: f0
    real(real64),    intent(in)    :: x1
    real(real64),    intent(in)    :: f1

    ! set boundary values, this%f is reversed due to chebyshev nodes going from 1 to -1
    this%x0            = x0
    this%x1            = x1
    this%f(max_degree) = f0
    this%f(0)          = f1
  end subroutine



  subroutine interval_buffer_init(this, n_intvl)
    !! initialize interval buffer
    class(interval_buffer), target, intent(out) :: this
    integer,                        intent(in)  :: n_intvl
      !! initial number of intervals

    integer :: capacity

    capacity = 8 * n_intvl

    this%input_index  = 0
    this%refine_index = 1
    this%n_intvl      = n_intvl
    this%n_refine     = 0

    allocate (this%buf(capacity,0:1))
    this%input  => this%buf(1:this%n_intvl, this%input_index)
    this%refine => this%buf(:,this%refine_index)
  end subroutine

  subroutine interval_buffer_reserve(this)
    !! reserve memory for refinement
    class(interval_buffer), target, intent(inout) :: this

    type(interval), allocatable :: tmp(:,:)

    if (2 * this%n_intvl > size(this%buf, 1)) then
      allocate (tmp(3 * this%n_intvl,0:1))
      tmp(1:this%n_intvl,this%input_index) = this%buf(1:this%n_intvl,this%input_index)
      call move_alloc(tmp, this%buf)
      this%input  => this%buf(1:this%n_intvl, this%input_index)
      this%refine => this%buf(:,this%refine_index)
    end if
  end subroutine

  subroutine interval_buffer_switch(this)
    !! switch input and refinement buffer
    class(interval_buffer), target, intent(inout) :: this

    this%input_index  = this%refine_index
    this%refine_index = 1 - this%input_index

    this%n_intvl  = this%n_refine
    this%n_refine = 0

    this%input  => this%buf(1:this%n_intvl, this%input_index)
    this%refine => this%buf(:,this%refine_index)
  end subroutine



  subroutine finished_buffer_init(this, capacity, max_degree)
    !! allocate memory for finished intervals buffer
    class(finished_buffer), intent(out) :: this
    integer,                intent(in)  :: capacity
      !! initial capacity (number of intervals)
    integer,                intent(in)  :: max_degree
      !! maximum polynomial degree

    this%n_intvl = 0
    this%n_coeff = 0

    allocate (this%x0(capacity), this%f_scale(capacity), this%use_log(capacity), &
      &       this%coeff_begin(capacity), this%coeff_end(capacity), this%coeff(capacity * (max_degree + 1)))
  end subroutine

  subroutine finished_buffer_push(this, intvl)
    !! add new entry to finished intervals buffer
    class(finished_buffer), intent(inout) :: this
    type(interval),         intent(in)    :: intvl
      !! interval data

    integer                   :: n_intvl, n_coeff, cap
    integer,      allocatable :: itmp(:)
    logical,      allocatable :: ltmp(:)
    real(real64), allocatable :: rtmp(:)

    n_intvl = this%n_intvl + 1
    n_coeff = this%n_coeff + 1 + intvl%degree

    ! resize buffers if necessary
    if (n_intvl > size(this%x0)) then
      ! growth factor of 1.5 < 1.618... is below golden ratio leading to better memory efficiency
      cap = n_intvl + n_intvl / 2

      allocate (rtmp(cap))
      rtmp(1:this%n_intvl) = this%x0(1:this%n_intvl)
      call move_alloc(rtmp, this%x0)

      allocate (rtmp(cap))
      rtmp(1:this%n_intvl) = this%f_scale(1:this%n_intvl)
      call move_alloc(rtmp, this%f_scale)

      allocate (ltmp(cap))
      ltmp(1:this%n_intvl) = this%use_log(1:this%n_intvl)
      call move_alloc(ltmp, this%use_log)

      allocate (itmp(cap))
      itmp(1:this%n_intvl) = this%coeff_begin(1:this%n_intvl)
      call move_alloc(itmp, this%coeff_begin)

      allocate (itmp(cap))
      itmp(1:this%n_intvl) = this%coeff_end(1:this%n_intvl)
      call move_alloc(itmp, this%coeff_end)
    end if
    if (n_coeff > size(this%coeff)) then
      ! growth factor of 1.5 < 1.618... is below golden ratio leading to better memory efficiency
      cap = n_coeff + n_coeff / 2

      allocate (rtmp(2 * n_coeff))
      rtmp(1:this%n_coeff) = this%coeff(1:this%n_coeff)
      call move_alloc(rtmp, this%coeff)
    end if

    this%x0(     n_intvl) = intvl%x0
    this%f_scale(n_intvl) = intvl%f_scale
    this%use_log(n_intvl) = intvl%use_log

    this%coeff_begin(n_intvl) = this%n_coeff + 1
    this%coeff_end(  n_intvl) = n_coeff
    if (intvl%use_log) then
      this%f_scale(n_intvl) = this%f_scale(n_intvl) * intvl%f_sgn
      this%coeff(this%n_coeff + 1:n_coeff) = intvl%log_cheb_coeff(0:intvl%degree)
    else
      this%coeff(this%n_coeff + 1:n_coeff) = intvl%cheb_coeff(0:intvl%degree)
    end if

    this%n_intvl = n_intvl
    this%n_coeff = n_coeff
  end subroutine



  subroutine sample_intervals(func, par, xcheb, intervals)
    !! evaluate function on sample points in intervals
    procedure(funct)              :: func
      !! function to interpolate
    real(real64),   intent(in)    :: par(:)
      !! addtional parameters passed to function
    real(real64),   intent(in)    :: xcheb(0:)
      !! chebyshev nodes
    type(interval), intent(inout) :: intervals(:)
      !! intervals to process

    integer :: i, j, n, max_degree

    n = size(intervals)
    max_degree = ubound(xcheb, 1)

    ! Evaluate function on chebyshev nodes inside of intervals. This might be expensive, so work in parallel.
    !$omp parallel do schedule(dynamic) default(none) &
    !$omp private(i,j) &
    !$omp shared(intervals,n,par,xcheb,max_degree)
    do i = 1, n
      associate (intvl => intervals(i))
        do j = 1, max_degree - 1
          intvl%f(j) = func(0.5 * ((intvl%x0 + intvl%x1) + (intvl%x1 - intvl%x0) * xcheb(j)), size(par), par)
        end do

        ! determine f_scale to keep function values on the order of 1
        intvl%f_scale = maxval(abs(intvl%f(0:max_degree)))
        if (intvl%f_scale == 0) intvl%f_scale = 1

        ! determine sign of function
        if (intvl%f(0) == 0) then
          intvl%f_sgn = 0
        else
          intvl%f_sgn = merge(1, -1, intvl%f(0) > 0)
          do j = 1, max_degree
            if (intvl%f(j) * intvl%f_sgn <= 0) then
              intvl%f_sgn = 0
              exit
            end if
          end do
        end if

        ! if same sign, we can try logarithmic interpolation as well
        if (intvl%f_sgn /= 0) then
          intvl%log_f(0:max_degree) = log(abs(intvl%f(0:max_degree)) / intvl%f_scale)
        end if
      end associate
    end do
    !$omp end parallel do
  end subroutine

  subroutine compute_fft(intervals, max_degree, fft)
    !! compute FFT for all intervals at once in batch mode
    type(interval),            intent(in)  :: intervals(:)
    integer,                   intent(in)  :: max_degree
      !! maximum polynomial degree
    real(real64), allocatable, intent(out) :: fft(:)
      !! output FFT coefficients for all intervals in flat array

    integer                        :: i, i0, i1, n_intvl, n_transform, status
    real(real64), allocatable      :: f(:)
    type(DFTI_DESCRIPTOR), pointer :: desc

    ! number of transformations: 1 per interval and an additional 1 for logarithmic interpolation if possible
    n_intvl     = size(intervals)
    n_transform = n_intvl + count(intervals%f_sgn /= 0)

    allocate (f(2 * max_degree * n_transform), fft(2 * (max_degree + 1) * n_transform))

    ! populate input buffer with function values and periodically extend them to get discrete cosine transform.
    i1 = 0
    do i = 1, n_intvl
      i0 = i1 + 1
      i1 = i1 + 2 * max_degree
      associate (fi => f(i0:i1))
        fi(           1:max_degree+1) = intervals(i)%f(0:max_degree) / intervals(i)%f_scale
        fi(max_degree+2:2*max_degree) = fi(max_degree:2:-1) ! periodical extension
      end associate

      if (intervals(i)%f_sgn /= 0) then
        i0 = i1 + 1
        i1 = i1 + 2 * max_degree
        associate (fi => f(i0:i1))
          fi(           1:max_degree+1) = intervals(i)%log_f(0:max_degree)
          fi(max_degree+2:2*max_degree) = fi(max_degree:2:-1) ! periodical extension
        end associate
      end if
    end do

    ! perform FFT in batch mode, size per FFT is 2 * max_degree
    status = DftiCreateDescriptor(desc, DFTI_DOUBLE, DFTI_REAL, 1, 2 * max_degree)
    if (status /= 0) error stop "create: "//DftiErrorMessage(status)

    ! work out-of-place (required by batch mode)
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    if (status /= 0) error stop "placement: "//DftiErrorMessage(status)

    ! input array is split into sub-arrays of size 2 * max_degree
    status = DftiSetValue(desc, DFTI_INPUT_DISTANCE, 2 * max_degree)
    if (status /= 0) error stop "input_distance: "//DftiErrorMessage(status)

    ! output array is split into sub-arrays of size 2 * (max_degree + 1), use half of that because of complex data size
    status = DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, max_degree + 1)
    if (status /= 0) error stop "output distance: "//DftiErrorMessage(status)

    ! set the number of transformations to compute
    status = DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, n_transform)
    if (status /= 0) error stop "transform: "//DftiErrorMessage(status)

    ! finish initialization of descriptor
    status = DftiCommitDescriptor(desc)
    if (status /= 0) error stop "commit: "//DftiErrorMessage(status)

    ! perform actual computation
    status = DftiComputeForward(desc, f, fft)
    if (status /= 0) error stop "compute: "//DftiErrorMessage(status)

    ! cleanup
    status = DftiFreeDescriptor(desc)
    if (status /= 0) error stop "free: "//DftiErrorMessage(status)
  end subroutine

  subroutine process_coeffs(max_degree, fft, rtol, intervals, finished, allow_refine, imprecise)
    !! process coefficients from FFT, store interval data in buffers or create new intervals for refinement
    integer,               intent(in)    :: max_degree
      !! maximum polynomial degree
    real(real64),          intent(in)    :: fft(:)
      !! FFT output coefficients for all intervals in conjugate-even storage
    real(real64),          intent(in)    :: rtol
      !! relative tolerance for degree reduction
    type(interval_buffer), intent(inout) :: intervals
      !! buffer for input and refinement intervals
    type(finished_buffer), intent(inout) :: finished
      !! buffer for finished intervals
    logical,               intent(in)    :: allow_refine
      !! allow refinement or not
    logical,               intent(out)   :: imprecise
      !! indicates that refinement is necessary, but not allowed
      !! only ever set if allow_refine is false, which can only happen on the last iteration of the main loop

    integer :: i, i0, i1
    logical :: refine

    ! extract Chebyshev coefficients
    imprecise  = .false.
    i1         = 0
    do i = 1, intervals%n_intvl
      associate (intvl => intervals%input(i))
        i0 = i1 + 1
        i1 = i1 + 2 * (max_degree + 1)

        ! the Chebyshev coefficients are given by the scaled real parts of the FFT coeffcients (imag should be zero)
        intvl%cheb_coeff(0:max_degree) = extract_chebyshev_coefficients(fft(i0:i1))

        ! determine degree by plateau detection
        intvl%degree = get_effective_degree(intvl%cheb_coeff(0:max_degree), rtol)

        ! check whether logarithmic interpolation is better
        intvl%use_log = .false.
        if (intvl%f_sgn /= 0) then
          i0 = i1 + 1
          i1 = i1 + 2 * (max_degree + 1)

          intvl%log_cheb_coeff(0:max_degree) = extract_chebyshev_coefficients(fft(i0:i1))
          intvl%log_degree = get_effective_degree(intvl%log_cheb_coeff(0:max_degree), rtol)

          if (intvl%log_degree < intvl%degree) then
            intvl%use_log = .true.
            intvl%degree = intvl%log_degree
          end if
        end if

        ! refinement is necessary if the coefficients do not vanish fast enough
        refine = (intvl%degree >= max_degree - 2)
        if (refine .and. .not. allow_refine) then
          refine    = .false.
          imprecise = .true.
        end if

        ! add interval data to finished interval buffer if no refinement required
        if (.not. refine) then
          call finished%push(intvl)

          ! skip refinement
          cycle
        end if
      end associate

      ! interval needs to be refined
      intervals%n_refine = intervals%n_refine + 2

      ! create 2 new intervals in other interval buffer, reuse function values on boundaries and midpoint
      ! values are given in reverse order:
      !   f(max_degree) is the leftmost, f(0) the rightmost, f(max_degree/2) the middle value
      associate (intvl  => intervals%input(i), &
        &        intvl1 => intervals%refine(intervals%n_refine-1), &
        &        intvl2 => intervals%refine(intervals%n_refine))
        call intvl1%set_boundaries(max_degree, intvl%x0, intvl%f(max_degree), 0.5 * (intvl%x0 + intvl%x1), intvl%f(max_degree / 2))
        call intvl2%set_boundaries(max_degree, intvl1%x1, intvl1%f(0), intvl%x1, intvl%f(0))
      end associate
    end do
  end subroutine

  pure function extract_chebyshev_coefficients(fft) result(coeff)
    !! extract Chebyshev coefficients from FFT output array
    real(real64), intent(in) :: fft(:)
      !! FFT coefficients in conjugate-even format (R0, 0, R1, I1, R2, I2, ...)
    real(real64)             :: coeff(0:(size(fft)-2)/2)
      !! coefficients from 0:N

    integer :: degree, fft_size

    fft_size = size(fft)
    degree   = (fft_size - 2) / 2

    ! the Chebyshev coefficients can be calculated from the scaled real parts of the FFT coefficients. The imaginary
    ! parts should be numerically zero and can be ignored.
    coeff         = fft(1:fft_size-1:2) / degree
    coeff(     0) = 0.5 * coeff(     0)
    coeff(degree) = 0.5 * coeff(degree)
  end function

  pure function get_effective_degree(coeff, rtol) result(degree)
    !! find the effective degree of a polynomial
    real(real64), intent(in) :: coeff(0:)
      !! coefficients, in ascending order
    real(real64), intent(in) :: rtol
      !! relative tolerance
    integer                  :: degree
      !! maximal degree after which coefficients can be ignored

    integer      :: max_degree
    real(real64) :: tol

    max_degree = ubound(coeff, 1)

    ! find the largest index (degree) after which all coefficients are below the tolerance
    tol = rtol * maxval(abs(coeff(0:max_degree)))
    do degree = max_degree, 1, -1
      if (abs(coeff(degree)) > tol) exit
    end do
  end function

end module
