! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

module integral_table_m

  use, intrinsic :: iso_c_binding,   only: c_double, c_int
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: omp_lib,         only: omp_get_max_threads, omp_get_thread_num

  use util_m, only: LOGGING_ERROR, LOGGING_INFO

  implicit none

  private
  public :: integrand, integral_table

  interface

    function integrand(x, n_par, par) result(f) bind(c)
      !! integrand for integral_table

      import :: c_double, c_int

      real(c_double), value, intent(in) :: x
        !! inner argument
      integer(c_int), value, intent(in) :: n_par
        !! number of additional parameters
      real(c_double),        intent(in) :: par(n_par)
        !! additional parameters
      real(c_double)                    :: f
        !! return integrand
    end function

  end interface

  type integral_table
    !! lookup table for definite integrals
    !! Q(a,b) = integral_{a}^{b} f(x) dx
    !! uses a tree structure to store definite integrals over many sub-intervals of varying sizes

    real(real64), allocatable :: par(:)
      !! additional parameters passed to integrand

    real(real64), allocatable :: x(:)
      !! x nodes
      !! minimum and maximum located at x(1) and x(2) (largest interval)
      !! rest of intervals can be accessed through interval midpoint indices
    real(real64), allocatable :: Q(:)
      !! Q(iQ) = integral_{x(i1)}^{x(i2)} f(t) dt
      !! largest interval located at iQ == 1
    integer,      allocatable :: midpoints(:)
      !! interval midpoint indices: x(midpoints(iQ)) = (x(i1) + x(i2)) / 2
  contains
    procedure :: init => integral_table_init
    procedure :: get  => integral_table_get

    procedure, private :: gauss => integral_table_gauss
  end type

  integer,      parameter :: NUM_GAUSS_LEGENDRE = 3
    !! number of gauss-legendre nodes used for low-level integration
  real(real64), parameter :: X_GAUSS_LEGENDRE(3) = [- sqrt(3.0_real64 / 5.0_real64), 0.0, sqrt(3.0_real64 / 5.0_real64) ]
    !! gauss-legendre nodes
  real(real64), parameter :: W_GAUSS_LEGENDRE(3) = [5.0_real64, 8.0_real64, 5.0_real64 ] / 9.0_real64
    !! gauss-legendre weights

contains

  subroutine integral_table_init(this, xmin, xmax, intg, par, rtol, atol, xtol, logging)
    class(integral_table), intent(out) :: this
    real(real64),          intent(in)  :: xmin
      !! minimal supported x
    real(real64),          intent(in)  :: xmax
      !! maximal supported x
    procedure(integrand)               :: intg
      !! integrand
    real(real64),          intent(in)  :: par(:)
      !! additional parameters passed to integrand
    real(real64),          intent(in)  :: rtol
      !! relative error tolerance
    real(real64),          intent(in)  :: atol
      !! absolute error tolerance
    real(real64),          intent(in)  :: xtol
      !! minimum interval size for x
    integer,               intent(in)  :: logging
      !! logging level

    integer, parameter :: MIN_LVL = 10

    integer                   :: lvl, max_lvl, i, i1, i2, i3, iQ, iQ0, thread_num, n_intervals, n_threads
    integer                   :: global_x_n, global_x_n0, global_Q_n, global_Q_n0, global_stack_n, global_midpoints_n
    integer                   :: thread_x_n, thread_Q_n, thread_stack_n
    integer,      allocatable :: stack_n(:), x_n(:), global_stack(:), thread_stack(:), global_midpoints(:), itmp(:)
    real(real64)              :: Q, err
    real(real64), allocatable :: global_x(:), global_Q(:), thread_x(:), thread_Q(:), tmp(:)

    ! save parameters
    this%par = par

    if(logging >= LOGGING_INFO) print "(A)", "Generating integral_table"

    max_lvl = max(ceiling(log((xmax - xmin) / xtol) / log(2.0)), MIN_LVL)

    ! init global memory
    allocate (global_x(128), global_Q(128), global_stack(128), global_midpoints(128))

    n_threads = omp_get_max_threads()
    allocate (stack_n(0:n_threads), x_n(0:n_threads), source = 0)

    ! 0-th level
    global_x_n          = 2
    global_x(1:2)       = [xmin, xmax]
    global_stack_n      = 3
    global_stack(1:3)   = [1, 2, 1]
    global_midpoints_n  = 1
    global_midpoints(1) = 0
    global_Q_n          = 1
    global_Q(1)         = this%gauss(intg, xmin, xmax)

    do lvl = 1, max_lvl
      n_intervals   = global_stack_n / 3
      if (n_intervals == 0) exit

      if(logging >= LOGGING_INFO) print "(I2,A,I2,A,I0)", lvl, " / ", max_lvl, ": ", n_intervals

      global_x_n0 = global_x_n
      global_Q_n0 = global_Q_n
      stack_n     = 0
      x_n         = 0

      !$omp parallel default(none) &
      !$omp private(Q,err,i,i1,i2,iQ,thread_num,thread_Q,thread_Q_n,thread_stack,thread_stack_n,thread_x,thread_x_n) &
      !$omp shared(this,atol,rtol,global_Q,global_Q_n,global_Q_n0,global_stack,global_stack_n,global_x,global_x_n, &
      !$omp        global_x_n0,global_midpoints,global_midpoints_n,lvl,n_intervals,n_threads,stack_n,x_n,tmp,itmp)

      thread_num = omp_get_thread_num() + 1

      ! init thread-local memory (overallocate)
      allocate (thread_x(n_intervals), thread_Q(2 * n_intervals), thread_stack(6 * n_intervals))
      thread_x_n = 0
      thread_Q_n = 0
      thread_stack_n = 0

      stack_n(thread_num) = 0
      x_n(thread_num) = 0

      !$omp do schedule(dynamic)
      do i = 1, n_intervals
        ! get interval from stack
        i1 = global_stack(i*3-2)
        i2 = global_stack(i*3-1)
        iQ = global_stack(i*3  )

        ! make room for one new point and two new intervals
        thread_x_n = thread_x_n + 1
        thread_Q_n = thread_Q_n + 2

        ! interval size and midpoint
        thread_x(thread_x_n) = 0.5 * (global_x(i1) + global_x(i2))

        ! left sub-interval
        thread_Q(thread_Q_n-1) = this%gauss(intg, global_x(i1), thread_x(thread_x_n))

        ! right sub-interval
        thread_Q(thread_Q_n) = this%gauss(intg, thread_x(thread_x_n), global_x(i2))

        ! refinement: compare sum over the two sub-intervals to previous value for complete interval
        Q   = thread_Q(thread_Q_n-1) + thread_Q(thread_Q_n)
        err = abs(Q - global_Q(iQ)) / (abs(Q) + atol)
        if ((lvl < MIN_LVL) .or. (err > rtol)) then
          ! add 2 new intervals to thread-local stack, add global_x_n0 to identify points that are in thread-local memory
          thread_stack_n = thread_stack_n + 6
          thread_stack(thread_stack_n-5:thread_stack_n) = [ &
            i1, global_x_n0 + thread_x_n, thread_Q_n-1,                   & ! left interval
            global_x_n0 + thread_x_n, i2, thread_Q_n                      & ! right interval
          ]
        end if

        ! update structure
        global_stack(i*3-2) = thread_num ! save thread number temporarily in global_stack
        global_midpoints(iQ) = global_x_n0 + thread_x_n
      end do
      !$omp end do nowait

      ! store number of generated stack intervals and points in global memory
      stack_n(thread_num) = thread_stack_n
      x_n(thread_num) = thread_x_n ! thread_Q_n == thread_x_n * 2

      !$omp barrier

      !$omp single

      ! accumulate numbers to get start indices in global memory
      do i = 2, n_threads
        stack_n(i) = stack_n(i-1) + stack_n(i)
        x_n(i) = x_n(i-1) + x_n(i)
      end do

      ! convert midpoints to global indices, use global_stack to get iQ and thread_num
      do i = 1, n_intervals
        iQ = global_stack(i*3)
        global_midpoints(iQ) = global_midpoints(iQ) + x_n(global_stack(i*3-2)-1)
      end do

      ! resize global memory
      global_x_n = global_x_n0 + x_n(n_threads)
      if (global_x_n > size(global_x)) then
        allocate (tmp(2 * global_x_n))
        tmp(1:global_x_n0) = global_x(1:global_x_n0)
        call move_alloc(tmp, global_x)
      end if
      global_Q_n = global_Q_n0 + 2 * x_n(n_threads)
      if (global_Q_n > size(global_Q)) then
        allocate (tmp(2 * global_Q_n))
        tmp(1:global_Q_n0) = global_Q(1:global_Q_n0)
        call move_alloc(tmp, global_Q)
      end if
      global_stack_n = stack_n(n_threads)
      if (global_stack_n > size(global_stack)) then
        deallocate (global_stack)
        allocate (global_stack(2 * global_stack_n))
      end if
      global_midpoints_n = global_midpoints_n + 2 * n_intervals
      if (global_midpoints_n > size(global_midpoints)) then
        allocate (itmp(2 * global_midpoints_n), source = 0)
        itmp(1:global_midpoints_n - 2 * n_intervals) = global_midpoints(1:global_midpoints_n - 2 * n_intervals)
        call move_alloc(itmp, global_midpoints)
      end if

      !$omp end single

      ! convert indices on thread-local stack to global indices
      do i = 1, thread_stack_n / 3
        if (thread_stack(3*i-2) > global_x_n0) thread_stack(3*i-2) = thread_stack(3*i-2) + x_n(thread_num-1)
        if (thread_stack(3*i-1) > global_x_n0) thread_stack(3*i-1) = thread_stack(3*i-1) + x_n(thread_num-1)
        thread_stack(3*i) = thread_stack(3*i) + global_Q_n0 + 2 * x_n(thread_num-1)
      end do

      ! copy thread-local values to global memory
      global_x((global_x_n0 +   x_n(thread_num-1) + 1):(global_x_n0 +   x_n(thread_num))) = thread_x(1:thread_x_n)
      global_Q((global_Q_n0 + 2*x_n(thread_num-1) + 1):(global_Q_n0 + 2*x_n(thread_num))) = thread_Q(1:thread_Q_n)
      global_stack((stack_n(thread_num-1) + 1):stack_n(thread_num)) = thread_stack(1:thread_stack_n)

      ! cleanup thread-local memory
      deallocate (thread_x, thread_Q, thread_stack)

      !$omp end parallel
    end do

    if ((logging >= LOGGING_ERROR) .and. (global_stack_n > 0)) print *, "Warning: some intervals are imprecise"

    ! delete Q values of upper levels
    global_stack_n = 1
    global_stack(1) = 1
    do while (global_stack_n > 0)
      iQ = global_stack(global_stack_n)
      global_stack_n = global_stack_n - 1

      ! check if leaf node
      i3 = global_midpoints(iQ)
      if (i3 == 0) cycle

      ! not a leaf node: clear value
      global_Q(iQ) = 0

      ! add children to stack
      global_stack_n = global_stack_n + 2
      if (global_stack_n > size(global_stack)) then
        allocate (itmp(2 * global_stack_n))
        itmp(1:global_stack_n-2) = global_stack(1:global_stack_n-2)
        call move_alloc(itmp, global_stack)
      end if
      global_stack(global_stack_n-1) = 2 * (i3 - 2) + 1
      global_stack(global_stack_n  ) = 2 * (i3 - 2)
    end do

    ! recalculate upper levels by summing over lower levels
    global_stack_n = 1
    global_stack(1) = 1
    do while (global_stack_n > 0)
      ! go to left-most leaf of subtree
      do
        iQ = global_stack(global_stack_n)
        i3 = global_midpoints(iQ)
        if (i3 == 0) exit ! found leaf

        ! add left sub-tree to stack
        global_stack_n = global_stack_n + 1
        if (global_stack_n > size(global_stack)) then
          allocate (itmp(2 * global_stack_n))
          itmp(1:global_stack_n-1) = global_stack(1:global_stack_n-1)
          call move_alloc(itmp, global_stack)
        end if
        global_stack(global_stack_n) = 2 * (i3 - 2)
      end do

      iQ0 = iQ
      do
        ! remove node from stack
        iQ = iQ0
        global_stack_n = global_stack_n - 1

        ! root node
        if (global_stack_n == 0) exit

        ! update Q value of direct parent
        iQ0 = global_stack(global_stack_n)
        global_Q(iQ0) = global_Q(iQ0) + global_Q(iQ)

        ! check if we come from left sub-tree
        i3 = global_midpoints(iQ0)
        if (iQ == 2 * (i3 - 2)) then
          ! add right sub-tree to stack and go back to top loop
          global_stack_n = global_stack_n + 1
          if (global_stack_n > size(global_stack)) then
            allocate (itmp(2 * global_stack_n))
            itmp(1:global_stack_n-1) = global_stack(1:global_stack_n-1)
            call move_alloc(itmp, global_stack)
          end if
          global_stack(global_stack_n) = 2 * (i3 - 2) + 1
          exit
        end if
      end do
    end do

    this%x         = global_x(1:global_x_n)
    this%Q         = global_Q(1:global_Q_n)
    this%midpoints = global_midpoints(1:global_midpoints_n)
  end subroutine

  recursive function integral_table_get(this, intg, a, b, rtol) result(Q)
    !! get value of integral
    class(integral_table), intent(in) :: this
    procedure(integrand)              :: intg
      !! integrand
    real(real64),          intent(in) :: a
      !! lower integration bound
    real(real64),          intent(in) :: b
      !! upper integration bound
    real(real64),          intent(in) :: rtol
      !! relative error tolerance
    real(real64)                      :: Q
      !! return value of integral

    integer      :: i1, i2, i3, j1, j2, j3, iQ

    if (b < a) then
      Q = - this%get(intg, b, a, rtol)
      return
    end if

    ! reset output
    Q = 0

    ! start at root node
    i1 = 1
    i2 = 2
    iQ = 1

    if ((a < this%x(i1)) .or. (b > this%x(i2))) then
      print "(A,ES25.16E3)", "a    = ", a
      print "(A,ES25.16E3)", "b    = ", b
      print "(A,ES25.16E3)", "xmin = ", this%x(i1)
      print "(A,ES25.16E3)", "xmax = ", this%x(i2)
      error stop "out of bounds"
    end if

    ! find last common ancestor
    do
      i3 = this%midpoints(iQ) ! x-index of interval midpoint

      if (i3 == 0) then
        ! both integration boundaries lie within the same unsplit interval => use gauss-legendre
        Q = this%gauss(intg, a, b)
        return
      end if

      if ((a <= this%x(i3)) .and. (b <= this%x(i3))) then
        ! both integration boundaries lie within the left sub-interval
        i2 = i3
        iQ = 2 * (i3 - 2) ! go to left sub-tree
      elseif ((a >= this%x(i3)) .and. (b >= this%x(i3))) then
        ! both integration boundaries lie within the right sub-interval
        i1 = i3
        iQ = 2 * (i3 - 2) + 1 ! go to right sub-tree
      else ! (a < this%x(i3)) .and. (b > this%x(i3))
        ! found last common ancestor, save i1, i2 and i3
        j1 = i1
        j2 = i2
        j3 = i3
        exit
      end if
    end do

    ! integrate from lower bound to midpoint
    i1 = j1
    i2 = j3
    iQ = 2 * (j3 - 2)
    do
      i3 = this%midpoints(iQ)

      if (abs(this%Q(iQ)) < rtol * abs(Q)) then
        ! optimization: stop going further down as Q will not be changed significantly within this sub-tree
        Q = Q + this%Q(iQ) * (this%x(i2) - a) / (this%x(i2) - this%x(i1)) ! approx. const. integrand
        exit
      end if

      if (i3 == 0) then
        ! unsplit interval => use gauss-legendre
        Q = Q + this%gauss(intg, a, this%x(i2))
        exit
      end if

      if (a <= this%x(i3)) then
        ! go to left sub-tree
        i2 = i3
        iQ = 2 * (i3 - 2)

        ! add complete integral over right sub-tree
        Q = Q + this%Q(iQ + 1)
      else
        ! go to right sub-tree
        i1 = i3
        iQ = 2 * (i3 - 2) + 1
      end if
    end do

    ! integrate from midpoint to upper bound
    i1 = j3
    i2 = j2
    iQ = 2 * (j3 - 2) + 1
    do
      i3 = this%midpoints(iQ)

      if (abs(this%Q(iQ)) < rtol * abs(Q)) then
        ! optimization: stop going further down as Q will not be changed significantly within this sub-tree
        Q = Q + this%Q(iQ) * (b - this%x(i1)) / (this%x(i2) - this%x(i1)) ! approx. const. integrand
        exit
      end if

      if (i3 == 0) then
        ! unsplit interval => use gauss-legendre
        Q = Q + this%gauss(intg, this%x(i1), b)
        exit
      end if

      if (b < this%x(i3)) then
        ! go to left sub-tree
        i2 = i3
        iQ = 2 * (i3 - 2)
      else
        ! go to right sub-tree
        i1 = i3
        iQ = 2 * (i3 - 2) + 1

        ! add complete integral over left sub-tree
        Q = Q + this%Q(iQ - 1)
      end if
    end do
  end function

  function integral_table_gauss(this, intg, a, b) result(Q)
    !! Gauss-Legendre quadrature
    class(integral_table), intent(in) :: this
    procedure(integrand)              :: intg
      !! integrand
    real(real64),          intent(in) :: a
      !! lower integration boundary
    real(real64),          intent(in) :: b
      !! upper integration boundary
    real(real64)                      :: Q
      !! return value of integral

    integer      :: i
    real(real64) :: x, w, f

    Q = 0
    do i = 1, NUM_GAUSS_LEGENDRE
      x = 0.5 * (a + b) + 0.5 * (b - a) * X_GAUSS_LEGENDRE(i)
      w =                 0.5 * (b - a) * W_GAUSS_LEGENDRE(i)
      Q = Q + w * intg(x, size(this%par), this%par)
    end do
  end function

end module
