! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

module distribution_table_m

  use, intrinsic :: iso_c_binding,   only: c_float128, c_int
  use, intrinsic :: iso_fortran_env, only: real64, real128
  use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_POSITIVE_INF
  use, intrinsic :: omp_lib,         only: omp_get_thread_num, omp_get_num_threads

  use quad_m, only: quad
  use util_m, only: bin_search, linspace, qsort

  implicit none

  private
  public :: distribution_table, density_of_states, distribution_density

  type distribution_table
    !! lookup table for cumulative distribution function

    integer :: kmax
      !! maximum derivative stored

    real(real64), allocatable :: eta(:)
      !! chemical potential
    real(real64), allocatable :: val(:,:)
      !! values and derivatives (0:kmax+1 x num_entries)
    logical,      allocatable :: lg(:,:)
      !! logarithmic interpolation? per interval (0:kmax+1 x num_entries)
  contains
    procedure :: init => distribution_table_init
    procedure :: gen  => distribution_table_gen
    procedure :: get  => distribution_table_get
    procedure :: inv  => distribution_table_inv
    procedure :: load => distribution_table_load
    procedure :: save => distribution_table_save
  end type

  interface
    function density_of_states(t) result(Z) bind(c)
      import :: real128

      !! get density of states
      real(real128), value, intent(in) :: t
        !! energy relative to band edge (in units of k_B T)
      real(real128)                    :: Z
        !! return density of states
    end function

    function distribution_density(u, k) result(f) bind(c)
      import :: c_float128, c_int

      !! k-th derivative of distribution density (e.g. fermi-dirac or maxwell-boltzmann)
      real(c_float128), value, intent(in) :: u
        !! energy relative to chemical potential/fermi level (in units of k_B T)
      integer(c_int),   value, intent(in) :: k
        !! k-th derivative (possible values from 0 to kmax)
      real(c_float128)                    :: f
        !! return k-th distribution density
    end function
  end interface

contains

  subroutine distribution_table_init(this, dos, t_min, t_max, dist, shift_eta, eta_min, eta_max, kmax)
    !! initialize distribution table
    class(distribution_table), intent(out) :: this
    procedure(density_of_states)           :: dos
    real(real64),              intent(in)  :: t_min
      !! minimal t for Z(t) with Z(t<t_min) = 0
    real(real64),              intent(in)  :: t_max
      !! maximal t for Z(t) with Z(t>t_max) = 0
    procedure(distribution_density)        :: dist
    logical,                   intent(in)  :: shift_eta
      !! shift integration by eta or not (might improve performance, should not change result)
    real(real64),              intent(in)  :: eta_min
      !! minimal supported eta
    real(real64),              intent(in)  :: eta_max
      !! maximal supported eta
    integer,                   intent(in)  :: kmax
      !! maximal derivative

    integer,      parameter :: N0 = 1024+1
    real(real64), parameter :: RTOL = 1e-13, ATOL = 1e-16

    integer              :: i, i1, i2, i3, j1, j2, j3, k, k1, k2, k3, ithread, nthreads, n, nstack, nveta, nvval, nvlg
    integer, allocatable :: nth(:), perm(:), stack(:), itmp(:)
    logical              :: status, lg1(0:kmax)
    logical, allocatable :: lg(:,:), vlg(:), ltmp(:)

    real(real64)              :: eta1, eta2, eta3, deta, val1(0:kmax), val2(0:kmax), err1(0:kmax), err2(0:kmax), sgn
    real(real64), allocatable :: eta(:), veta(:), val(:,:), vval(:), rtmp(:)

    print "(A)", "Generating table..."

    this%kmax = kmax

    ! initial coarse grid
    allocate (eta(N0), val(0:kmax+1,N0), lg(0:kmax,N0))
    eta = linspace(eta_min, eta_max, N0)

    !$omp parallel default(none) &
    !$omp private(i,i1,i2,i3,j1,j2,j3,k,k1,k2,k3,ithread,lg1,eta1,eta2,eta3,deta,val1,val2,err1,err2,sgn,stack,nstack,vlg,nvlg,veta,nveta,vval,nvval,itmp,ltmp,rtmp) &
    !$omp shared(this,t_min,t_max,shift_eta,kmax,nthreads,n,nth,eta,val,lg)

    ithread = omp_get_thread_num() + 1

    ! allocate number of table entries per thread
    !$omp single
    nthreads = omp_get_num_threads()
    allocate (nth(0:nthreads + 1), source = 0)
    !$omp end single

    ! generate entries for coarse grid
    !$omp do schedule(dynamic)
    do i = 1, N0
      call this%gen(dos, t_min, t_max, dist, shift_eta, eta(i), val(:,i))
      lg(:,i) = .false.

      ! print "(I4,A,I4)", i, " / ", N0
    end do
    !$omp end do

    ! copy coarse grid to thread-local memory
    nvlg   = N0 * (kmax + 1)
    nveta  = N0
    nvval  = N0 * (kmax + 2)
    allocate (vlg(4 * nvlg), veta(4 * nveta), vval(4 * nvval))
    vlg( 1:nvlg)  = reshape(lg, [nvlg])
    veta(1:nveta) = eta
    vval(1:nvval) = reshape(val, [nvval])

    ! interval stack with initial size of 16
    nstack = 0
    allocate (stack(16))

    ! refinement
    !$omp do schedule(dynamic) reduction(.or.:lg)
    do i = 1, N0 - 1
      ! add single interval to stack
      nstack   = 2
      stack(1) = i
      stack(2) = i + 1

      ! refine intervals on stack
      do while (nstack > 0)
        ! get interval from stack
        i1     = stack(nstack - 1)
        i2     = stack(nstack    )
        nstack = nstack - 2

        eta1 = veta(i1)
        eta2 = veta(i2)
        eta3 = 0.5 * (eta1 + eta2)
        deta = eta2 - eta1

        ! save midpoint
        nveta = nveta + 1
        if (nveta > size(veta)) then
          allocate (rtmp(2 * nveta))
          rtmp(1:size(veta)) = veta
          call move_alloc(rtmp, veta)
        end if
        veta(nveta) = eta3
        i3 = nveta

        ! data indices
        j1 = (i1 - 1) * (kmax + 2) + 1
        j2 = (i2 - 1) * (kmax + 2) + 1
        j3 = (i3 - 1) * (kmax + 2) + 1
        k1 = i1 * (kmax + 2)
        k2 = i2 * (kmax + 2)
        k3 = i3 * (kmax + 2)

        ! make room for new point
        if (k3 > size(vval)) then
          allocate (rtmp(2 * k3))
          rtmp(1:size(vval)) = vval
          call move_alloc(rtmp, vval)
        end if
        nvval = k3

        ! generate data for new point
        call this%gen(dos, t_min, t_max, dist, shift_eta, eta3, vval(j3:k3))

        ! direct interpolation from existing data
        val1(0:kmax) = 0.5 * (vval(j1:k1-1) + vval(j2:k2-1)) + 0.125 * deta * (vval(j1+1:k1) - vval(j2+1:k2))
        err1 = abs(val1 - vval(j3:k3-1)) / (abs(vval(j3:k3-1)) + ATOL)

        ! logarithmic interpolation from existing data
        lg1 = .false.
        do k = 0, kmax
          ! values must have the same sign and not be zero
          if (vval(j1+k) * vval(j2+k) <= 0) cycle
          sgn = sign(1.0, vval(j1+k))

          ! logarithmic interpolation
          val2(k) = sgn * sqrt(vval(j1+k) * vval(j2+k)) * exp(0.125 * deta * (vval(j1+k+1) / vval(j1+k) - vval(j2+k+1) / vval(j2+k)))

          ! decide if direct or logarithmic interpolation is better
          err2(k) = abs(val2(k) - vval(j3+k)) / (abs(vval(j3+k)) + ATOL)
          if (err2(k) < err1(k)) lg1(k) = .true.
        end do

        ! lg for interval i3 --- i2
        nvlg = nvlg + size(lg1)
        if (nvlg > size(vlg)) then
          allocate (ltmp(2 * nvlg))
          ltmp(1:size(vlg)) = vlg
          call move_alloc(ltmp, vlg)
        end if
        vlg(nvlg-size(lg1)+1:nvlg) = lg1

        ! check if refinement is necessary
        if (any(merge(err2, err1, lg1) > RTOL)) then
          ! add 2 new intervals to stack
          nstack = nstack + 4
          if (nstack > size(stack)) then
            allocate (itmp(2 * nstack))
            itmp(1:size(stack)) = stack
          end if
          stack(nstack - 3) = i1
          stack(nstack - 2) = i3
          stack(nstack - 1) = i3
          stack(nstack    ) = i2
        else
          ! lg for interval i1 -- i3
          vlg(((i1 - 1) * (kmax + 1) + 1):(i1 * (kmax + 1))) = lg1
        end if
      end do

      ! store lg for coarse grid
      lg(:,i) = vlg(((i - 1) * (kmax + 1) + 1):(i * (kmax + 1)))

      ! print "(I4,A,I4)", i, " / ", N0 - 1
    end do
    !$omp end do nowait

    ! save number of points generated by this thread
    nth(ithread + 1) = nveta - N0

    !$omp barrier

    !$omp single

    ! count elements
    nth(0) = 1
    nth(1) = N0
    do i = 1, nthreads + 1
      nth(i) = nth(i) + nth(i-1)
    end do
    n = nth(nthreads + 1) - 1

    ! allocate global memory + fill in coarse grid
    allocate (this%eta(n), this%val(0:kmax+1,n), this%lg(0:kmax,n))
    this%eta(         1:N0) = eta
    this%val(0:kmax+1,1:N0) = val
    this%lg( 0:kmax,  1:N0) = lg

    !$omp end single

    ! copy local values to global memory
    j1 = N0 * (kmax + 2) + 1
    k1 = nveta * (kmax + 2)
    this%eta(  nth(ithread):nth(ithread+1)-1) = veta(N0+1:nveta)
    this%val(:,nth(ithread):nth(ithread+1)-1) = reshape(vval(j1:k1), [kmax+2, nveta - N0])
    this%lg( :,nth(ithread):nth(ithread+1)-1) = reshape( vlg(j1-N0:nvlg), [kmax+1, nveta - N0])

    ! cleanup thread-local memory
    deallocate (stack, vlg, veta, vval)

    !$omp end parallel

    ! sort by eta
    allocate (perm(size(this%eta)))
    call qsort(this%eta, perm)
    this%val = this%val(:,perm)
    this%lg  = this%lg( :,perm)
  end subroutine

  subroutine distribution_table_gen(this, dos, t_min, t_max, dist, shift_eta, eta, val)
    !! generate single entry
    class(distribution_table), intent(in)  :: this
    procedure(density_of_states)           :: dos
    real(real64),              intent(in)  :: t_min
      !! minimal t for Z(t) with Z(t<t_min) = 0
    real(real64),              intent(in)  :: t_max
      !! maximal t for Z(t) with Z(t>t_max) = 0
    procedure(distribution_density)        :: dist
    logical,                   intent(in)  :: shift_eta
      !! shift integration by eta or not
    real(real64),              intent(in)  :: eta
      !! chemical potential
    real(real64),              intent(out) :: val(0:)
      !! output value + derivatives

    integer       :: k
    real(real128) :: t_min16, t_max16, p16(0), dFda16, dFdb16, dFdp16(0), val16

    t_min16 = t_min
    t_max16 = t_max

    do k = 0, this%kmax + 1
      if (shift_eta) then
        call quad(integrand, t_min16 - eta, t_max16 - eta, p16, val16, dFda16, dFdb16, dFdp16, rtol = 1e-15_16, max_levels = 18)
      else
        call quad(integrand, t_min16, t_max16, p16, val16, dFda16, dFdb16, dFdp16, rtol = 1e-15_16, max_levels = 18)
      end if
      val(k) = real(val16)
    end do

  contains

    subroutine integrand(t, p, g, dgdt, dgdp)
      real(real128), intent(in)  :: t
      real(real128), intent(in)  :: p(:)
      real(real128), intent(out) :: g
      real(real128), intent(out) :: dgdt
      real(real128), intent(out) :: dgdp(:)

      dgdt = 0

      if (shift_eta) then
        g = dos(t + eta) * dist(t, k)
      else
        g = dos(t) * dist(t - eta, k)
      end if
      if (mod(k, 2) == 1) g = -g
    end subroutine

  end subroutine

  subroutine distribution_table_get(this, eta, k, val, dvaldeta)
    !! get k-th derivative of cumulative distribution function from table
    class(distribution_table), intent(in)  :: this
    real(real64),              intent(in)  :: eta
      !! chemical potential
    integer,                   intent(in)  :: k
      !! derivative, must be between 0 and kmax
    real(real64),              intent(out) :: val
      !! output k-th derivative of cumulative distribution function
    real(real64),              intent(out) :: dvaldeta
      !! output derivative of val wrt eta

    integer      :: i
    real(real64) :: eta1, eta2, deta, val1, dval1, val2, dval2, t, h00, h10, h01, h11, g00, g10, g01, g11

    ! find interval
    i = bin_search(this%eta, eta)
    if (i == ubound(this%eta, 1)) i = i - 1

    eta1 = this%eta(i)
    eta2 = this%eta(i+1)
    deta = eta2 - eta1
    t    = (eta - eta1) / deta

    h00 = (1 + 2 * t) * (1 - t)**2
    h10 = t * (1 - t)**2
    h01 = t**2 * (3 - 2 * t)
    h11 = t**2 * (t - 1)

    g00 = 6 * t * (t - 1)
    g10 = (3 * t - 1) * (t - 1)
    g01 = - 6 * t * (t - 1)
    g11 = t * (3 * t - 2)

    val1  = this%val(k,  i  )
    dval1 = this%val(k+1,i  )
    val2  = this%val(k,  i+1)
    dval2 = this%val(k+1,i+1)

    if (this%lg(k,i)) then
      ! logarithmic interpolation
      val      =  h00 * log(abs(val1)) + h10 * deta * dval1 / val1 + h01 * log(abs(val2)) + h11 * deta * dval2 / val2
      dvaldeta = (g00 * log(abs(val1)) + g10 * deta * dval1 / val1 + g01 * log(abs(val2)) + g11 * deta * dval2 / val2) / deta

      val      = sign(1.0, val1) * exp(val)
      dvaldeta = val * dvaldeta
    else
      ! direct interpolation
      val      =  h00 * val1 + h10 * deta * dval1 + h01 * val2 + h11 * deta * dval2
      dvaldeta = (g00 * val1 + g10 * deta * dval1 + g01 * val2 + g11 * deta * dval2) / deta
    end if
  end subroutine

  subroutine distribution_table_inv(this, F, eta, detadF)
    !! get inverse of distribution function (only for k = 0)
    class(distribution_table), intent(in)  :: this
    real(real64),              intent(in)  :: F
      !! value of distribution function
    real(real64),              intent(out) :: eta
      !! output corresponding eta
    real(real64),              intent(out) :: detadF
    !! output derivative of eta wrt F

    real(real64), parameter :: ATOL   = 5e-13
    integer,      parameter :: MAX_IT = 10

    integer      :: i, it
    real(real64) :: eta1, eta2, deta, F1, dF1, F2, dF2, t, t_min, t_max, t_old, res, dresdt, dt, err

    if ((F < this%val(0,1)) .or. (F > this%val(0,size(this%eta)))) then
      print "(A,ES25.16E3)", "F   = ", F
      print "(A,ES25.16E3)", "min = ", this%val(0,1)
      print "(A,ES25.16E3)", "max = ", this%val(0,size(this%eta))
      error stop "F is out of range"
    end if

    ! find interval
    i = bin_search(this%val(0,:), F)

    eta1 = this%eta(i)
    eta2 = this%eta(i+1)
    deta = eta2 - eta1

    F1   = this%val(0,i  )
    dF1  = this%val(1,i  )
    F2   = this%val(0,i+1)
    dF2  = this%val(1,i+1)

    if (i == ubound(this%val, 2)) i = i - 1

    ! start value
    if (this%lg(0,i)) then
      t = log(F / F1) / log(F2 / F1)
    else
      t = (F - F1) / (F2 - F1)
    end if

    ! bounds
    t_min = 0
    t_max = 1

    ! Newton iteration to get t
    err = huge(1.0)
    it  = 0
    do while (err > ATOL)
      it = it + 1
      if (it > MAX_IT) then
        print "(A,ES25.16E3)", "F = ", F
        error stop "Newton did not converge"
      end if

      ! evaluate resiudal and get Newton update
      call residual(t, res, dresdt)
      dt = - res / dresdt
      err = abs(dt)

      ! update bounds
      if (dt > 0) then
        t_min = t
      else
        t_max = t
      end if

      ! update solution
      t_old = t
      t     = t + dt

      ! bisection
      if ((t < t_min) .or. (t > t_max) .or. ((t_old == t_min) .and. (t == t_max))) then
        t = 0.5 * (t_min + t_max)
        err = min(err, t_max - t_min)
      end if
    end do

    ! get eta and derivative
    call residual(t, res, dresdt)
    eta    = eta1 + deta * t
    detadF = deta / dresdt

  contains

    subroutine residual(t, res, dresdt)
      real(real64), intent(in)  :: t
      real(real64), intent(out) :: res
      real(real64), intent(out) :: dresdt

      real(real64) :: h00, h10, h01, h11, g00, g10, g01, g11

      h00 = (1 + 2 * t) * (1 - t)**2
      h10 = t * (1 - t)**2
      h01 = t**2 * (3 - 2 * t)
      h11 = t**2 * (t - 1)

      g00 = 6 * t * (t - 1)
      g10 = (3 * t - 1) * (t - 1)
      g01 = - 6 * t * (t - 1)
      g11 = t * (3 * t - 2)

      if (this%lg(0,i)) then
        ! logarithmic interpolation
        res    = h00 * log(F1) + h10 * deta * dF1 / F1 + h01 * log(F2) + h11 * deta * dF2 / F2
        dresdt = g00 * log(F1) + g10 * deta * dF1 / F1 + g01 * log(F2) + g11 * deta * dF2 / F2

        res    = exp(res)
        dresdt = res * dresdt
        res    = res - F
      else
        ! direct interpolation
        res    = h00 * F1 + h10 * deta * dF1 + h01 * F2 + h11 * deta * dF2 - F
        dresdt = g00 * F1 + g10 * deta * dF1 + g01 * F2 + g11 * deta * dF2
      end if
    end subroutine

  end subroutine

  subroutine distribution_table_load(this, fname, status)
    !! load distribution table from file
    class(distribution_table), intent(out) :: this
    character(*),              intent(in)  :: fname
      !! filename
    logical,                   intent(out) :: status
      !! success (true) or fail (false)

    integer :: funit, num_entries

    status = .false.

    ! check if file exists
    inquire (file = fname, exist = status)
    if (.not. status) return

    ! load values
    open (newunit = funit, file = fname, status = "old", action = "read", form = "unformatted")
    read (funit) this%kmax
    read (funit) num_entries
    allocate (this%eta(num_entries), this%val(0:this%kmax+1,num_entries), this%lg(0:this%kmax,num_entries))
    read (funit) this%eta
    read (funit) this%val
    read (funit) this%lg
    close (funit)

    status = .true.
  end subroutine

  subroutine distribution_table_save(this, fname)
    !! save distribution table to file
    class(distribution_table), intent(in) :: this
    character(*),              intent(in) :: fname
      !! filename

    integer :: funit

    open (newunit = funit, file = fname, status = "replace", action = "write", form = "unformatted")
    write (funit) this%kmax
    write (funit) size(this%eta)
    write (funit) this%eta
    write (funit) this%val
    write (funit) this%lg
    close (funit)
  end subroutine

end module
