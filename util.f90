! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

module util_m

  use ieee_arithmetic
  use iso_c_binding

  implicit none

  private
  public int2str, cstrlen, c2fstring, expm1, linspace, bin_search, qsort

  interface expm1
    module procedure :: expm1, expm1_16
  end interface

contains

  function int2str(i, fmt) result(str)
    !! convert integer to string
    integer,                intent(in)  :: i
    character(*), optional, intent(in)  :: fmt
    character(:),           allocatable :: str

    character(256) :: tmp

    if (present(fmt)) then
      write(tmp, fmt) i
    else
      write(tmp, "(I0)") i
    end if
    str = trim(tmp)
  end function

  pure function cstrlen(cstr) result(len)
    !! get length of c string
    character(1), intent(in) :: cstr(*)
    integer                  :: len

    ! go through cstr until null character is found
    len = 1
    do while (cstr(len) /= c_null_char)
      len = len + 1
    end do

    ! discard trailing null character
    len = len - 1
  end function

  pure function c2fstring(cstr) result(fstr)
    !! convert c string to fortran string
    character(1), intent(in) :: cstr(*)
    character(cstrlen(cstr)) :: fstr

    fstr = transfer(cstr(1:len(fstr)), fstr)
  end function

  elemental function expm1(x) result(e)
    !! exp(x) - 1; accurate even for x close to 0
    real, intent(in) :: x
    real             :: e

    if (ieee_class(x) == IEEE_POSITIVE_INF) then
      e = x
    else
      e = exp(x)

      if (e == 1.0) then
        e = x
      else if (e - 1.0 == - 1.0) then
        e = -1
      else if (ieee_is_finite(e)) then
        e = (e - 1.0) * x / log(e)
      end if
    end if
  end function

  elemental function expm1_16(x) result(e)
    !! quadprecision exp(x) - 1; accurate even for x close to 0
    real(kind=16), intent(in) :: x
    real(kind=16)             :: e

    if (ieee_class(x) == IEEE_POSITIVE_INF) then
      e = x
    else
      e = exp(x)

      if (e == 1.0) then
        e = x
      else if (e - 1.0 == - 1.0) then
        e = -1
      else if (ieee_is_finite(e)) then
        e = (e - 1.0) * x / log(e)
      end if
    end if
  end function

  pure function linspace(x0, x1, nx) result(x)
    !! create array of linear spaced values
    real,    intent(in) :: x0
      !! start value of x
    real,    intent(in) :: x1
      !! end value of x
    integer, intent(in) :: nx
      !! number of values
    real                :: x(nx)
      !! return array x

    integer :: i
    real    :: dx

    if (nx < 1) return
    x(nx) = x1
    if (nx == 1) return

    ! spacing between values
    dx = (x1 - x0) / (nx - 1)

    x(1) = x0
    do i = 2, nx - 1
      x(i) = x0 + (i - 1) * dx
    end do
  end function

  function bin_search(x, x1) result(i)
    !! find index of nearest value in sorted array
    real, intent(in) :: x(:)
      !! sorted array
    real, intent(in) :: x1
      !! value to find
    integer          :: i
      !! return array index

    integer :: i0, i1

    ! starting interval is the whole array
    i0 = 1
    i1 = size(x)

    ! test if x1 is outside of interval [x(1), x(end)]
    if (x1 <= x(i0)) then
      i = i0
      return
    end if
    if (x1 >= x(i1)) then
      i = i1
      return
    end if

    ! binary search
    do while (i1 > (i0 + 1))
      i = (i1 + i0) / 2
      if      (x(i) < x1) then
        i0 = i
      else if (x(i) > x1) then
        i1 = i
      else
        return
      end if
    end do

    ! pick smaller index
    i = i0
  end function

  subroutine qsort(x, perm)
    real,    intent(inout) :: x(:)
    integer, intent(out)   :: perm(:)

    integer :: i

    ! init permutation array
    perm = [(i, i = 1, size(x))]

    ! start recursive quick sort
    call qsort_rec(1, size(x))

  contains

    recursive subroutine qsort_rec(i0, i1)
      integer, intent(in) :: i0
      integer, intent(in) :: i1

      integer :: mid, l, r
      real    :: pivot

      ! 0 or 1 element
      if (i0 >= i1) return

      ! quicksort with optimal sorting networks and insertion sort as base cases
      select case (i1 - i0)
        case (1) ! 2 elements
          call swap_cmp(i0,   i1  )
        case (2) ! 3 elements
          call swap_cmp(i0+1, i1  )
          call swap_cmp(i0  , i1  )
          call swap_cmp(i0  , i0+1)
        case (3) ! 4 elements
          call swap_cmp(i0  , i0+1)
          call swap_cmp(i1-1, i1  )
          call swap_cmp(i0  , i1-1)
          call swap_cmp(i0+1, i1  )
          call swap_cmp(i0+1, i1-1)
        case (4) ! 5 elements
          call swap_cmp(i0  , i0+1)
          call swap_cmp(i1-1, i1  )
          call swap_cmp(i0+2, i1  )
          call swap_cmp(i0+2, i1-1)
          call swap_cmp(i0  , i1-1)
          call swap_cmp(i0  , i0+2)
          call swap_cmp(i0+1, i1  )
          call swap_cmp(i0+1, i1-1)
          call swap_cmp(i0+1, i0+2)
        case (5:15) ! 6 - 16 elements : insertion sort
          do r = i0 + 1, i1
            pivot = x(r)
            mid = perm(r)
            do l = r - 1, i0, -1
              if (x(l) <= pivot) exit
              x(l+1) = x(l)
              perm(l+1) = perm(l)
            end do
            x(l+1) = pivot
            perm(l+1) = mid
          end do
        case default ! > 16 elements : quick sort
          ! set pivot to median of three (sort them in the process)
          mid = (i0 + i1) / 2
          call swap_cmp(mid, i1 )
          call swap_cmp(i0,  i1 )
          call swap_cmp(i0,  mid)
          pivot = x(mid)

          ! sort elements left or right according to pivot
          l = i0
          r = i1
          do while (.true.)
            l = l + 1
            do while (x(l) < pivot)
              l = l + 1
            end do
            r = r - 1
            do while (x(r) > pivot)
              r = r - 1
            end do
            if (l >= r) exit
            call swap(l, r)
          end do

          ! recursion on two sub-arrays
          call qsort_rec(i0,  r )
          call qsort_rec(r+1, i1)
      end select
    end subroutine

    subroutine swap_cmp(i0, i1)
      integer, intent(in) :: i0
      integer, intent(in) :: i1

      if (x(i1) < x(i0)) call swap(i0, i1)
    end subroutine

    subroutine swap(i0, i1)
      integer, intent(in) :: i0
      integer, intent(in) :: i1

      integer :: itmp
      real    :: tmp

      tmp      = x(i0)
      x(i0)    = x(i1)
      x(i1)    = tmp

      itmp     = perm(i0)
      perm(i0) = perm(i1)
      perm(i1) = itmp
    end subroutine

  end subroutine

end module
