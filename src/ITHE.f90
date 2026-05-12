! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

module ITHE_m
  !! C/Fortran interface module, see include/ITHE.h for documentation

  use, intrinsic :: iso_c_binding,   only: c_bool, c_char, c_double, c_f_pointer, c_f_procpointer, c_float128, c_funptr, &
    &                                      c_int, c_loc, c_null_ptr, c_ptr
  use, intrinsic :: iso_fortran_env, only: real64, real128

  use generalized_SG_m, only: distribution, inv_distribution, int_distribution, get_edge_current
  use cheb_table_m,     only: cheb_table, funct
  use integral_table_m, only: integral_table, integrand
  use quad_m,           only: quad, integrand_64, integrand_128, integrand_derivatives_64, integrand_derivatives_128
  use util_m,           only: c2fstring

  implicit none

contains

  subroutine ITHE_quad_64(intg, a, b, n_par, par, I, rtol, err, min_levels, max_levels, n_calls) bind(c, name = "ITHE_quad_64")
    type(c_funptr), value, intent(in)  :: intg
    real(c_double), value, intent(in)  :: a
    real(c_double), value, intent(in)  :: b
    integer(c_int), value, intent(in)  :: n_par
    real(c_double),        intent(in)  :: par(n_par)
    real(c_double),        intent(out) :: I
    real(c_double), value, intent(in)  :: rtol
    real(c_double),        intent(out) :: err
    integer(c_int), value, intent(in)  :: min_levels
    integer(c_int), value, intent(in)  :: max_levels
    integer(c_int),        intent(out) :: n_calls

    procedure(integrand_64), pointer :: fintg

    call c_f_procpointer(intg, fintg)

    call quad(fintg, a, b, par, I, rtol = rtol, err = err, min_levels = min_levels, max_levels = max_levels, n_calls = n_calls)
  end subroutine

  subroutine ITHE_quad_128(intg, a, b, n_par, par, I, rtol, err, min_levels, max_levels, n_calls) bind(c, name = "ITHE_quad_128")
    type(c_funptr),   value, intent(in)  :: intg
    real(c_float128), value, intent(in)  :: a
    real(c_float128), value, intent(in)  :: b
    integer(c_int),   value, intent(in)  :: n_par
    real(c_float128),        intent(in)  :: par(n_par)
    real(c_float128),        intent(out) :: I
    real(c_float128), value, intent(in)  :: rtol
    real(c_float128),        intent(out) :: err
    integer(c_int),   value, intent(in)  :: min_levels
    integer(c_int),   value, intent(in)  :: max_levels
    integer(c_int),          intent(out) :: n_calls

    procedure(integrand_128), pointer :: fintg

    call c_f_procpointer(intg, fintg)

    call quad(fintg, a, b, par, I, rtol = rtol, err = err, min_levels = min_levels, max_levels = max_levels, n_calls = n_calls)
  end subroutine

  subroutine ITHE_quad_derivatives_64(intg, a, b, n_par, par, I, dI_da, dI_db, dI_dpar, rtol, err, min_levels, max_levels, n_calls) bind(c, name = "ITHE_quad_derivatives_64")
    type(c_funptr), value, intent(in)  :: intg
    real(c_double), value, intent(in)  :: a
    real(c_double), value, intent(in)  :: b
    integer(c_int), value, intent(in)  :: n_par
    real(c_double),        intent(in)  :: par(n_par)
    real(c_double),        intent(out) :: I
    real(c_double),        intent(out) :: dI_da
    real(c_double),        intent(out) :: dI_db
    real(c_double),        intent(out) :: dI_dpar(n_par)
    real(c_double), value, intent(in)  :: rtol
    real(c_double),        intent(out) :: err
    integer(c_int), value, intent(in)  :: min_levels
    integer(c_int), value, intent(in)  :: max_levels
    integer(c_int),        intent(out) :: n_calls

    procedure(integrand_derivatives_64), pointer :: fintg

    call c_f_procpointer(intg, fintg)

    call quad(fintg, a, b, par, I, dI_da, dI_db, dI_dpar, rtol = rtol, err = err, min_levels = min_levels, max_levels = max_levels, n_calls = n_calls)
  end subroutine

  subroutine ITHE_quad_derivatives_128(intg, a, b, n_par, par, I, dI_da, dI_db, dI_dpar, rtol, err, min_levels, max_levels, n_calls) bind(c, name = "ITHE_quad_derivatives_128")
    type(c_funptr),   value, intent(in)  :: intg
    real(c_float128), value, intent(in)  :: a
    real(c_float128), value, intent(in)  :: b
    integer(c_int),   value, intent(in)  :: n_par
    real(c_float128),        intent(in)  :: par(n_par)
    real(c_float128),        intent(out) :: I
    real(c_float128),        intent(out) :: dI_da
    real(c_float128),        intent(out) :: dI_db
    real(c_float128),        intent(out) :: dI_dpar(n_par)
    real(c_float128), value, intent(in)  :: rtol
    real(c_float128),        intent(out) :: err
    integer(c_int),   value, intent(in)  :: min_levels
    integer(c_int),   value, intent(in)  :: max_levels
    integer(c_int),          intent(out) :: n_calls

    procedure(integrand_derivatives_128), pointer :: fintg

    call c_f_procpointer(intg, fintg)

    call quad(fintg, a, b, par, I, dI_da, dI_db, dI_dpar, rtol = rtol, err = err, min_levels = min_levels, max_levels = max_levels, n_calls = n_calls)
  end subroutine

  function ITHE_integral_table_init(xmin, xmax, intg, n_par, par, rtol, atol, xtol, logging) result(tab) bind(c, name = "ITHE_integral_table_init")
    real(c_double), value, intent(in) :: xmin
    real(c_double), value, intent(in) :: xmax
    type(c_funptr), value, intent(in) :: intg
    integer(c_int), value, intent(in) :: n_par
    real(c_double),        intent(in) :: par(n_par)
    real(c_double), value, intent(in) :: rtol
    real(c_double), value, intent(in) :: atol
    real(c_double), value, intent(in) :: xtol
    integer(c_int), value, intent(in) :: logging
    type(c_ptr)                       :: tab

    type(integral_table), pointer :: ftab
    procedure(integrand), pointer :: fintg

    call c_f_procpointer(intg, fintg)
    allocate (ftab)
    call ftab%init(xmin, xmax, fintg, par, rtol, atol, xtol, logging)
    tab = c_loc(ftab)
  end function

  subroutine ITHE_integral_table_destruct(tab) bind(c, name = "ITHE_integral_table_destruct")
    type(c_ptr), value, intent(in) :: tab

    type(integral_table), pointer :: ftab

    call c_f_pointer(tab, ftab)
    deallocate (ftab)
  end subroutine

  function ITHE_integral_table_get(tab, intg, a, b, rtol) result(Q) bind(c, name = "ITHE_integral_table_get")
    type(c_ptr),    value, intent(in) :: tab
    type(c_funptr), value, intent(in) :: intg
    real(c_double), value, intent(in) :: a
    real(c_double), value, intent(in) :: b
    real(c_double), value, intent(in) :: rtol
    real(c_double)                    :: Q

    type(integral_table), pointer :: ftab
    procedure(integrand), pointer :: fintg

    call c_f_pointer(tab, ftab)
    call c_f_procpointer(intg, fintg)
    Q = ftab%get(fintg, a, b, rtol)
  end function

  function ITHE_cheb_table_init(func, n_par, par, xmin, xmax, dx_min, dx_max, rtol, max_degree, logging) result(tab) bind(c, name = "ITHE_cheb_table_init")
    type(c_funptr), value, intent(in) :: func
    integer(c_int), value, intent(in) :: n_par
    real(c_double),        intent(in) :: par(n_par)
    real(c_double), value, intent(in) :: xmin
    real(c_double), value, intent(in) :: xmax
    real(c_double), value, intent(in) :: dx_min
    real(c_double), value, intent(in) :: dx_max
    real(c_double), value, intent(in) :: rtol
    integer(c_int), value, intent(in) :: max_degree
    integer(c_int), value, intent(in) :: logging
    type(c_ptr)                       :: tab

    type(cheb_table), pointer :: ftab
    procedure(funct), pointer :: ffunc

    call c_f_procpointer(func, ffunc)
    allocate (ftab)
    call ftab%init(ffunc, par, xmin, xmax, dx_min, dx_max, rtol, max_degree, logging)
    tab = c_loc(ftab)
  end function

  subroutine ITHE_cheb_table_destruct(tab) bind(c, name = "ITHE_cheb_table_destruct")
    type(c_ptr), value, intent(in) :: tab

    type(cheb_table), pointer :: ftab

    call c_f_pointer(tab, ftab)
    deallocate (ftab)
  end subroutine

  function ITHE_cheb_table_eval(tab, x) result(f) bind(c, name = "ITHE_cheb_table_eval")
    type(c_ptr),    value, intent(in) :: tab
    real(c_double), value, intent(in) :: x
    real(c_double)                    :: f

    type(cheb_table), pointer :: ftab

    call c_f_pointer(tab, ftab)
    f = ftab%eval(x)
  end function

  subroutine ITHE_get_edge_current(dist, inv_dist, int_dist, FL, FR, dpot, j, dj_dFL, dj_dFR, dj_ddpot, logging, status) bind(c, name = "ITHE_get_edge_current")
    type(c_funptr), value, intent(in)  :: dist
    type(c_funptr), value, intent(in)  :: inv_dist
    type(c_funptr), value, intent(in)  :: int_dist
    real(c_double), value, intent(in)  :: FL
    real(c_double), value, intent(in)  :: FR
    real(c_double), value, intent(in)  :: dpot
    real(c_double),        intent(out) :: j
    real(c_double),        intent(out) :: dj_dFL
    real(c_double),        intent(out) :: dj_dFR
    real(c_double),        intent(out) :: dj_ddpot
    integer(c_int), value, intent(in)  :: logging
    integer(c_int),        intent(out) :: status

    integer                               :: fstatus
    real(real64)                          :: F(2), dj_dF(2)
    procedure(distribution),      pointer :: fdist
    procedure(inv_distribution),  pointer :: finv_dist
    procedure(int_distribution),  pointer :: fint_dist

    call c_f_procpointer(    dist,     fdist)
    call c_f_procpointer(inv_dist, finv_dist)
    call c_f_procpointer(int_dist, fint_dist)

    F = [FL, FR]
    call get_edge_current(fdist, finv_dist, fint_dist, F, dpot, j, dj_dF, dj_ddpot, int(logging), fstatus)
    dj_dFL = dj_dF(1)
    dj_dFR = dj_dF(2)
    status = int(fstatus, c_int)
  end subroutine

end module
