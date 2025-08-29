! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

module C_interface_m

  use, intrinsic :: iso_c_binding,   only: c_bool, c_char, c_double, c_f_pointer, c_f_procpointer, c_funptr, c_int, &
    &                                      c_loc, c_null_ptr, c_ptr
  use, intrinsic :: iso_fortran_env, only: real64

  use generalized_SG_m,     only: distribution, inv_distribution, get_current_tab, get_current
  use distribution_table_m, only: distribution_table, density_of_states, distribution_density
  use util_m,               only: c2fstring

  implicit none

contains

  !> Initialize lookup table for cumulative distribution function and its derivatives
  !>
  !> @param dos        Function pointer to the density of states (quad-precision)
  !> @param t_min      Lower integration bound (can be -INFINITY)
  !> @param t_max      Upper integration bound (can be +INFINITY)
  !> @param dist       Function pointer to the k-th derivative of the distribution function (quad-precision)
  !> @param shift_eta  If true, shift integration by eta to potentially improve quadrature performance
  !> @param eta_min    Minimum supported chemical potential
  !> @param eta_max    Maximum supported chemical potential
  !> @param kmax       Maximum supported derivative of the integral with respect to eta
  !> @return           Pointer to the newly created table
  function ITHE_distribution_table_init(dos, t_min, t_max, dist, shift_eta, eta_min, eta_max, kmax) result(tab) bind(c, name = "ITHE_distribution_table_init")
    type(c_funptr),  value, intent(in) :: dos
    real(c_double),  value, intent(in) :: t_min
    real(c_double),  value, intent(in) :: t_max
    type(c_funptr),  value, intent(in) :: dist
    logical(c_bool), value, intent(in) :: shift_eta
    real(c_double),  value, intent(in) :: eta_min
    real(c_double),  value, intent(in) :: eta_max
    integer(c_int),  value, intent(in) :: kmax
    type(c_ptr)                        :: tab

    logical                                  :: fshift_eta
    type(distribution_table),        pointer :: ftab
    procedure(density_of_states),    pointer :: fdos
    procedure(distribution_density), pointer :: fdist

    call c_f_procpointer(dos, fdos)
    call c_f_procpointer(dist, fdist)
    fshift_eta = logical(shift_eta)

    allocate (ftab)
    call ftab%init(fdos, t_min, t_max, fdist, fshift_eta, eta_min, eta_max, kmax)
    tab = c_loc(ftab)
  end function

  !> Destruct lookup table and free memory
  !>
  !> @param tab  Pointer to lookup table to be destroyed
  subroutine ITHE_distribution_table_destruct(tab) bind(c, name = "ITHE_distribution_table_destruct")
    type(c_ptr), value, intent(in) :: tab

    type(distribution_table), pointer :: ftab

    call c_f_pointer(tab, ftab)
    deallocate (ftab)
  end subroutine

  !> Lookup specific value of cumulative distribution function or its derivatives
  !>
  !> @param tab       Pointer to lookup table
  !> @param eta       Chemical potential
  !> @param k         Derivative order, 0 <= k <= kmax
  !> @param val       Output parameter for the computed value
  !> @param dvaldeta  Output parameter for the derivative of val with respect to eta
  subroutine ITHE_distribution_table_get(tab, eta, k, val, dvaldeta) bind(c, name = "ITHE_distribution_table_get")
    type(c_ptr),    value, intent(in)  :: tab
    real(c_double), value, intent(in)  :: eta
    integer(c_int), value, intent(in)  :: k
    real(c_double),        intent(out) :: val
    real(c_double),        intent(out) :: dvaldeta

    type(distribution_table), pointer :: ftab

    call c_f_pointer(tab, ftab)
    call ftab%get(eta, k, val, dvaldeta)
  end subroutine

  !> Get chemical potential from value of cumulative distribution function (inverse lookup)
  !>
  !> @param tab     Pointer to lookup table
  !> @param F       Value of cumulative distribution function
  !> @param eta     Output parameter for the chemical potential
  !> @param detadF  Output parameter for the derivative of eta with respect to F
  subroutine ITHE_distribution_table_inv(tab, F, eta, detadF) bind(c, name = "ITHE_distribution_table_inv")
    type(c_ptr),    value, intent(in)  :: tab
    real(c_double), value, intent(in)  :: F
    real(c_double),        intent(out) :: eta
    real(c_double),        intent(out) :: detadF

    type(distribution_table), pointer :: ftab

    call c_f_pointer(tab, ftab)
    call ftab%inv(F, eta, detadF)
  end subroutine

  !> Load lookup table from file
  !>
  !> @param fname  Filename to load the table from
  !> @return       Pointer to the newly loaded table
  function ITHE_distribution_table_load(fname) result(tab) bind(c, name = "ITHE_distribution_table_load")
    character(len=1,kind=c_char), intent(in) :: fname(*)
    type(c_ptr)                              :: tab

    logical                           :: status
    type(distribution_table), pointer :: ftab

    allocate (ftab)
    call ftab%load(c2fstring(fname), status)
    if (.not. status) then
      deallocate (ftab)
      tab = c_null_ptr
    else
      tab = c_loc(ftab)
    end if
  end function

  !> Save lookup table to file
  !>
  !> @param tab    Pointer to lookup table
  !> @param fname  Filename to save the table to
  subroutine ITHE_distribution_table_save(tab, fname) bind(c, name = "ITHE_distribution_table_save")
    type(c_ptr), value,           intent(in) :: tab
    character(len=1,kind=c_char), intent(in) :: fname(*)

    type(distribution_table), pointer :: ftab

    call c_f_pointer(tab, ftab)
    call ftab%save(c2fstring(fname))
  end subroutine

  !> Get normalized edge current density using lookup table for cumulative distribution function
  !>
  !> @param tab      Pointer to lookup table
  !> @param FL       Value of F on left-hand end of edge
  !> @param FR       Value of F on right-hand end of edge
  !> @param dpot     Normalized potential drop along edge
  !> @param j        Output parameter for normalized current density
  !> @param djdFL    Output parameter for derivative of j with respect to FL
  !> @param djdFR    Output parameter for derivative of j with respect to FR
  !> @param djddpot  Output parameter for derivative of j with respect to dpot
  subroutine ITHE_get_current_tab(tab, FL, FR, dpot, j, djdFL, djdFR, djddpot) bind(c, name = "ITHE_get_current_tab")
    type(c_ptr),    value, intent(in)  :: tab
    real(c_double), value, intent(in)  :: FL
    real(c_double), value, intent(in)  :: FR
    real(c_double), value, intent(in)  :: dpot
    real(c_double),        intent(out) :: j
    real(c_double),        intent(out) :: djdFL
    real(c_double),        intent(out) :: djdFR
    real(c_double),        intent(out) :: djddpot

    real(real64)                      :: n(2), djdn(2)
    type(distribution_table), pointer :: ftab

    call c_f_pointer(tab, ftab)
    n = [FL, FR]
    call get_current_tab(ftab, n, dpot, j, djdn, djddpot)
    djdFL = djdn(1)
    djdFR = djdn(2)
  end subroutine

  !> Get normalized edge current density using explicit routines for calculating
  !> the cumulative distribution function, its derivatives and its inverse
  !>
  !> @param dist     Function pointer to k-th derivative of cumulative distribution function
  !> @param idist    Function pointer to inverse of cumulative distribution function
  !> @param FL       Value of F on left-hand end of edge
  !> @param FR       Value of F on right-hand end of edge
  !> @param dpot     Normalized potential drop along edge
  !> @param j        Output parameter for normalized current density
  !> @param djdFL    Output parameter for derivative of j with respect to FL
  !> @param djdFR    Output parameter for derivative of j with respect to FR
  !> @param djddpot  Output parameter for derivative of j with respect to dpot
  subroutine ITHE_get_current(dist, idist, FL, FR, dpot, j, djdFL, djdFR, djddpot) bind(c, name = "ITHE_get_current")
    type(c_funptr), value, intent(in)  :: dist
    type(c_funptr), value, intent(in)  :: idist
    real(c_double), value, intent(in)  :: FL
    real(c_double), value, intent(in)  :: FR
    real(c_double), value, intent(in)  :: dpot
    real(c_double),        intent(out) :: j
    real(c_double),        intent(out) :: djdFL
    real(c_double),        intent(out) :: djdFR
    real(c_double),        intent(out) :: djddpot

    real(real64)                         :: n(2), djdn(2)
    procedure(distribution),     pointer :: fdist
    procedure(inv_distribution), pointer :: fidist

    call c_f_procpointer( dist,  fdist)
    call c_f_procpointer(idist, fidist)

    n = [FL, FR]
    call get_current(fdist, fidist, n, dpot, j, djdn, djddpot)
    djdFL = djdn(1)
    djdFR = djdn(2)
  end subroutine

end module
