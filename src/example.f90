! SPDX-License-Identifier: MIT
! This file is licensed under the MIT License.
! See https://opensource.org/licenses/MIT

program example

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, IEEE_NEGATIVE_INF, IEEE_POSITIVE_INF
  use, intrinsic :: iso_c_binding,   only: c_double, c_int
  use, intrinsic :: iso_fortran_env, only: real64, real128

  use distribution_m, only: init_fermi_dirac, fermi_dirac_1h, inv_fermi_dirac_1h, int_fermi_dirac_1h, &
    &                       maxwell_boltzmann, inv_maxwell_boltzmann, int_maxwell_boltzmann, &
    &                       init_band_tailing, band_tailing, inv_band_tailing, int_band_tailing, &
    &                       init_gauss_fermi, gauss_fermi, inv_gauss_fermi, int_gauss_fermi
  use generalized_SG_m, only: get_edge_current
  use util_m,           only: ber, LOGGING_ERROR, LOGGING_INFO

  implicit none

  integer      :: status
  real(real64) :: F(2), dpot, j, j_SG, dj_dF(2), dj_ddpot

  F    = [1.0_real64, 42.0_real64]
  dpot = 7.0_real64

  print "(A)", "Maxwell-Boltzmann"
  j_SG = scharfetter_gummel(F, dpot)
  call get_edge_current(maxwell_boltzmann, inv_maxwell_boltzmann, int_maxwell_boltzmann, F, dpot, j, dj_dF, dj_ddpot, LOGGING_INFO, status)
  print "(A,2ES25.16E3)", "  FL, FR   = ", F
  print "(A,ES25.16E3)",  "  dpot     = ", dpot
  print "(A,ES25.16E3)",  "  j_SG     = ", j_SG
  print "(A,2ES25.16E3)", "  j, err   = ", j, abs(j - j_SG) / abs(j_SG)
  print "(A,2ES25.16E3)", "  dj_dF    = ", dj_dF
  print "(A,ES25.16E3)",  "  dj_ddpot = ", dj_ddpot
  print *
  print *
  print *

  print "(A)", "Fermi-Dirac 1/2"
  call init_fermi_dirac(-100.0_real64, 200.0_real64, LOGGING_INFO)
  call get_edge_current(fermi_dirac_1h, inv_fermi_dirac_1h, int_fermi_dirac_1h, F, dpot, j, dj_dF, dj_ddpot, LOGGING_INFO, status)
  print "(A,2ES25.16E3)", "  FL, FR   = ", F
  print "(A,ES25.16E3)",  "  dpot     = ", dpot
  print "(A,2ES25.16E3)", "  j        = ", j
  print "(A,2ES25.16E3)", "  dj_dF    = ", dj_dF
  print "(A,ES25.16E3)",  "  dj_ddpot = ", dj_ddpot
  print *
  print *
  print *

  print "(A)", "Fermi-Dirac with band-tailing"
  call init_band_tailing(2.5_real64, -100.0_real64, 200.0_real64, LOGGING_INFO)
  call get_edge_current(band_tailing, inv_band_tailing, int_band_tailing, F, dpot, j, dj_dF, dj_ddpot, LOGGING_INFO, status)
  print "(A,2ES25.16E3)", "  FL, FR   = ", F
  print "(A,ES25.16E3)",  "  dpot     = ", dpot
  print "(A,2ES25.16E3)", "  j        = ", j
  print "(A,2ES25.16E3)", "  dj_dF    = ", dj_dF
  print "(A,ES25.16E3)",  "  dj_ddpot = ", dj_ddpot
  print *
  print *
  print *

  print "(A)", "Gauss-Fermi"
  call init_gauss_fermi(4.0_real64, -40.0_real64, 40.0_real64, LOGGING_INFO)
  call get_edge_current(gauss_fermi, inv_gauss_fermi, int_gauss_fermi, F, dpot, j, dj_dF, dj_ddpot, LOGGING_INFO, status)
  print "(A,2ES25.16E3)", "  FL, FR   = ", F
  print "(A,ES25.16E3)",  "  dpot     = ", dpot
  print "(A,2ES25.16E3)", "  j        = ", j
  print "(A,2ES25.16E3)", "  dj_dF    = ", dj_dF
  print "(A,ES25.16E3)",  "  dj_ddpot = ", dj_ddpot

contains

  function scharfetter_gummel(F, dpot) result(j)
    real(real64), intent(in) :: F(2)
    real(real64), intent(in) :: dpot
    real(real64)             :: j

    j = ber(-dpot) * F(1) - ber(dpot) * F(2)
  end function

end program
