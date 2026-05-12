// SPDX-License-Identifier: MIT
// This file is licensed under the MIT License.
// See https://opensource.org/licenses/MIT

#include <math.h>
#include <quadmath.h>
#include <stdbool.h>
#include <stdio.h>

#include "ITHE.h"
#include "distribution.h"

static double ber(double x) {
  if (x < -37.0) return -x;
  if (fabs(x) < 1e-9) return 1.0 - 0.5 * x;
  if (fabs(x) < 1e-4) return 1.0 + x * (-0.5 + x / 12.0);
  if (x < 712.5) return 0.5 * x * exp(-0.5 * x) / sinh(0.5 * x);
  return 0.0;
}

static double scharfetter_gummel(double FL, double FR, double dpot) {
  return ber(-dpot) * FL - ber(dpot) * FR;
}

int main(void) {
  int status;
  double FL = 1.0;
  double FR = 42.0;
  double dpot = 7.0;
  double j_SG, j, dj_dFL, dj_dFR, dj_ddpot;

  ITHE_init_fermi_dirac(-100.0, 200.0, LOGGING_INFO);
  ITHE_init_band_tailing(2.5, -100.0, 200.0, LOGGING_INFO);
  ITHE_init_gauss_fermi(4.0, -40.0, 40.0, LOGGING_INFO);

  printf("Maxwell-Boltzmann\n");
  j_SG = scharfetter_gummel(FL, FR, dpot);
  ITHE_get_edge_current(ITHE_maxwell_boltzmann, ITHE_inv_maxwell_boltzmann,  ITHE_int_maxwell_boltzmann,
                        FL, FR,  dpot,  &j,  &dj_dFL, &dj_dFR,  &dj_ddpot, LOGGING_INFO, &status);
  printf("  FL, FR   = %25.16e %25.16e\n", FL, FR);
  printf("  dpot     = %25.16e\n", dpot);
  printf("  j_SG     = %25.16e\n", j_SG);
  printf("  j, err   = %25.16e %25.16e\n", j, fabs(j - j_SG) / fabs(j_SG));
  printf("  dj_dF    = %25.16e %25.16e\n", dj_dFL, dj_dFR);
  printf("  dj_ddpot = %25.16e\n\n\n\n", dj_ddpot);

  printf("Fermi-Dirac 1/2\n");
  ITHE_get_edge_current(ITHE_fermi_dirac_1h, ITHE_inv_fermi_dirac_1h, ITHE_int_fermi_dirac_1h,
                        FL, FR, dpot, &j, &dj_dFL, &dj_dFR, &dj_ddpot, LOGGING_INFO, &status);
  printf("  FL, FR   = %25.16e %25.16e\n", FL, FR);
  printf("  dpot     = %25.16e\n", dpot);
  printf("  j        = %25.16e\n", j);
  printf("  dj_dF    = %25.16e %25.16e\n", dj_dFL, dj_dFR);
  printf("  dj_ddpot = %25.16e\n\n\n\n", dj_ddpot);

  printf("Fermi-Dirac with band-tailing\n");
  ITHE_get_edge_current(ITHE_band_tailing, ITHE_inv_band_tailing, ITHE_int_band_tailing,
                        FL, FR, dpot, &j, &dj_dFL, &dj_dFR, &dj_ddpot, LOGGING_INFO, &status);
  printf("  FL, FR   = %25.16e %25.16e\n", FL, FR);
  printf("  dpot     = %25.16e\n", dpot);
  printf("  j        = %25.16e\n", j);
  printf("  dj_dF    = %25.16e %25.16e\n", dj_dFL, dj_dFR);
  printf("  dj_ddpot = %25.16e\n\n\n\n", dj_ddpot);

  printf("Gauss-Fermi\n");
  ITHE_get_edge_current(ITHE_gauss_fermi, ITHE_inv_gauss_fermi, ITHE_int_gauss_fermi,
                        FL, FR, dpot, &j, &dj_dFL, &dj_dFR, &dj_ddpot, LOGGING_INFO, &status);
  printf("  FL, FR   = %25.16e %25.16e\n", FL, FR);
  printf("  dpot     = %25.16e\n", dpot);
  printf("  j        = %25.16e\n", j);
  printf("  dj_dF    = %25.16e %25.16e\n", dj_dFL, dj_dFR);
  printf("  dj_ddpot = %25.16e\n", dj_ddpot);

  return 0;
}


