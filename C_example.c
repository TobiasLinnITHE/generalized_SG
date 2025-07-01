// SPDX-License-Identifier: MIT
// This file is licensed under the MIT License.
// See https://opensource.org/licenses/MIT

#include <math.h>
#include <quadmath.h>
#include <stdbool.h>
#include <stdio.h>

#include "ITHE.h"

/**
 * Square root density of states function.
 *
 * @param t Energy parameter
 * @return  Density of states value (2*sqrt(t/π) for t >= 0, 0 otherwise)
 */
__float128 sqrt_dos(__float128 t) {
  if (t < 0) {
    return 0;
  } else {
    return 2 * sqrtq(t / M_PIq);
  }
}

/**
 * Gaussian density of states function.
 *
 * @param t Energy parameter
 * @return  Density of states value (Gaussian distribution with sigma=4)
 */
__float128 gaussian_dos(__float128 t) {
  const __float128 SIGMA = 4;
  return 1.0 / (sqrtq(2 * M_PIq) * SIGMA) * expq(- t*t / (2 * SIGMA * SIGMA));
}

/**
 * Fermi-Dirac distribution density function and its derivatives.
 *
 * @param u Energy parameter divided by temperature
 * @param k Derivative order (0-4)
 * @return  k-th derivative of the Fermi-Dirac distribution
 */
__float128 fermi_dirac_density(__float128 u, int k) {
  __float128 e, f, f0;

  e = expq(u);
  f0 = 1 / (1 + e);
  f  = 0;

  switch (k) {
  case 0:
    f = f0;
    break;
  case 1:
    if (finiteq(e)) f = (- e * f0) * f0;
    break;
  case 2:
    if (finiteq(e)) f = ((e * f0) * (expm1q(u) * f0)) * f0;
    break;
  case 3:
    if (finiteq(e)) f = (e * f0) * (- 1 + 6 * (e * f0) * (1 - e * f0)) * f0;
    break;
  case 4:
    if (finiteq(e)) f = (- (e * f0) + 14 * (e * f0)*(e * f0) - 36 * powq(e * f0, 3) + 24 * powq(e * f0, 4)) * f0;
    break;
  }

  return f;
}

/**
 * Maxwell-Boltzmann distribution function and its derivative.
 *
 * @param eta    Chemical potential
 * @param k      Derivative order (ignored, always computes 0th and 1st derivatives)
 * @param F      Output parameter for the distribution function value
 * @param dFdeta Output parameter for the derivative with respect to eta
 */
void maxwell_boltzmann(double eta, int k, double *F, double *dFdeta) {
  *F      = exp(eta);
  *dFdeta = *F;
}

/**
 * Inverse of the Maxwell-Boltzmann distribution function.
 *
 * @param F      Value of the distribution function
 * @param eta    Output parameter for the chemical potential
 * @param detadF Output parameter for the derivative of eta with respect to F
 */
void inv_maxwell_boltzmann(double F, double *eta, double *detadF) {
  *eta = log(F);
  *detadF = 1.0 / F;
}

/**
 * Bernoulli function used in the Scharfetter-Gummel scheme.
 *
 * @param x Input parameter
 * @return  Bernoulli function value B(x)
 */
double ber(double x) {
  if (x < -37) {
    return - x;
  } else if (abs(x) < 1e-9) {
    return 1.0 - 0.5 * x;
  } else if (abs(x) < 1e-4) {
    return 1.0 + x * (-0.5 + x / 12.0);
  } else if (x < 712.5) {
    return 0.5 * x * exp(-0.5 * x) / sinh(0.5 * x);
  } else {
    return 0.0;
  }
}

/**
 * Scharfetter-Gummel flux approximation.
 *
 * @param FL   Value on left-hand end of edge
 * @param FR   Value on right-hand end of edge
 * @param dpot Normalized potential drop along edge
 * @return     Approximated flux
 */
double scharfetter_gummel(double FL, double FR, double dpot) {
  return ber(-dpot) * FL - ber(dpot) * FR;
}

/**
 * Main function demonstrating the use of ITHE functions with different distribution types.
 *
 * @param argc Command line argument count (unused)
 * @param argv Command line arguments (unused)
 * @return     Exit code (0 for success)
 */
int main(int argc, char **argv) {
  double FL, FR, dpot, j, j_SG, djdFL, djdFR, djddpot;

  FL   = 1.0;
  FR   = 42.0;
  dpot = 7.0;

  j_SG = scharfetter_gummel(FL, FR, dpot);
  ITHE_get_current(maxwell_boltzmann, inv_maxwell_boltzmann, FL, FR, dpot, &j, &djdFL, &djdFR, &djddpot);
  printf("Maxwell-Boltzmann\n");
  printf("FLR     = %24.16e %24.16e\n", FL, FR);
  printf("dpot    = %24.16e\n", dpot);
  printf("j_SG    = %24.16e\n", j_SG);
  printf("j, err  = %24.16e %24.16e\n", j, fabs(j - j_SG) / fabs(j_SG));
  printf("djdFLR  = %24.16e %24.16e\n", djdFL, djdFR);
  printf("djddpot = %24.16e\n\n", djddpot);

  void *tab = ITHE_distribution_table_load("tab_FD.bin");
  if (!tab) {
    tab = ITHE_distribution_table_init(sqrt_dos, 0.0, INFINITY, fermi_dirac_density, false, -50.0, 50.0, 3);
    ITHE_distribution_table_save(tab, "tab_FD.bin");
  }
  ITHE_get_current_tab(tab, FL, FR, dpot, &j, &djdFL, &djdFR, &djddpot);
  printf("Fermi-Dirac\n");
  printf("FLR     = %24.16e %24.16e\n", FL, FR);
  printf("dpot    = %24.16e\n", dpot);
  printf("j       = %24.16e\n", j);
  printf("djdFLR  = %24.16e %24.16e\n", djdFL, djdFR);
  printf("djddpot = %24.16e\n\n", djddpot);

  ITHE_distribution_table_destruct(tab);

  FL   = 0.2;
  FR   = 0.8;
  dpot = 10.0;

  tab = ITHE_distribution_table_load("tab_GF.bin");
  if (!tab) {
    tab = ITHE_distribution_table_init(gaussian_dos, -INFINITY, INFINITY, fermi_dirac_density, true, -50.0, 50.0, 3);
    ITHE_distribution_table_save(tab, "tab_GF.bin");
  }
  ITHE_get_current_tab(tab, FL, FR, dpot, &j, &djdFL, &djdFR, &djddpot);
  printf("Gauss-Fermi\n");
  printf("FL, FR  = %24.16e %24.16e\n", FL, FR);
  printf("dpot    = %24.16e\n", dpot);
  printf("j       = %24.16e\n", j);
  printf("djdFLR  = %24.16e %24.16e\n", djdFL, djdFR);
  printf("djddpot = %24.16e\n\n", djddpot);

  return 0;
}
