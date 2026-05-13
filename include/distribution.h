// SPDX-License-Identifier: MIT
// This file is licensed under the MIT License.
// See https://opensource.org/licenses/MIT

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Evaluate the k-th derivative of the Fermi-Dirac distribution density.
 *
 * @param u Distribution argument (typically energy offset).
 * @param k Derivative order in [0, 5].
 * @return Distribution density derivative value.
 */
extern __float128 ITHE_fermi_dirac_density(__float128 u, int k);

/**
 * Evaluate the k-th derivative of the Maxwell-Boltzmann distribution.
 *
 * @param eta Normalized chemical potential.
 * @param k Derivative order
 * @return Distribution value.
 */
extern double ITHE_maxwell_boltzmann(double eta, int k);

/**
 * Evaluate the inverse of the Maxwell-Boltzmann distribution.
 *
 * @param F Distribution value.
 * @return Normalized chemical potential.
 */
extern double ITHE_inv_maxwell_boltzmann(double F);

/**
 * Evaluate the integral of the Maxwell-Boltzmann distribution raised to k.
 *
 * @param eta1 Lower integration bound.
 * @param eta2 Upper integration bound.
 * @param k Integrand exponent.
 * @return Integral value.
 */
extern double ITHE_int_maxwell_boltzmann(double eta1, double eta2, int k);

/**
 * Initialize lookup tables for Fermi-Dirac 1/2 integrals.
 *
 * @param eta_min Minimum supported eta.
 * @param eta_max Maximum supported eta.
 * @param logging Logging level.
 */
extern void ITHE_init_fermi_dirac(double eta_min, double eta_max, int logging);

/**
 * Evaluate the k-th derivative of the Fermi-Dirac integral of order 1/2.
 *
 * @param eta Normalized chemical potential.
 * @param k Derivative order in [0, 5].
 * @return Distribution value.
 */
extern double ITHE_fermi_dirac_1h(double eta, int k);

/**
 * Evaluate the inverse Fermi-Dirac integral of order 1/2.
 *
 * @param F Distribution value.
 * @return Normalized chemical potential.
 */
extern double ITHE_inv_fermi_dirac_1h(double F);

/**
 * Evaluate the integral of the Fermi-Dirac 1/2 distribution raised to k.
 *
 * @param eta1 Lower integration bound.
 * @param eta2 Upper integration bound.
 * @param k Integrand exponent.
 * @return Integral value.
 */
extern double ITHE_int_fermi_dirac_1h(double eta1, double eta2, int k);

/**
 * Initialize band-tailing distribution tables.
 *
 * @param t0 Band-tailing cutoff energy.
 * @param eta_min Minimum supported eta.
 * @param eta_max Maximum supported eta.
 * @param logging Logging level.
 */
extern void ITHE_init_band_tailing(double t0, double eta_min, double eta_max, int logging);

/**
 * Evaluate the k-th derivative of the band-tailing distribution.
 *
 * @param eta Normalized chemical potential.
 * @param k Derivative order in [0, 5].
 * @return Distribution value.
 */
extern double ITHE_band_tailing(double eta, int k);

/**
 * Evaluate the inverse of the band-tailing distribution.
 *
 * @param F Distribution value.
 * @return Normalized chemical potential.
 */
extern double ITHE_inv_band_tailing(double F);

/**
 * Evaluate the integral of the band-tailing distribution raised to k.
 *
 * @param eta1 Lower integration bound.
 * @param eta2 Upper integration bound.
 * @param k Integrand exponent.
 * @return Integral value.
 */
extern double ITHE_int_band_tailing(double eta1, double eta2, int k);

/**
 * Initialize Gauss-Fermi distribution tables.
 *
 * @param sigma Gaussian width parameter.
 * @param eta_min Minimum supported eta.
 * @param eta_max Maximum supported eta.
 * @param logging Logging level.
 */
extern void ITHE_init_gauss_fermi(double sigma, double eta_min, double eta_max, int logging);

/**
 * Evaluate the k-th derivative of the Gauss-Fermi distribution.
 *
 * @param eta Normalized chemical potential.
 * @param k Derivative order in [0, 5].
 * @return Distribution value.
 */
extern double ITHE_gauss_fermi(double eta, int k);

/**
 * Evaluate the inverse Gauss-Fermi distribution.
 *
 * @param F Distribution value.
 * @return Normalized chemical potential.
 */
extern double ITHE_inv_gauss_fermi(double F);

/**
 * Evaluate the integral of the Gauss-Fermi distribution raised to k.
 *
 * @param eta1 Lower integration bound.
 * @param eta2 Upper integration bound.
 * @param k Integrand exponent.
 * @return Integral value.
 */
extern double ITHE_int_gauss_fermi(double eta1, double eta2, int k);

#ifdef __cplusplus
}
#endif

#endif
