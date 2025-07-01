// SPDX-License-Identifier: MIT
// This file is licensed under the MIT License.
// See https://opensource.org/licenses/MIT

#ifndef ITHE_H
#define ITHE_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Initialize lookup table for cumulative distribution function and its derivatives.
 *
 * @param dos         Function pointer to the density of states (quad-precision).
 * @param t_min       Lower integration bound (can be -INFINITY).
 * @param t_max       Upper integration bound (can be +INFINITY).
 * @param dist        Function pointer to the k-th derivative of the distribution function (quad-precision).
 * @param shift_eta   If true, shift integration by eta to potentially improve quadrature performance.
 * @param eta_min     Minimum supported chemical potential.
 * @param eta_max     Maximum supported chemical potential.
 * @param kmax        Maximum supported derivative of the integral with respect to eta.
 * @return            Pointer to the newly created table.
 */
extern void* ITHE_distribution_table_init(
  __float128 (*dos)(__float128 t),
  double t_min,
  double t_max,
  __float128 (*dist)(__float128 u, int k),
  bool shift_eta,
  double eta_min,
  double eta_max,
  int kmax
);

/**
 * Destruct lookup table and free memory.
 *
 * @param tab Pointer to lookup table to be destroyed.
 */
extern void ITHE_distribution_table_destruct(
  void *tab
);

/**
 * Lookup specific value of cumulative distribution function or its derivatives.
 *
 * @param tab       Pointer to lookup table.
 * @param eta       Chemical potential.
 * @param k         Derivative order, 0 <= k <= kmax.
 * @param val       Output parameter for the computed value.
 * @param dvaldeta  Output parameter for the derivative of val with respect to eta.
 */
extern void ITHE_distribution_table_get(
  void *tab,
  double eta,
  int k,
  double *val,
  double *dvaldeta
);

/**
 * Get chemical potential from value of cumulative distribution function (inverse lookup).
 *
 * @param tab     Pointer to lookup table.
 * @param F       Value of cumulative distribution function.
 * @param eta     Output parameter for the chemical potential.
 * @param detadF  Output parameter for the derivative of eta with respect to F.
 */
extern void ITHE_distribution_table_inv(
  void *tab,
  double F,
  double *eta,
  double *detadF
);

/**
 * Load lookup table from file.
 *
 * @param fname  Filename to load the table from.
 * @return       Pointer to the newly loaded table.
 */
extern void* ITHE_distribution_table_load(
  const char *fname
);

/**
 * Save lookup table to file.
 *
 * @param tab    Pointer to lookup table.
 * @param fname  Filename to save the table to.
 */
extern void ITHE_distribution_table_save(
  void *tab,
  const char *fname
);

/**
 * Get normalized edge current density using lookup table for cumulative distribution function.
 *
 * @param tab      Pointer to lookup table.
 * @param FL       Value of F on left-hand end of edge.
 * @param FR       Value of F on right-hand end of edge.
 * @param dpot     Normalized potential drop along edge.
 * @param j        Output parameter for normalized current density.
 * @param djdFL    Output parameter for derivative of j with respect to FL.
 * @param djdFR    Output parameter for derivative of j with respect to FR.
 * @param djddpot  Output parameter for derivative of j with respect to dpot.
 */
extern void ITHE_get_current_tab(
  void *tab,
  double FL,
  double FR,
  double dpot,
  double *j,
  double *djdFL,
  double *djdFR,
  double *djddpot
);

/**
 * Get normalized edge current density using explicit routines for calculating
 * the cumulative distribution function, its derivatives and its inverse.
 *
 * @param dist     Function pointer to k-th derivative of cumulative distribution function.
 * @param idist    Function pointer to inverse of cumulative distribution function.
 * @param FL       Value of F on left-hand end of edge.
 * @param FR       Value of F on right-hand end of edge.
 * @param dpot     Normalized potential drop along edge.
 * @param j        Output parameter for normalized current density.
 * @param djdFL    Output parameter for derivative of j with respect to FL.
 * @param djdFR    Output parameter for derivative of j with respect to FR.
 * @param djddpot  Output parameter for derivative of j with respect to dpot.
 */
extern void ITHE_get_current(
  void (*dist)(double eta, int k, double *F, double *dFdeta),
  void (*idist)(double F, double *eta, double *detadF),
  double FL,
  double FR,
  double dpot,
  double *j,
  double *djdFL,
  double *djdFR,
  double *djddpot
);

#ifdef __cplusplus
}
#endif

#endif
