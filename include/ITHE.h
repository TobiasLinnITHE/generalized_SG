// SPDX-License-Identifier: MIT
// This file is licensed under the MIT License.
// See https://opensource.org/licenses/MIT

#ifndef ITHE_H
#define ITHE_H

#include <quadmath.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  STATUS_SUCCESS          = 0,
  STATUS_INVALID_INPUT    = 1,
  STATUS_ITERATION_FAILED = 2,
  STATUS_UNKNOWN_ERROR    = 3
};

enum {
  LOGGING_SILENT = 0,
  LOGGING_ERROR  = 1,
  LOGGING_INFO   = 2,
  LOGGING_DEBUG  = 3
};

typedef double (*ITHE_integrand_64)(double x, int n_par, const double *par);
typedef __float128 (*ITHE_integrand_128)(__float128 x, int n_par, const __float128 *par);
typedef void (*ITHE_integrand_derivatives_64)(double x, int n_par, const double *par, double *f, double *df_dx, double *df_dpar);
typedef void (*ITHE_integrand_derivatives_128)(__float128 x, int n_par, const __float128 *par, __float128 *f, __float128 *df_dx, __float128 *df_dpar);

typedef double (*ITHE_distribution)(double eta, int k);
typedef void (*ITHE_inv_distribution)(double F, double *eta, double *deta_dF);
typedef double (*ITHE_int_distribution)(double eta1, double eta2, int k);

/**
 * Adaptive quadrature using double exponential integration in double precision.
 *
 * Wraps quad_m:quad for an integrand provided as a C function pointer.
 * The routine estimates relative integration error and refines adaptively.
 *
 * Recommended settings for typical real64 workloads:
 * - rtol = 1e-13 (good accuracy/cost balance; default in Fortran backend)
 * - min_levels = 2 (ensures at least two refinement levels)
 * - max_levels = 16 (usually sufficient while keeping cost bounded)
 *
 * @param intg        C function pointer to integrand: double f(double x, int n_par, const double *par)
 * @param a           Lower integration bound
 * @param b           Upper integration bound
 * @param n_par       Number of parameters in par
 * @param par         Parameter array forwarded to integrand
 * @param I           Output integral estimate
 * @param rtol        Relative tolerance target (recommendation: 1e-13)
 * @param err         Output error overestimate
 * @param min_levels  Minimum refinement levels (recommendation: 2)
 * @param max_levels  Maximum refinement levels (recommendation: 16)
 * @param n_calls     Output number of integrand evaluations
 */
extern void ITHE_quad_64(
  ITHE_integrand_64 intg,
  double a,
  double b,
  int n_par,
  const double *par,
  double *I,
  double rtol,
  double *err,
  int min_levels,
  int max_levels,
  int *n_calls
);

/**
 * Adaptive quadrature using double exponential integration in quadruple precision.
 *
 * Wraps quad_m:quad for an integrand provided as a C function pointer.
 * The routine estimates relative integration error and refines adaptively.
 *
 * Recommended settings for typical real128 workloads:
 * - rtol = 1e-13 (conservative, robust default)
 * - min_levels = 2 (ensures at least two refinement levels)
 * - max_levels = 16 (usually sufficient while keeping cost bounded)
 *
 * @param intg        C function pointer to integrand: __float128 f(__float128 x, int n_par, const __float128 *par)
 * @param a           Lower integration bound
 * @param b           Upper integration bound
 * @param n_par       Number of parameters in par
 * @param par         Parameter array forwarded to integrand
 * @param I           Output integral estimate
 * @param rtol        Relative tolerance target (recommendation: 1e-13)
 * @param err         Output error overestimate
 * @param min_levels  Minimum refinement levels (recommendation: 2)
 * @param max_levels  Maximum refinement levels (recommendation: 16)
 * @param n_calls     Output number of integrand evaluations
 */
extern void ITHE_quad_128(
  ITHE_integrand_128 intg,
  __float128 a,
  __float128 b,
  int n_par,
  const __float128 *par,
  __float128 *I,
  __float128 rtol,
  __float128 *err,
  int min_levels,
  int max_levels,
  int *n_calls
);

/**
 * Adaptive quadrature using double exponential integration with derivative outputs in double precision.
 *
 * Wraps quad_m:quad with an integrand callback that also provides
 * derivatives with respect to the integration variable and parameters.
 *
 * Recommended settings for typical real64 workloads:
 * - rtol = 1e-13 (good accuracy/cost balance; default in Fortran backend)
 * - min_levels = 2 (ensures at least two refinement levels)
 * - max_levels = 16 (usually sufficient while keeping cost bounded)
 *
 * @param intg        C function pointer to derivative integrand
 *                    void f(double x, int n_par, const double *par, double *f, double *df_dx, double *df_dpar)
 * @param a           Lower integration bound
 * @param b           Upper integration bound
 * @param n_par       Number of parameters in par
 * @param par         Parameter array forwarded to integrand
 * @param I           Output integral estimate
 * @param dI_da       Output derivative of integral with respect to a
 * @param dI_db       Output derivative of integral with respect to b
 * @param dI_dpar     Output derivative vector with respect to par
 * @param rtol        Relative tolerance target (recommendation: 1e-13)
 * @param err         Output error overestimate
 * @param min_levels  Minimum refinement levels (recommendation: 2)
 * @param max_levels  Maximum refinement levels (recommendation: 16)
 * @param n_calls     Output number of integrand evaluations
 */
extern void ITHE_quad_derivatives_64(
  ITHE_integrand_derivatives_64 intg,
  double a,
  double b,
  int n_par,
  const double *par,
  double *I,
  double *dI_da,
  double *dI_db,
  double *dI_dpar,
  double rtol,
  double *err,
  int min_levels,
  int max_levels,
  int *n_calls
);

/**
 * Adaptive quadrature using double exponential integration with derivative outputs in quadruple precision.
 *
 * Wraps quad_m:quad with an integrand callback that also provides
 * derivatives with respect to the integration variable and parameters.
 *
 * Recommended settings for typical real128 workloads:
 * - rtol = 1e-13 (conservative, robust default)
 * - min_levels = 2 (ensures at least two refinement levels)
 * - max_levels = 16 (usually sufficient while keeping cost bounded)
 *
 * @param intg        C function pointer to derivative integrand
 *                    void f(__float128 x, int n_par, const __float128 *par, __float128 *f, __float128 *df_dx, __float128 *df_dpar)
 * @param a           Lower integration bound
 * @param b           Upper integration bound
 * @param n_par       Number of parameters in par
 * @param par         Parameter array forwarded to integrand
 * @param I           Output integral estimate
 * @param dI_da       Output derivative of integral with respect to a
 * @param dI_db       Output derivative of integral with respect to b
 * @param dI_dpar     Output derivative vector with respect to par
 * @param rtol        Relative tolerance target (recommendation: 1e-13)
 * @param err         Output error overestimate
 * @param min_levels  Minimum refinement levels (recommendation: 2)
 * @param max_levels  Maximum refinement levels (recommendation: 16)
 * @param n_calls     Output number of integrand evaluations
 */
extern void ITHE_quad_derivatives_128(
  ITHE_integrand_derivatives_128 intg,
  __float128 a,
  __float128 b,
  int n_par,
  const __float128 *par,
  __float128 *I,
  __float128 *dI_da,
  __float128 *dI_db,
  __float128 *dI_dpar,
  __float128 rtol,
  __float128 *err,
  int min_levels,
  int max_levels,
  int *n_calls
);

/**
 * Initialize an integral lookup table.
 *
 * Builds an integral_table_m:integral_table for fast repeated evaluation of
 * Q(a,b) = integral_a^b f(x) dx on [xmin, xmax] and returns an opaque handle.
 *
 * The returned handle owns allocated memory and must be released with
 * ITHE_integral_table_destruct.
 *
 * @param xmin     Minimum supported x value for subsequent lookups
 * @param xmax     Maximum supported x value for subsequent lookups
 * @param intg     C function pointer to integrand callback: double f(double x, int n_par, const double *par)
 * @param n_par    Number of entries in par
 * @param par      Parameter array passed unchanged to intg
 * @param rtol     Relative tolerance used during table construction
 * @param atol     Absolute tolerance used during table construction
 * @param xtol     Minimum x-interval size used for adaptive splitting
 * @param logging  Logging level forwarded to the Fortran table builder
 * @return         Opaque pointer handle to the initialized integral table
 */
extern void *ITHE_integral_table_init(
  double xmin,
  double xmax,
  ITHE_integrand_64 intg,
  int n_par,
  const double *par,
  double rtol,
  double atol,
  double xtol,
  int logging
);

/**
 * Destruct integral table and free memory.
 *
 * @param tab  Opaque pointer to integral table to be destroyed
 */
extern void ITHE_integral_table_destruct(
  void *tab
);

/**
 * Evaluate a definite integral using a precomputed integral table.
 *
 * Computes Q(a,b) = integral_a^b f(x) dx using the table referenced by tab.
 * For unresolved sub-intervals, the provided callback is used for local
 * quadrature refinement.
 *
 * @param tab    Opaque pointer returned by ITHE_integral_table_init
 * @param intg   C function pointer to integrand callback: double f(double x, int n_par, const double *par)
 * @param a      Lower integration bound
 * @param b      Upper integration bound
 * @param rtol   Relative tolerance used during lookup/refinement
 * @return       Integral value Q(a,b)
 */
extern double ITHE_integral_table_get(
  void *tab,
  ITHE_integrand_64 intg,
  double a,
  double b,
  double rtol
);

/**
 * Initialize a Chebyshev interpolation table.
 *
 * Builds a cheb_table_m:cheb_table for fast repeated evaluation of an
 * expensive function over [xmin, xmax].
 *
 * The returned handle owns allocated memory and must be released with
 * ITHE_cheb_table_destruct.
 *
 * @param func        C function pointer to callback: double f(double x, int n_par, const double *par)
 * @param n_par       Number of entries in par
 * @param par         Parameter array passed unchanged to func
 * @param xmin        Minimum supported x value
 * @param xmax        Maximum supported x value
 * @param dx_min      Minimum interval size for adaptive refinement
 * @param dx_max      Maximum interval size for initial partitioning
 * @param rtol        Relative interpolation tolerance
 * @param max_degree  Maximum Chebyshev polynomial degree per interval
 * @param logging     Logging level forwarded to the Fortran table builder
 * @return            Opaque pointer handle to the initialized Chebyshev table
 */
extern void *ITHE_cheb_table_init(
  ITHE_integrand_64 func,
  int n_par,
  const double *par,
  double xmin,
  double xmax,
  double dx_min,
  double dx_max,
  double rtol,
  int max_degree,
  int logging
);

/**
 * Destruct Chebyshev table and free memory.
 *
 * @param tab  Opaque pointer to Chebyshev table to be destroyed
 */
extern void ITHE_cheb_table_destruct(
  void *tab
);

/**
 * Evaluate a function value from a precomputed Chebyshev table.
 *
 * @param tab  Opaque pointer returned by ITHE_cheb_table_init
 * @param x    Evaluation point (must satisfy xmin <= x <= xmax)
 * @return     Interpolated function value at x
 */
extern double ITHE_cheb_table_eval(
  void *tab,
  double x
);

/**
 * Get normalized edge current density.
 *
 * Wrapper around generalized_SG_m:get_edge_current.
 * Uses explicit callbacks for distribution, inverse distribution, and
 * integrated distribution.
 *
 * @param dist      Function pointer to k-th derivative of cumulative distribution function: double dist(double eta, int k)
 * @param inv_dist  Function pointer to inverse of cumulative distribution function: void inv_dist(double F, double *eta, double *deta_dF)
 * @param int_dist  Function pointer to integral of cumulative distribution function: double int_dist(double eta1, double eta2, int k)
 * @param FL        Value of F on left-hand end of edge
 * @param FR        Value of F on right-hand end of edge
 * @param dpot      Normalized potential drop along edge
 * @param j         Output parameter for normalized current density
 * @param dj_dFL    Output parameter for derivative of j with respect to FL
 * @param dj_dFR    Output parameter for derivative of j with respect to FR
 * @param dj_ddpot  Output parameter for derivative of j with respect to dpot
 * @param logging   Logging level
 * @param status    Output status code
 */
extern void ITHE_get_edge_current(
  ITHE_distribution dist,
  ITHE_inv_distribution inv_dist,
  ITHE_int_distribution int_dist,
  double FL,
  double FR,
  double dpot,
  double *j,
  double *dj_dFL,
  double *dj_dFR,
  double *dj_ddpot,
  int logging,
  int *status
);

#ifdef __cplusplus
}
#endif

#endif
