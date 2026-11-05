// What the compiled expected information by simulation and by integration
// reads of a family, and what every family's piece shares.
//
// For each family of the scalar registry, a header pt_expect_<family>.h
// defines in namespace d7x four inline functions, which reproduce the R
// methods distrib_rng(), distrib_pdf(log = FALSE) and distrib_quantile()
// (lower tail, not log) at ONE parameter setting, value for value:
//
//   int  <fam>_rng_ok(const double* th);
//       0 where <fam>_rng() reproduces distrib_rng(distrib, nsim, th) at
//       this setting, nonzero otherwise; it draws nothing, since the caller
//       hands the whole computation back to R before any draw.
//   void <fam>_rng(const double* th, int n, double* out);
//       the n draws of distrib_rng(distrib, n, th), consuming R's random
//       stream in R's order; the caller brackets it with GetRNGstate() and
//       PutRNGstate().
//   void <fam>_pdf(const double* y, int n, const double* th, int vb, double* out);
//       distrib_pdf(distrib, y, th, log = FALSE) at n points.
//   int  <fam>_quantile(const double* p, int np, const double* th, int vb, double* out);
//       distrib_quantile(distrib, p, th) at np probabilities; 0, or nonzero
//       where R signals an error or the value cannot be reproduced.
//
// th holds the parameters in the family's order followed by its constants
// (the size of a binomial or a beta-binomial). vb is the value of a test
// that R evaluates over a whole vector (weibull3, betabinom2); the other
// families ignore it.
//
// The arithmetic is R's: x^2 is x * x and any other power R_POW(); a sum
// that R takes with sum(), rowSums() or mean() is a long double sum; no
// expression is simplified.

#ifndef D7_PT_EXPECT_COMMON_H
#define D7_PT_EXPECT_COMMON_H

// Rcpp first: Rmath is then reached through the R:: namespace (R::rnorm,
// R::qgamma, R::dt, R::lgammafn, ...), and unif_rand(), norm_rand() and
// exp_rand() are global
#include <Rcpp.h>
#include <cmath>
#include <cfloat>
#include <vector>
#include <numericals7/n7_expect.h>

namespace d7x {

using n7::r_pow;
using n7::r_fmod;
using n7::r_pmax;

// R's x^y
inline double R_POW(double x, double y) { return n7::r_pow(x, y); }

// R's ifelse(test, yes, no) on scalars, with both branches evaluated by the
// caller
inline double ifelse(bool t, double a, double b) { return t ? a : b; }

// R's vectorized exp(), expm1(), log1p(), sin(), tan() and cos() are not
// the plain libm calls: arithmetic.c of R 4.6.0 (Rexp() ... Rcos(), lines
// 1300-1356) returns 1 + x, x or 1 - x^2/2 near zero. Every call that R
// makes at the R level goes through these; code transcribed from a compiled
// kernel keeps std::.
inline double r_exp(double x) { return (std::fabs(x) <= std::sqrt(DBL_EPSILON)) ? 1. + x : std::exp(x); }
inline double r_expm1(double x) { return (std::fabs(x) <= DBL_EPSILON) ? x : std::expm1(x); }
inline double r_log1p(double x) { return (std::fabs(x) <= DBL_EPSILON) ? x : std::log1p(x); }
inline double r_sin(double x) { return (std::fabs(x) <= std::sqrt(3. * DBL_EPSILON)) ? x : std::sin(x); }
inline double r_tan(double x) { return (std::fabs(x) <= std::sqrt(1.5 * DBL_EPSILON)) ? x : std::tan(x); }
inline double r_cos(double x) {
    if (std::fabs(x) < std::sqrt(std::sqrt(12. * DBL_EPSILON))) return (1. - x * x * 0.5);
    return std::cos(x);
}

// R's sign()
inline double r_sign(double x) {
    if (ISNAN(x)) return x;
    return (double) ((x > 0) - (x < 0));
}

}  // namespace d7x

#endif
