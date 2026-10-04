#ifndef D7_PT_GPD_H
#define D7_PT_GPD_H

#include <Rcpp.h>
#include <cmath>

// The generalized Pareto in (sigma, xi). With z = y/sigma, t = 1 + xi z and
// u = xi z,
//   log f = -log sigma - log t - W,  W = log(t)/xi = z phi(u),
//   phi(u) = log(1 + u)/u,           S(y) = exp(-W).
// Every component below is a closed form derived offline with sympy
// (stabilita/gen_gpd.py), one exported function per order and per surface,
// and checked numerically in the tests and against exact values in
// stabilita/gpd_measure.R. The forms contain no division by xi: the
// derivatives of W that name z are elementary in 1/t, and the pure xi
// derivatives are z^(j+1) phi^(j)(u), computed by gpd_phi without
// cancellation at xi = 0.

// a*b = p + e exactly (Dekker's split, no fused multiply-add)
static inline double gpd_two_prod_err(double a, double b, double p) {
  const double sp = 134217729.0;
  double c = sp * a, ah = c - (c - a), al = a - ah;
  c = sp * b;
  double bh = c - (c - b), bl = b - bh;
  return ((ah * bh - p) + ah * bl + al * bh) + al * bl;
}

// t = 1 + xi y/sigma = (sigma + xi y)/sigma, with sigma + xi y formed from
// the exact product: near the upper end of the support, where t is a small
// difference, it keeps its relative accuracy
static inline double gpd_t(double y, double sigma, double xi) {
  const double p = xi * y;
  if (!(std::fabs(p) < 1e290)) return 1.0 + p / sigma;
  const double e = gpd_two_prod_err(xi, y, p);
  const double s = sigma + p, bb = s - sigma;
  const double err = (sigma - (s - bb)) + (p - bb);
  return (s + (err + e)) / sigma;
}

// log t: from u where u is small, from t (exact to rounding) elsewhere
static inline double gpd_logt(double u, double t) {
  return std::fabs(u) < 0.5 ? std::log1p(u) : std::log(t);
}

// phi^(j)(u) = (-1)^j j! u^-(j+1) sum_{i>j} v^i/i, v = u/t. Where |v| <= 0.75
// the tail is summed as v^(j+1) sum_m v^m/(m + j + 1), so that
// phi^(j) = (-1)^j j! t^-(j+1) sum_m v^m/(m + j + 1), with no division by u;
// elsewhere it is log t minus the head of the series.
static inline double gpd_phi(int j, double u, double v, double t, double lt) {
  double fact = 1.0;
  for (int i = 2; i <= j; i++) fact *= i;
  const double sg = (j % 2) ? -fact : fact;
  if (std::fabs(v) <= 0.75) {
    double s = 0.0, vm = 1.0;
    for (int m = 0; m < 400; m++) {
      const double term = vm / (m + j + 1.0);
      s += term;
      if (std::fabs(term) < 1e-17 * std::fabs(s)) break;
      vm *= v;
    }
    double tp = t;
    for (int i = 0; i < j; i++) tp *= t;
    return sg * s / tp;
  }
  double head = 0.0, vi = 1.0;
  for (int i = 1; i <= j; i++) {
    vi *= v;
    head += vi / i;
  }
  double up = u;
  for (int i = 0; i < j; i++) up *= u;
  return sg * (lt - head) / up;
}

static inline double gpd_ipow(double x, int n) {
  double r = 1.0;
  for (int j = 0; j < n; j++) r *= x;
  return r;
}

// The components the scalar registry reads, one function per component and
// order, called by the vector kernels (gpd.cpp, dexpected_kernels.cpp) and
// by the registry (d7_ccallable.cpp). The caller checks the support with
// gpd_in_support() and passes z = y/sigma, T1 = 1/t, is = 1/sigma and the
// series values F_j = gpd_phi(j, ...).

namespace d7 {

inline bool gpd_in_support(double y, double t) {
  return !(y < 0.0 || !(t > 0.0));
}

inline double gpd_score_sigma(double z, double x, double T1, double is) {
  const double w0 = T1*z;
  return (w0*x + w0 - 1) * gpd_ipow(is, 1);
}

inline double gpd_score_xi(double z, double T1, double F1) {
  return (-z*(F1*z + T1));
}

inline double gpd_hess_sigma_sigma(double z, double x, double T1, double is) {
  const double w0 = T1*z;
  const double w1 = 2*w0;
  const double w2 = std::pow(z, 2);
  const double w3 = std::pow(T1, 2);
  const double w4 = w2*w3;
  return (-w1*x - w1 + w4*std::pow(x, 2) + w4*x + 1) * gpd_ipow(is, 2);
}

inline double gpd_hess_xi_xi(double z, double T1, double F2) {
  const double w2 = std::pow(z, 2);
  const double w3 = std::pow(T1, 2);
  return (w2*(-F2*z + w3));
}

// the expected information exists only for xi > -1/2, NA elsewhere
inline double gpd_expected_sigma_sigma(double sv, double x) {
  const double is = 1.0 / sv;
  const double w0 = 1.0/(2*x + 1);
  return x > -1.0/2.0 ? (-w0) * gpd_ipow(is, 2) : NA_REAL;
}

inline double gpd_expected_xi_xi(double x) {
  const double w0 = 1.0/(2*x + 1);
  const double w1 = w0/(x + 1);
  return x > -1.0/2.0 ? (-2*w1) : NA_REAL;
}

// With d = 1 + 2 xi, E_ss = -1/(d s^2) and E_xx = -2F, F = 1/(d (1 + xi)),
// F' = -(3 + 4 xi) F^2
inline double gpd_dexpected_sigma_sigma_sigma(double s, double x) {
  if (x <= -0.5) return NA_REAL;
  double d = 1.0 + 2.0 * x;
  double is = 1.0 / s, is2 = is * is, is3 = is2 * is;
  double id = 1.0 / d;
  return 2.0 * id * is3;
}

inline double gpd_dexpected_xi_xi_xi(double x) {
  if (x <= -0.5) return NA_REAL;
  double d = 1.0 + 2.0 * x, F = 1.0 / (d * (1.0 + x)), c = 3.0 + 4.0 * x;
  double F1 = -c * F * F;
  return -2.0 * F1;
}

inline void gpd_score_curv(int k, double y, const double* th, double* out) {
  const double sv = th[0], x = th[1];
  const double t = gpd_t(y, sv, x);
  if (!gpd_in_support(y, t)) {
    out[0] = R_NaN;
    out[1] = R_NaN;
    return;
  }
  const double z = y / sv, u = x * z;
  const double T1 = 1.0 / t;
  const double is = 1.0 / sv;
  if (k == 0) {
    out[0] = gpd_score_sigma(z, x, T1, is);
    out[1] = gpd_hess_sigma_sigma(z, x, T1, is);
  } else {
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    out[0] = gpd_score_xi(z, T1, gpd_phi(1, u, v, t, LT));
    out[1] = gpd_hess_xi_xi(z, T1, gpd_phi(2, u, v, t, LT));
  }
}

inline void gpd_info_dinfo(int k, double y, const double* th, double* out) {
  const double sv = th[0], x = th[1];
  if (k == 0) {
    out[0] = gpd_expected_sigma_sigma(sv, x);
    out[1] = gpd_dexpected_sigma_sigma_sigma(sv, x);
  } else {
    out[0] = gpd_expected_xi_xi(x);
    out[1] = gpd_dexpected_xi_xi_xi(x);
  }
}

} // namespace d7

#endif
