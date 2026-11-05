#ifndef D7_PT_GPD_H
#define D7_PT_GPD_H

#include <Rcpp.h>
#include <cmath>
#include "pt_parts.h"

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

// Off-diagonal entry of the observed Hessian, as gpd_hessian_cpp computes it.
inline double gpd_hess_sigma_xi(double z, double x, double T1, double is) {
  const double w0 = T1*z;
  return (w0*(-w0*x - w0 + 1)) * gpd_ipow(is, 1);
}

// Off-diagonal entry of the expected information, as gpd_expected_hessian_cpp
// computes it (NA where it does not exist).
inline double gpd_expected_sigma_xi(double sv, double x) {
  const double is = 1.0 / sv;
  const double w0 = 1.0/(2*x + 1);
  const double w1 = w0/(x + 1);
  return x > -1.0/2.0 ? (-w1) * gpd_ipow(is, 1) : NA_REAL;
}

// The derivatives of the expected information entries that
// gpd_dexpected1_cpp computes inline: d E_ss / d xi, d E_sx / d sigma and
// d E_sx / d xi, with d = 1 + 2 xi and F = 1/(d (1 + xi)).
inline double gpd_dexpected_sigma_sigma_xi(double s, double x) {
  double d = 1.0 + 2.0 * x;
  double is = 1.0 / s, is2 = is * is;
  double id = 1.0 / d, id2 = id * id;
  return 2.0 * id2 * is2;
}

inline double gpd_dexpected_sigma_xi_sigma(double s, double x) {
  double d = 1.0 + 2.0 * x, F = 1.0 / (d * (1.0 + x));
  double is = 1.0 / s, is2 = is * is;
  return F * is2;
}

inline double gpd_dexpected_sigma_xi_xi(double s, double x) {
  double d = 1.0 + 2.0 * x, F = 1.0 / (d * (1.0 + x)), c = 3.0 + 4.0 * x;
  double F1 = -c * F * F;
  double is = 1.0 / s;
  return -F1 * is;
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

// The routers of the joint registry, one per order: every component at one
// observation, on the parameter scale (sigma, xi). The Hessian and the
// expected information fill a full P x P matrix and the derivative of the
// expected information a P x P x P array, column-major, out[a + P b + P^2 c]
// being the derivative of E[l_ab] in parameter c. Outside the support the
// observed components are NaN, as in gpd_score_curv(). The expected
// information is NA where xi <= -1/2, and so is its derivative.
inline void gpd_joint_grad(double y, const double* th, double* out) {
  const double sv = th[0], x = th[1];
  const double t = gpd_t(y, sv, x);
  if (!gpd_in_support(y, t)) {
    out[0] = R_NaN;
    out[1] = R_NaN;
    return;
  }
  const double z = y / sv, u = x * z;
  const double T1 = 1.0 / t;
  const double LT = gpd_logt(u, t);
  const double v = u / t;
  const double F1 = gpd_phi(1, u, v, t, LT);
  const double is = 1.0 / sv;
  out[0] = gpd_score_sigma(z, x, T1, is);
  out[1] = gpd_score_xi(z, T1, F1);
}

inline void gpd_joint_hess(double y, const double* th, double* out) {
  const double sv = th[0], x = th[1];
  const double t = gpd_t(y, sv, x);
  if (!gpd_in_support(y, t)) {
    out[0] = R_NaN;
    out[3] = R_NaN;
    out[1] = out[2] = R_NaN;
    return;
  }
  const double z = y / sv, u = x * z;
  const double T1 = 1.0 / t;
  const double LT = gpd_logt(u, t);
  const double v = u / t;
  const double F2 = gpd_phi(2, u, v, t, LT);
  const double is = 1.0 / sv;
  out[0] = gpd_hess_sigma_sigma(z, x, T1, is);
  out[3] = gpd_hess_xi_xi(z, T1, F2);
  out[1] = out[2] = gpd_hess_sigma_xi(z, x, T1, is);
}

inline void gpd_joint_info(double y, const double* th, double* out) {
  const double sv = th[0], x = th[1];
  out[0] = gpd_expected_sigma_sigma(sv, x);
  out[3] = gpd_expected_xi_xi(x);
  out[1] = out[2] = gpd_expected_sigma_xi(sv, x);
}

inline void gpd_joint_dinfo(double y, const double* th, double* out) {
  const double sv = th[0], x = th[1];
  for (int j = 0; j < 8; ++j) out[j] = 0.0;
  if (x <= -0.5) {
    for (int j = 0; j < 8; ++j) out[j] = NA_REAL;
    return;
  }
  out[0] = gpd_dexpected_sigma_sigma_sigma(sv, x);
  out[4] = gpd_dexpected_sigma_sigma_xi(sv, x);
  out[1] = out[2] = gpd_dexpected_sigma_xi_sigma(sv, x);
  out[3] = 0.0;
  out[5] = out[6] = gpd_dexpected_sigma_xi_xi(sv, x);
  out[7] = gpd_dexpected_xi_xi_xi(x);
}

// The parts route (pt_parts.h). The buffer, with the lowest level that fills
// each entry in brackets. Outside the support every observed order is NaN, as
// in gpd_joint_grad(); the costly entries are then NaN as well.
//   [0] level
//   [1] y                           (NaN for the expected parts)
//   [2] sigma
//   [3] xi
//   [4] t = gpd_t(y, sigma, xi)     (always)
//   [5] 1 inside the support, 0 outside (always; 1 for the expected parts)
//   [6] z = y / sigma               (always)
//   [7] u = xi z                    (always)
//   [8] T1 = 1 / t                  (always)
//   [9] is = 1 / sigma              (always)
//   [10] v = u / t                  (always)
//   [11] LT = gpd_logt(u, t)        (1)
//   [12] F1 = gpd_phi(1, ...)       (1)
//   [13] F2 = gpd_phi(2, ...)       (2)
//   [14] F3 = gpd_phi(3, ...)       (3)
//   [15] F4 = gpd_phi(4, ...)       (4)
// The expected parts need no special function: the expected entries are
// arithmetic in sigma and xi, read from b[2], b[3] and b[9].

inline void gpd_parts(double y, const double* th, int level, double* b) {
  const double sv = th[0], x = th[1];
  const double t = gpd_t(y, sv, x);
  const bool in = gpd_in_support(y, t);
  const double z = y / sv, u = x * z, v = u / t;
  const double T1 = 1.0 / t, is = 1.0 / sv;
  b[0] = level; b[1] = y; b[2] = sv; b[3] = x; b[4] = t;
  b[5] = in ? 1.0 : 0.0;
  b[6] = z; b[7] = u; b[8] = T1; b[9] = is; b[10] = v;
  b[11] = (level >= 1 && in) ? gpd_logt(u, t) : R_NaN;
  b[12] = (level >= 1 && in) ? gpd_phi(1, u, v, t, b[11]) : R_NaN;
  b[13] = (level >= 2 && in) ? gpd_phi(2, u, v, t, b[11]) : R_NaN;
  b[14] = (level >= 3 && in) ? gpd_phi(3, u, v, t, b[11]) : R_NaN;
  b[15] = (level >= 4 && in) ? gpd_phi(4, u, v, t, b[11]) : R_NaN;
}

inline void gpd_parts_expected(const double* th, int level, double* b) {
  const double sv = th[0], x = th[1];
  b[0] = level; b[1] = R_NaN; b[2] = sv; b[3] = x; b[4] = R_NaN;
  b[5] = 1.0;
  b[6] = R_NaN; b[7] = R_NaN; b[8] = R_NaN; b[9] = 1.0 / sv; b[10] = R_NaN;
  for (int j = 11; j <= 15; ++j) b[j] = R_NaN;
}

inline void gpd_parts_grad(const double* b, double* out) {
  if (b[5] == 0.0) {
    out[0] = R_NaN;
    out[1] = R_NaN;
    return;
  }
  const double z = b[6], x = b[3], T1 = b[8], is = b[9], F1 = b[12];
  out[0] = gpd_score_sigma(z, x, T1, is);
  out[1] = gpd_score_xi(z, T1, F1);
}

inline void gpd_parts_hess(const double* b, double* out) {
  if (b[5] == 0.0) {
    out[0] = R_NaN;
    out[3] = R_NaN;
    out[1] = out[2] = R_NaN;
    return;
  }
  const double z = b[6], x = b[3], T1 = b[8], is = b[9], F2 = b[13];
  out[0] = gpd_hess_sigma_sigma(z, x, T1, is);
  out[3] = gpd_hess_xi_xi(z, T1, F2);
  out[1] = out[2] = gpd_hess_sigma_xi(z, x, T1, is);
}

// from the expected parts
inline void gpd_parts_info(const double* b, double* out) {
  const double sv = b[2], x = b[3];
  out[0] = gpd_expected_sigma_sigma(sv, x);
  out[3] = gpd_expected_xi_xi(x);
  out[1] = out[2] = gpd_expected_sigma_xi(sv, x);
}

// the unique components of the third derivative in the order of
// gpd_deriv3_cpp: sigma sigma sigma, sigma sigma xi, sigma xi xi, xi xi xi
inline void gpd_parts_d3(const double* b, double* out) {
  if (b[5] == 0.0) {
    for (int j = 0; j < 8; ++j) out[j] = R_NaN;
    return;
  }
  const double z = b[6], x = b[3], T1 = b[8], is = b[9], F3 = b[14];
  const double w0 = T1*z;
  const double w1 = 3*w0;
  const double w2 = std::pow(T1, 2)*std::pow(z, 2);
  const double w3 = 3*w2;
  const double w4 = std::pow(x, 2);
  const double w5 = std::pow(T1, 3);
  const double w6 = std::pow(z, 3);
  const double w7 = w5*w6;
  const double w8 = 2*w2;
  double uq[4];
  uq[0] = (2*w1*x + 2*w1 - 2*w3*w4 - 2*w3*x + 2*w4*w7 + 2*w7*std::pow(x, 3) - 2) * gpd_ipow(is, 3);
  uq[1] = (w0*(4*T1*x*z + 3*T1*z - w4*w8 - w8*x - 2)) * gpd_ipow(is, 2);
  uq[2] = (w8*(w0*x + w0 - 1)) * gpd_ipow(is, 1);
  uq[3] = (-w6*(F3*z + 2*w5));
  sym_fill(2, 3, uq, out);
}

// the unique components of the fourth derivative, in the order of
// gpd_deriv4_cpp
inline void gpd_parts_d4(const double* b, double* out) {
  if (b[5] == 0.0) {
    for (int j = 0; j < 16; ++j) out[j] = R_NaN;
    return;
  }
  const double z = b[6], x = b[3], T1 = b[8], is = b[9], F4 = b[15];
  const double w0 = T1*z;
  const double w1 = 4*w0;
  const double w2 = -w1;
  const double w3 = std::pow(T1, 2);
  const double w4 = std::pow(z, 2);
  const double w5 = w3*w4;
  const double w6 = 6*w5;
  const double w7 = std::pow(x, 3);
  const double w8 = std::pow(T1, 4);
  const double w9 = std::pow(z, 4);
  const double w10 = w8*w9;
  const double w11 = std::pow(x, 2);
  const double w12 = std::pow(T1, 3)*std::pow(z, 3);
  const double w13 = 4*w12;
  const double w14 = 6*w0;
  const double w15 = w0*x;
  const double w16 = 3*w12;
  double uq[5];
  uq[0] = (-6*w1*x + 6*w10*w7 + 6*w10*std::pow(x, 4) - 6*w11*w13 + 6*w11*w6 - 6*w13*w7 + 6*w2 + 6*w6*x + 6) * gpd_ipow(is, 4);
  uq[1] = (2*w0*(-w11*w16 + 9*w11*w3*w4 - w14 - 9*w15 - w16*w7 + 8*w3*w4*x + 3)) * gpd_ipow(is, 3);
  uq[2] = (2*w5*(3*w11*w5 - w14*x + w2 + 3*w5*x + 3)) * gpd_ipow(is, 2);
  uq[3] = (6*w12*(-w0 - w15 + 1)) * gpd_ipow(is, 1);
  uq[4] = (w9*(-F4*z + 6*w8));
  sym_fill(2, 4, uq, out);
}

} // namespace d7

#endif
