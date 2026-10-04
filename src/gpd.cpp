#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
using namespace Rcpp;

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

// [[Rcpp::export]]
NumericVector gpd_logpdf_cpp(NumericVector y, NumericVector sigma,
                             NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector out(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) { out[i] = R_NegInf; return; }
    const double z = yv / sv, u = x * z;
    const double lt = gpd_logt(u, t);
    out[i] = -std::log(sv) - lt - z * gpd_phi(0, u, u / t, t, lt);
  });
  return out;
}

// order 1 of log f in (sigma, xi)
// [[Rcpp::export]]
List gpd_gradient_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma(n);
  NumericVector o_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) {
      o_sigma[i] = R_NaN;
      o_xi[i] = R_NaN;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F1 = gpd_phi(1, u, v, t, LT);
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = T1*z;
    o_sigma[i] = (w0*x + w0 - 1) * gpd_ipow(is, 1);
    o_xi[i] = (-z*(F1*z + T1));
  });
  return List::create(Named("sigma") = o_sigma, Named("xi") = o_xi);
}

// order 2 of log f in (sigma, xi)
// [[Rcpp::export]]
List gpd_hessian_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma(n);
  NumericVector o_xi_xi(n);
  NumericVector o_sigma_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) {
      o_sigma_sigma[i] = R_NaN;
      o_xi_xi[i] = R_NaN;
      o_sigma_xi[i] = R_NaN;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F2 = gpd_phi(2, u, v, t, LT);
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = T1*z;
    const double w1 = 2*w0;
    const double w2 = std::pow(z, 2);
    const double w3 = std::pow(T1, 2);
    const double w4 = w2*w3;
    o_sigma_sigma[i] = (-w1*x - w1 + w4*std::pow(x, 2) + w4*x + 1) * gpd_ipow(is, 2);
    o_xi_xi[i] = (w2*(-F2*z + w3));
    o_sigma_xi[i] = (w0*(-w0*x - w0 + 1)) * gpd_ipow(is, 1);
  });
  return List::create(Named("sigma_sigma") = o_sigma_sigma, Named("xi_xi") = o_xi_xi, Named("sigma_xi") = o_sigma_xi);
}

// order 2, expected; NA where it does not exist
// [[Rcpp::export]]
List gpd_expected_hessian_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma(n);
  NumericVector o_xi_xi(n);
  NumericVector o_sigma_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x];
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = 1.0/(2*x + 1);
    const double w1 = w0/(x + 1);
    o_sigma_sigma[i] = x > -1.0/2.0 ? (-w0) * gpd_ipow(is, 2) : NA_REAL;
    o_xi_xi[i] = x > -1.0/2.0 ? (-2*w1) : NA_REAL;
    o_sigma_xi[i] = x > -1.0/2.0 ? (-w1) * gpd_ipow(is, 1) : NA_REAL;
  });
  return List::create(Named("sigma_sigma") = o_sigma_sigma, Named("xi_xi") = o_xi_xi, Named("sigma_xi") = o_sigma_xi);
}

// order 3 of log f in (sigma, xi)
// [[Rcpp::export]]
List gpd_deriv3_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_xi(n);
  NumericVector o_sigma_xi_xi(n);
  NumericVector o_xi_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) {
      o_sigma_sigma_sigma[i] = R_NaN;
      o_sigma_sigma_xi[i] = R_NaN;
      o_sigma_xi_xi[i] = R_NaN;
      o_xi_xi_xi[i] = R_NaN;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F3 = gpd_phi(3, u, v, t, LT);
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = T1*z;
    const double w1 = 3*w0;
    const double w2 = std::pow(T1, 2)*std::pow(z, 2);
    const double w3 = 3*w2;
    const double w4 = std::pow(x, 2);
    const double w5 = std::pow(T1, 3);
    const double w6 = std::pow(z, 3);
    const double w7 = w5*w6;
    const double w8 = 2*w2;
    o_sigma_sigma_sigma[i] = (2*w1*x + 2*w1 - 2*w3*w4 - 2*w3*x + 2*w4*w7 + 2*w7*std::pow(x, 3) - 2) * gpd_ipow(is, 3);
    o_sigma_sigma_xi[i] = (w0*(4*T1*x*z + 3*T1*z - w4*w8 - w8*x - 2)) * gpd_ipow(is, 2);
    o_sigma_xi_xi[i] = (w8*(w0*x + w0 - 1)) * gpd_ipow(is, 1);
    o_xi_xi_xi[i] = (-w6*(F3*z + 2*w5));
  });
  return List::create(Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_xi") = o_sigma_sigma_xi, Named("sigma_xi_xi") = o_sigma_xi_xi, Named("xi_xi_xi") = o_xi_xi_xi);
}

// order 3, expected; NA where it does not exist
// [[Rcpp::export]]
List gpd_deriv3_expected_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_xi(n);
  NumericVector o_sigma_xi_xi(n);
  NumericVector o_xi_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x];
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = 1.0/(3*x + 1);
    const double w1 = 4*w0;
    const double w2 = 1.0/(2*x + 1);
    const double w3 = w0*w2/(x + 1);
    o_sigma_sigma_sigma[i] = x > -1.0/3.0 ? (w1) * gpd_ipow(is, 3) : NA_REAL;
    o_sigma_sigma_xi[i] = x > -1.0/3.0 ? (w1*w2) * gpd_ipow(is, 2) : NA_REAL;
    o_sigma_xi_xi[i] = x > -1.0/3.0 ? (8*w3) * gpd_ipow(is, 1) : NA_REAL;
    o_xi_xi_xi[i] = x > -1.0/3.0 ? (24*w3) : NA_REAL;
  });
  return List::create(Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_xi") = o_sigma_sigma_xi, Named("sigma_xi_xi") = o_sigma_xi_xi, Named("xi_xi_xi") = o_xi_xi_xi);
}

// order 4 of log f in (sigma, xi)
// [[Rcpp::export]]
List gpd_deriv4_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_xi(n);
  NumericVector o_sigma_sigma_xi_xi(n);
  NumericVector o_sigma_xi_xi_xi(n);
  NumericVector o_xi_xi_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) {
      o_sigma_sigma_sigma_sigma[i] = R_NaN;
      o_sigma_sigma_sigma_xi[i] = R_NaN;
      o_sigma_sigma_xi_xi[i] = R_NaN;
      o_sigma_xi_xi_xi[i] = R_NaN;
      o_xi_xi_xi_xi[i] = R_NaN;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F4 = gpd_phi(4, u, v, t, LT);
    const double is = 1.0 / sv;
    (void) is;
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
    o_sigma_sigma_sigma_sigma[i] = (-6*w1*x + 6*w10*w7 + 6*w10*std::pow(x, 4) - 6*w11*w13 + 6*w11*w6 - 6*w13*w7 + 6*w2 + 6*w6*x + 6) * gpd_ipow(is, 4);
    o_sigma_sigma_sigma_xi[i] = (2*w0*(-w11*w16 + 9*w11*w3*w4 - w14 - 9*w15 - w16*w7 + 8*w3*w4*x + 3)) * gpd_ipow(is, 3);
    o_sigma_sigma_xi_xi[i] = (2*w5*(3*w11*w5 - w14*x + w2 + 3*w5*x + 3)) * gpd_ipow(is, 2);
    o_sigma_xi_xi_xi[i] = (6*w12*(-w0 - w15 + 1)) * gpd_ipow(is, 1);
    o_xi_xi_xi_xi[i] = (w9*(-F4*z + 6*w8));
  });
  return List::create(Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_xi") = o_sigma_sigma_sigma_xi, Named("sigma_sigma_xi_xi") = o_sigma_sigma_xi_xi, Named("sigma_xi_xi_xi") = o_sigma_xi_xi_xi, Named("xi_xi_xi_xi") = o_xi_xi_xi_xi);
}

// order 4, expected; NA where it does not exist
// [[Rcpp::export]]
List gpd_deriv4_expected_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_xi(n);
  NumericVector o_sigma_sigma_xi_xi(n);
  NumericVector o_sigma_xi_xi_xi(n);
  NumericVector o_xi_xi_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x];
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = 1.0/(4*x + 1);
    const double w1 = 18*w0;
    const double w2 = 1.0/(3*x + 1);
    const double w3 = w0*w2/(2*x + 1);
    const double w4 = w3/(x + 1);
    o_sigma_sigma_sigma_sigma[i] = x > -1.0/4.0 ? (-w1) * gpd_ipow(is, 4) : NA_REAL;
    o_sigma_sigma_sigma_xi[i] = x > -1.0/4.0 ? (-w1*w2) * gpd_ipow(is, 3) : NA_REAL;
    o_sigma_sigma_xi_xi[i] = x > -1.0/4.0 ? (-36*w3) * gpd_ipow(is, 2) : NA_REAL;
    o_sigma_xi_xi_xi[i] = x > -1.0/4.0 ? (-108*w4) * gpd_ipow(is, 1) : NA_REAL;
    o_xi_xi_xi_xi[i] = x > -1.0/4.0 ? (-432*w4) : NA_REAL;
  });
  return List::create(Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_xi") = o_sigma_sigma_sigma_xi, Named("sigma_sigma_xi_xi") = o_sigma_sigma_xi_xi, Named("sigma_xi_xi_xi") = o_sigma_xi_xi_xi, Named("xi_xi_xi_xi") = o_xi_xi_xi_xi);
}

// order 5 of log f in (sigma, xi)
// [[Rcpp::export]]
List gpd_deriv5_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_sigma_xi(n);
  NumericVector o_sigma_sigma_sigma_xi_xi(n);
  NumericVector o_sigma_sigma_xi_xi_xi(n);
  NumericVector o_sigma_xi_xi_xi_xi(n);
  NumericVector o_xi_xi_xi_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) {
      o_sigma_sigma_sigma_sigma_sigma[i] = R_NaN;
      o_sigma_sigma_sigma_sigma_xi[i] = R_NaN;
      o_sigma_sigma_sigma_xi_xi[i] = R_NaN;
      o_sigma_sigma_xi_xi_xi[i] = R_NaN;
      o_sigma_xi_xi_xi_xi[i] = R_NaN;
      o_xi_xi_xi_xi_xi[i] = R_NaN;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F5 = gpd_phi(5, u, v, t, LT);
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = T1*z;
    const double w1 = 5*w0;
    const double w2 = std::pow(T1, 2)*std::pow(z, 2);
    const double w3 = 10*w2;
    const double w4 = std::pow(x, 4);
    const double w5 = std::pow(T1, 5);
    const double w6 = std::pow(z, 5);
    const double w7 = w5*w6;
    const double w8 = std::pow(x, 2);
    const double w9 = std::pow(T1, 3);
    const double w10 = std::pow(z, 3);
    const double w11 = w10*w9;
    const double w12 = 10*w11;
    const double w13 = std::pow(x, 3);
    const double w14 = std::pow(T1, 4)*std::pow(z, 4);
    const double w15 = 5*w14;
    const double w16 = w2*x;
    const double w17 = w2*w8;
    const double w18 = 4*w14;
    const double w19 = w0*x;
    o_sigma_sigma_sigma_sigma_sigma[i] = (24*w1*x + 24*w1 + 24*w12*w13 + 24*w12*w8 - 24*w13*w15 - 24*w15*w4 - 24*w3*w8 - 24*w3*x + 24*w4*w7 + 24*w7*std::pow(x, 5) - 24) * gpd_ipow(is, 5);
    o_sigma_sigma_sigma_sigma_xi[i] = (6*w0*(16*T1*x*z + 10*T1*z + 16*w10*w13*w9 + 15*w10*w8*w9 - w13*w18 - 20*w16 - 24*w17 - w18*w4 - 4)) * gpd_ipow(is, 4);
    o_sigma_sigma_sigma_xi_xi[i] = (4*w2*(10*w0 + 6*w11*w13 + 6*w11*w8 - 15*w16 - 18*w17 + 18*w19 - 6)) * gpd_ipow(is, 3);
    o_sigma_sigma_xi_xi_xi[i] = (6*w11*(8*T1*x*z + 5*T1*z - 4*w16 - 4*w17 - 4)) * gpd_ipow(is, 2);
    o_sigma_xi_xi_xi_xi[i] = (24*w14*(w0 + w19 - 1)) * gpd_ipow(is, 1);
    o_xi_xi_xi_xi_xi[i] = (-w6*(F5*z + 24*w5));
  });
  return List::create(Named("sigma_sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_sigma_xi") = o_sigma_sigma_sigma_sigma_xi, Named("sigma_sigma_sigma_xi_xi") = o_sigma_sigma_sigma_xi_xi, Named("sigma_sigma_xi_xi_xi") = o_sigma_sigma_xi_xi_xi, Named("sigma_xi_xi_xi_xi") = o_sigma_xi_xi_xi_xi, Named("xi_xi_xi_xi_xi") = o_xi_xi_xi_xi_xi);
}

// order 1 of log f in y
// [[Rcpp::export]]
List gpd_dy1_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_y(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) {
      o_y[i] = R_NaN;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double is = 1.0 / sv;
    (void) is;
    o_y[i] = (-T1*(x + 1)) * gpd_ipow(is, 1);
  });
  return List::create(Named("y") = o_y);
}

// order 2 of log f in y
// [[Rcpp::export]]
List gpd_dy2_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_y(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) {
      o_y[i] = R_NaN;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double is = 1.0 / sv;
    (void) is;
    o_y[i] = (std::pow(T1, 2)*x*(x + 1)) * gpd_ipow(is, 2);
  });
  return List::create(Named("y") = o_y);
}

// order 3 of log f in y
// [[Rcpp::export]]
List gpd_dy3_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_y(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) {
      o_y[i] = R_NaN;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double is = 1.0 / sv;
    (void) is;
    o_y[i] = (-2*std::pow(T1, 3)*std::pow(x, 2)*(x + 1)) * gpd_ipow(is, 3);
  });
  return List::create(Named("y") = o_y);
}

// order 4 of log f in y
// [[Rcpp::export]]
List gpd_dy4_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_y(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) {
      o_y[i] = R_NaN;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double is = 1.0 / sv;
    (void) is;
    o_y[i] = (6*std::pow(T1, 4)*std::pow(x, 3)*(x + 1)) * gpd_ipow(is, 4);
  });
  return List::create(Named("y") = o_y);
}

// 1 in y and 1 in (sigma, xi)
// [[Rcpp::export]]
List gpd_cross_y_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma(n);
  NumericVector o_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) {
      o_sigma[i] = R_NaN;
      o_xi[i] = R_NaN;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double is = 1.0 / sv;
    (void) is;
    if (u > 0.0) {
      // R = u/t and 1/t, both in (0, 1), in the basis R^a t^-b
      const double R = u / t;
      (void) R;
      const double w0 = std::pow(T1, 2);
      o_sigma[i] = (w0*(x + 1)) * gpd_ipow(is, 2);
      o_xi[i] = (w0*(z - 1)) * gpd_ipow(is, 1);
    } else {
      const double w0 = T1*z;
      const double w1 = w0*x - 1;
      o_sigma[i] = (T1*(-w0*std::pow(x, 2) - w1 + x)) * gpd_ipow(is, 2);
      o_xi[i] = (T1*(w0 + w1)) * gpd_ipow(is, 1);
    }
  });
  return List::create(Named("sigma") = o_sigma, Named("xi") = o_xi);
}

// 2 in y and 1 in (sigma, xi)
// [[Rcpp::export]]
List gpd_cross2_y_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma(n);
  NumericVector o_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) {
      o_sigma[i] = R_NaN;
      o_xi[i] = R_NaN;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double is = 1.0 / sv;
    (void) is;
    if (u > 0.0) {
      // R = u/t and 1/t, both in (0, 1), in the basis R^a t^-b
      const double R = u / t;
      (void) R;
      const double w0 = 2*x;
      o_sigma[i] = (-std::pow(T1, 3)*w0*(x + 1)) * gpd_ipow(is, 3);
      o_xi[i] = (std::pow(T1, 2)*(-R + T1*(w0 + 1))) * gpd_ipow(is, 2);
    } else {
      const double w0 = std::pow(T1, 2);
      const double w1 = T1*z;
      const double w2 = w1*std::pow(x, 2);
      const double w3 = 2*x;
      o_sigma[i] = (w0*w3*(w1*x + w2 - x - 1)) * gpd_ipow(is, 3);
      o_xi[i] = (w0*(-w1*w3 - 2*w2 + 2*x + 1)) * gpd_ipow(is, 2);
    }
  });
  return List::create(Named("sigma") = o_sigma, Named("xi") = o_xi);
}

// 1 in y and 2 in (sigma, xi)
// [[Rcpp::export]]
List gpd_grad_y_hess_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma(n);
  NumericVector o_xi_xi(n);
  NumericVector o_sigma_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) {
      o_sigma_sigma[i] = R_NaN;
      o_xi_xi[i] = R_NaN;
      o_sigma_xi[i] = R_NaN;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double is = 1.0 / sv;
    (void) is;
    if (u > 0.0) {
      // R = u/t and 1/t, both in (0, 1), in the basis R^a t^-b
      const double R = u / t;
      (void) R;
      const double w0 = std::pow(T1, 3);
      const double w1 = 2*z;
      o_sigma_sigma[i] = (-2*w0*(x + 1)) * gpd_ipow(is, 3);
      o_xi_xi[i] = (-w0*w1*(z - 1)) * gpd_ipow(is, 1);
      o_sigma_xi[i] = (-std::pow(T1, 2)*(R + T1*(w1 - 1))) * gpd_ipow(is, 2);
    } else {
      const double w0 = std::pow(x, 2);
      const double w1 = std::pow(T1, 2);
      const double w2 = w1*std::pow(z, 2);
      const double w3 = w0*w2;
      const double w4 = T1*z;
      const double w5 = w4*x;
      o_sigma_sigma[i] = (2*T1*(2*T1*w0*z + 2*T1*x*z - w2*std::pow(x, 3) - w3 - x - 1)) * gpd_ipow(is, 3);
      o_xi_xi[i] = (2*w1*z*(-w4 - w5 + 1)) * gpd_ipow(is, 1);
      o_sigma_xi[i] = (T1*(2*w2*x + 2*w3 - 2*w4 - 3*w5 + 1)) * gpd_ipow(is, 2);
    }
  });
  return List::create(Named("sigma_sigma") = o_sigma_sigma, Named("xi_xi") = o_xi_xi, Named("sigma_xi") = o_sigma_xi);
}

// 2 in y and 2 in (sigma, xi)
// [[Rcpp::export]]
List gpd_hess_y_hess_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma(n);
  NumericVector o_xi_xi(n);
  NumericVector o_sigma_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv < 0.0 || !(t > 0.0)) {
      o_sigma_sigma[i] = R_NaN;
      o_xi_xi[i] = R_NaN;
      o_sigma_xi[i] = R_NaN;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double is = 1.0 / sv;
    (void) is;
    if (u > 0.0) {
      // R = u/t and 1/t, both in (0, 1), in the basis R^a t^-b
      const double R = u / t;
      (void) R;
      const double w0 = 2*std::pow(T1, 3);
      o_sigma_sigma[i] = (6*std::pow(T1, 4)*x*(x + 1)) * gpd_ipow(is, 4);
      o_xi_xi[i] = (w0*(R*(z - 2) - T1*(2*z - 1))) * gpd_ipow(is, 2);
      o_sigma_xi[i] = (w0*(R*(x + 2) - T1*(2*x + 1))) * gpd_ipow(is, 3);
    } else {
      const double w0 = std::pow(T1, 2);
      const double w1 = 2*x;
      const double w2 = T1*z;
      const double w3 = std::pow(x, 2);
      const double w4 = 2*w2;
      const double w5 = w0*std::pow(z, 2);
      const double w6 = w3*w5;
      const double w7 = w5*std::pow(x, 3);
      const double w8 = -4*w2*x + 3*w6 + 1;
      const double w9 = 2*w0;
      o_sigma_sigma[i] = (6*w0*x*(-w1*w2 - w3*w4 + w6 + w7 + x + 1)) * gpd_ipow(is, 4);
      o_xi_xi[i] = (w9*(-w4 + 3*w5*x + w8)) * gpd_ipow(is, 2);
      o_sigma_xi[i] = (w9*(5*T1*w3*z - w1 - 3*w7 - w8)) * gpd_ipow(is, 3);
    }
  });
  return List::create(Named("sigma_sigma") = o_sigma_sigma, Named("xi_xi") = o_xi_xi, Named("sigma_xi") = o_sigma_xi);
}

// order 1 of log F in (sigma, xi)
// [[Rcpp::export]]
List gpd_grad_cdf_lower_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma(n);
  NumericVector o_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv <= 0.0 || !(t > 0.0)) {
      // below the support F = 0, past its end F = 1: constant
      o_sigma[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      o_xi[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F1 = gpd_phi(1, u, v, t, LT);
    const double W = z * gpd_phi(0, u, u / t, t, LT);
    const double q = 1.0 / std::expm1(W);
    const double is = 1.0 / sv;
    (void) is;
    o_sigma[i] = (-T1*q*z) * gpd_ipow(is, 1);
    o_xi[i] = (F1*q*std::pow(z, 2));
  });
  return List::create(Named("sigma") = o_sigma, Named("xi") = o_xi);
}

// order 1 of log S in (sigma, xi)
// [[Rcpp::export]]
List gpd_grad_cdf_upper_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma(n);
  NumericVector o_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv <= 0.0 || !(t > 0.0)) {
      // below the support F = 0, past its end F = 1: constant
      o_sigma[i] = 0.0;
      o_xi[i] = 0.0;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F1 = gpd_phi(1, u, v, t, LT);
    const double is = 1.0 / sv;
    (void) is;
    o_sigma[i] = (T1*z) * gpd_ipow(is, 1);
    o_xi[i] = (-F1*std::pow(z, 2));
  });
  return List::create(Named("sigma") = o_sigma, Named("xi") = o_xi);
}

// order 1 of S = 1 - F in (sigma, xi)
// [[Rcpp::export]]
List gpd_grad_surv_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma(n);
  NumericVector o_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv <= 0.0 || !(t > 0.0)) {
      // below the support F = 0, past its end F = 1: constant
      o_sigma[i] = 0.0;
      o_xi[i] = 0.0;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F1 = gpd_phi(1, u, v, t, LT);
    const double W = z * gpd_phi(0, u, u / t, t, LT);
    const double Sv = std::exp(-W);
    const double is = 1.0 / sv;
    (void) is;
    o_sigma[i] = (Sv*T1*z) * gpd_ipow(is, 1);
    o_xi[i] = (-F1*Sv*std::pow(z, 2));
  });
  return List::create(Named("sigma") = o_sigma, Named("xi") = o_xi);
}

// order 2 of log F in (sigma, xi)
// [[Rcpp::export]]
List gpd_hess_cdf_lower_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma(n);
  NumericVector o_sigma_xi(n);
  NumericVector o_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv <= 0.0 || !(t > 0.0)) {
      // below the support F = 0, past its end F = 1: constant
      o_sigma_sigma[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      o_sigma_xi[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      o_xi_xi[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F1 = gpd_phi(1, u, v, t, LT);
    const double F2 = gpd_phi(2, u, v, t, LT);
    const double W = z * gpd_phi(0, u, u / t, t, LT);
    const double q = 1.0 / std::expm1(W);
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = T1*z;
    const double w1 = q*w0;
    const double w2 = F1*z;
    const double w3 = std::pow(F1, 2)*z;
    o_sigma_sigma[i] = (w1*(-w0*x - w0 - w1 + 2)) * gpd_ipow(is, 2);
    o_sigma_xi[i] = (T1*q*std::pow(z, 2)*(T1 + q*w2 + w2)) * gpd_ipow(is, 1);
    o_xi_xi[i] = (q*std::pow(z, 3)*(F2 - q*w3 - w3));
  });
  return List::create(Named("sigma_sigma") = o_sigma_sigma, Named("sigma_xi") = o_sigma_xi, Named("xi_xi") = o_xi_xi);
}

// order 2 of log S in (sigma, xi)
// [[Rcpp::export]]
List gpd_hess_cdf_upper_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma(n);
  NumericVector o_sigma_xi(n);
  NumericVector o_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv <= 0.0 || !(t > 0.0)) {
      // below the support F = 0, past its end F = 1: constant
      o_sigma_sigma[i] = 0.0;
      o_sigma_xi[i] = 0.0;
      o_xi_xi[i] = 0.0;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F2 = gpd_phi(2, u, v, t, LT);
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = T1*z;
    o_sigma_sigma[i] = (w0*(w0*x - 2)) * gpd_ipow(is, 2);
    o_sigma_xi[i] = (-std::pow(T1, 2)*std::pow(z, 2)) * gpd_ipow(is, 1);
    o_xi_xi[i] = (-F2*std::pow(z, 3));
  });
  return List::create(Named("sigma_sigma") = o_sigma_sigma, Named("sigma_xi") = o_sigma_xi, Named("xi_xi") = o_xi_xi);
}

// order 2 of S = 1 - F in (sigma, xi)
// [[Rcpp::export]]
List gpd_hess_surv_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma(n);
  NumericVector o_sigma_xi(n);
  NumericVector o_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv <= 0.0 || !(t > 0.0)) {
      // below the support F = 0, past its end F = 1: constant
      o_sigma_sigma[i] = 0.0;
      o_sigma_xi[i] = 0.0;
      o_xi_xi[i] = 0.0;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F1 = gpd_phi(1, u, v, t, LT);
    const double F2 = gpd_phi(2, u, v, t, LT);
    const double W = z * gpd_phi(0, u, u / t, t, LT);
    const double Sv = std::exp(-W);
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = T1*z;
    o_sigma_sigma[i] = (Sv*w0*(w0*x + w0 - 2)) * gpd_ipow(is, 2);
    o_sigma_xi[i] = (-Sv*T1*std::pow(z, 2)*(F1*z + T1)) * gpd_ipow(is, 1);
    o_xi_xi[i] = (Sv*std::pow(z, 3)*(std::pow(F1, 2)*z - F2));
  });
  return List::create(Named("sigma_sigma") = o_sigma_sigma, Named("sigma_xi") = o_sigma_xi, Named("xi_xi") = o_xi_xi);
}

// order 3 of log F in (sigma, xi)
// [[Rcpp::export]]
List gpd_deriv3_cdf_lower_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_xi(n);
  NumericVector o_sigma_xi_xi(n);
  NumericVector o_xi_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv <= 0.0 || !(t > 0.0)) {
      // below the support F = 0, past its end F = 1: constant
      o_sigma_sigma_sigma[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      o_sigma_sigma_xi[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      o_sigma_xi_xi[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      o_xi_xi_xi[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F1 = gpd_phi(1, u, v, t, LT);
    const double F2 = gpd_phi(2, u, v, t, LT);
    const double F3 = gpd_phi(3, u, v, t, LT);
    const double W = z * gpd_phi(0, u, u / t, t, LT);
    const double q = 1.0 / std::expm1(W);
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = std::pow(T1, 2);
    const double w1 = std::pow(z, 2);
    const double w2 = w0*w1;
    const double w3 = 3*w2;
    const double w4 = q*w3;
    const double w5 = std::pow(q, 2);
    const double w6 = 2*w2;
    const double w7 = T1*z;
    const double w8 = q*w7;
    const double w9 = 3*T1;
    const double w10 = 2*F1;
    const double w11 = w10*z;
    const double w12 = 2*w0;
    const double w13 = w12*z;
    const double w14 = T1*w1;
    const double w15 = F1*w14;
    const double w16 = w15*x;
    const double w17 = std::pow(F1, 2)*w1;
    const double w18 = 3*q;
    const double w19 = 2*w5;
    const double w20 = F2*z;
    const double w21 = 3*F1;
    const double w22 = std::pow(F1, 3)*w1;
    o_sigma_sigma_sigma[i] = (w8*(6*T1*q*z + 6*T1*x*z + 6*T1*z - w2 - w3*x - w4*x - w4 - w5*w6 - w6*std::pow(x, 2) - 6)) * gpd_ipow(is, 3);
    o_sigma_sigma_xi[i] = (q*w14*(F1*q*w1*w9 - q*w11 + q*w13 + q*w16 + w10*w14*w5 - w11 + w13*x + w13 + w15 + w16 - w9)) * gpd_ipow(is, 2);
    o_sigma_xi_xi[i] = (T1*q*std::pow(z, 3)*(F2*q*z + F2*z - w10*w7 - w10*w8 - w12 - w17*w18 - w17*w19 - w17)) * gpd_ipow(is, 1);
    o_xi_xi_xi[i] = (q*std::pow(z, 4)*(F3 - q*w20*w21 + w18*w22 + w19*w22 - w20*w21 + w22));
  });
  return List::create(Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_xi") = o_sigma_sigma_xi, Named("sigma_xi_xi") = o_sigma_xi_xi, Named("xi_xi_xi") = o_xi_xi_xi);
}

// order 3 of log S in (sigma, xi)
// [[Rcpp::export]]
List gpd_deriv3_cdf_upper_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_xi(n);
  NumericVector o_sigma_xi_xi(n);
  NumericVector o_xi_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv <= 0.0 || !(t > 0.0)) {
      // below the support F = 0, past its end F = 1: constant
      o_sigma_sigma_sigma[i] = 0.0;
      o_sigma_sigma_xi[i] = 0.0;
      o_sigma_xi_xi[i] = 0.0;
      o_xi_xi_xi[i] = 0.0;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F3 = gpd_phi(3, u, v, t, LT);
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = T1*z;
    const double w1 = w0*x;
    const double w2 = std::pow(T1, 2)*std::pow(z, 2);
    o_sigma_sigma_sigma[i] = (2*w0*(-3*w1 + w2*std::pow(x, 2) + 3)) * gpd_ipow(is, 3);
    o_sigma_sigma_xi[i] = (w2*(3 - 2*w1)) * gpd_ipow(is, 2);
    o_sigma_xi_xi[i] = (2*std::pow(T1, 3)*std::pow(z, 3)) * gpd_ipow(is, 1);
    o_xi_xi_xi[i] = (-F3*std::pow(z, 4));
  });
  return List::create(Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_xi") = o_sigma_sigma_xi, Named("sigma_xi_xi") = o_sigma_xi_xi, Named("xi_xi_xi") = o_xi_xi_xi);
}

// order 3 of S = 1 - F in (sigma, xi)
// [[Rcpp::export]]
List gpd_deriv3_surv_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_xi(n);
  NumericVector o_sigma_xi_xi(n);
  NumericVector o_xi_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv <= 0.0 || !(t > 0.0)) {
      // below the support F = 0, past its end F = 1: constant
      o_sigma_sigma_sigma[i] = 0.0;
      o_sigma_sigma_xi[i] = 0.0;
      o_sigma_xi_xi[i] = 0.0;
      o_xi_xi_xi[i] = 0.0;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F1 = gpd_phi(1, u, v, t, LT);
    const double F2 = gpd_phi(2, u, v, t, LT);
    const double F3 = gpd_phi(3, u, v, t, LT);
    const double W = z * gpd_phi(0, u, u / t, t, LT);
    const double Sv = std::exp(-W);
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = T1*z;
    const double w1 = 6*w0;
    const double w2 = std::pow(T1, 2);
    const double w3 = std::pow(z, 2);
    const double w4 = w2*w3;
    const double w5 = 2*w2;
    const double w6 = w5*z;
    const double w7 = T1*w3;
    const double w8 = F1*w7;
    o_sigma_sigma_sigma[i] = (Sv*w0*(-w1*x - w1 + 2*w4*std::pow(x, 2) + 3*w4*x + w4 + 6)) * gpd_ipow(is, 3);
    o_sigma_sigma_xi[i] = (Sv*w7*(2*F1*z + 3*T1 - w6*x - w6 - w8*x - w8)) * gpd_ipow(is, 2);
    o_sigma_xi_xi[i] = (Sv*T1*std::pow(z, 3)*(std::pow(F1, 2)*w3 + 2*F1*w0 - F2*z + w5)) * gpd_ipow(is, 1);
    o_xi_xi_xi[i] = (Sv*std::pow(z, 4)*(-std::pow(F1, 3)*w3 + 3*F1*F2*z - F3));
  });
  return List::create(Named("sigma_sigma_sigma") = o_sigma_sigma_sigma, Named("sigma_sigma_xi") = o_sigma_sigma_xi, Named("sigma_xi_xi") = o_sigma_xi_xi, Named("xi_xi_xi") = o_xi_xi_xi);
}

// order 4 of log F in (sigma, xi)
// [[Rcpp::export]]
List gpd_deriv4_cdf_lower_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_xi(n);
  NumericVector o_sigma_sigma_xi_xi(n);
  NumericVector o_sigma_xi_xi_xi(n);
  NumericVector o_xi_xi_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv <= 0.0 || !(t > 0.0)) {
      // below the support F = 0, past its end F = 1: constant
      o_sigma_sigma_sigma_sigma[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      o_sigma_sigma_sigma_xi[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      o_sigma_sigma_xi_xi[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      o_sigma_xi_xi_xi[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      o_xi_xi_xi_xi[i] = (yv <= 0.0 && !(t <= 0.0)) ? R_NaN : 0.0;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F1 = gpd_phi(1, u, v, t, LT);
    const double F2 = gpd_phi(2, u, v, t, LT);
    const double F3 = gpd_phi(3, u, v, t, LT);
    const double F4 = gpd_phi(4, u, v, t, LT);
    const double W = z * gpd_phi(0, u, u / t, t, LT);
    const double q = 1.0 / std::expm1(W);
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = T1*z;
    const double w1 = 36*w0;
    const double w2 = std::pow(T1, 3);
    const double w3 = std::pow(z, 3);
    const double w4 = w2*w3;
    const double w5 = std::pow(T1, 2);
    const double w6 = std::pow(z, 2);
    const double w7 = 7*q;
    const double w8 = 6*w4;
    const double w9 = std::pow(q, 2);
    const double w10 = std::pow(x, 2);
    const double w11 = 12*w9;
    const double w12 = w11*w4;
    const double w13 = std::pow(q, 3);
    const double w14 = 11*w10*w4;
    const double w15 = q*w0;
    const double w16 = 12*T1;
    const double w17 = 6*F1;
    const double w18 = w17*z;
    const double w19 = w5*z;
    const double w20 = 15*w19;
    const double w21 = T1*w6;
    const double w22 = w17*w21;
    const double w23 = w2*w6;
    const double w24 = w3*w5;
    const double w25 = F1*w24;
    const double w26 = q*w21;
    const double w27 = w22*x;
    const double w28 = 9*w23;
    const double w29 = q*w28;
    const double w30 = F1*w6;
    const double w31 = w16*w9;
    const double w32 = w25*x;
    const double w33 = 6*w2;
    const double w34 = w33*w6;
    const double w35 = 9*q;
    const double w36 = w17*w24;
    const double w37 = 2*w10*w25;
    const double w38 = 2*F2*z;
    const double w39 = w33*z;
    const double w40 = std::pow(F1, 2);
    const double w41 = T1*w3;
    const double w42 = w40*w41;
    const double w43 = w5*w6;
    const double w44 = 4*F1;
    const double w45 = w43*w44;
    const double w46 = w42*x;
    const double w47 = w45*x;
    const double w48 = 6*w13;
    const double w49 = F3*z;
    const double w50 = q*w49;
    const double w51 = 3*F2;
    const double w52 = std::pow(F1, 3)*w3;
    const double w53 = w17*w19;
    const double w54 = w40*w6;
    const double w55 = T1*w54;
    const double w56 = 3*std::pow(F2, 2)*z;
    const double w57 = std::pow(F1, 4)*w3;
    o_sigma_sigma_sigma_sigma[i] = (w15*(-q*w1 - q*w14 - 18*q*w4*x + 36*q*w5*w6*x + 36*q*w5*w6 - w1*x - w1 + 24*w10*w5*w6 - w12*x - w12 - w13*w8 - w14 - w4*w7 - w4 + 24*w5*w6*w9 + 36*w5*w6*x + 12*w5*w6 - w8*std::pow(x, 3) - w8*x + 24)) * gpd_ipow(is, 4);
    o_sigma_sigma_sigma_xi[i] = (w26*(-18*F1*w26 + q*w18 - q*w20 - q*w27 + q*w37 + w10*w34 + w11*w25 + w13*w36 + w16 + w18 - 16*w19*x - w20 - w22 + 3*w23 + w25*w7 + w25 - w27 + w28*x + w29*x + w29 - w30*w31 + w32*w35 + 3*w32 + w34*w9 + w36*w9*x + w37)) * gpd_ipow(is, 3);
    o_sigma_sigma_xi_xi[i] = (q*w41*(6*F1*T1*q*z + 6*F1*T1*z - 12*F1*q*w43 + F2*T1*q*w6*x + 3*F2*T1*q*w6 + 2*F2*T1*w6*w9 + F2*T1*w6*x + F2*T1*w6 - q*w38 - q*w39 + 6*q*w40*w6 - 3*q*w46 - q*w47 - w3*w31*w40 - 8*w30*w5*w9 - w38 - w39*x - w39 + 4*w40*w6*w9 + 2*w40*w6 - w42*w48 - w42*w7 - w42 - w45 - 2*w46*w9 - w46 - w47 + 8*w5)) * gpd_ipow(is, 2);
    o_sigma_xi_xi_xi[i] = (T1*q*std::pow(z, 4)*(-F2*w17*w6*w9 - F2*w30*w35 + 6*T1*w54*w9 + q*w53 - w0*w51 + w11*w52 - w15*w51 - w30*w51 + w33 + w35*w55 + w48*w52 + w49 + w50 + w52*w7 + w52 + w53 + 3*w55)) * gpd_ipow(is, 1);
    o_xi_xi_xi_xi[i] = (q*std::pow(z, 5)*(18*F2*q*w40*w6 + 12*F2*w40*w6*w9 + 6*F2*w40*w6 + F4 - q*w56 - w11*w57 - w44*w49 - w44*w50 - w48*w57 - w56 - w57*w7 - w57));
  });
  return List::create(Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_xi") = o_sigma_sigma_sigma_xi, Named("sigma_sigma_xi_xi") = o_sigma_sigma_xi_xi, Named("sigma_xi_xi_xi") = o_sigma_xi_xi_xi, Named("xi_xi_xi_xi") = o_xi_xi_xi_xi);
}

// order 4 of log S in (sigma, xi)
// [[Rcpp::export]]
List gpd_deriv4_cdf_upper_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_xi(n);
  NumericVector o_sigma_sigma_xi_xi(n);
  NumericVector o_sigma_xi_xi_xi(n);
  NumericVector o_xi_xi_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv <= 0.0 || !(t > 0.0)) {
      // below the support F = 0, past its end F = 1: constant
      o_sigma_sigma_sigma_sigma[i] = 0.0;
      o_sigma_sigma_sigma_xi[i] = 0.0;
      o_sigma_sigma_xi_xi[i] = 0.0;
      o_sigma_xi_xi_xi[i] = 0.0;
      o_xi_xi_xi_xi[i] = 0.0;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F4 = gpd_phi(4, u, v, t, LT);
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = 6*T1*z;
    const double w1 = std::pow(T1, 3)*std::pow(z, 3);
    const double w2 = std::pow(T1, 2)*std::pow(z, 2);
    const double w3 = w2*std::pow(x, 2);
    o_sigma_sigma_sigma_sigma[i] = (w0*(w0*x + w1*std::pow(x, 3) - 4*w3 - 4)) * gpd_ipow(is, 4);
    o_sigma_sigma_sigma_xi[i] = (2*w2*(8*T1*x*z - 3*w3 - 6)) * gpd_ipow(is, 3);
    o_sigma_sigma_xi_xi[i] = (2*w1*(3*T1*x*z - 4)) * gpd_ipow(is, 2);
    o_sigma_xi_xi_xi[i] = (-6*std::pow(T1, 4)*std::pow(z, 4)) * gpd_ipow(is, 1);
    o_xi_xi_xi_xi[i] = (-F4*std::pow(z, 5));
  });
  return List::create(Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_xi") = o_sigma_sigma_sigma_xi, Named("sigma_sigma_xi_xi") = o_sigma_sigma_xi_xi, Named("sigma_xi_xi_xi") = o_sigma_xi_xi_xi, Named("xi_xi_xi_xi") = o_xi_xi_xi_xi);
}

// order 4 of S = 1 - F in (sigma, xi)
// [[Rcpp::export]]
List gpd_deriv4_surv_cpp(NumericVector y, NumericVector sigma, NumericVector xi, int threads = 1) {
  const int n = y.size();
  const int n_s = sigma.size(), n_x = xi.size();
  NumericVector o_sigma_sigma_sigma_sigma(n);
  NumericVector o_sigma_sigma_sigma_xi(n);
  NumericVector o_sigma_sigma_xi_xi(n);
  NumericVector o_sigma_xi_xi_xi(n);
  NumericVector o_xi_xi_xi_xi(n);
  d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
    const double sv = sigma[i % n_s], x = xi[i % n_x], yv = y[i];
    const double t = gpd_t(yv, sv, x);
    if (yv <= 0.0 || !(t > 0.0)) {
      // below the support F = 0, past its end F = 1: constant
      o_sigma_sigma_sigma_sigma[i] = 0.0;
      o_sigma_sigma_sigma_xi[i] = 0.0;
      o_sigma_sigma_xi_xi[i] = 0.0;
      o_sigma_xi_xi_xi[i] = 0.0;
      o_xi_xi_xi_xi[i] = 0.0;
      return;
    }
    const double z = yv / sv, u = x * z;
    (void) z; (void) u;
    const double T1 = 1.0 / t;
    const double LT = gpd_logt(u, t);
    const double v = u / t;
    const double F1 = gpd_phi(1, u, v, t, LT);
    const double F2 = gpd_phi(2, u, v, t, LT);
    const double F3 = gpd_phi(3, u, v, t, LT);
    const double F4 = gpd_phi(4, u, v, t, LT);
    const double W = z * gpd_phi(0, u, u / t, t, LT);
    const double Sv = std::exp(-W);
    const double is = 1.0 / sv;
    (void) is;
    const double w0 = T1*z;
    const double w1 = 36*w0;
    const double w2 = std::pow(T1, 3);
    const double w3 = std::pow(z, 3);
    const double w4 = w2*w3;
    const double w5 = std::pow(T1, 2);
    const double w6 = std::pow(z, 2);
    const double w7 = w5*w6;
    const double w8 = w7*x;
    const double w9 = 6*w4;
    const double w10 = std::pow(x, 2);
    const double w11 = 6*F1;
    const double w12 = w11*z;
    const double w13 = w2*w6;
    const double w14 = F1*w3*w5;
    const double w15 = 6*w2;
    const double w16 = T1*w6;
    const double w17 = w15*z;
    const double w18 = F2*w16;
    const double w19 = std::pow(F1, 2);
    const double w20 = w19*w6;
    const double w21 = T1*w3;
    const double w22 = w19*w21;
    const double w23 = 4*F1;
    const double w24 = F3*z;
    o_sigma_sigma_sigma_sigma[i] = (Sv*w0*(w1*x + w1 + 11*w10*w4 - 24*w10*w7 + w4 - 12*w7 - 36*w8 + w9*std::pow(x, 3) + w9*x - 24)) * gpd_ipow(is, 4);
    o_sigma_sigma_sigma_xi[i] = (Sv*w16*(6*F1*T1*w6*x + 6*F1*T1*w6 - 12*T1 - 2*w10*w14 - w10*w15*w6 - w12 - 9*w13*x - 3*w13 - 3*w14*x - w14 + 16*w5*x*z + 15*w5*z)) * gpd_ipow(is, 3);
    o_sigma_sigma_xi_xi[i] = (Sv*w21*(2*F2*z - w0*w11 + w17*x + w17 - w18*x - w18 - 2*w20 + w22*x + w22 + w23*w7 + w23*w8 - 8*w5)) * gpd_ipow(is, 2);
    o_sigma_xi_xi_xi[i] = (Sv*T1*std::pow(z, 4)*(-std::pow(F1, 3)*w3 + 3*F1*F2*w6 + 3*F2*T1*z - 3*T1*w20 - w12*w5 - w15 - w24)) * gpd_ipow(is, 1);
    o_xi_xi_xi_xi[i] = (Sv*std::pow(z, 5)*(std::pow(F1, 4)*w3 + 3*std::pow(F2, 2)*z - 6*F2*w20 - F4 + w23*w24));
  });
  return List::create(Named("sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma, Named("sigma_sigma_sigma_xi") = o_sigma_sigma_sigma_xi, Named("sigma_sigma_xi_xi") = o_sigma_sigma_xi_xi, Named("sigma_xi_xi_xi") = o_sigma_xi_xi_xi, Named("xi_xi_xi_xi") = o_xi_xi_xi_xi);
}
