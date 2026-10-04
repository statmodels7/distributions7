#include <Rcpp.h>
#include <cmath>
#include <cfloat>
#include "d7_par.h"
using namespace Rcpp;

// The generalized gamma in Stacy's form, scale a, shape d and power p. With
// k = d/p and U = p log(y/a), the log-density is
//   log p - log y + k U - e^U - lgamma(k),
// and e^U = (y/a)^p is Gamma(k, 1). Every component below is a closed form
// derived offline with sympy (stabilita/gen_gengamma1.py), one exported
// function per order and per surface: a polynomial in U, e^U and a/y whose
// coefficients depend on (d, p) alone, times a power of 1/a. The polygamma
// functions are taken at k + 1, where the poles in 1/k that the derivatives
// combine have cancelled symbolically. Checked numerically in the tests and
// against exact values in stabilita/gengamma1_measure.R.

// log(y/a), also where y/a overflows or underflows
static inline double gg1_logratio(double y, double a) {
  const double r = y / a;
  if (r > DBL_MIN && r < DBL_MAX) return std::log(r);
  return std::log(y) - std::log(a);
}

// a^-m (a/y)^l from ia = 1/a and Yi = a/y, and through their logarithms
// where the direct product overflows or underflows although the factor
// itself need not
static inline double gg1_scale(double ia, double Yi, int m, int l,
                               double lia, double lYi) {
  double r = 1.0;
  for (int j = 0; j < m; j++) r *= ia;
  for (int j = 0; j < l; j++) r *= Yi;
  if (std::isfinite(r) && r != 0.0) return r;
  return std::exp(m * lia + l * lYi);
}

// order 1 of log f in (a, d, p)
// [[Rcpp::export]]
List gengamma1_gradient_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a(n);
  NumericVector o_d(n);
  NumericVector o_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[6]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    const double PB0 = R::digamma(k1);
    double* K = P.K;
    const double c0 = 1.0/p;
    const double c1 = PB0*d;
    K[0] = p;
    K[1] = -d;
    K[2] = c0;
    K[3] = -c0*(c1 - p)/d;
    K[4] = -c0;
    K[5] = c1/std::pow(p, 2);
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double UU = pv * gg1_logratio(y[i], av);
    const double TT = std::exp(UU);
    const double lia = -std::log(av);
    const double ia = 1.0 / av;
    const double S_1_0 = gg1_scale(ia, 1.0, 1, 0, lia, 0.0);
    o_a[i] = S_1_0*(K[0]*TT + K[1]);
    o_d[i] = K[2]*UU + K[3];
    o_p[i] = K[4]*TT*UU + K[5];
  });
  return List::create(Named("a") = o_a, Named("d") = o_d, Named("p") = o_p);
}

// order 2 of log f in (a, d, p)
// [[Rcpp::export]]
List gengamma1_hessian_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a_a(n);
  NumericVector o_d_d(n);
  NumericVector o_p_p(n);
  NumericVector o_a_d(n);
  NumericVector o_a_p(n);
  NumericVector o_d_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[9]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    const double PB0 = R::digamma(k1);
    const double PB1 = R::trigamma(k1);
    double* K = P.K;
    const double c0 = std::pow(d, 2);
    const double c1 = std::pow(p, 2);
    const double c2 = 1.0/c1;
    const double c3 = PB1*d;
    const double c4 = PB0*p;
    K[0] = -p*(p + 1);
    K[1] = d;
    K[2] = -c2*(PB1*c0 + c1)/c0;
    K[3] = -c2;
    K[4] = -d*(c3 + 2*c4)/std::pow(p, 4);
    K[5] = -1;
    K[6] = 1;
    K[7] = 1;
    K[8] = (c3 + c4)/std::pow(p, 3);
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double UU = pv * gg1_logratio(y[i], av);
    const double TT = std::exp(UU);
    const double lia = -std::log(av);
    const double ia = 1.0 / av;
    const double S_1_0 = gg1_scale(ia, 1.0, 1, 0, lia, 0.0);
    const double S_2_0 = gg1_scale(ia, 1.0, 2, 0, lia, 0.0);
    o_a_a[i] = S_2_0*(K[0]*TT + K[1]);
    o_d_d[i] = K[2];
    o_p_p[i] = K[3]*TT*std::pow(UU, 2) + K[4];
    o_a_d[i] = K[5]*S_1_0;
    o_a_p[i] = S_1_0*TT*(K[6]*UU + K[7]);
    o_d_p[i] = K[8];
  });
  return List::create(Named("a_a") = o_a_a, Named("d_d") = o_d_d, Named("p_p") = o_p_p, Named("a_d") = o_a_d, Named("a_p") = o_a_p, Named("d_p") = o_d_p);
}

// order 2, expected
// [[Rcpp::export]]
List gengamma1_expected_hessian_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a_a(n);
  NumericVector o_d_d(n);
  NumericVector o_p_p(n);
  NumericVector o_a_d(n);
  NumericVector o_a_p(n);
  NumericVector o_d_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[6]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    const double PB0 = R::digamma(k1);
    const double PB1 = R::trigamma(k1);
    double* K = P.K;
    const double c0 = std::pow(d, 2);
    const double c1 = std::pow(p, 2);
    const double c2 = PB1*d;
    const double c3 = PB0*p;
    K[0] = -d*p;
    K[1] = -(PB1*c0 + c1)/(c0*c1);
    K[2] = -d*(std::pow(PB0, 2)*p + PB1*p + c2 + 2*c3)/std::pow(p, 4);
    K[3] = -1;
    K[4] = d*(PB0 + 1)/p;
    K[5] = (c2 + c3)/std::pow(p, 3);
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double lia = -std::log(av);
    const double ia = 1.0 / av;
    const double S_1_0 = gg1_scale(ia, 1.0, 1, 0, lia, 0.0);
    const double S_2_0 = gg1_scale(ia, 1.0, 2, 0, lia, 0.0);
    o_a_a[i] = K[0]*S_2_0;
    o_d_d[i] = K[1];
    o_p_p[i] = K[2];
    o_a_d[i] = K[3]*S_1_0;
    o_a_p[i] = K[4]*S_1_0;
    o_d_p[i] = K[5];
  });
  return List::create(Named("a_a") = o_a_a, Named("d_d") = o_d_d, Named("p_p") = o_p_p, Named("a_d") = o_a_d, Named("a_p") = o_a_p, Named("d_p") = o_d_p);
}

// order 3 of log f in (a, d, p)
// [[Rcpp::export]]
List gengamma1_deriv3_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a_a_a(n);
  NumericVector o_a_a_d(n);
  NumericVector o_a_a_p(n);
  NumericVector o_a_d_d(n);
  NumericVector o_a_d_p(n);
  NumericVector o_a_p_p(n);
  NumericVector o_d_d_d(n);
  NumericVector o_d_d_p(n);
  NumericVector o_d_p_p(n);
  NumericVector o_p_p_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[12]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    const double PB0 = R::digamma(k1);
    const double PB1 = R::trigamma(k1);
    const double PB2 = R::psigamma(k1, 2.0);
    double* K = P.K;
    const double c0 = std::pow(p, 2);
    const double c1 = 2*p;
    const double c2 = 1.0/p;
    const double c3 = std::pow(d, 3);
    const double c4 = std::pow(p, 3);
    const double c5 = 1.0/c4;
    const double c6 = PB2*std::pow(d, 2);
    const double c7 = PB1*d*p;
    const double c8 = PB0*c0;
    K[0] = p*(c0 + 3*p + 2);
    K[1] = -2*d;
    K[2] = 1;
    K[3] = -p - 1;
    K[4] = -c1 - 1;
    K[5] = c2;
    K[6] = 2*c2;
    K[7] = -c5*(PB2*c3 - 2*c4)/c3;
    K[8] = (PB1*c1 + PB2*d)/std::pow(p, 4);
    K[9] = -(c6 + 4*c7 + 2*c8)/std::pow(p, 5);
    K[10] = -c5;
    K[11] = d*(c6 + 6*c7 + 6*c8)/std::pow(p, 6);
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double UU = pv * gg1_logratio(y[i], av);
    const double TT = std::exp(UU);
    const double lia = -std::log(av);
    const double ia = 1.0 / av;
    const double S_1_0 = gg1_scale(ia, 1.0, 1, 0, lia, 0.0);
    const double S_2_0 = gg1_scale(ia, 1.0, 2, 0, lia, 0.0);
    const double S_3_0 = gg1_scale(ia, 1.0, 3, 0, lia, 0.0);
    o_a_a_a[i] = S_3_0*(K[0]*TT + K[1]);
    o_a_a_d[i] = K[2]*S_2_0;
    o_a_a_p[i] = S_2_0*TT*(K[3]*UU + K[4]);
    o_a_d_d[i] = 0.0;
    o_a_d_p[i] = 0.0;
    o_a_p_p[i] = S_1_0*TT*UU*(K[5]*UU + K[6]);
    o_d_d_d[i] = K[7];
    o_d_d_p[i] = K[8];
    o_d_p_p[i] = K[9];
    o_p_p_p[i] = K[10]*TT*std::pow(UU, 3) + K[11];
  });
  return List::create(Named("a_a_a") = o_a_a_a, Named("a_a_d") = o_a_a_d, Named("a_a_p") = o_a_a_p, Named("a_d_d") = o_a_d_d, Named("a_d_p") = o_a_d_p, Named("a_p_p") = o_a_p_p, Named("d_d_d") = o_d_d_d, Named("d_d_p") = o_d_d_p, Named("d_p_p") = o_d_p_p, Named("p_p_p") = o_p_p_p);
}

// order 3, expected
// [[Rcpp::export]]
List gengamma1_deriv3_expected_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a_a_a(n);
  NumericVector o_a_a_d(n);
  NumericVector o_a_a_p(n);
  NumericVector o_a_d_d(n);
  NumericVector o_a_d_p(n);
  NumericVector o_a_p_p(n);
  NumericVector o_d_d_d(n);
  NumericVector o_d_d_p(n);
  NumericVector o_d_p_p(n);
  NumericVector o_p_p_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[8]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    const double PB0 = R::digamma(k1);
    const double PB1 = R::trigamma(k1);
    const double PB2 = R::psigamma(k1, 2.0);
    double* K = P.K;
    const double c0 = d*p;
    const double c1 = 2*p;
    const double c2 = std::pow(p, 2);
    const double c3 = 2*PB0;
    const double c4 = std::pow(d, 3);
    const double c5 = std::pow(p, 3);
    const double c6 = PB2*std::pow(d, 2);
    const double c7 = PB1*c0;
    const double c8 = PB0*c2;
    K[0] = c0*(p + 3);
    K[1] = 1;
    K[2] = -d*(PB0*p + PB0 + c1 + 1)/p;
    K[3] = d*(std::pow(PB0, 2) + PB1 + c3)/c2;
    K[4] = -(PB2*c4 - 2*c5)/(c4*c5);
    K[5] = (PB1*c1 + PB2*d)/std::pow(p, 4);
    K[6] = -(c2*c3 + c6 + 4*c7)/std::pow(p, 5);
    K[7] = -d*(std::pow(PB0, 3)*c2 + 3*PB1*c8 + PB2*c2 - c6 - 6*c7 - 6*c8)/std::pow(p, 6);
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double lia = -std::log(av);
    const double ia = 1.0 / av;
    const double S_1_0 = gg1_scale(ia, 1.0, 1, 0, lia, 0.0);
    const double S_2_0 = gg1_scale(ia, 1.0, 2, 0, lia, 0.0);
    const double S_3_0 = gg1_scale(ia, 1.0, 3, 0, lia, 0.0);
    o_a_a_a[i] = K[0]*S_3_0;
    o_a_a_d[i] = K[1]*S_2_0;
    o_a_a_p[i] = K[2]*S_2_0;
    o_a_d_d[i] = 0.0;
    o_a_d_p[i] = 0.0;
    o_a_p_p[i] = K[3]*S_1_0;
    o_d_d_d[i] = K[4];
    o_d_d_p[i] = K[5];
    o_d_p_p[i] = K[6];
    o_p_p_p[i] = K[7];
  });
  return List::create(Named("a_a_a") = o_a_a_a, Named("a_a_d") = o_a_a_d, Named("a_a_p") = o_a_a_p, Named("a_d_d") = o_a_d_d, Named("a_d_p") = o_a_d_p, Named("a_p_p") = o_a_p_p, Named("d_d_d") = o_d_d_d, Named("d_d_p") = o_d_d_p, Named("d_p_p") = o_d_p_p, Named("p_p_p") = o_p_p_p);
}

// order 4 of log f in (a, d, p)
// [[Rcpp::export]]
List gengamma1_deriv4_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a_a_a_a(n);
  NumericVector o_a_a_a_d(n);
  NumericVector o_a_a_a_p(n);
  NumericVector o_a_a_d_d(n);
  NumericVector o_a_a_d_p(n);
  NumericVector o_a_a_p_p(n);
  NumericVector o_a_d_d_d(n);
  NumericVector o_a_d_d_p(n);
  NumericVector o_a_d_p_p(n);
  NumericVector o_a_p_p_p(n);
  NumericVector o_d_d_d_d(n);
  NumericVector o_d_d_d_p(n);
  NumericVector o_d_d_p_p(n);
  NumericVector o_d_p_p_p(n);
  NumericVector o_p_p_p_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[16]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    const double PB0 = R::digamma(k1);
    const double PB1 = R::trigamma(k1);
    const double PB2 = R::psigamma(k1, 2.0);
    const double PB3 = R::psigamma(k1, 3.0);
    double* K = P.K;
    const double c0 = std::pow(p, 3);
    const double c1 = std::pow(p, 2);
    const double c2 = 6*c1;
    const double c3 = 6*d;
    const double c4 = 3*p;
    const double c5 = 1.0/p;
    const double c6 = 1.0/c1;
    const double c7 = std::pow(d, 4);
    const double c8 = std::pow(p, 4);
    const double c9 = 1.0/c8;
    const double c10 = std::pow(d, 2);
    const double c11 = PB2*p;
    const double c12 = PB3*std::pow(d, 3);
    const double c13 = PB0*c0;
    const double c14 = PB1*c1*d;
    const double c15 = c10*c11;
    K[0] = -p*(c0 + c2 + 11*p + 6);
    K[1] = c3;
    K[2] = -2;
    K[3] = c1 + c4 + 2;
    K[4] = 3*c1 + 6*p + 2;
    K[5] = -c5*(p + 1);
    K[6] = -2*c5*(2*p + 1);
    K[7] = -2;
    K[8] = c6;
    K[9] = 3*c6;
    K[10] = -c9*(PB3*c7 + 6*c8)/c7;
    K[11] = (PB2*c4 + PB3*d)/std::pow(p, 5);
    K[12] = -(PB1*c2 + PB3*c10 + c11*c3)/std::pow(p, 6);
    K[13] = (c12 + 6*c13 + 18*c14 + 9*c15)/std::pow(p, 7);
    K[14] = -c9;
    K[15] = -d*(c12 + 24*c13 + 36*c14 + 12*c15)/std::pow(p, 8);
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double UU = pv * gg1_logratio(y[i], av);
    const double TT = std::exp(UU);
    const double lia = -std::log(av);
    const double ia = 1.0 / av;
    const double S_1_0 = gg1_scale(ia, 1.0, 1, 0, lia, 0.0);
    const double S_2_0 = gg1_scale(ia, 1.0, 2, 0, lia, 0.0);
    const double S_3_0 = gg1_scale(ia, 1.0, 3, 0, lia, 0.0);
    const double S_4_0 = gg1_scale(ia, 1.0, 4, 0, lia, 0.0);
    const double u0 = std::pow(UU, 2);
    o_a_a_a_a[i] = S_4_0*(K[0]*TT + K[1]);
    o_a_a_a_d[i] = K[2]*S_3_0;
    o_a_a_a_p[i] = S_3_0*TT*(K[3]*UU + K[4]);
    o_a_a_d_d[i] = 0.0;
    o_a_a_d_p[i] = 0.0;
    o_a_a_p_p[i] = S_2_0*TT*(K[5]*u0 + K[6]*UU + K[7]);
    o_a_d_d_d[i] = 0.0;
    o_a_d_d_p[i] = 0.0;
    o_a_d_p_p[i] = 0.0;
    o_a_p_p_p[i] = S_1_0*TT*u0*(K[8]*UU + K[9]);
    o_d_d_d_d[i] = K[10];
    o_d_d_d_p[i] = K[11];
    o_d_d_p_p[i] = K[12];
    o_d_p_p_p[i] = K[13];
    o_p_p_p_p[i] = K[14]*TT*std::pow(UU, 4) + K[15];
  });
  return List::create(Named("a_a_a_a") = o_a_a_a_a, Named("a_a_a_d") = o_a_a_a_d, Named("a_a_a_p") = o_a_a_a_p, Named("a_a_d_d") = o_a_a_d_d, Named("a_a_d_p") = o_a_a_d_p, Named("a_a_p_p") = o_a_a_p_p, Named("a_d_d_d") = o_a_d_d_d, Named("a_d_d_p") = o_a_d_d_p, Named("a_d_p_p") = o_a_d_p_p, Named("a_p_p_p") = o_a_p_p_p, Named("d_d_d_d") = o_d_d_d_d, Named("d_d_d_p") = o_d_d_d_p, Named("d_d_p_p") = o_d_d_p_p, Named("d_p_p_p") = o_d_p_p_p, Named("p_p_p_p") = o_p_p_p_p);
}

// order 4, expected
// [[Rcpp::export]]
List gengamma1_deriv4_expected_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a_a_a_a(n);
  NumericVector o_a_a_a_d(n);
  NumericVector o_a_a_a_p(n);
  NumericVector o_a_a_d_d(n);
  NumericVector o_a_a_d_p(n);
  NumericVector o_a_a_p_p(n);
  NumericVector o_a_d_d_d(n);
  NumericVector o_a_d_d_p(n);
  NumericVector o_a_d_p_p(n);
  NumericVector o_a_p_p_p(n);
  NumericVector o_d_d_d_d(n);
  NumericVector o_d_d_d_p(n);
  NumericVector o_d_d_p_p(n);
  NumericVector o_d_p_p_p(n);
  NumericVector o_p_p_p_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[10]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    const double PB0 = R::digamma(k1);
    const double PB1 = R::trigamma(k1);
    const double PB2 = R::psigamma(k1, 2.0);
    const double PB3 = R::psigamma(k1, 3.0);
    double* K = P.K;
    const double c0 = 6*p;
    const double c1 = std::pow(p, 2);
    const double c2 = 2*PB0;
    const double c3 = PB0*p;
    const double c4 = std::pow(PB0, 2);
    const double c5 = std::pow(p, 3);
    const double c6 = 3*PB1;
    const double c7 = std::pow(d, 4);
    const double c8 = std::pow(p, 4);
    const double c9 = PB2*p;
    const double c10 = std::pow(d, 2);
    const double c11 = PB1*c1;
    const double c12 = PB3*std::pow(d, 3);
    const double c13 = PB0*c5;
    const double c14 = c11*d;
    const double c15 = c10*c9;
    K[0] = -d*p*(c0 + c1 + 11);
    K[1] = -2;
    K[2] = d*(PB0*c1 + c0 + 3*c1 + c2 + 3*c3 + 2)/p;
    K[3] = -d*(PB1*p + PB1 + c2 + 4*c3 + c4*p + c4 + 2*p)/c1;
    K[4] = d*(std::pow(PB0, 3) + PB0*c6 + PB2 + 3*c4 + c6)/c5;
    K[5] = -(PB3*c7 + 6*c8)/(c7*c8);
    K[6] = (PB3*d + 3*c9)/std::pow(p, 5);
    K[7] = -(PB2*c0*d + PB3*c10 + 6*c11)/std::pow(p, 6);
    K[8] = (c12 + 6*c13 + 18*c14 + 9*c15)/std::pow(p, 7);
    K[9] = -d*(std::pow(PB0, 4)*c5 + 3*std::pow(PB1, 2)*c5 + 6*PB1*c4*c5 + 4*PB2*c13 + PB3*c5 + c12 + 24*c13 + 36*c14 + 12*c15)/std::pow(p, 8);
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double lia = -std::log(av);
    const double ia = 1.0 / av;
    const double S_1_0 = gg1_scale(ia, 1.0, 1, 0, lia, 0.0);
    const double S_2_0 = gg1_scale(ia, 1.0, 2, 0, lia, 0.0);
    const double S_3_0 = gg1_scale(ia, 1.0, 3, 0, lia, 0.0);
    const double S_4_0 = gg1_scale(ia, 1.0, 4, 0, lia, 0.0);
    o_a_a_a_a[i] = K[0]*S_4_0;
    o_a_a_a_d[i] = K[1]*S_3_0;
    o_a_a_a_p[i] = K[2]*S_3_0;
    o_a_a_d_d[i] = 0.0;
    o_a_a_d_p[i] = 0.0;
    o_a_a_p_p[i] = K[3]*S_2_0;
    o_a_d_d_d[i] = 0.0;
    o_a_d_d_p[i] = 0.0;
    o_a_d_p_p[i] = 0.0;
    o_a_p_p_p[i] = K[4]*S_1_0;
    o_d_d_d_d[i] = K[5];
    o_d_d_d_p[i] = K[6];
    o_d_d_p_p[i] = K[7];
    o_d_p_p_p[i] = K[8];
    o_p_p_p_p[i] = K[9];
  });
  return List::create(Named("a_a_a_a") = o_a_a_a_a, Named("a_a_a_d") = o_a_a_a_d, Named("a_a_a_p") = o_a_a_a_p, Named("a_a_d_d") = o_a_a_d_d, Named("a_a_d_p") = o_a_a_d_p, Named("a_a_p_p") = o_a_a_p_p, Named("a_d_d_d") = o_a_d_d_d, Named("a_d_d_p") = o_a_d_d_p, Named("a_d_p_p") = o_a_d_p_p, Named("a_p_p_p") = o_a_p_p_p, Named("d_d_d_d") = o_d_d_d_d, Named("d_d_d_p") = o_d_d_d_p, Named("d_d_p_p") = o_d_d_p_p, Named("d_p_p_p") = o_d_p_p_p, Named("p_p_p_p") = o_p_p_p_p);
}

// order 5 of log f in (a, d, p)
// [[Rcpp::export]]
List gengamma1_deriv5_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a_a_a_a_a(n);
  NumericVector o_a_a_a_a_d(n);
  NumericVector o_a_a_a_a_p(n);
  NumericVector o_a_a_a_d_d(n);
  NumericVector o_a_a_a_d_p(n);
  NumericVector o_a_a_a_p_p(n);
  NumericVector o_a_a_d_d_d(n);
  NumericVector o_a_a_d_d_p(n);
  NumericVector o_a_a_d_p_p(n);
  NumericVector o_a_a_p_p_p(n);
  NumericVector o_a_d_d_d_d(n);
  NumericVector o_a_d_d_d_p(n);
  NumericVector o_a_d_d_p_p(n);
  NumericVector o_a_d_p_p_p(n);
  NumericVector o_a_p_p_p_p(n);
  NumericVector o_d_d_d_d_d(n);
  NumericVector o_d_d_d_d_p(n);
  NumericVector o_d_d_d_p_p(n);
  NumericVector o_d_d_p_p_p(n);
  NumericVector o_d_p_p_p_p(n);
  NumericVector o_p_p_p_p_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[20]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    const double PB0 = R::digamma(k1);
    const double PB1 = R::trigamma(k1);
    const double PB2 = R::psigamma(k1, 2.0);
    const double PB3 = R::psigamma(k1, 3.0);
    const double PB4 = R::psigamma(k1, 4.0);
    double* K = P.K;
    const double c0 = std::pow(p, 4);
    const double c1 = std::pow(p, 2);
    const double c2 = std::pow(p, 3);
    const double c3 = 11*p;
    const double c4 = 1.0/p;
    const double c5 = p + 1;
    const double c6 = 1.0/c1;
    const double c7 = 1.0/c2;
    const double c8 = std::pow(d, 5);
    const double c9 = std::pow(p, 5);
    const double c10 = 1.0/c9;
    const double c11 = PB3*p;
    const double c12 = std::pow(d, 2);
    const double c13 = PB2*c1;
    const double c14 = std::pow(d, 3);
    const double c15 = PB1*c2;
    const double c16 = PB4*std::pow(d, 4);
    const double c17 = PB0*c0;
    const double c18 = c15*d;
    const double c19 = c11*c14;
    const double c20 = c12*c13;
    K[0] = p*(c0 + 35*c1 + 10*c2 + 50*p + 24);
    K[1] = -24*d;
    K[2] = 6;
    K[3] = -6*c1 - c2 - c3 - 6;
    K[4] = -18*c1 - 4*c2 - 2*c3 - 6;
    K[5] = c4*(c1 + 3*p + 2);
    K[6] = 2*c4*(3*c1 + 6*p + 2);
    K[7] = 6*c5;
    K[8] = -c5*c6;
    K[9] = -3*c6*(2*p + 1);
    K[10] = -6*c4;
    K[11] = c7;
    K[12] = 4*c7;
    K[13] = -c10*(PB4*c8 - 24*c9)/c8;
    K[14] = (PB4*d + 4*c11)/std::pow(p, 6);
    K[15] = -(PB4*c12 + 8*c11*d + 12*c13)/std::pow(p, 7);
    K[16] = (PB4*c14 + 12*c11*c12 + 36*c13*d + 24*c15)/std::pow(p, 8);
    K[17] = -(c16 + 24*c17 + 96*c18 + 16*c19 + 72*c20)/std::pow(p, 9);
    K[18] = -c10;
    K[19] = d*(c16 + 120*c17 + 240*c18 + 20*c19 + 120*c20)/std::pow(p, 10);
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double UU = pv * gg1_logratio(y[i], av);
    const double TT = std::exp(UU);
    const double lia = -std::log(av);
    const double ia = 1.0 / av;
    const double S_1_0 = gg1_scale(ia, 1.0, 1, 0, lia, 0.0);
    const double S_2_0 = gg1_scale(ia, 1.0, 2, 0, lia, 0.0);
    const double S_3_0 = gg1_scale(ia, 1.0, 3, 0, lia, 0.0);
    const double S_4_0 = gg1_scale(ia, 1.0, 4, 0, lia, 0.0);
    const double S_5_0 = gg1_scale(ia, 1.0, 5, 0, lia, 0.0);
    const double u0 = std::pow(UU, 2);
    o_a_a_a_a_a[i] = S_5_0*(K[0]*TT + K[1]);
    o_a_a_a_a_d[i] = K[2]*S_4_0;
    o_a_a_a_a_p[i] = S_4_0*TT*(K[3]*UU + K[4]);
    o_a_a_a_d_d[i] = 0.0;
    o_a_a_a_d_p[i] = 0.0;
    o_a_a_a_p_p[i] = S_3_0*TT*(K[5]*u0 + K[6]*UU + K[7]);
    o_a_a_d_d_d[i] = 0.0;
    o_a_a_d_d_p[i] = 0.0;
    o_a_a_d_p_p[i] = 0.0;
    o_a_a_p_p_p[i] = S_2_0*TT*UU*(K[10] + K[8]*u0 + K[9]*UU);
    o_a_d_d_d_d[i] = 0.0;
    o_a_d_d_d_p[i] = 0.0;
    o_a_d_d_p_p[i] = 0.0;
    o_a_d_p_p_p[i] = 0.0;
    o_a_p_p_p_p[i] = S_1_0*TT*std::pow(UU, 3)*(K[11]*UU + K[12]);
    o_d_d_d_d_d[i] = K[13];
    o_d_d_d_d_p[i] = K[14];
    o_d_d_d_p_p[i] = K[15];
    o_d_d_p_p_p[i] = K[16];
    o_d_p_p_p_p[i] = K[17];
    o_p_p_p_p_p[i] = K[18]*TT*std::pow(UU, 5) + K[19];
  });
  return List::create(Named("a_a_a_a_a") = o_a_a_a_a_a, Named("a_a_a_a_d") = o_a_a_a_a_d, Named("a_a_a_a_p") = o_a_a_a_a_p, Named("a_a_a_d_d") = o_a_a_a_d_d, Named("a_a_a_d_p") = o_a_a_a_d_p, Named("a_a_a_p_p") = o_a_a_a_p_p, Named("a_a_d_d_d") = o_a_a_d_d_d, Named("a_a_d_d_p") = o_a_a_d_d_p, Named("a_a_d_p_p") = o_a_a_d_p_p, Named("a_a_p_p_p") = o_a_a_p_p_p, Named("a_d_d_d_d") = o_a_d_d_d_d, Named("a_d_d_d_p") = o_a_d_d_d_p, Named("a_d_d_p_p") = o_a_d_d_p_p, Named("a_d_p_p_p") = o_a_d_p_p_p, Named("a_p_p_p_p") = o_a_p_p_p_p, Named("d_d_d_d_d") = o_d_d_d_d_d, Named("d_d_d_d_p") = o_d_d_d_d_p, Named("d_d_d_p_p") = o_d_d_d_p_p, Named("d_d_p_p_p") = o_d_d_p_p_p, Named("d_p_p_p_p") = o_d_p_p_p_p, Named("p_p_p_p_p") = o_p_p_p_p_p);
}

// first derivatives of the expected Hessian
// [[Rcpp::export]]
List gengamma1_dexpected1_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a_a_a(n);
  NumericVector o_a_a_d(n);
  NumericVector o_a_a_p(n);
  NumericVector o_d_d_a(n);
  NumericVector o_d_d_d(n);
  NumericVector o_d_d_p(n);
  NumericVector o_p_p_a(n);
  NumericVector o_p_p_d(n);
  NumericVector o_p_p_p(n);
  NumericVector o_a_d_a(n);
  NumericVector o_a_d_d(n);
  NumericVector o_a_d_p(n);
  NumericVector o_a_p_a(n);
  NumericVector o_a_p_d(n);
  NumericVector o_a_p_p(n);
  NumericVector o_d_p_a(n);
  NumericVector o_d_p_d(n);
  NumericVector o_d_p_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[13]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    const double PB0 = R::digamma(k1);
    const double PB1 = R::trigamma(k1);
    const double PB2 = R::psigamma(k1, 2.0);
    double* K = P.K;
    const double c0 = 2*p;
    const double c1 = std::pow(d, 3);
    const double c2 = std::pow(p, 3);
    const double c3 = 1.0/c2;
    const double c4 = PB2*d;
    const double c5 = (PB1*c0 + c4)/std::pow(p, 4);
    const double c6 = std::pow(p, -5);
    const double c7 = std::pow(p, 2);
    const double c8 = PB1*c7;
    const double c9 = std::pow(PB0, 2)*c7;
    const double c10 = PB2*std::pow(d, 2);
    const double c11 = PB1*d;
    const double c12 = c11*p;
    const double c13 = PB0*c7;
    const double c14 = c10 + 4*c12 + 2*c13;
    const double c15 = PB0*p;
    const double c16 = 2*c11*c15 + c4*p;
    const double c17 = c11 + c15 + p;
    K[0] = c0*d;
    K[1] = -p;
    K[2] = -d;
    K[3] = -c3*(PB2*c1 - 2*c2)/c1;
    K[4] = c5;
    K[5] = -c6*(c14 + c16 + c8 + c9);
    K[6] = d*(c10 + 6*c12 + 6*c13 + c16 + 3*c8 + 3*c9)/std::pow(p, 6);
    K[7] = 1;
    K[8] = -d*(PB0 + 1)/p;
    K[9] = c17/c7;
    K[10] = -c17*c3*d;
    K[11] = c5;
    K[12] = -c14*c6;
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double lia = -std::log(av);
    const double ia = 1.0 / av;
    const double S_1_0 = gg1_scale(ia, 1.0, 1, 0, lia, 0.0);
    const double S_2_0 = gg1_scale(ia, 1.0, 2, 0, lia, 0.0);
    const double S_3_0 = gg1_scale(ia, 1.0, 3, 0, lia, 0.0);
    o_a_a_a[i] = K[0]*S_3_0;
    o_a_a_d[i] = K[1]*S_2_0;
    o_a_a_p[i] = K[2]*S_2_0;
    o_d_d_a[i] = 0.0;
    o_d_d_d[i] = K[3];
    o_d_d_p[i] = K[4];
    o_p_p_a[i] = 0.0;
    o_p_p_d[i] = K[5];
    o_p_p_p[i] = K[6];
    o_a_d_a[i] = K[7]*S_2_0;
    o_a_d_d[i] = 0.0;
    o_a_d_p[i] = 0.0;
    o_a_p_a[i] = K[8]*S_2_0;
    o_a_p_d[i] = K[9]*S_1_0;
    o_a_p_p[i] = K[10]*S_1_0;
    o_d_p_a[i] = 0.0;
    o_d_p_d[i] = K[11];
    o_d_p_p[i] = K[12];
  });
  return List::create(Named("a_a_a") = o_a_a_a, Named("a_a_d") = o_a_a_d, Named("a_a_p") = o_a_a_p, Named("d_d_a") = o_d_d_a, Named("d_d_d") = o_d_d_d, Named("d_d_p") = o_d_d_p, Named("p_p_a") = o_p_p_a, Named("p_p_d") = o_p_p_d, Named("p_p_p") = o_p_p_p, Named("a_d_a") = o_a_d_a, Named("a_d_d") = o_a_d_d, Named("a_d_p") = o_a_d_p, Named("a_p_a") = o_a_p_a, Named("a_p_d") = o_a_p_d, Named("a_p_p") = o_a_p_p, Named("d_p_a") = o_d_p_a, Named("d_p_d") = o_d_p_d, Named("d_p_p") = o_d_p_p);
}

// second derivatives of the expected Hessian
// [[Rcpp::export]]
List gengamma1_dexpected2_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a_a_a_a(n);
  NumericVector o_a_a_d_d(n);
  NumericVector o_a_a_p_p(n);
  NumericVector o_a_a_a_d(n);
  NumericVector o_a_a_a_p(n);
  NumericVector o_a_a_d_p(n);
  NumericVector o_d_d_a_a(n);
  NumericVector o_d_d_d_d(n);
  NumericVector o_d_d_p_p(n);
  NumericVector o_d_d_a_d(n);
  NumericVector o_d_d_a_p(n);
  NumericVector o_d_d_d_p(n);
  NumericVector o_p_p_a_a(n);
  NumericVector o_p_p_d_d(n);
  NumericVector o_p_p_p_p(n);
  NumericVector o_p_p_a_d(n);
  NumericVector o_p_p_a_p(n);
  NumericVector o_p_p_d_p(n);
  NumericVector o_a_d_a_a(n);
  NumericVector o_a_d_d_d(n);
  NumericVector o_a_d_p_p(n);
  NumericVector o_a_d_a_d(n);
  NumericVector o_a_d_a_p(n);
  NumericVector o_a_d_d_p(n);
  NumericVector o_a_p_a_a(n);
  NumericVector o_a_p_d_d(n);
  NumericVector o_a_p_p_p(n);
  NumericVector o_a_p_a_d(n);
  NumericVector o_a_p_a_p(n);
  NumericVector o_a_p_d_p(n);
  NumericVector o_d_p_a_a(n);
  NumericVector o_d_p_d_d(n);
  NumericVector o_d_p_p_p(n);
  NumericVector o_d_p_a_d(n);
  NumericVector o_d_p_a_p(n);
  NumericVector o_d_p_d_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[20]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    const double PB0 = R::digamma(k1);
    const double PB1 = R::trigamma(k1);
    const double PB2 = R::psigamma(k1, 2.0);
    const double PB3 = R::psigamma(k1, 3.0);
    double* K = P.K;
    const double c0 = 6*p;
    const double c1 = 2*p;
    const double c2 = 2*d;
    const double c3 = std::pow(d, 4);
    const double c4 = std::pow(p, 4);
    const double c5 = 1.0/c4;
    const double c6 = std::pow(p, -6);
    const double c7 = std::pow(d, 2);
    const double c8 = PB3*c7;
    const double c9 = PB2*d;
    const double c10 = std::pow(p, 2);
    const double c11 = 6*PB1*c10 + c0*c9 + c8;
    const double c12 = -c11*c6;
    const double c13 = std::pow(p, -5);
    const double c14 = PB3*d;
    const double c15 = 3*p;
    const double c16 = c13*(PB2*c15 + c14);
    const double c17 = 2*c10;
    const double c18 = std::pow(PB1, 2)*c1;
    const double c19 = PB0*c10;
    const double c20 = PB0*c1;
    const double c21 = PB3*std::pow(d, 3);
    const double c22 = std::pow(p, 3);
    const double c23 = 12*c22;
    const double c24 = std::pow(PB0, 2);
    const double c25 = PB0*c22;
    const double c26 = c10*c9;
    const double c27 = PB2*c7;
    const double c28 = c27*p;
    const double c29 = PB1*d;
    const double c30 = c10*c29;
    const double c31 = c19*c29;
    const double c32 = c18*c7 + c20*c27 + c8*p;
    const double c33 = std::pow(p, -7);
    const double c34 = 3*c22;
    const double c35 = c21 + 6*c25 + 9*c28 + 18*c30;
    const double c36 = 1.0/c22;
    const double c37 = PB0*p + c29 + p;
    K[0] = -c0*d;
    K[1] = c1;
    K[2] = c2;
    K[3] = -1;
    K[4] = -c5*(PB3*c3 + 6*c4)/c3;
    K[5] = c12;
    K[6] = c16;
    K[7] = -c6*(4*PB1*c19 + PB2*c17 + c11 + c14*p + c18*d + c20*c9);
    K[8] = -d*(PB1*c23 + c21 + c23*c24 + 24*c25 + 8*c26 + 12*c28 + 36*c30 + 16*c31 + c32)/std::pow(p, 8);
    K[9] = c33*(PB1*c34 + c24*c34 + 5*c26 + 10*c31 + c32 + c35);
    K[10] = -2;
    K[11] = c2*(PB0 + 1)/p;
    K[12] = c36*(PB1*c1 + c9);
    K[13] = c13*d*(PB0*c17 + c17 + c27 + 4*c29*p);
    K[14] = -c37/c10;
    K[15] = c36*c37*d;
    K[16] = -c5*(c10 + c15*c29 + c19 + c27);
    K[17] = c16;
    K[18] = c33*c35;
    K[19] = c12;
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double lia = -std::log(av);
    const double ia = 1.0 / av;
    const double S_1_0 = gg1_scale(ia, 1.0, 1, 0, lia, 0.0);
    const double S_2_0 = gg1_scale(ia, 1.0, 2, 0, lia, 0.0);
    const double S_3_0 = gg1_scale(ia, 1.0, 3, 0, lia, 0.0);
    const double S_4_0 = gg1_scale(ia, 1.0, 4, 0, lia, 0.0);
    o_a_a_a_a[i] = K[0]*S_4_0;
    o_a_a_d_d[i] = 0.0;
    o_a_a_p_p[i] = 0.0;
    o_a_a_a_d[i] = K[1]*S_3_0;
    o_a_a_a_p[i] = K[2]*S_3_0;
    o_a_a_d_p[i] = K[3]*S_2_0;
    o_d_d_a_a[i] = 0.0;
    o_d_d_d_d[i] = K[4];
    o_d_d_p_p[i] = K[5];
    o_d_d_a_d[i] = 0.0;
    o_d_d_a_p[i] = 0.0;
    o_d_d_d_p[i] = K[6];
    o_p_p_a_a[i] = 0.0;
    o_p_p_d_d[i] = K[7];
    o_p_p_p_p[i] = K[8];
    o_p_p_a_d[i] = 0.0;
    o_p_p_a_p[i] = 0.0;
    o_p_p_d_p[i] = K[9];
    o_a_d_a_a[i] = K[10]*S_3_0;
    o_a_d_d_d[i] = 0.0;
    o_a_d_p_p[i] = 0.0;
    o_a_d_a_d[i] = 0.0;
    o_a_d_a_p[i] = 0.0;
    o_a_d_d_p[i] = 0.0;
    o_a_p_a_a[i] = K[11]*S_3_0;
    o_a_p_d_d[i] = K[12]*S_1_0;
    o_a_p_p_p[i] = K[13]*S_1_0;
    o_a_p_a_d[i] = K[14]*S_2_0;
    o_a_p_a_p[i] = K[15]*S_2_0;
    o_a_p_d_p[i] = K[16]*S_1_0;
    o_d_p_a_a[i] = 0.0;
    o_d_p_d_d[i] = K[17];
    o_d_p_p_p[i] = K[18];
    o_d_p_a_d[i] = 0.0;
    o_d_p_a_p[i] = 0.0;
    o_d_p_d_p[i] = K[19];
  });
  return List::create(Named("a_a_a_a") = o_a_a_a_a, Named("a_a_d_d") = o_a_a_d_d, Named("a_a_p_p") = o_a_a_p_p, Named("a_a_a_d") = o_a_a_a_d, Named("a_a_a_p") = o_a_a_a_p, Named("a_a_d_p") = o_a_a_d_p, Named("d_d_a_a") = o_d_d_a_a, Named("d_d_d_d") = o_d_d_d_d, Named("d_d_p_p") = o_d_d_p_p, Named("d_d_a_d") = o_d_d_a_d, Named("d_d_a_p") = o_d_d_a_p, Named("d_d_d_p") = o_d_d_d_p, Named("p_p_a_a") = o_p_p_a_a, Named("p_p_d_d") = o_p_p_d_d, Named("p_p_p_p") = o_p_p_p_p, Named("p_p_a_d") = o_p_p_a_d, Named("p_p_a_p") = o_p_p_a_p, Named("p_p_d_p") = o_p_p_d_p, Named("a_d_a_a") = o_a_d_a_a, Named("a_d_d_d") = o_a_d_d_d, Named("a_d_p_p") = o_a_d_p_p, Named("a_d_a_d") = o_a_d_a_d, Named("a_d_a_p") = o_a_d_a_p, Named("a_d_d_p") = o_a_d_d_p, Named("a_p_a_a") = o_a_p_a_a, Named("a_p_d_d") = o_a_p_d_d, Named("a_p_p_p") = o_a_p_p_p, Named("a_p_a_d") = o_a_p_a_d, Named("a_p_a_p") = o_a_p_a_p, Named("a_p_d_p") = o_a_p_d_p, Named("d_p_a_a") = o_d_p_a_a, Named("d_p_d_d") = o_d_p_d_d, Named("d_p_p_p") = o_d_p_p_p, Named("d_p_a_d") = o_d_p_a_d, Named("d_p_a_p") = o_d_p_a_p, Named("d_p_d_p") = o_d_p_d_p);
}

// 1 in y and 1 in (a, d, p)
// [[Rcpp::export]]
List gengamma1_cross_y_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a(n);
  NumericVector o_d(n);
  NumericVector o_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[4]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    double* K = P.K;
    K[0] = std::pow(p, 2);
    K[1] = 1;
    K[2] = -1;
    K[3] = -1;
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double UU = pv * gg1_logratio(y[i], av);
    const double TT = std::exp(UU);
    const double lia = -std::log(av);
    const double lYi = std::log(av) - std::log(y[i]);
    const double Yi = av / y[i];
    const double ia = 1.0 / av;
    const double S_1_1 = gg1_scale(ia, Yi, 1, 1, lia, lYi);
    const double S_2_1 = gg1_scale(ia, Yi, 2, 1, lia, lYi);
    o_a[i] = K[0]*S_2_1*TT;
    o_d[i] = K[1]*S_1_1;
    o_p[i] = S_1_1*TT*(K[2]*UU + K[3]);
  });
  return List::create(Named("a") = o_a, Named("d") = o_d, Named("p") = o_p);
}

// 2 in y and 1 in (a, d, p)
// [[Rcpp::export]]
List gengamma1_cross2_y_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a(n);
  NumericVector o_d(n);
  NumericVector o_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[4]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    double* K = P.K;
    const double c0 = p - 1;
    K[0] = c0*std::pow(p, 2);
    K[1] = -1;
    K[2] = -c0;
    K[3] = 1 - 2*p;
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double UU = pv * gg1_logratio(y[i], av);
    const double TT = std::exp(UU);
    const double lia = -std::log(av);
    const double lYi = std::log(av) - std::log(y[i]);
    const double Yi = av / y[i];
    const double ia = 1.0 / av;
    const double S_2_2 = gg1_scale(ia, Yi, 2, 2, lia, lYi);
    const double S_3_2 = gg1_scale(ia, Yi, 3, 2, lia, lYi);
    o_a[i] = K[0]*S_3_2*TT;
    o_d[i] = K[1]*S_2_2;
    o_p[i] = S_2_2*TT*(K[2]*UU + K[3]);
  });
  return List::create(Named("a") = o_a, Named("d") = o_d, Named("p") = o_p);
}

// 1 in y and 2 in (a, d, p)
// [[Rcpp::export]]
List gengamma1_grad_y_hess_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a_a(n);
  NumericVector o_d_d(n);
  NumericVector o_p_p(n);
  NumericVector o_a_d(n);
  NumericVector o_a_p(n);
  NumericVector o_d_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[5]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    double* K = P.K;
    const double c0 = 1.0/p;
    K[0] = -std::pow(p, 2)*(p + 1);
    K[1] = -c0;
    K[2] = -2*c0;
    K[3] = p;
    K[4] = 2*p;
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double UU = pv * gg1_logratio(y[i], av);
    const double TT = std::exp(UU);
    const double lia = -std::log(av);
    const double lYi = std::log(av) - std::log(y[i]);
    const double Yi = av / y[i];
    const double ia = 1.0 / av;
    const double S_1_1 = gg1_scale(ia, Yi, 1, 1, lia, lYi);
    const double S_2_1 = gg1_scale(ia, Yi, 2, 1, lia, lYi);
    const double S_3_1 = gg1_scale(ia, Yi, 3, 1, lia, lYi);
    o_a_a[i] = K[0]*S_3_1*TT;
    o_d_d[i] = 0.0;
    o_p_p[i] = S_1_1*TT*UU*(K[1]*UU + K[2]);
    o_a_d[i] = 0.0;
    o_a_p[i] = S_2_1*TT*(K[3]*UU + K[4]);
    o_d_p[i] = 0.0;
  });
  return List::create(Named("a_a") = o_a_a, Named("d_d") = o_d_d, Named("p_p") = o_p_p, Named("a_d") = o_a_d, Named("a_p") = o_a_p, Named("d_p") = o_d_p);
}

// 2 in y and 2 in (a, d, p)
// [[Rcpp::export]]
List gengamma1_hess_y_hess_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_a_a(n);
  NumericVector o_d_d(n);
  NumericVector o_p_p(n);
  NumericVector o_a_d(n);
  NumericVector o_a_p(n);
  NumericVector o_d_p(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[6]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    double* K = P.K;
    const double c0 = std::pow(p, 2);
    const double c1 = 1.0/p;
    const double c2 = p - 1;
    K[0] = c0*(1 - c0);
    K[1] = -c1*c2;
    K[2] = -2*c1*(2*p - 1);
    K[3] = -2;
    K[4] = c2*p;
    K[5] = p*(3*p - 2);
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double UU = pv * gg1_logratio(y[i], av);
    const double TT = std::exp(UU);
    const double lia = -std::log(av);
    const double lYi = std::log(av) - std::log(y[i]);
    const double Yi = av / y[i];
    const double ia = 1.0 / av;
    const double S_2_2 = gg1_scale(ia, Yi, 2, 2, lia, lYi);
    const double S_3_2 = gg1_scale(ia, Yi, 3, 2, lia, lYi);
    const double S_4_2 = gg1_scale(ia, Yi, 4, 2, lia, lYi);
    o_a_a[i] = K[0]*S_4_2*TT;
    o_d_d[i] = 0.0;
    o_p_p[i] = S_2_2*TT*(K[1]*std::pow(UU, 2) + K[2]*UU + K[3]);
    o_a_d[i] = 0.0;
    o_a_p[i] = S_3_2*TT*(K[4]*UU + K[5]);
    o_d_p[i] = 0.0;
  });
  return List::create(Named("a_a") = o_a_a, Named("d_d") = o_d_d, Named("p_p") = o_p_p, Named("a_d") = o_a_d, Named("a_p") = o_a_p, Named("d_p") = o_d_p);
}

// order 1 of log f in y
// [[Rcpp::export]]
List gengamma1_dy1_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_y(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[2]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    double* K = P.K;
    K[0] = -p;
    K[1] = d - 1;
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double UU = pv * gg1_logratio(y[i], av);
    const double TT = std::exp(UU);
    const double lia = -std::log(av);
    const double lYi = std::log(av) - std::log(y[i]);
    const double Yi = av / y[i];
    const double ia = 1.0 / av;
    const double S_1_1 = gg1_scale(ia, Yi, 1, 1, lia, lYi);
    o_y[i] = S_1_1*(K[0]*TT + K[1]);
  });
  return List::create(Named("y") = o_y);
}

// order 2 of log f in y
// [[Rcpp::export]]
List gengamma1_dy2_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_y(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[2]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    double* K = P.K;
    K[0] = p*(1 - p);
    K[1] = 1 - d;
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double UU = pv * gg1_logratio(y[i], av);
    const double TT = std::exp(UU);
    const double lia = -std::log(av);
    const double lYi = std::log(av) - std::log(y[i]);
    const double Yi = av / y[i];
    const double ia = 1.0 / av;
    const double S_2_2 = gg1_scale(ia, Yi, 2, 2, lia, lYi);
    o_y[i] = S_2_2*(K[0]*TT + K[1]);
  });
  return List::create(Named("y") = o_y);
}

// order 3 of log f in y
// [[Rcpp::export]]
List gengamma1_dy3_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_y(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[2]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    double* K = P.K;
    K[0] = p*(-std::pow(p, 2) + 3*p - 2);
    K[1] = 2*(d - 1);
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double UU = pv * gg1_logratio(y[i], av);
    const double TT = std::exp(UU);
    const double lia = -std::log(av);
    const double lYi = std::log(av) - std::log(y[i]);
    const double Yi = av / y[i];
    const double ia = 1.0 / av;
    const double S_3_3 = gg1_scale(ia, Yi, 3, 3, lia, lYi);
    o_y[i] = S_3_3*(K[0]*TT + K[1]);
  });
  return List::create(Named("y") = o_y);
}

// order 4 of log f in y
// [[Rcpp::export]]
List gengamma1_dy4_cpp(NumericVector y, NumericVector a_, NumericVector d_, NumericVector p_, int threads = 1) {
  const int n = y.size();
  const int n_a = a_.size(), n_d = d_.size(), n_p = p_.size();
  NumericVector o_y(n);
  // the coefficients depend on (d, p) alone: once when they are scalars
  struct Prm { double K[2]; };
  auto make = [&](std::size_t i, Prm& P) {
    const double d = d_[i % n_d], p = p_[i % n_p];
    const double k1 = d / p + 1.0;
    (void) k1; (void) d; (void) p;
    double* K = P.K;
    K[0] = p*(-std::pow(p, 3) + 6*std::pow(p, 2) - 11*p + 6);
    K[1] = 6*(1 - d);
  };
  const bool scalar = n_d == 1 && n_p == 1;
  Prm P0;
  if (scalar) make(0, P0);
  d7::par_for(n, threads, scalar ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    Prm Pi;
    if (!scalar) make(i, Pi);
    const double* K = scalar ? P0.K : Pi.K;
    (void) K;
    const double av = a_[i % n_a], pv = p_[i % n_p];
    (void) pv;
    const double UU = pv * gg1_logratio(y[i], av);
    const double TT = std::exp(UU);
    const double lia = -std::log(av);
    const double lYi = std::log(av) - std::log(y[i]);
    const double Yi = av / y[i];
    const double ia = 1.0 / av;
    const double S_4_4 = gg1_scale(ia, Yi, 4, 4, lia, lYi);
    o_y[i] = S_4_4*(K[0]*TT + K[1]);
  });
  return List::create(Named("y") = o_y);
}
