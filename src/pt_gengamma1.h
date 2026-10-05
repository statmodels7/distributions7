#ifndef D7_PT_GENGAMMA1_H
#define D7_PT_GENGAMMA1_H

#include <Rcpp.h>
#include <cmath>
#include <cfloat>

// The generalized gamma in Stacy's form, scale a, shape d and power p: one
// function per component and order for the quantities the scalar registry
// reads, written from the coefficient tables of the generated kernels in
// gengamma1.cpp (which records the derivation), each with its entry's own
// expression, and called by those kernels and by the registry
// (d7_ccallable.cpp). With k1 = d/p + 1 and U = p log(y/a), the caller
// passes the polygammas PB_j = psigamma(k1, j), U, T = e^U and the scale
// factors S_m = a^-m from gg1_scale().

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

namespace d7 {

inline double gengamma1_score_a(double d, double p, double TT, double S1) {
  const double K0 = p;
  const double K1 = -d;
  return S1*(K0*TT + K1);
}

inline double gengamma1_score_d(double d, double p, double PB0, double UU) {
  const double c0 = 1.0/p;
  const double c1 = PB0*d;
  const double K2 = c0;
  const double K3 = -c0*(c1 - p)/d;
  return K2*UU + K3;
}

inline double gengamma1_score_p(double d, double p, double PB0, double UU,
                                double TT) {
  const double c0 = 1.0/p;
  const double c1 = PB0*d;
  const double K4 = -c0;
  const double K5 = c1/std::pow(p, 2);
  return K4*TT*UU + K5;
}

inline double gengamma1_hess_a_a(double d, double p, double TT, double S2) {
  const double K0 = -p*(p + 1);
  const double K1 = d;
  return S2*(K0*TT + K1);
}

inline double gengamma1_hess_d_d(double d, double p, double PB1) {
  const double c0 = std::pow(d, 2);
  const double c1 = std::pow(p, 2);
  const double c2 = 1.0/c1;
  return -c2*(PB1*c0 + c1)/c0;
}

inline double gengamma1_hess_p_p(double d, double p, double PB0, double PB1,
                                 double UU, double TT) {
  const double c1 = std::pow(p, 2);
  const double c2 = 1.0/c1;
  const double c3 = PB1*d;
  const double c4 = PB0*p;
  const double K3 = -c2;
  const double K4 = -d*(c3 + 2*c4)/std::pow(p, 4);
  return K3*TT*std::pow(UU, 2) + K4;
}

inline double gengamma1_expected_a_a(double d, double p, double S2) {
  const double K0 = -d*p;
  return K0*S2;
}

inline double gengamma1_expected_d_d(double d, double p, double PB1) {
  const double c0 = std::pow(d, 2);
  const double c1 = std::pow(p, 2);
  return -(PB1*c0 + c1)/(c0*c1);
}

inline double gengamma1_expected_p_p(double d, double p, double PB0,
                                     double PB1) {
  const double c2 = PB1*d;
  const double c3 = PB0*p;
  return -d*(std::pow(PB0, 2)*p + PB1*p + c2 + 2*c3)/std::pow(p, 4);
}

inline double gengamma1_dexpected_a_a_a(double d, double p, double S3) {
  const double c0 = 2*p;
  const double K0 = c0*d;
  return K0*S3;
}

inline double gengamma1_dexpected_d_d_d(double d, double p, double PB2) {
  const double c1 = std::pow(d, 3);
  const double c2 = std::pow(p, 3);
  const double c3 = 1.0/c2;
  return -c3*(PB2*c1 - 2*c2)/c1;
}

inline double gengamma1_dexpected_p_p_p(double d, double p, double PB0,
                                        double PB1, double PB2) {
  const double c4 = PB2*d;
  const double c7 = std::pow(p, 2);
  const double c8 = PB1*c7;
  const double c9 = std::pow(PB0, 2)*c7;
  const double c10 = PB2*std::pow(d, 2);
  const double c11 = PB1*d;
  const double c12 = c11*p;
  const double c13 = PB0*c7;
  const double c15 = PB0*p;
  const double c16 = 2*c11*c15 + c4*p;
  return d*(c10 + 6*c12 + 6*c13 + c16 + 3*c8 + 3*c9)/std::pow(p, 6);
}

inline void gengamma1_score_curv(int k, double y, const double* th,
                                 double* out) {
  const double a = th[0], d = th[1], p = th[2];
  const double UU = p * gg1_logratio(y, a);
  const double TT = std::exp(UU);
  if (k == 0) {
    const double lia = -std::log(a);
    const double ia = 1.0 / a;
    out[0] = gengamma1_score_a(d, p, TT, gg1_scale(ia, 1.0, 1, 0, lia, 0.0));
    out[1] = gengamma1_hess_a_a(d, p, TT, gg1_scale(ia, 1.0, 2, 0, lia, 0.0));
    return;
  }
  const double k1 = d / p + 1.0;
  const double PB0 = R::digamma(k1);
  if (k == 1) {
    out[0] = gengamma1_score_d(d, p, PB0, UU);
    out[1] = gengamma1_hess_d_d(d, p, R::trigamma(k1));
  } else {
    out[0] = gengamma1_score_p(d, p, PB0, UU, TT);
    out[1] = gengamma1_hess_p_p(d, p, PB0, R::trigamma(k1), UU, TT);
  }
}

inline void gengamma1_info_dinfo(int k, double y, const double* th,
                                 double* out) {
  const double a = th[0], d = th[1], p = th[2];
  if (k == 0) {
    const double lia = -std::log(a);
    const double ia = 1.0 / a;
    out[0] = gengamma1_expected_a_a(d, p, gg1_scale(ia, 1.0, 2, 0, lia, 0.0));
    out[1] = gengamma1_dexpected_a_a_a(d, p,
                                       gg1_scale(ia, 1.0, 3, 0, lia, 0.0));
    return;
  }
  const double k1 = d / p + 1.0;
  const double PB1 = R::trigamma(k1), PB2 = R::psigamma(k1, 2.0);
  if (k == 1) {
    out[0] = gengamma1_expected_d_d(d, p, PB1);
    out[1] = gengamma1_dexpected_d_d_d(d, p, PB2);
  } else {
    const double PB0 = R::digamma(k1);
    out[0] = gengamma1_expected_p_p(d, p, PB0, PB1);
    out[1] = gengamma1_dexpected_p_p_p(d, p, PB0, PB1, PB2);
  }
}

} // namespace d7

#endif
