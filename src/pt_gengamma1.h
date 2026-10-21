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

// The components of the joint registry that the scalar routers do not read.
// Each one reproduces an expression of gengamma1.cpp (the make() lambdas of
// gengamma1_hessian_cpp, gengamma1_expected_hessian_cpp and
// gengamma1_dexpected1_cpp, their constants K included), with the same
// operations in the same order. The arguments are the quantities the kernels
// use: d, p, the polygammas PB0 = digamma(k1), PB1 = trigamma(k1) and
// PB2 = psigamma(k1, 2) at k1 = d/p + 1, UU, TT and the scale factors S_m.

// order 2, the Hessian entries a_d, a_p and d_p
inline double gengamma1_hess_a_d(double S1) {
  const double K5 = -1;
  return K5*S1;
}

inline double gengamma1_hess_a_p(double UU, double TT, double S1) {
  const double K6 = 1, K7 = 1;
  return S1*TT*(K6*UU + K7);
}

inline double gengamma1_hess_d_p(double d, double p, double PB0, double PB1) {
  const double c3 = PB1*d;
  const double c4 = PB0*p;
  return (c3 + c4)/std::pow(p, 3);
}

// order 2, expected: the entries a_d, a_p and d_p
inline double gengamma1_expected_a_d(double S1) {
  const double K3 = -1;
  return K3*S1;
}

inline double gengamma1_expected_a_p(double d, double p, double PB0, double S1) {
  const double K4 = d*(PB0 + 1)/p;
  return K4*S1;
}

inline double gengamma1_expected_d_p(double d, double p, double PB0, double PB1) {
  const double c2 = PB1*d;
  const double c3 = PB0*p;
  return (c2 + c3)/std::pow(p, 3);
}

// first derivatives of the expected Hessian, named pair then differentiating
// parameter; the entries of the kernel that vanish are written as 0.0 in the
// router. The kernel's d_p_d and d_d_p are one and the same expression (K[11]
// and K[4]), so one function serves both.
inline double gengamma1_dexpected_a_a_d(double p, double S2) {
  const double K1 = -p;
  return K1*S2;
}

inline double gengamma1_dexpected_a_a_p(double d, double S2) {
  const double K2 = -d;
  return K2*S2;
}

inline double gengamma1_dexpected_d_d_p(double d, double p, double PB1,
                                        double PB2) {
  const double c0 = 2*p;
  const double c4 = PB2*d;
  return (PB1*c0 + c4)/std::pow(p, 4);
}

inline double gengamma1_dexpected_p_p_d(double d, double p, double PB0,
                                        double PB1, double PB2) {
  const double c4 = PB2*d;
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
  return -c6*(c14 + c16 + c8 + c9);
}

inline double gengamma1_dexpected_a_d_a(double S2) {
  const double K7 = 1;
  return K7*S2;
}

inline double gengamma1_dexpected_a_p_a(double d, double p, double PB0,
                                        double S2) {
  const double K8 = -d*(PB0 + 1)/p;
  return K8*S2;
}

inline double gengamma1_dexpected_a_p_d(double d, double p, double PB0,
                                        double PB1, double S1) {
  const double c7 = std::pow(p, 2);
  const double c11 = PB1*d;
  const double c15 = PB0*p;
  const double c17 = c11 + c15 + p;
  const double K9 = c17/c7;
  return K9*S1;
}

inline double gengamma1_dexpected_a_p_p(double d, double p, double PB0,
                                        double PB1, double S1) {
  const double c2 = std::pow(p, 3);
  const double c3 = 1.0/c2;
  const double c11 = PB1*d;
  const double c15 = PB0*p;
  const double c17 = c11 + c15 + p;
  const double K10 = -c17*c3*d;
  return K10*S1;
}

inline double gengamma1_dexpected_d_p_p(double d, double p, double PB0,
                                        double PB1, double PB2) {
  const double c6 = std::pow(p, -5);
  const double c7 = std::pow(p, 2);
  const double c10 = PB2*std::pow(d, 2);
  const double c11 = PB1*d;
  const double c12 = c11*p;
  const double c13 = PB0*c7;
  const double c14 = c10 + 4*c12 + 2*c13;
  const double K12 = -c14*c6;
  return K12;
}

// The routers of the joint registry: the components of the log-density at one
// observation on the parameter scale, P = 3 and the parameters in the order
// a, d, p. The Hessian and the expected information fill a full 3 x 3 matrix
// (out[a + 3 b], both triangles) and the derivative of the expected
// information a 3 x 3 x 3 array (out[a + 3 b + 9 c] is the derivative of
// E[l_ab] in parameter c, symmetric in a and b).
inline void gengamma1_joint_grad(double y, const double* th, double* out) {
  const double a = th[0], d = th[1], p = th[2];
  const double UU = p * gg1_logratio(y, a);
  const double TT = std::exp(UU);
  const double lia = -std::log(a);
  const double ia = 1.0 / a;
  const double k1 = d / p + 1.0;
  const double PB0 = R::digamma(k1);
  out[0] = gengamma1_score_a(d, p, TT, gg1_scale(ia, 1.0, 1, 0, lia, 0.0));
  out[1] = gengamma1_score_d(d, p, PB0, UU);
  out[2] = gengamma1_score_p(d, p, PB0, UU, TT);
}

inline void gengamma1_joint_hess(double y, const double* th, double* out) {
  const double a = th[0], d = th[1], p = th[2];
  const double UU = p * gg1_logratio(y, a);
  const double TT = std::exp(UU);
  const double lia = -std::log(a);
  const double ia = 1.0 / a;
  const double S1 = gg1_scale(ia, 1.0, 1, 0, lia, 0.0);
  const double S2 = gg1_scale(ia, 1.0, 2, 0, lia, 0.0);
  const double k1 = d / p + 1.0;
  const double PB0 = R::digamma(k1);
  const double PB1 = R::trigamma(k1);
  out[0] = gengamma1_hess_a_a(d, p, TT, S2);
  out[4] = gengamma1_hess_d_d(d, p, PB1);
  out[8] = gengamma1_hess_p_p(d, p, PB0, PB1, UU, TT);
  out[1] = out[3] = gengamma1_hess_a_d(S1);
  out[2] = out[6] = gengamma1_hess_a_p(UU, TT, S1);
  out[5] = out[7] = gengamma1_hess_d_p(d, p, PB0, PB1);
}

inline void gengamma1_joint_info(double y, const double* th, double* out) {
  const double a = th[0], d = th[1], p = th[2];
  const double lia = -std::log(a);
  const double ia = 1.0 / a;
  const double S1 = gg1_scale(ia, 1.0, 1, 0, lia, 0.0);
  const double S2 = gg1_scale(ia, 1.0, 2, 0, lia, 0.0);
  const double k1 = d / p + 1.0;
  const double PB0 = R::digamma(k1);
  const double PB1 = R::trigamma(k1);
  out[0] = gengamma1_expected_a_a(d, p, S2);
  out[4] = gengamma1_expected_d_d(d, p, PB1);
  out[8] = gengamma1_expected_p_p(d, p, PB0, PB1);
  out[1] = out[3] = gengamma1_expected_a_d(S1);
  out[2] = out[6] = gengamma1_expected_a_p(d, p, PB0, S1);
  out[5] = out[7] = gengamma1_expected_d_p(d, p, PB0, PB1);
}

inline void gengamma1_joint_dinfo(double y, const double* th, double* out) {
  const double a = th[0], d = th[1], p = th[2];
  const double lia = -std::log(a);
  const double ia = 1.0 / a;
  const double S1 = gg1_scale(ia, 1.0, 1, 0, lia, 0.0);
  const double S2 = gg1_scale(ia, 1.0, 2, 0, lia, 0.0);
  const double S3 = gg1_scale(ia, 1.0, 3, 0, lia, 0.0);
  const double k1 = d / p + 1.0;
  const double PB0 = R::digamma(k1);
  const double PB1 = R::trigamma(k1);
  const double PB2 = R::psigamma(k1, 2.0);
  for (int j = 0; j < 27; ++j) out[j] = 0.0;
  // E[l_aa]: in a, in d, in p
  out[0] = gengamma1_dexpected_a_a_a(d, p, S3);
  out[9] = gengamma1_dexpected_a_a_d(p, S2);
  out[18] = gengamma1_dexpected_a_a_p(d, S2);
  // E[l_dd]: in a (zero), in d, in p
  out[4] = 0.0;
  out[13] = gengamma1_dexpected_d_d_d(d, p, PB2);
  out[22] = gengamma1_dexpected_d_d_p(d, p, PB1, PB2);
  // E[l_pp]: in a (zero), in d, in p
  out[8] = 0.0;
  out[17] = gengamma1_dexpected_p_p_d(d, p, PB0, PB1, PB2);
  out[26] = gengamma1_dexpected_p_p_p(d, p, PB0, PB1, PB2);
  // E[l_ad]: in a, in d (zero), in p (zero)
  out[3] = out[1] = gengamma1_dexpected_a_d_a(S2);
  out[12] = out[10] = 0.0;
  out[21] = out[19] = 0.0;
  // E[l_ap]: in a, in d, in p
  out[6] = out[2] = gengamma1_dexpected_a_p_a(d, p, PB0, S2);
  out[15] = out[11] = gengamma1_dexpected_a_p_d(d, p, PB0, PB1, S1);
  out[24] = out[20] = gengamma1_dexpected_a_p_p(d, p, PB0, PB1, S1);
  // E[l_dp]: in a (zero), in d, in p
  out[7] = out[5] = 0.0;
  out[16] = out[14] = gengamma1_dexpected_d_d_p(d, p, PB1, PB2);
  out[25] = out[23] = gengamma1_dexpected_d_p_p(d, p, PB0, PB1, PB2);
}

} // namespace d7

#endif
