#ifndef D7_PT_WEIBULL3_H
#define D7_PT_WEIBULL3_H

#include <Rcpp.h>
#include <cmath>
#include "pt_parts.h"

// The Weibull in its mean m and shape s: one function per component and
// order for the quantities the scalar registry reads, lifted from the
// generated kernels in weibull3.cpp (which records the derivation) and
// called by them and by the registry (d7_ccallable.cpp). The quantities of
// the parameters alone, with x1 = 1 + 1/s, LG = lgamma(x1) and
// P_k = psigamma(x1, k), are formed once per parameter value by
// weibull3_<kernel>_par(); the caller passes that struct, C = log(y/m),
// eC = exp(C s), eL = exp(LG s) and E = eC eL.

namespace d7 {

struct Weibull3GradientPar { double LG, P0, t1, t2; };

inline Weibull3GradientPar weibull3_gradient_par(double m, double s) {
  const double x1 = 1.0 + 1.0 / s;
  const double LG = R::lgammafn(x1);
  const double P0 = R::psigamma(x1, 0.0);
  const double t1 = 1.0/s;
  const double t2 = P0*t1;
  Weibull3GradientPar P;
  P.LG = LG;
  P.P0 = P0;
  P.t1 = t1;
  P.t2 = t2;
  return P;
}

struct Weibull3HessianPar { double LG, P0, P1, t2, t3, t6, t7, t8; };

inline Weibull3HessianPar weibull3_hessian_par(double m, double s) {
  const double x1 = 1.0 + 1.0 / s;
  const double LG = R::lgammafn(x1);
  const double P0 = R::psigamma(x1, 0.0);
  const double P1 = R::psigamma(x1, 1.0);
  const double t2 = LG*s;
  const double t3 = std::exp(t2);
  const double t6 = std::pow(s, -2);
  const double t7 = std::pow(s, -3);
  const double t8 = 1.0/s;
  Weibull3HessianPar P;
  P.LG = LG;
  P.P0 = P0;
  P.P1 = P1;
  P.t2 = t2;
  P.t3 = t3;
  P.t6 = t6;
  P.t7 = t7;
  P.t8 = t8;
  return P;
}

struct Weibull3ExpectedPar { double LG, P0, t0; };

inline Weibull3ExpectedPar weibull3_expected_par(double m, double s) {
  const double x1 = 1.0 + 1.0 / s;
  const double LG = R::lgammafn(x1);
  const double P0 = R::psigamma(x1, 0.0);
  const double t0 = std::pow(s, 2);
  Weibull3ExpectedPar P;
  P.LG = LG;
  P.P0 = P0;
  P.t0 = t0;
  return P;
}

struct Weibull3DexpectedPar { double LG, P0, P1, t0, t1, t2; };

inline Weibull3DexpectedPar weibull3_dexpected_par(double m, double s) {
  const double x1 = 1.0 + 1.0 / s;
  const double LG = R::lgammafn(x1);
  const double P0 = R::psigamma(x1, 0.0);
  const double P1 = R::psigamma(x1, 1.0);
  const double t0 = std::pow(s, 2);
  const double t1 = std::pow(m, -2);
  const double t2 = P1/s;
  Weibull3DexpectedPar P;
  P.LG = LG;
  P.P0 = P0;
  P.P1 = P1;
  P.t0 = t0;
  P.t1 = t1;
  P.t2 = t2;
  return P;
}

inline double weibull3_score_mean(const Weibull3GradientPar& P, double m, double s, double E) {
  return s*(E - 1)/m;
}

inline double weibull3_score_sigma(const Weibull3GradientPar& P, double C, double E) {
  const double LG = P.LG;
  const double t1 = P.t1;
  const double t2 = P.t2;
  return -C*E + C - LG*E + LG + E*t2 + t1 - t2;
}

inline double weibull3_hess_mean_mean(const Weibull3HessianPar& P, double m, double s, double eC) {
  const double t3 = P.t3;
  const double t4 = eC*t3;
  const double t5 = t4 - 1;
  return s*(-s*t4 - t5)/std::pow(m, 2);
}

inline double weibull3_hess_sigma_sigma(const Weibull3HessianPar& P, double C, double eC) {
  const double LG = P.LG;
  const double P0 = P.P0;
  const double P1 = P.P1;
  const double t3 = P.t3;
  const double t6 = P.t6;
  const double t7 = P.t7;
  const double t8 = P.t8;
  const double t4 = eC*t3;
  return -std::pow(C, 2)*t4 - 2*C*LG*t4 + 2*C*P0*eC*t3*t8 - std::pow(LG, 2)*t4 + 2*LG*P0*eC*t3*t8 - std::pow(P0, 2)*t4*t6 - P1*t4*t7 + P1*t7 - t6;
}

inline double weibull3_expected_mean_mean(const Weibull3ExpectedPar& P, double m) {
  const double t0 = P.t0;
  return -1.0*t0/std::pow(m, 2);
}

inline double weibull3_expected_sigma_sigma(const Weibull3ExpectedPar& P) {
  const double P0 = P.P0;
  const double t0 = P.t0;
  return (-1.0*std::pow(P0, 2) + 0.84556867019693428*P0 - 1.8236806608528794)/t0;
}

inline double weibull3_dexpected_mean_mean_mean(const Weibull3DexpectedPar& P, double m) {
  const double t0 = P.t0;
  return 2.0*t0/std::pow(m, 3);
}

inline double weibull3_dexpected_sigma_sigma_sigma(const Weibull3DexpectedPar& P, double s) {
  const double P0 = P.P0;
  const double t2 = P.t2;
  return (2.0*std::pow(P0, 2) + 2.0*P0*t2 - 1.6911373403938686*P0 - 0.84556867019693428*t2 + 3.6473613217057588)/std::pow(s, 3);
}

// Off-diagonal entries of the observed Hessian, the expected information and
// its derivatives, as the kernels of weibull3.cpp compute them inline. The
// caller passes eC = exp(C s) with C = log(y/m) for the observed entry.
inline double weibull3_hess_mean_sigma(const Weibull3HessianPar& P, double C,
                                       double s, double eC, double m) {
  const double t0 = C*s;
  const double t2 = P.t2;
  const double P0 = P.P0;
  const double t3 = P.t3;
  const double t4 = eC*t3;
  const double t5 = t4 - 1;
  return (-P0*t4 + t0*t4 + t2*t4 + t5)/m;
}

inline double weibull3_expected_mean_sigma(const Weibull3ExpectedPar& P,
                                           double m) {
  return (0.42278433509846714 - 1.0*P.P0)/m;
}

// d E[mean mean] / d sigma
inline double weibull3_dexpected_mean_mean_sigma(const Weibull3DexpectedPar& P,
                                                 double s) {
  return -2.0*s*P.t1;
}

// d E[mean sigma] / d mean
inline double weibull3_dexpected_mean_sigma_mean(const Weibull3DexpectedPar& P) {
  return P.t1*(1.0*P.P0 - 0.42278433509846714);
}

// d E[mean sigma] / d sigma
inline double weibull3_dexpected_mean_sigma_sigma(const Weibull3DexpectedPar& P,
                                                  double m) {
  return 1.0*P.P1/(m*P.t0);
}

inline void weibull3_score_curv(int k, double y, const double* th,
                                double* out) {
  const double m = th[0], s = th[1];
  const Weibull3GradientPar Pg = weibull3_gradient_par(m, s);
  const Weibull3HessianPar Ph = weibull3_hessian_par(m, s);
  const double C = std::log(y / m);
  const double eC = std::exp(C*s);
  const double E = eC*std::exp(Pg.LG*s);
  if (k == 0) {
    out[0] = weibull3_score_mean(Pg, m, s, E);
    out[1] = weibull3_hess_mean_mean(Ph, m, s, eC);
  } else {
    out[0] = weibull3_score_sigma(Pg, C, E);
    out[1] = weibull3_hess_sigma_sigma(Ph, C, eC);
  }
}

inline void weibull3_info_dinfo(int k, double y, const double* th,
                                double* out) {
  const double m = th[0], s = th[1];
  const Weibull3ExpectedPar Pe = weibull3_expected_par(m, s);
  const Weibull3DexpectedPar Pd = weibull3_dexpected_par(m, s);
  if (k == 0) {
    out[0] = weibull3_expected_mean_mean(Pe, m);
    out[1] = weibull3_dexpected_mean_mean_mean(Pd, m);
  } else {
    out[0] = weibull3_expected_sigma_sigma(Pe);
    out[1] = weibull3_dexpected_sigma_sigma_sigma(Pd, s);
  }
}

// The routers of the joint registry, one per order: every component at one
// observation, on the parameter scale (mean, sigma). The Hessian and the
// expected information fill a full P x P matrix and the derivative of the
// expected information a P x P x P array, column-major, out[a + P b + P^2 c]
// being the derivative of E[l_ab] in parameter c.
inline void weibull3_joint_grad(double y, const double* th, double* out) {
  const double m = th[0], s = th[1];
  const Weibull3GradientPar Pg = weibull3_gradient_par(m, s);
  const double C = std::log(y / m);
  const double eC = std::exp(C*s);
  const double E = eC*std::exp(Pg.LG*s);
  out[0] = weibull3_score_mean(Pg, m, s, E);
  out[1] = weibull3_score_sigma(Pg, C, E);
}

inline void weibull3_joint_hess(double y, const double* th, double* out) {
  const double m = th[0], s = th[1];
  const Weibull3HessianPar Ph = weibull3_hessian_par(m, s);
  const double C = std::log(y / m);
  const double eC = std::exp(C*s);
  out[0] = weibull3_hess_mean_mean(Ph, m, s, eC);
  out[3] = weibull3_hess_sigma_sigma(Ph, C, eC);
  out[1] = out[2] = weibull3_hess_mean_sigma(Ph, C, s, eC, m);
}

inline void weibull3_joint_info(double y, const double* th, double* out) {
  const double m = th[0], s = th[1];
  const Weibull3ExpectedPar Pe = weibull3_expected_par(m, s);
  out[0] = weibull3_expected_mean_mean(Pe, m);
  out[3] = weibull3_expected_sigma_sigma(Pe);
  out[1] = out[2] = weibull3_expected_mean_sigma(Pe, m);
}

inline void weibull3_joint_dinfo(double y, const double* th, double* out) {
  const double m = th[0], s = th[1];
  const Weibull3DexpectedPar Pd = weibull3_dexpected_par(m, s);
  for (int j = 0; j < 8; ++j) out[j] = 0.0;
  out[0] = weibull3_dexpected_mean_mean_mean(Pd, m);
  out[4] = weibull3_dexpected_mean_mean_sigma(Pd, s);
  out[1] = out[2] = weibull3_dexpected_mean_sigma_mean(Pd);
  out[3] = 0.0;
  out[5] = out[6] = weibull3_dexpected_mean_sigma_sigma(Pd, m);
  out[7] = weibull3_dexpected_sigma_sigma_sigma(Pd, s);
}

// The parts route (pt_parts.h). The buffer, with the lowest level that fills
// each entry in brackets. x1 = 1 + 1/s, and P_k = psigamma(x1, k).
//   [0] level
//   [1] y                         (NaN for the expected parts)
//   [2] m = mean
//   [3] s = sigma
//   [4] LG = lgamma(x1)           (1; NaN for the expected parts, which do not read it)
//   [5] P0                        (1; 2 for the expected parts)
//   [6] P1                        (2; 3 for the expected parts)
//   [7] P2                        (3)
//   [8] P3                        (4)
//   [9] t3 = exp(LG s)            (1)
//   [10] C = log(y / m)           (1)
//   [11] eC = exp(C s)            (1)
//   [12] E = eC t3                (1)
// The entries of the quantities that only the expected information reads at
// level 2 (P0) and its derivatives (P1, level 3) follow the same rule.

inline void weibull3_parts(double y, const double* th, int level, double* b) {
    const double m = th[0], s = th[1];
    const double x1 = 1.0 + 1.0 / s;
    b[0] = level; b[1] = y; b[2] = m; b[3] = s;
    b[4] = level >= 1 ? R::lgammafn(x1) : R_NaN;
    b[5] = level >= 1 ? R::psigamma(x1, 0.0) : R_NaN;
    b[6] = level >= 2 ? R::psigamma(x1, 1.0) : R_NaN;
    b[7] = level >= 3 ? R::psigamma(x1, 2.0) : R_NaN;
    b[8] = level >= 4 ? R::psigamma(x1, 3.0) : R_NaN;
    b[9] = level >= 1 ? std::exp(b[4] * s) : R_NaN;
    b[10] = level >= 1 ? std::log(y / m) : R_NaN;
    b[11] = level >= 1 ? std::exp(b[10] * s) : R_NaN;
    b[12] = level >= 1 ? b[11] * b[9] : R_NaN;
}

inline void weibull3_parts_expected(const double* th, int level, double* b) {
    const double m = th[0], s = th[1];
    const double x1 = 1.0 + 1.0 / s;
    b[0] = level; b[1] = R_NaN; b[2] = m; b[3] = s;
    b[4] = R_NaN;
    b[5] = level >= 2 ? R::psigamma(x1, 0.0) : R_NaN;
    b[6] = level >= 3 ? R::psigamma(x1, 1.0) : R_NaN;
    b[7] = R_NaN; b[8] = R_NaN;
    b[9] = R_NaN; b[10] = R_NaN; b[11] = R_NaN; b[12] = R_NaN;
}

inline void weibull3_parts_grad(const double* b, double* out) {
    const double m = b[2], s = b[3], LG = b[4], P0 = b[5], C = b[10], E = b[12];
    out[0] = s*(E - 1)/m;
    const double t1 = 1.0/s;
    const double t2 = P0*t1;
    out[1] = -C*E + C - LG*E + LG + E*t2 + t1 - t2;
}

inline void weibull3_parts_hess(const double* b, double* out) {
    const double m = b[2], s = b[3], LG = b[4], P0 = b[5], P1 = b[6];
    const double t3 = b[9], C = b[10], eC = b[11];
    const double t4 = eC*t3;
    const double t5 = t4 - 1;
    out[0] = s*(-s*t4 - t5)/std::pow(m, 2);
    const double t6 = std::pow(s, -2), t7 = std::pow(s, -3), t8 = 1.0/s;
    out[3] = -std::pow(C, 2)*t4 - 2*C*LG*t4 + 2*C*P0*eC*t3*t8 - std::pow(LG, 2)*t4 + 2*LG*P0*eC*t3*t8 - std::pow(P0, 2)*t4*t6 - P1*t4*t7 + P1*t7 - t6;
    const double t0 = C*s, t2 = LG*s;
    out[1] = out[2] = (-P0*t4 + t0*t4 + t2*t4 + t5)/m;
}

// from the expected parts
inline void weibull3_parts_info(const double* b, double* out) {
    const double m = b[2], s = b[3], P0 = b[5];
    const double t0 = std::pow(s, 2);
    out[0] = -1.0*t0/std::pow(m, 2);
    out[3] = (-1.0*std::pow(P0, 2) + 0.84556867019693428*P0 - 1.8236806608528794)/t0;
    out[1] = out[2] = (0.42278433509846714 - 1.0*P0)/m;
}

// the unique components of the third derivative in the order of
// weibull3_deriv3_cpp: mean mean mean, mean mean sigma, mean sigma sigma,
// sigma sigma sigma
inline void weibull3_parts_d3(const double* b, double* out) {
    const double m = b[2], s = b[3], LG = b[4], P0 = b[5], P1 = b[6], P2 = b[7];
    const double t3 = b[9], C = b[10], t1 = b[11];
    const double t2 = LG*s;
    const double t7 = std::pow(s, 2);
    const double t10 = 2*LG;
    const double t12 = std::pow(LG, 2);
    const double t13 = 1.0/t7;
    const double t14 = 1.0/s;
    const double t15 = std::pow(P0, 2);
    const double t16 = std::pow(s, -3);
    const double t17 = std::pow(s, -5);
    const double t18 = std::pow(s, -4);
    const double t19 = 3*P1;
    const double t23 = t13*t15;
    const double t0 = C*s;
    const double t4 = t1*t3;
    const double t5 = 2*t4;
    const double t6 = 3*t4;
    const double t8 = t4*t7;
    const double t9 = 2*C;
    const double t11 = std::pow(C, 2);
    const double t20 = C*t6;
    const double t21 = LG*t6;
    const double t22 = t16*t19*t4;
    double uq[4];
    uq[0] = s*(s*t6 + t5 + t8 - 2)/std::pow(m, 3);
    uq[1] = (-C*t8 - LG*t8 + P0*s*t1*t3 + P0*t1*t3 - s*t5 - t0*t4 - t2*t4 - t4 + 1)/std::pow(m, 2);
    uq[2] = t4*(-P0*t10 - 2*P0*t14 - P0*t9 + P1*t13 + s*t11 + s*t12 + t0*t10 + t10 + t14*t15 + t9)/m;
    uq[3] = -std::pow(C, 3)*t4 + 6*C*LG*P0*t1*t14*t3 - C*t22 - std::pow(LG, 3)*t4 - LG*t22 + std::pow(P0, 3)*t1*t16*t3 + 3*P0*P1*t1*t18*t3 + 3*P0*t1*t11*t14*t3 + 3*P0*t1*t12*t14*t3 + 3*P1*t1*t18*t3 + P2*t1*t17*t3 - P2*t17 - t11*t21 - t12*t20 + 2*t16 - t18*t19 - t20*t23 - t21*t23;
    sym_fill(2, 3, uq, out);
}

// the unique components of the fourth derivative, in the order of
// weibull3_deriv4_cpp
inline void weibull3_parts_d4(const double* b, double* out) {
    const double m = b[2], s = b[3], LG = b[4], P0 = b[5], P1 = b[6], P2 = b[7], P3 = b[8];
    const double t3 = b[9], C = b[10], t1 = b[11];
    const double t2 = LG*s;
    const double t7 = std::pow(s, 3);
    const double t9 = std::pow(s, 2);
    const double t11 = 3*P0;
    const double t15 = 2*LG;
    const double t16 = std::pow(P0, 2);
    const double t18 = std::pow(LG, 2);
    const double t19 = 1.0/t9;
    const double t20 = P1*t19;
    const double t21 = 1.0/s;
    const double t22 = t16*t21;
    const double t26 = 3*t18;
    const double t28 = std::pow(LG, 3);
    const double t29 = std::pow(s, -4);
    const double t30 = P0*t21;
    const double t32 = 3*LG;
    const double t33 = 1.0/t7;
    const double t34 = P1*t33;
    const double t35 = t16*t19;
    const double t36 = std::pow(P0, 3);
    const double t37 = std::pow(s, -5);
    const double t38 = std::pow(s, -6);
    const double t39 = std::pow(s, -7);
    const double t41 = P1*t37;
    const double t43 = P2*t38;
    const double t0 = C*s;
    const double t4 = t1*t3;
    const double t5 = 6*t4;
    const double t6 = s*t4;
    const double t8 = t4*t7;
    const double t10 = 2*t4;
    const double t12 = t4*t9;
    const double t13 = 3*t12;
    const double t14 = 2*C;
    const double t17 = std::pow(C, 2);
    const double t23 = 6*C;
    const double t24 = LG*t23;
    const double t25 = 3*t17;
    const double t27 = std::pow(C, 3);
    const double t31 = 3*C;
    const double t40 = 4*t4;
    const double t42 = 12*t4*t41;
    const double t44 = t17*t5;
    const double t45 = 12*C*LG*t4;
    const double t46 = t18*t5;
    double uq[5];
    uq[0] = s*(-t5*t9 - t5 - 11*t6 - t8 + 6)/std::pow(m, 4);
    uq[1] = (C*t13 + C*t8 + LG*t13 + LG*t8 - P0*t10 - P0*t12 + s*t5 + t0*t10 + t10*t2 + t10 - t11*t6 + t13 - 2)/std::pow(m, 3);
    uq[2] = t4*(2*C*P0*s + 2*C*P0 + 2*LG*P0*s + 2*LG*P0 - LG*t14*t9 + 2*P0*t21 + 4*P0 - P1*t21 - s*t17 - s*t18 - t0*t15 - 4*t0 - t14 - t15 - t16 - t17*t9 - t18*t9 - 4*t2 - t20 - t22 - 2)/std::pow(m, 2);
    uq[3] = t4*(-6*LG*t30 - P0*t24 - P0*t25 - P0*t26 - P2*t29 + s*t27 + s*t28 + t0*t26 - t11*t34 - t19*t36 + t2*t25 + t20*t31 + t20*t32 + t22*t31 + t22*t32 - t23*t30 + t24 + t25 + t26 + 3*t35)/m;
    uq[4] = -std::pow(C, 4)*t4 + 12*C*P0*P1*t1*t29*t3 + 12*C*P0*t1*t18*t21*t3 + 12*C*P1*t1*t29*t3 + 4*C*P2*t1*t3*t37 + 4*C*t1*t3*t33*t36 - C*t28*t40 - std::pow(LG, 4)*t4 + 12*LG*P0*P1*t1*t29*t3 + 12*LG*P0*t1*t17*t21*t3 + 12*LG*P1*t1*t29*t3 + 4*LG*P2*t1*t3*t37 + 4*LG*t1*t3*t33*t36 - LG*t27*t40 - std::pow(P0, 4)*t29*t4 + 4*P0*t1*t21*t27*t3 + 4*P0*t1*t21*t28*t3 - 4*P0*t4*t43 - P0*t42 - 3*std::pow(P1, 2)*t38*t4 + 12*P1*t37 + 8*P2*t38 - P3*t39*t4 + P3*t39 - t16*t41*t5 - t18*t44 - 6*t29 - t34*t44 - t34*t45 - t34*t46 - t35*t44 - t35*t45 - t35*t46 - 8*t4*t43 - t42;
    sym_fill(2, 4, uq, out);
}

} // namespace d7

#endif
