#ifndef D7_PT_WEIBULL3_H
#define D7_PT_WEIBULL3_H

#include <Rcpp.h>
#include <cmath>

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

} // namespace d7

#endif
