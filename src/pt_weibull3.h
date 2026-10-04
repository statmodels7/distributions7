#ifndef D7_PT_WEIBULL3_H
#define D7_PT_WEIBULL3_H

#include <Rcpp.h>
#include <cmath>

// The Weibull in its mean m and shape s: one function per component and
// order for the quantities the scalar registry reads, lifted from the
// generated kernels in weibull3.cpp (which records the derivation) and
// called by them and by the registry (d7_ccallable.cpp). With
// x1 = 1 + 1/s and C = log(y/m) the caller passes LG = lgamma(x1),
// P_k = psigamma(x1, k), eC = exp(C s), eL = exp(LG s) and E = eC eL.

namespace d7 {

inline double weibull3_score_mean(double m, double s, double E) {
  return s*(E - 1)/m;
}

inline double weibull3_score_sigma(double s, double C, double LG, double P0, double E) {
  const double t1 = 1.0/s;
  const double t2 = P0*t1;
  return -C*E + C - LG*E + LG + E*t2 + t1 - t2;
}

inline double weibull3_hess_mean_mean(double m, double s, double eC, double eL) {
  const double t4 = eC*eL;
  const double t5 = t4 - 1;
  return s*(-s*t4 - t5)/std::pow(m, 2);
}

inline double weibull3_hess_sigma_sigma(double s, double C, double LG, double P0, double P1, double eC, double eL) {
  const double t6 = std::pow(s, -2);
  const double t7 = std::pow(s, -3);
  const double t8 = 1.0/s;
  const double t4 = eC*eL;
  return -std::pow(C, 2)*t4 - 2*C*LG*t4 + 2*C*P0*eC*eL*t8 - std::pow(LG, 2)*t4 + 2*LG*P0*eC*eL*t8 - std::pow(P0, 2)*t4*t6 - P1*t4*t7 + P1*t7 - t6;
}

inline double weibull3_expected_mean_mean(double m, double s) {
  const double t0 = std::pow(s, 2);
  return -1.0*t0/std::pow(m, 2);
}

inline double weibull3_expected_sigma_sigma(double s, double P0) {
  const double t0 = std::pow(s, 2);
  return (-1.0*std::pow(P0, 2) + 0.84556867019693428*P0 - 1.8236806608528794)/t0;
}

inline double weibull3_dexpected_mean_mean_mean(double m, double s) {
  const double t0 = std::pow(s, 2);
  return 2.0*t0/std::pow(m, 3);
}

inline double weibull3_dexpected_sigma_sigma_sigma(double s, double P0, double P1) {
  const double t2 = P1/s;
  return (2.0*std::pow(P0, 2) + 2.0*P0*t2 - 1.6911373403938686*P0 - 0.84556867019693428*t2 + 3.6473613217057588)/std::pow(s, 3);
}

inline void weibull3_score_curv(int k, double y, const double* th,
                                double* out) {
  const double m = th[0], s = th[1];
  const double x1 = 1.0 + 1.0 / s;
  const double LG = R::lgammafn(x1);
  const double C = std::log(y / m);
  const double eC = std::exp(C*s), eL = std::exp(LG*s);
  if (k == 0) {
    out[0] = weibull3_score_mean(m, s, eC*eL);
    out[1] = weibull3_hess_mean_mean(m, s, eC, eL);
  } else {
    const double P0 = R::psigamma(x1, 0.0), P1 = R::psigamma(x1, 1.0);
    out[0] = weibull3_score_sigma(s, C, LG, P0, eC*eL);
    out[1] = weibull3_hess_sigma_sigma(s, C, LG, P0, P1, eC, eL);
  }
}

inline void weibull3_info_dinfo(int k, double y, const double* th,
                                double* out) {
  const double m = th[0], s = th[1];
  if (k == 0) {
    out[0] = weibull3_expected_mean_mean(m, s);
    out[1] = weibull3_dexpected_mean_mean_mean(m, s);
  } else {
    const double x1 = 1.0 + 1.0 / s;
    const double P0 = R::psigamma(x1, 0.0), P1 = R::psigamma(x1, 1.0);
    out[0] = weibull3_expected_sigma_sigma(s, P0);
    out[1] = weibull3_dexpected_sigma_sigma_sigma(s, P0, P1);
  }
}

} // namespace d7

#endif
