#ifndef D7_PT_LOGNORMAL1_H
#define D7_PT_LOGNORMAL1_H

#include <cmath>
#include "pt_gaussian2.h"
#include "pt_parts.h"

// Lognormal in the mean and the variance of log y: one function per
// component and order for the quantities the scalar registry reads, called
// by the vector kernels (lognormal.cpp) and by the registry
// (d7_ccallable.cpp). With res = log y - mu; the caller passes log y. The
// derivatives of the expected information are the gaussian's by its
// variance (pt_gaussian2.h).

namespace d7 {

inline double lognormal1_score_mu(double res, double s2) {
    return res / s2;
}

inline double lognormal1_score_sigma2(double res, double s2) {
    return (res * res - s2) / (2.0 * s2 * s2);
}

inline double lognormal1_hess_mu_mu(double s2) {
    return -1.0 / s2;
}

inline double lognormal1_hess_sigma2_sigma2(double res, double s2) {
    double s4 = s2 * s2;
    return 0.5 / s4 - (res * res) / (s4 * s2);
}

inline double lognormal1_expected_mu_mu(double s2) {
    return -1.0 / s2;
}

inline double lognormal1_expected_sigma2_sigma2(double s2) {
    return -0.5 / (s2 * s2);
}

inline void lognormal1_score_curv(int k, double y, const double* th,
                                  double* out) {
    double s2 = th[1];
    double res = std::log(y) - th[0];
    if (k == 0) {
        out[0] = lognormal1_score_mu(res, s2);
        out[1] = lognormal1_hess_mu_mu(s2);
    } else {
        out[0] = lognormal1_score_sigma2(res, s2);
        out[1] = lognormal1_hess_sigma2_sigma2(res, s2);
    }
}

inline void lognormal1_info_dinfo(int k, double y, const double* th,
                                  double* out) {
    double s2 = th[1];
    if (k == 0) {
        out[0] = lognormal1_expected_mu_mu(s2);
        out[1] = gaussian2_dexpected_mu_mu_mu();
    } else {
        out[0] = lognormal1_expected_sigma2_sigma2(s2);
        out[1] = gaussian2_dexpected_sigma2_sigma2_sigma2(1.0 / s2);
    }
}

// Mixed component of the observed Hessian, -res / s4, with s4 = s2^2 as
// lognormal_hessian_cpp forms it.
inline double lognormal1_hess_mu_sigma2(double res, double s2) {
    double s4 = s2 * s2;
    return -res / s4;
}

// The routers of the joint registry, one per order (see pt_gaussian2.h). The
// derivatives of the expected information are the gaussian's, as in the
// registry of the scalar routers above.
inline void lognormal1_joint_grad(double y, const double* th, double* out) {
    double s2 = th[1];
    double res = std::log(y) - th[0];
    out[0] = lognormal1_score_mu(res, s2);
    out[1] = lognormal1_score_sigma2(res, s2);
}

inline void lognormal1_joint_hess(double y, const double* th, double* out) {
    double s2 = th[1];
    double res = std::log(y) - th[0];
    out[0] = lognormal1_hess_mu_mu(s2);
    out[3] = lognormal1_hess_sigma2_sigma2(res, s2);
    out[1] = out[2] = lognormal1_hess_mu_sigma2(res, s2);
}

inline void lognormal1_joint_info(double y, const double* th, double* out) {
    double s2 = th[1];
    out[0] = lognormal1_expected_mu_mu(s2);
    out[3] = lognormal1_expected_sigma2_sigma2(s2);
    out[1] = out[2] = 0.0;
}

inline void lognormal1_joint_dinfo(double y, const double* th, double* out) {
    double inv = 1.0 / th[1];
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[0] = gaussian2_dexpected_mu_mu_mu();
    out[4] = gaussian2_dexpected_mu_mu_sigma2(inv);
    out[7] = gaussian2_dexpected_sigma2_sigma2_sigma2(inv);
}

// The parts route (pt_parts.h). The one costly quantity is log y, which every
// order reads; the buffer holds it with the inputs:
//   [0] level                       always
//   [1] y                           observed parts; NaN in the expected parts
//   [2] mu                          always
//   [3] sigma2                      always
//   [4] log y                       level >= 1; NaN in the expected parts
// Orders one and two form res = log y - mu, as the joint routers do; orders
// three and four form r = mu - log y, as lognormal_hd.cpp does
// (lognormal_deriv3_cpp, lognormal_deriv4_cpp). The expected Hessian reads
// only mu and sigma2.
inline void lognormal1_parts(double y, const double* th, int level, double* b) {
    b[0] = level; b[1] = y; b[2] = th[0]; b[3] = th[1];
    b[4] = level >= 1 ? std::log(y) : R_NaN;
}

inline void lognormal1_parts_expected(const double* th, int level, double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = th[1]; b[4] = R_NaN;
}

inline void lognormal1_parts_grad(const double* b, double* out) {
    const double res = b[4] - b[2];
    out[0] = lognormal1_score_mu(res, b[3]);
    out[1] = lognormal1_score_sigma2(res, b[3]);
}

inline void lognormal1_parts_hess(const double* b, double* out) {
    const double res = b[4] - b[2];
    out[0] = lognormal1_hess_mu_mu(b[3]);
    out[3] = lognormal1_hess_sigma2_sigma2(res, b[3]);
    out[1] = out[2] = lognormal1_hess_mu_sigma2(res, b[3]);
}

// from the expected parts
inline void lognormal1_parts_info(const double* b, double* out) {
    out[0] = lognormal1_expected_mu_mu(b[3]);
    out[3] = lognormal1_expected_sigma2_sigma2(b[3]);
    out[1] = out[2] = 0.0;
}

// third order: lognormal_deriv3_cpp, unique components in the sorted order
// mu mu mu, mu mu sigma2, mu sigma2 sigma2, sigma2 sigma2 sigma2
inline void lognormal1_parts_d3(const double* b, double* out) {
    const double s2 = b[3];
    const double r = b[2] - b[4];
    const double s2_2 = s2 * s2, s2_3 = s2_2 * s2, s2_4 = s2_2 * s2_2;
    double u[4];
    u[0] = 0.0;
    u[1] = 1.0 / s2_2;
    u[2] = -2.0 * r / s2_3;
    u[3] = (-s2 + 3.0 * r * r) / s2_4;
    sym_fill(2, 3, u, out);
}

// fourth order: lognormal_deriv4_cpp, unique components in the sorted order
// mu mu mu mu, mu mu mu sigma2, mu mu sigma2 sigma2, mu sigma2 sigma2 sigma2,
// sigma2 sigma2 sigma2 sigma2
inline void lognormal1_parts_d4(const double* b, double* out) {
    const double s2 = b[3];
    const double r = b[2] - b[4];
    const double s2_3 = s2 * s2 * s2, s2_4 = s2_3 * s2, s2_5 = s2_4 * s2;
    double u[5];
    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = -2.0 / s2_3;
    u[3] = 6.0 * r / s2_4;
    u[4] = 3.0 * (s2 - 4.0 * r * r) / s2_5;
    sym_fill(2, 4, u, out);
}

} // namespace d7

#endif
