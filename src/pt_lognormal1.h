#ifndef D7_PT_LOGNORMAL1_H
#define D7_PT_LOGNORMAL1_H

#include <cmath>
#include "pt_gaussian2.h"

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

} // namespace d7

#endif
