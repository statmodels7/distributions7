#ifndef D7_PT_GUMBEL_H
#define D7_PT_GUMBEL_H

#include <Rcpp.h>
#include <cmath>
#include "pt_constants.h"

// The Gumbel (maximum) in location and scale, with z = (y - mu)/sigma and
// w = exp(-z): one function per component and order for the quantities the
// scalar registry reads, called by the vector kernels (gumbel.cpp,
// dexpected_kernels.cpp) and by the registry (d7_ccallable.cpp). The caller
// passes z and w.

namespace d7 {

inline double gumbel_score_mu(double w, double s) {
    return (1.0 - w) / s;
}

inline double gumbel_score_sigma(double z, double w, double s) {
    return (z * (1.0 - w) - 1.0) / s;
}

inline double gumbel_hess_mu_mu(double w, double s) {
    double s2 = s * s;
    return -w / s2;
}

inline double gumbel_hess_sigma_sigma(double z, double w, double s) {
    double s2 = s * s;
    return (1.0 - 2.0 * z + 2.0 * z * w - z * z * w) / s2;
}

// E_mm = -1/s^2, E_ss = -((1 - gamma)^2 + pi^2/6)/s^2, E_ms = (1 - gamma)/s^2
inline double gumbel_expected_mu_mu(double s) {
    double s2 = s * s;
    return -1.0 / s2;
}

inline double gumbel_expected_sigma_sigma(double s) {
    double s2 = s * s;
    return -kGumbelInfo / s2;
}

// nothing moves with the location
inline double gumbel_dexpected_mu_mu_mu() {
    return 0.0;
}

inline double gumbel_dexpected_sigma_sigma_sigma(double s) {
    const double k = -kGumbelInfo;
    double u = 1.0 / s, u3 = u * u * u;
    return -2.0 * k * u3;
}

inline void gumbel_score_curv(int k, double y, const double* th,
                              double* out) {
    double s = th[1];
    double z = (y - th[0]) / s, w = std::exp(-z);
    if (k == 0) {
        out[0] = gumbel_score_mu(w, s);
        out[1] = gumbel_hess_mu_mu(w, s);
    } else {
        out[0] = gumbel_score_sigma(z, w, s);
        out[1] = gumbel_hess_sigma_sigma(z, w, s);
    }
}

inline void gumbel_info_dinfo(int k, double y, const double* th,
                              double* out) {
    double s = th[1];
    if (k == 0) {
        out[0] = gumbel_expected_mu_mu(s);
        out[1] = gumbel_dexpected_mu_mu_mu();
    } else {
        out[0] = gumbel_expected_sigma_sigma(s);
        out[1] = gumbel_dexpected_sigma_sigma_sigma(s);
    }
}

// The mixed components, written as gumbel.cpp writes them. The expected
// mu-sigma entry is (1 - gamma)/sigma^2, positive; its derivative in sigma
// is -2 k u^3 with k = 1 - gamma, and the derivative of E_mm in sigma is
// -2 k u^3 with k = -1.
inline double gumbel_hess_mu_sigma(double z, double w, double s) {
    return -(1.0 - w + z * w) / (s * s);
}

inline double gumbel_expected_mu_sigma(double s) {
    return (1.0 - kEulerGamma) / (s * s);
}

inline double gumbel_dexpected_mu_mu_sigma(double s) {
    const double k = -1.0;
    double u = 1.0 / s, u3 = u * u * u;
    return -2.0 * k * u3;
}

inline double gumbel_dexpected_mu_sigma_sigma(double s) {
    const double k = 1.0 - kEulerGamma;
    double u = 1.0 / s, u3 = u * u * u;
    return -2.0 * k * u3;
}

// The routers of the joint registry, one per order, in (mu, sigma); the
// layout is the one pt_cauchy.h describes.
inline void gumbel_joint_grad(double y, const double* th, double* out) {
    double s = th[1];
    double z = (y - th[0]) / s, w = std::exp(-z);
    out[0] = gumbel_score_mu(w, s);
    out[1] = gumbel_score_sigma(z, w, s);
}

inline void gumbel_joint_hess(double y, const double* th, double* out) {
    double s = th[1];
    double z = (y - th[0]) / s, w = std::exp(-z);
    out[0] = gumbel_hess_mu_mu(w, s);
    out[3] = gumbel_hess_sigma_sigma(z, w, s);
    out[1] = out[2] = gumbel_hess_mu_sigma(z, w, s);
}

inline void gumbel_joint_info(double y, const double* th, double* out) {
    double s = th[1];
    out[0] = gumbel_expected_mu_mu(s);
    out[3] = gumbel_expected_sigma_sigma(s);
    out[1] = out[2] = gumbel_expected_mu_sigma(s);
}

inline void gumbel_joint_dinfo(double y, const double* th, double* out) {
    double s = th[1];
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[4] = gumbel_dexpected_mu_mu_sigma(s);
    out[7] = gumbel_dexpected_sigma_sigma_sigma(s);
    out[5] = out[6] = gumbel_dexpected_mu_sigma_sigma(s);
}

} // namespace d7

#endif
