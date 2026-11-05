#ifndef D7_PT_GUMBEL_H
#define D7_PT_GUMBEL_H

#include <Rcpp.h>
#include <cmath>
#include "pt_constants.h"
#include "pt_parts.h"

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

// The parts route (pt_parts.h). The buffer, with z = (y - mu)/sigma and
// w = exp(-z) as the joint routers form them:
//
//   [0] level
//   [1] y                           NaN in the expected buffer
//   [2] mu
//   [3] sigma
//   [4] z = (y - mu)/sigma          level 1
//   [5] w = exp(-z)                 level 1
//
// The expected information reads sigma only.
inline void gumbel_parts(double y, const double* th, int level, double* b) {
    const double m = th[0], s = th[1];
    b[0] = level; b[1] = y; b[2] = m; b[3] = s;
    b[4] = level >= 1 ? (y - m) / s : R_NaN;
    b[5] = level >= 1 ? std::exp(-b[4]) : R_NaN;
}

inline void gumbel_parts_expected(const double* th, int level, double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = th[1];
    b[4] = R_NaN; b[5] = R_NaN;
}

inline void gumbel_parts_grad(const double* b, double* out) {
    out[0] = gumbel_score_mu(b[5], b[3]);
    out[1] = gumbel_score_sigma(b[4], b[5], b[3]);
}

inline void gumbel_parts_hess(const double* b, double* out) {
    out[0] = gumbel_hess_mu_mu(b[5], b[3]);
    out[3] = gumbel_hess_sigma_sigma(b[4], b[5], b[3]);
    out[1] = out[2] = gumbel_hess_mu_sigma(b[4], b[5], b[3]);
}

// from the expected parts
inline void gumbel_parts_info(const double* b, double* out) {
    out[0] = gumbel_expected_mu_mu(b[3]);
    out[3] = gumbel_expected_sigma_sigma(b[3]);
    out[1] = out[2] = gumbel_expected_mu_sigma(b[3]);
}

// third order: gumbel_deriv3_cpp (gumbel_hd.cpp), unique components in the
// sorted order mu mu mu, mu mu sigma, mu sigma sigma, sigma sigma sigma
inline void gumbel_parts_d3(const double* b, double* out) {
    const double z = b[4], w = b[5], s = b[3];
    const double s3 = s * s * s;
    const double z2 = z * z;
    double u[4];
    u[0] = -w / s3;
    u[1] = w * (2.0 - z) / s3;
    u[2] = (2.0 - 2.0 * w + 4.0 * z * w - z2 * w) / s3;
    u[3] = (-2.0 + 6.0 * z - 6.0 * z * w + 6.0 * z2 * w - z2 * z * w) / s3;
    sym_fill(2, 3, u, out);
}

// fourth order: gumbel_deriv4_cpp (gumbel_hd.cpp), unique components in the
// sorted order mu mu mu mu, mu mu mu sigma, mu mu sigma sigma, mu sigma sigma
// sigma, sigma sigma sigma sigma
inline void gumbel_parts_d4(const double* b, double* out) {
    const double z = b[4], w = b[5], s = b[3];
    const double s4 = s * s * s * s;
    const double z2 = z * z, z3 = z2 * z, z4 = z2 * z2;
    double u[5];
    u[0] = -w / s4;
    u[1] = w * (3.0 - z) / s4;
    u[2] = -w * (6.0 - 6.0 * z + z2) / s4;
    u[3] = (-6.0 + 6.0 * w - 18.0 * z * w + 9.0 * z2 * w - z3 * w) / s4;
    u[4] = (6.0 - 24.0 * z + 24.0 * z * w - 36.0 * z2 * w + 12.0 * z3 * w
            - z4 * w) / s4;
    sym_fill(2, 4, u, out);
}

} // namespace d7

#endif
