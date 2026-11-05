#ifndef D7_PT_WEIBULL1_H
#define D7_PT_WEIBULL1_H

#include <Rcpp.h>
#include <cmath>
#include "pt_constants.h"
#include "pt_parts.h"

// The Weibull in its scale mu and shape sigma, with z = y/mu, lz = log z
// and u = z^sigma = exp(sigma lz): one function per component and order
// for the quantities the scalar registry reads, called by the vector
// kernels (weibull1.cpp, dexpected_kernels.cpp) and by the registry
// (d7_ccallable.cpp). The caller passes lz and u.

namespace d7 {

inline double weibull1_score_mu(double u, double m, double s) {
    return s * (u - 1.0) / m;
}

inline double weibull1_score_sigma(double lz, double u, double s) {
    return 1.0 / s + (1.0 - u) * lz;
}

inline double weibull1_hess_mu_mu(double u, double m, double s) {
    return s * (1.0 - (1.0 + s) * u) / (m * m);
}

inline double weibull1_hess_sigma_sigma(double lz, double u, double s) {
    return -1.0 / (s * s) - u * lz * lz;
}

// E_mm = -s^2/m^2, E_ss = -((1 - gamma)^2 + pi^2/6)/s^2, E_ms = (1 - gamma)/m
inline double weibull1_expected_mu_mu(double m, double s) {
    return -(s * s) / (m * m);
}

inline double weibull1_expected_sigma_sigma(double s) {
    return -kGumbelInfo / (s * s);
}

inline double weibull1_dexpected_mu_mu_mu(double m, double s) {
    double im = 1.0 / m, im2 = im * im, im3 = im2 * im;
    return 2.0 * s * s * im3;
}

inline double weibull1_dexpected_sigma_sigma_sigma(double s) {
    double is = 1.0 / s, is3 = is * is * is;
    return 2.0 * kGumbelInfo * is3;
}

// The off-diagonal components, written from the expressions of the vector
// kernels (weibull1.cpp, dexpected_kernels.cpp) with the same operations in
// the same order. lz = log(y/mu) and u = exp(sigma lz), as in the routers.
// H_ms = (u - 1 + sigma u lz)/mu
inline double weibull1_hess_mu_sigma(double lz, double u, double m, double s) {
    return (u - 1.0 + s * u * lz) / m;
}

// E_ms = (1 - gamma)/mu
inline double weibull1_expected_mu_sigma(double m) {
    return (1.0 - d7::kEulerGamma) / m;
}

// d E_mm / d sigma = -2 s/mu^2
inline double weibull1_dexpected_mu_mu_sigma(double m, double s) {
    double im = 1.0 / m, im2 = im * im;
    return -2.0 * s * im2;
}

// d E_ms / d mu = -(1 - gamma)/mu^2
inline double weibull1_dexpected_mu_sigma_mu(double m) {
    double im = 1.0 / m, im2 = im * im;
    const double k = 1.0 - d7::kEulerGamma;
    return -k * im2;
}

inline void weibull1_score_curv(int k, double y, const double* th,
                                double* out) {
    double m = th[0], s = th[1];
    double lz = std::log(y / m), u = std::exp(s * lz);
    if (k == 0) {
        out[0] = weibull1_score_mu(u, m, s);
        out[1] = weibull1_hess_mu_mu(u, m, s);
    } else {
        out[0] = weibull1_score_sigma(lz, u, s);
        out[1] = weibull1_hess_sigma_sigma(lz, u, s);
    }
}

inline void weibull1_info_dinfo(int k, double y, const double* th,
                                double* out) {
    double m = th[0], s = th[1];
    if (k == 0) {
        out[0] = weibull1_expected_mu_mu(m, s);
        out[1] = weibull1_dexpected_mu_mu_mu(m, s);
    } else {
        out[0] = weibull1_expected_sigma_sigma(s);
        out[1] = weibull1_dexpected_sigma_sigma_sigma(s);
    }
}

// The routers of the joint registry, one per order: every component at one
// observation, on the parameter scale, with th = (mu, sigma). The Hessian and
// the expected information fill a full 2 x 2 matrix and the derivative of the
// expected information a 2 x 2 x 2 array, column-major, out[a + 2 b + 4 c]
// being the derivative of E[l_ab] in parameter c. Every entry is written.
inline void weibull1_joint_grad(double y, const double* th, double* out) {
    double m = th[0], s = th[1];
    double lz = std::log(y / m), u = std::exp(s * lz);
    out[0] = weibull1_score_mu(u, m, s);
    out[1] = weibull1_score_sigma(lz, u, s);
}

inline void weibull1_joint_hess(double y, const double* th, double* out) {
    double m = th[0], s = th[1];
    double lz = std::log(y / m), u = std::exp(s * lz);
    out[0] = weibull1_hess_mu_mu(u, m, s);
    out[3] = weibull1_hess_sigma_sigma(lz, u, s);
    out[1] = out[2] = weibull1_hess_mu_sigma(lz, u, m, s);
}

inline void weibull1_joint_info(double y, const double* th, double* out) {
    double m = th[0], s = th[1];
    out[0] = weibull1_expected_mu_mu(m, s);
    out[3] = weibull1_expected_sigma_sigma(s);
    out[1] = out[2] = weibull1_expected_mu_sigma(m);
}

inline void weibull1_joint_dinfo(double y, const double* th, double* out) {
    double m = th[0], s = th[1];
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[0] = weibull1_dexpected_mu_mu_mu(m, s);
    out[4] = weibull1_dexpected_mu_mu_sigma(m, s);
    out[1] = out[2] = weibull1_dexpected_mu_sigma_mu(m);
    out[7] = weibull1_dexpected_sigma_sigma_sigma(s);
}

// The parts route (pt_parts.h). The buffer, with the lowest level that fills
// each entry in brackets:
//   [0] level
//   [1] y                        (NaN for the expected parts)
//   [2] mu
//   [3] sigma
//   [4] lz = log(y / mu)         (1)
//   [5] u = exp(sigma lz)        (1)
// The expected parts need no special function: the expected information
// reads the constants kGumbelInfo and kEulerGamma, so [4] and [5] stay NaN.

inline void weibull1_parts(double y, const double* th, int level, double* b) {
    const double m = th[0], s = th[1];
    b[0] = level; b[1] = y; b[2] = m; b[3] = s;
    b[4] = level >= 1 ? std::log(y / m) : R_NaN;
    b[5] = level >= 1 ? std::exp(s * b[4]) : R_NaN;
}

inline void weibull1_parts_expected(const double* th, int level, double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = th[1];
    b[4] = R_NaN; b[5] = R_NaN;
}

inline void weibull1_parts_grad(const double* b, double* out) {
    out[0] = weibull1_score_mu(b[5], b[2], b[3]);
    out[1] = weibull1_score_sigma(b[4], b[5], b[3]);
}

inline void weibull1_parts_hess(const double* b, double* out) {
    out[0] = weibull1_hess_mu_mu(b[5], b[2], b[3]);
    out[3] = weibull1_hess_sigma_sigma(b[4], b[5], b[3]);
    out[1] = out[2] = weibull1_hess_mu_sigma(b[4], b[5], b[2], b[3]);
}

// from the expected parts
inline void weibull1_parts_info(const double* b, double* out) {
    out[0] = weibull1_expected_mu_mu(b[2], b[3]);
    out[3] = weibull1_expected_sigma_sigma(b[3]);
    out[1] = out[2] = weibull1_expected_mu_sigma(b[2]);
}

// the unique components of the third derivative in the order of
// weibull_deriv3_cpp: mu mu mu, mu mu sigma, mu sigma sigma, sigma sigma sigma
inline void weibull1_parts_d3(const double* b, double* out) {
    const double m = b[2], s = b[3], L = b[4], u = b[5];
    double m2 = m * m, m3 = m2 * m, s3 = s * s * s;
    double uq[4];
    uq[0] = s * (-2.0 + (1.0 + s) * (2.0 + s) * u) / m3;
    uq[1] = (1.0 - (1.0 + 2.0 * s) * u - s * (1.0 + s) * L * u) / m2;
    uq[2] = L * u * (2.0 + s * L) / m;
    uq[3] = 2.0 / s3 - u * L * L * L;
    sym_fill(2, 3, uq, out);
}

// the unique components of the fourth derivative, in the order of
// weibull_deriv4_cpp
inline void weibull1_parts_d4(const double* b, double* out) {
    const double m = b[2], s = b[3], L = b[4], u = b[5];
    double m2 = m * m, m3 = m2 * m, m4 = m2 * m2, s4 = s * s * s * s;
    double L2 = L * L;
    double uq[5];
    uq[0] = s * (6.0 - (1.0 + s) * (2.0 + s) * (3.0 + s) * u) / m4;
    uq[1] = (-2.0 + (2.0 + 6.0 * s + 3.0 * s * s) * u
             + s * (1.0 + s) * (2.0 + s) * L * u) / m3;
    uq[2] = -u * (2.0 + 2.0 * (1.0 + 2.0 * s) * L
                  + s * (1.0 + s) * L2) / m2;
    uq[3] = L2 * u * (3.0 + s * L) / m;
    uq[4] = -6.0 / s4 - u * L2 * L2;
    sym_fill(2, 4, uq, out);
}

} // namespace d7

#endif
