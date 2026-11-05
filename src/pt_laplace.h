#ifndef D7_PT_LAPLACE_H
#define D7_PT_LAPLACE_H

#include <Rcpp.h>
#include <cmath>
#include "pt_parts.h"

// The Laplace in location mu and scale b, l = -log(2b) - |y - mu|/b: one
// function per component and order for the quantities the scalar registry
// reads, called by the vector kernels (laplace.cpp) and by the registry
// (d7_ccallable.cpp). With r = y - mu the caller passes sgn(r), which is 0
// at r = 0, one point of the subdifferential, and |r|.

namespace d7 {

inline double laplace_sign(double r) {
    return (r > 0.0) - (r < 0.0);
}

inline double laplace_score_mu(double sg, double b) {
    return sg / b;
}

inline double laplace_score_sigma(double ar, double b) {
    return (ar / b - 1.0) / b;
}

// zero away from the kink, where the density is not differentiable twice
inline double laplace_hess_mu_mu() {
    return 0.0;
}

// (b - 2|r|)/b^3 written as (1 - 2|r|/b)/b^2 with 1/b formed once, so that
// a tiny scale does not underflow b^3 where the value is finite
inline double laplace_hess_sigma_sigma(double ar, double b) {
    const double ib = 1.0 / b;
    return (1.0 - 2.0 * ar * ib) * ib * ib;
}

inline double laplace_expected_mu_mu(double b) {
    const double ib = 1.0 / b;
    return -ib * ib;
}

inline double laplace_expected_sigma_sigma(double b) {
    const double ib = 1.0 / b;
    return -ib * ib;
}

// E_mm = E_bb = -1/b^2 do not depend on mu
inline double laplace_dexpected_mu_mu_mu() {
    return 0.0;
}

inline double laplace_dexpected_sigma_sigma_sigma(double b) {
    const double ib = 1.0 / b;
    return 2.0 * ib * ib * ib;
}

inline void laplace_score_curv(int k, double y, const double* th,
                               double* out) {
    double r = y - th[0], b = th[1];
    if (k == 0) {
        out[0] = laplace_score_mu(laplace_sign(r), b);
        out[1] = laplace_hess_mu_mu();
    } else {
        double ar = std::fabs(r);
        out[0] = laplace_score_sigma(ar, b);
        out[1] = laplace_hess_sigma_sigma(ar, b);
    }
}

inline void laplace_info_dinfo(int k, double y, const double* th,
                               double* out) {
    double b = th[1];
    if (k == 0) {
        out[0] = laplace_expected_mu_mu(b);
        out[1] = laplace_dexpected_mu_mu_mu();
    } else {
        out[0] = laplace_expected_sigma_sigma(b);
        out[1] = laplace_dexpected_sigma_sigma_sigma(b);
    }
}

// The mixed components, written as laplace.cpp writes them: the Hessian's
// -sgn(r)/b^2 and the derivative of E_mm = -1/b^2 in b, 2/b^3.
inline double laplace_hess_mu_sigma(double sg, double b) {
    return -sg * (1.0 / b) * (1.0 / b);
}

inline double laplace_dexpected_mu_mu_sigma(double b) {
    return 2.0 / (b * b * b);
}

// The routers of the joint registry, one per order, in (mu, sigma); the
// layout is the one pt_cauchy.h describes. The Laplace kernels take the
// sign of the residual and its modulus.
inline void laplace_joint_grad(double y, const double* th, double* out) {
    double r = y - th[0], b = th[1];
    out[0] = laplace_score_mu(laplace_sign(r), b);
    out[1] = laplace_score_sigma(std::fabs(r), b);
}

inline void laplace_joint_hess(double y, const double* th, double* out) {
    double r = y - th[0], b = th[1];
    out[0] = laplace_hess_mu_mu();
    out[3] = laplace_hess_sigma_sigma(std::fabs(r), b);
    out[1] = out[2] = laplace_hess_mu_sigma(laplace_sign(r), b);
}

inline void laplace_joint_info(double y, const double* th, double* out) {
    double b = th[1];
    out[0] = laplace_expected_mu_mu(b);
    out[3] = laplace_expected_sigma_sigma(b);
    out[1] = out[2] = 0.0;
}

// E_ms = 0 for every sigma, so both of its derivatives are zero
inline void laplace_joint_dinfo(double y, const double* th, double* out) {
    double b = th[1];
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[4] = laplace_dexpected_mu_mu_sigma(b);
    out[7] = laplace_dexpected_sigma_sigma_sigma(b);
}

// The parts route (pt_parts.h). The buffer, with r = y - mu, the sign of r
// and its modulus as the joint routers form them, and 1/b, which the orders
// three and four read:
//
//   [0] level
//   [1] y                           NaN in the expected buffer
//   [2] mu
//   [3] b
//   [4] 1/b                         always
//   [5] sgn(r), 0 at r = 0          level 1
//   [6] |r|                         level 1
inline void laplace_parts(double y, const double* th, int level, double* b) {
    const double m = th[0], bb = th[1];
    const double r = y - m;
    b[0] = level; b[1] = y; b[2] = m; b[3] = bb; b[4] = 1.0 / bb;
    b[5] = level >= 1 ? laplace_sign(r) : R_NaN;
    b[6] = level >= 1 ? std::fabs(r) : R_NaN;
}

inline void laplace_parts_expected(const double* th, int level, double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = th[1];
    b[4] = 1.0 / th[1];
    b[5] = R_NaN; b[6] = R_NaN;
}

inline void laplace_parts_grad(const double* b, double* out) {
    out[0] = laplace_score_mu(b[5], b[3]);
    out[1] = laplace_score_sigma(b[6], b[3]);
}

inline void laplace_parts_hess(const double* b, double* out) {
    out[0] = laplace_hess_mu_mu();
    out[3] = laplace_hess_sigma_sigma(b[6], b[3]);
    out[1] = out[2] = laplace_hess_mu_sigma(b[5], b[3]);
}

// from the expected parts
inline void laplace_parts_info(const double* b, double* out) {
    out[0] = laplace_expected_mu_mu(b[3]);
    out[3] = laplace_expected_sigma_sigma(b[3]);
    out[1] = out[2] = 0.0;
}

// third order: laplace_deriv3_cpp (laplace_hd.cpp), unique components in the
// sorted order mu mu mu, mu mu sigma, mu sigma sigma, sigma sigma sigma
inline void laplace_parts_d3(const double* b, double* out) {
    const double s = b[5], ib = b[4];
    const double ib3 = ib * ib * ib;
    const double u = b[6] * ib;
    double uu[4];
    uu[0] = 0.0;
    uu[1] = 0.0;
    uu[2] = 2.0 * s * ib3;
    uu[3] = (6.0 * u - 2.0) * ib3;
    sym_fill(2, 3, uu, out);
}

// fourth order: laplace_deriv4_cpp (laplace_hd.cpp), unique components in
// the sorted order mu mu mu mu, mu mu mu sigma, mu mu sigma sigma, mu sigma
// sigma sigma, sigma sigma sigma sigma
inline void laplace_parts_d4(const double* b, double* out) {
    const double s = b[5], ib = b[4];
    const double ib2 = ib * ib, ib4 = ib2 * ib2;
    const double u = b[6] * ib;
    double uu[5];
    uu[0] = 0.0;
    uu[1] = 0.0;
    uu[2] = 0.0;
    uu[3] = -6.0 * s * ib4;
    uu[4] = (6.0 - 24.0 * u) * ib4;
    sym_fill(2, 4, uu, out);
}

} // namespace d7

#endif
