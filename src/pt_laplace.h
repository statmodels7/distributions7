#ifndef D7_PT_LAPLACE_H
#define D7_PT_LAPLACE_H

#include <cmath>

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

} // namespace d7

#endif
