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

} // namespace d7

#endif
