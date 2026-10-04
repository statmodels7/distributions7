#ifndef D7_PT_GAUSSIAN1_H
#define D7_PT_GAUSSIAN1_H

#include <Rcpp.h>

// Gaussian in the mean and the standard deviation: one function per
// component and order for the quantities the scalar registry reads (the
// score, the diagonal of the Hessian, the diagonal of the expected
// information and its derivative in its own parameter). The vector kernels
// in gaussian.cpp and the registry in d7_ccallable.cpp both call these, so
// the two routes run the same arithmetic. Every component is written in
// z = (y - mu)/sigma and inv = 1/sigma, never in a positive power of the
// scale; gaussian.cpp records why.

namespace d7 {

inline double gaussian1_score_mu(double z, double inv) {
    return z * inv;
}

inline double gaussian1_score_sigma(double z, double inv) {
    return (z * z - 1.0) * inv;
}

inline double gaussian1_hess_mu_mu(double inv) {
    double inv2 = inv * inv;
    return -inv2;
}

inline double gaussian1_hess_sigma_sigma(double z, double inv) {
    double inv2 = inv * inv;
    return (1.0 - 3.0 * z * z) * inv2;
}

inline double gaussian1_expected_mu_mu(double inv) {
    double inv2 = inv * inv;
    return -inv2;
}

inline double gaussian1_expected_sigma_sigma(double inv) {
    double inv2 = inv * inv;
    return -2.0 * inv2;
}

// E_mm = -1/sigma^2 does not depend on mu
inline double gaussian1_dexpected_mu_mu_mu() {
    return 0.0;
}

inline double gaussian1_dexpected_sigma_sigma_sigma(double inv) {
    double inv3 = inv * inv * inv;
    return 4.0 * inv3;
}

// The routers of the scalar registry: k is the 0-based parameter index,
// th the parameter vector at one observation, out[0] and out[1] the two
// quantities on the parameter scale.
inline void gaussian1_score_curv(int k, double y, const double* th,
                                 double* out) {
    double inv = 1.0 / th[1];
    double z = (y - th[0]) / th[1];
    if (k == 0) {
        out[0] = gaussian1_score_mu(z, inv);
        out[1] = gaussian1_hess_mu_mu(inv);
    } else {
        out[0] = gaussian1_score_sigma(z, inv);
        out[1] = gaussian1_hess_sigma_sigma(z, inv);
    }
}

inline void gaussian1_info_dinfo(int k, double y, const double* th,
                                 double* out) {
    double inv = 1.0 / th[1];
    if (k == 0) {
        out[0] = gaussian1_expected_mu_mu(inv);
        out[1] = gaussian1_dexpected_mu_mu_mu();
    } else {
        out[0] = gaussian1_expected_sigma_sigma(inv);
        out[1] = gaussian1_dexpected_sigma_sigma_sigma(inv);
    }
}

} // namespace d7

#endif
