#ifndef D7_PT_GAUSSIAN1_H
#define D7_PT_GAUSSIAN1_H

#include <Rcpp.h>

// Gaussian in the mean and the standard deviation: one function per
// component and order for the quantities the registry in d7_ccallable.cpp
// reads (the score, the Hessian, the expected information and its first
// derivative in every parameter). The vector kernels in gaussian.cpp call
// the components of the scalar registry (the score, the diagonal of the
// Hessian and of the expected information, and the latter's derivative in
// its own parameter), so the two routes run the same arithmetic there. The
// off-diagonal components that only the joint registry reads are written
// from the kernels' expressions, which the kernels keep inline beside the
// intermediates they share, as student_t1.cpp does for the reason
// recorded in pt_student_t1.h; test-joint-registry.R holds the two to each
// other. Every
// component is written in z = (y - mu)/sigma and inv = 1/sigma, never in a
// positive power of the scale; gaussian.cpp records why.

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

inline double gaussian1_hess_mu_sigma(double z, double inv) {
    double inv2 = inv * inv;
    return -2.0 * z * inv2;
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

inline double gaussian1_dexpected_mu_mu_sigma(double inv) {
    double inv3 = inv * inv * inv;
    return 2.0 * inv3;
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

// The routers of the joint registry, one per order: every component at one
// observation, on the parameter scale. The Hessian and the expected
// information fill a full P x P matrix and the derivative of the expected
// information a P x P x P array, column-major, out[a + P b + P^2 c] being
// the derivative of E[l_ab] in parameter c.
inline void gaussian1_joint_grad(double y, const double* th, double* out) {
    double inv = 1.0 / th[1];
    double z = (y - th[0]) / th[1];
    out[0] = gaussian1_score_mu(z, inv);
    out[1] = gaussian1_score_sigma(z, inv);
}

inline void gaussian1_joint_hess(double y, const double* th, double* out) {
    double inv = 1.0 / th[1];
    double z = (y - th[0]) / th[1];
    out[0] = gaussian1_hess_mu_mu(inv);
    out[3] = gaussian1_hess_sigma_sigma(z, inv);
    out[1] = out[2] = gaussian1_hess_mu_sigma(z, inv);
}

inline void gaussian1_joint_info(double y, const double* th, double* out) {
    double inv = 1.0 / th[1];
    out[0] = gaussian1_expected_mu_mu(inv);
    out[3] = gaussian1_expected_sigma_sigma(inv);
    out[1] = out[2] = 0.0;
}

inline void gaussian1_joint_dinfo(double y, const double* th, double* out) {
    double inv = 1.0 / th[1];
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[0] = gaussian1_dexpected_mu_mu_mu();
    out[4] = gaussian1_dexpected_mu_mu_sigma(inv);
    out[7] = gaussian1_dexpected_sigma_sigma_sigma(inv);
}

} // namespace d7

#endif
