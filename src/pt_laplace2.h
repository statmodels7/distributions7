#ifndef D7_PT_LAPLACE2_H
#define D7_PT_LAPLACE2_H

#include <cmath>
#include "pt_laplace.h"

// The Laplace in location mu and rate lambda, l = log(lambda/2) -
// lambda |y - mu|: one function per component and order for the quantities
// the scalar registry reads, called by the vector kernels (laplace.cpp) and
// by the registry (d7_ccallable.cpp). With r = y - mu the caller passes
// sgn(r) (laplace_sign) and |r|.

namespace d7 {

inline double laplace2_score_mu(double sg, double lam) {
    return lam * sg;
}

inline double laplace2_score_lambda(double ar, double lam) {
    return 1.0 / lam - ar;
}

// zero away from the kink
inline double laplace2_hess_mu_mu() {
    return 0.0;
}

inline double laplace2_hess_lambda_lambda(double lam) {
    return -1.0 / (lam * lam);
}

inline double laplace2_expected_mu_mu(double lam) {
    return -lam * lam;
}

inline double laplace2_expected_lambda_lambda(double lam) {
    return -1.0 / (lam * lam);
}

inline double laplace2_dexpected_mu_mu_mu() {
    return 0.0;
}

inline double laplace2_dexpected_lambda_lambda_lambda(double lam) {
    return 2.0 / (lam * lam * lam);
}

inline void laplace2_score_curv(int k, double y, const double* th,
                                double* out) {
    double r = y - th[0], lam = th[1];
    if (k == 0) {
        out[0] = laplace2_score_mu(laplace_sign(r), lam);
        out[1] = laplace2_hess_mu_mu();
    } else {
        out[0] = laplace2_score_lambda(std::fabs(r), lam);
        out[1] = laplace2_hess_lambda_lambda(lam);
    }
}

inline void laplace2_info_dinfo(int k, double y, const double* th,
                                double* out) {
    double lam = th[1];
    if (k == 0) {
        out[0] = laplace2_expected_mu_mu(lam);
        out[1] = laplace2_dexpected_mu_mu_mu();
    } else {
        out[0] = laplace2_expected_lambda_lambda(lam);
        out[1] = laplace2_dexpected_lambda_lambda_lambda(lam);
    }
}

// The mixed components, written as laplace.cpp writes them: the Hessian's
// mu-lambda entry is the sign of the residual itself, and the derivative of
// E_mm = -lambda^2 in lambda is -2 lambda.
inline double laplace2_hess_mu_lambda(double sg) {
    return sg;
}

inline double laplace2_dexpected_mu_mu_lambda(double lam) {
    return -2.0 * lam;
}

// The routers of the joint registry, one per order, in (mu, lambda); the
// layout is the one pt_cauchy.h describes.
inline void laplace2_joint_grad(double y, const double* th, double* out) {
    double r = y - th[0], lam = th[1];
    out[0] = laplace2_score_mu(laplace_sign(r), lam);
    out[1] = laplace2_score_lambda(std::fabs(r), lam);
}

inline void laplace2_joint_hess(double y, const double* th, double* out) {
    double r = y - th[0], lam = th[1];
    out[0] = laplace2_hess_mu_mu();
    out[3] = laplace2_hess_lambda_lambda(lam);
    out[1] = out[2] = laplace2_hess_mu_lambda(laplace_sign(r));
}

inline void laplace2_joint_info(double y, const double* th, double* out) {
    double lam = th[1];
    out[0] = laplace2_expected_mu_mu(lam);
    out[3] = laplace2_expected_lambda_lambda(lam);
    out[1] = out[2] = 0.0;
}

// E_ml = 0 for every rate, so both of its derivatives are zero
inline void laplace2_joint_dinfo(double y, const double* th, double* out) {
    double lam = th[1];
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[4] = laplace2_dexpected_mu_mu_lambda(lam);
    out[7] = laplace2_dexpected_lambda_lambda_lambda(lam);
}

} // namespace d7

#endif
