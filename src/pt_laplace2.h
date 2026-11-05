#ifndef D7_PT_LAPLACE2_H
#define D7_PT_LAPLACE2_H

#include <Rcpp.h>
#include <cmath>
#include "pt_laplace.h"
#include "pt_parts.h"

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

// The parts route (pt_parts.h). The buffer, with r = y - mu, the sign of r
// and its modulus as the joint routers form them:
//
//   [0] level
//   [1] y                           NaN in the expected buffer
//   [2] mu
//   [3] lambda
//   [4] sgn(r), 0 at r = 0          level 1
//   [5] |r|                         level 1
inline void laplace2_parts(double y, const double* th, int level, double* b) {
    const double m = th[0], lam = th[1];
    const double r = y - m;
    b[0] = level; b[1] = y; b[2] = m; b[3] = lam;
    b[4] = level >= 1 ? laplace_sign(r) : R_NaN;
    b[5] = level >= 1 ? std::fabs(r) : R_NaN;
}

inline void laplace2_parts_expected(const double* th, int level, double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = th[1];
    b[4] = R_NaN; b[5] = R_NaN;
}

inline void laplace2_parts_grad(const double* b, double* out) {
    out[0] = laplace2_score_mu(b[4], b[3]);
    out[1] = laplace2_score_lambda(b[5], b[3]);
}

inline void laplace2_parts_hess(const double* b, double* out) {
    out[0] = laplace2_hess_mu_mu();
    out[3] = laplace2_hess_lambda_lambda(b[3]);
    out[1] = out[2] = laplace2_hess_mu_lambda(b[4]);
}

// from the expected parts
inline void laplace2_parts_info(const double* b, double* out) {
    out[0] = laplace2_expected_mu_mu(b[3]);
    out[3] = laplace2_expected_lambda_lambda(b[3]);
    out[1] = out[2] = 0.0;
}

// third order: no compiled kernel, the R closed form of
// distrib_deriv3.Laplace2Distrib written in C++. Every component with mu is
// zero, and lambda^3 carries 2/lambda^3.
inline void laplace2_parts_d3(const double* b, double* out) {
    const double lam = b[3];
    const double l2 = lam * lam;
    double u[4];
    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[3] = 2.0 / (l2 * lam);
    sym_fill(2, 3, u, out);
}

// fourth order: no compiled kernel, the R closed form of
// distrib_deriv4.Laplace2Distrib, -6/lambda^4 in the pure rate component
inline void laplace2_parts_d4(const double* b, double* out) {
    const double lam = b[3];
    const double l2 = lam * lam;
    double u[5];
    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[3] = 0.0;
    u[4] = -6.0 / (l2 * l2);
    sym_fill(2, 4, u, out);
}

} // namespace d7

#endif
