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

} // namespace d7

#endif
