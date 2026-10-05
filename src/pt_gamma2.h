#ifndef D7_PT_GAMMA2_H
#define D7_PT_GAMMA2_H

#include <Rcpp.h>
#include <cmath>
#include "psi_diff.h"

// Gamma in the mean and the variance, with shape a = mu^2/v and rate
// lambda = mu/v: one function per component and order for the quantities
// the scalar registry reads, called by the vector kernels (gamma.cpp,
// dexpected_kernels.cpp) and by the registry (d7_ccallable.cpp). The caller
// passes dg = psi(a), tg = psi'(a), log lambda and log y; the derivatives
// of the expected information read r1 = 1/a - psi'(a) and
// r2 = -1/a^2 - psi''(a) (psi_diff.h), through g = a psi'(a) - 1 = -a r1.

namespace d7 {

inline double gamma2_score_mu(double y, double m, double s2, double dg,
                              double log_lambda, double log_y) {
    return (-2.0 * m * dg + 2.0 * m * log_lambda + m + 2.0 * m * log_y - y) / s2;
}

inline double gamma2_score_sigma2(double y, double m, double s2, double dg,
                                  double log_lambda, double log_y) {
    return -(m * (-m * dg + m + m * (log_lambda + log_y) - y)) / (s2 * s2);
}

inline double gamma2_hess_mu_mu(double m, double s2, double dg, double tg,
                                double log_lambda, double log_y) {
    return (-(4.0 * m * m * tg) / s2 - 2.0 * dg + 2.0 * log_lambda +
            2.0 * log_y + 3.0) / s2;
}

inline double gamma2_hess_sigma2_sigma2(double y, double m, double s2,
                                        double dg, double tg,
                                        double log_lambda, double log_y) {
    return -(m * (2.0 * m * s2 * dg + m * m * m * tg +
                  s2 * (-2.0 * m * log_lambda - 3.0 * m - 2.0 * m * log_y +
                        2.0 * y))) / (s2 * s2 * s2 * s2);
}

inline double gamma2_expected_mu_mu(double m, double s2, double tg) {
    return (3.0 * s2 - 4.0 * m * m * tg) / (s2 * s2);
}

inline double gamma2_expected_sigma2_sigma2(double m, double s2, double tg) {
    return -(m * m * (m * m * tg - s2)) / (s2 * s2 * s2 * s2);
}

inline double gamma2_dexpected_mu_mu_mu(double m, double v, double r1,
                                        double r2) {
    double a = m * m / v;
    double g1 = -r1 - a * r2;
    double iv = 1.0 / v, iv2 = iv * iv;
    return -8.0 * m * g1 * iv2;
}

inline double gamma2_dexpected_sigma2_sigma2_sigma2(double m, double v,
                                                    double r1, double r2) {
    double a = m * m / v;
    double g = -a * r1, g1 = -r1 - a * r2;
    double iv = 1.0 / v, iv2 = iv * iv, iv3 = iv2 * iv;
    return a * (3.0 * g + a * g1) * iv3;
}

inline void gamma2_score_curv(int k, double y, const double* th,
                              double* out) {
    double m = th[0], s2 = th[1];
    double alpha = m * m / s2, lambda = m / s2;
    double dg = R::digamma(alpha), tg = R::trigamma(alpha);
    double log_lambda = std::log(lambda), log_y = std::log(y);
    if (k == 0) {
        out[0] = gamma2_score_mu(y, m, s2, dg, log_lambda, log_y);
        out[1] = gamma2_hess_mu_mu(m, s2, dg, tg, log_lambda, log_y);
    } else {
        out[0] = gamma2_score_sigma2(y, m, s2, dg, log_lambda, log_y);
        out[1] = gamma2_hess_sigma2_sigma2(y, m, s2, dg, tg, log_lambda,
                                           log_y);
    }
}

inline void gamma2_info_dinfo(int k, double y, const double* th,
                              double* out) {
    double m = th[0], s2 = th[1];
    double a = m * m / s2;
    double tg = R::trigamma(a);
    double r1 = d7::psi1_rest(a), r2 = d7::psi2_rest(a);
    if (k == 0) {
        out[0] = gamma2_expected_mu_mu(m, s2, tg);
        out[1] = gamma2_dexpected_mu_mu_mu(m, s2, r1, r2);
    } else {
        out[0] = gamma2_expected_sigma2_sigma2(m, s2, tg);
        out[1] = gamma2_dexpected_sigma2_sigma2_sigma2(m, s2, r1, r2);
    }
}

} // namespace d7

#endif
