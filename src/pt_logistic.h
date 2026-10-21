#ifndef D7_PT_LOGISTIC_H
#define D7_PT_LOGISTIC_H

#include <Rcpp.h>
#include <cmath>

// Logistic in location and scale: one function per component and order for
// the quantities the scalar registry reads, called by the vector kernels
// (logistic.cpp, dexpected_kernels.cpp) and by the registry
// (d7_ccallable.cpp). With res = y - mu, the caller passes
// t = tanh(res/(2 sigma)).

namespace d7 {

inline double logistic_score_mu(double t, double s) {
    return t / s;
}

inline double logistic_score_sigma(double res, double t, double s) {
    return (res * t - s) / (s * s);
}

// sech^2(z/2) = 4e/(1 + e)^2 with e = exp(-|z|): 1 - tanh^2(z/2) cancels
// in the tails, where the value is about 4 exp(-|z|)
inline double logistic_sech2_half(double res, double s) {
    double e = std::exp(-std::fabs(res / s));
    return 4.0 * e / ((1.0 + e) * (1.0 + e));
}

inline double logistic_hess_mu_mu(double res, double s) {
    double s2 = s * s;
    return -logistic_sech2_half(res, s) / (2.0 * s2);
}

inline double logistic_hess_sigma_sigma(double res, double t, double s) {
    double s2 = s * s;
    double z_half = 0.5 * res / s;
    double sech2_z = logistic_sech2_half(res, s);
    return (1.0 - 4.0 * z_half * t - 2.0 * z_half * z_half * sech2_z) / s2;
}

inline double logistic_expected_mu_mu(double s) {
    double s2 = s * s;
    return -1.0 / (3.0 * s2);
}

inline double logistic_expected_sigma_sigma(double s) {
    double s2 = s * s;
    return -(3.0 + M_PI * M_PI) / (9.0 * s2);
}

// E_mm and E_ss are k/sigma^2; nothing moves with mu
inline double logistic_dexpected_mu_mu_mu() {
    return 0.0;
}

inline double logistic_dexpected_sigma_sigma_sigma(double s) {
    const double k = -(3.0 + M_PI * M_PI) / 9.0;
    double u = 1.0 / s, u3 = u * u * u;
    return -2.0 * k * u3;
}

inline void logistic_score_curv(int k, double y, const double* th,
                                double* out) {
    double res = y - th[0], s = th[1];
    double t = std::tanh(0.5 * res / s);
    if (k == 0) {
        out[0] = logistic_score_mu(t, s);
        out[1] = logistic_hess_mu_mu(res, s);
    } else {
        out[0] = logistic_score_sigma(res, t, s);
        out[1] = logistic_hess_sigma_sigma(res, t, s);
    }
}

inline void logistic_info_dinfo(int k, double y, const double* th,
                                double* out) {
    double s = th[1];
    if (k == 0) {
        out[0] = logistic_expected_mu_mu(s);
        out[1] = logistic_dexpected_mu_mu_mu();
    } else {
        out[0] = logistic_expected_sigma_sigma(s);
        out[1] = logistic_dexpected_sigma_sigma_sigma(s);
    }
}

// The mixed Hessian component, written as logistic.cpp writes it: sech^2 is
// formed as 1 - t^2 from the tanh the caller passes, not through
// logistic_sech2_half.
inline double logistic_hess_mu_sigma(double res, double t, double s) {
    double s2 = s * s;
    double z_half = 0.5 * res / s;
    double sech2_z = 1.0 - t * t;
    return -(t + z_half * sech2_z) / s2;
}

// The derivative of E_mm = -1/(3 sigma^2) in sigma, -2 k u^3 with k = -1/3
inline double logistic_dexpected_mu_mu_sigma(double s) {
    const double k = -1.0 / 3.0;
    double u = 1.0 / s, u3 = u * u * u;
    return -2.0 * k * u3;
}

// The routers of the joint registry, one per order, in (mu, sigma); the
// layout is the one pt_cauchy.h describes.
inline void logistic_joint_grad(double y, const double* th, double* out) {
    double res = y - th[0], s = th[1];
    double t = std::tanh(0.5 * res / s);
    out[0] = logistic_score_mu(t, s);
    out[1] = logistic_score_sigma(res, t, s);
}

inline void logistic_joint_hess(double y, const double* th, double* out) {
    double res = y - th[0], s = th[1];
    double t = std::tanh(0.5 * res / s);
    out[0] = logistic_hess_mu_mu(res, s);
    out[3] = logistic_hess_sigma_sigma(res, t, s);
    out[1] = out[2] = logistic_hess_mu_sigma(res, t, s);
}

inline void logistic_joint_info(double y, const double* th, double* out) {
    double s = th[1];
    out[0] = logistic_expected_mu_mu(s);
    out[3] = logistic_expected_sigma_sigma(s);
    out[1] = out[2] = 0.0;
}

// E_ms = 0 for every sigma, so both of its derivatives are zero
inline void logistic_joint_dinfo(double y, const double* th, double* out) {
    double s = th[1];
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[4] = logistic_dexpected_mu_mu_sigma(s);
    out[7] = logistic_dexpected_sigma_sigma_sigma(s);
}

} // namespace d7

#endif
