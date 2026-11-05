#ifndef D7_PT_LOGISTIC_H
#define D7_PT_LOGISTIC_H

#include <Rcpp.h>
#include <cmath>
#include "pt_parts.h"

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

// The parts route (pt_parts.h). Three costly quantities are shared, each
// computed once per observation: tanh(res/(2 sigma)) for the orders one and
// two, sech^2 of the same argument for the second order (an exponential), and
// the sigmoid pair of z = res/sigma for the orders three and four, which
// logistic_hd.cpp computes on its own stable branch. Each order keeps its own
// route, as the joint routers and the vector kernels do.
//   [0] level                       always
//   [1] y                           observed parts; NaN in the expected parts
//   [2] mu                          observed parts; NaN in the expected parts
//   [3] sigma                       always
//   [4] res = y - mu                observed parts; NaN in the expected parts
//   [5] t = tanh(res/(2 sigma))     level >= 1; NaN in the expected parts
//   [6] sech^2(res/(2 sigma))       level >= 2; NaN in the expected parts
//   [7] sigmoid(z), z = res/sigma   level >= 3; NaN in the expected parts
//   [8] 1 - sigmoid(z)              level >= 3; NaN in the expected parts

// the sigmoid pair of logistic_hd.cpp (sigmoid_pair), copied under this prefix
inline void logistic_parts_sigmoid_pair(double z, double& t, double& u) {
    if (z >= 0.0) {
        double e = std::exp(-z);
        t = 1.0 / (1.0 + e);
        u = e / (1.0 + e);
    } else {
        double e = std::exp(z);
        t = e / (1.0 + e);
        u = 1.0 / (1.0 + e);
    }
}

// logistic_hess_mu_mu() with sech^2 given
inline double logistic_parts_hess_mu_mu(double sech2, double s) {
    double s2 = s * s;
    return -sech2 / (2.0 * s2);
}

// logistic_hess_sigma_sigma() with sech^2 given
inline double logistic_parts_hess_sigma_sigma(double res, double t,
                                              double sech2_z, double s) {
    double s2 = s * s;
    double z_half = 0.5 * res / s;
    return (1.0 - 4.0 * z_half * t - 2.0 * z_half * z_half * sech2_z) / s2;
}

inline void logistic_parts(double y, const double* th, int level, double* b) {
    const double m = th[0], s = th[1];
    const double res = y - m;
    b[0] = level; b[1] = y; b[2] = m; b[3] = s; b[4] = res;
    if (level >= 1) b[5] = std::tanh(0.5 * res / s); else b[5] = R_NaN;
    if (level >= 2) b[6] = logistic_sech2_half(res, s); else b[6] = R_NaN;
    if (level >= 3) {
        double t, u;
        logistic_parts_sigmoid_pair(res / s, t, u);
        b[7] = t; b[8] = u;
    } else {
        b[7] = R_NaN; b[8] = R_NaN;
    }
}

inline void logistic_parts_expected(const double* th, int level, double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = R_NaN; b[3] = th[1];
    for (int j = 4; j <= 8; ++j) b[j] = R_NaN;
}

// from the observed parts (the joint grad)
inline void logistic_parts_grad(const double* b, double* out) {
    const double res = b[4], t = b[5], s = b[3];
    out[0] = logistic_score_mu(t, s);
    out[1] = logistic_score_sigma(res, t, s);
}

inline void logistic_parts_hess(const double* b, double* out) {
    const double res = b[4], t = b[5], sech2 = b[6], s = b[3];
    out[0] = logistic_parts_hess_mu_mu(sech2, s);
    out[3] = logistic_parts_hess_sigma_sigma(res, t, sech2, s);
    out[1] = out[2] = logistic_hess_mu_sigma(res, t, s);
}

// from the expected parts
inline void logistic_parts_info(const double* b, double* out) {
    out[0] = logistic_expected_mu_mu(b[3]);
    out[3] = logistic_expected_sigma_sigma(b[3]);
    out[1] = out[2] = 0.0;
}

// third order: logistic_deriv3_cpp, unique components in the sorted order
// mu mu mu, mu mu sigma, mu sigma sigma, sigma sigma sigma
inline void logistic_parts_d3(const double* b, double* out) {
    const double s = b[3];
    const double z = b[4] / s, z2 = z * z, z3 = z2 * z;
    const double t = b[7], u = b[8];
    const double tu = t * u;
    const double g1 = u - t;
    const double g2 = -2.0 * tu;
    const double g3 = g2 * g1;
    const double s3 = s * s * s;
    double w[4];
    w[0] = -g3 / s3;
    w[1] = -(2.0 * g2 + z * g3) / s3;
    w[2] = -(2.0 * g1 + 4.0 * z * g2 + z2 * g3) / s3;
    w[3] = -(2.0 + 6.0 * z * g1 + 6.0 * z2 * g2 + z3 * g3) / s3;
    sym_fill(2, 3, w, out);
}

// fourth order: logistic_deriv4_cpp, unique components in the sorted order
// mu mu mu mu, mu mu mu sigma, mu mu sigma sigma, mu sigma sigma sigma,
// sigma sigma sigma sigma
inline void logistic_parts_d4(const double* b, double* out) {
    const double s = b[3];
    const double s2 = s * s, s4 = s2 * s2;
    const double z = b[4] / s, z2 = z * z, z3 = z2 * z, z4 = z2 * z2;
    const double t = b[7], u = b[8];
    const double tu = t * u;
    const double g1 = u - t;
    const double g2 = -2.0 * tu;
    const double g3 = g2 * g1;
    const double g4 = g2 * (1.0 - 6.0 * tu);
    double w[5];
    w[0] = g4 / s4;
    w[1] = (3.0 * g3 + z * g4) / s4;
    w[2] = (6.0 * g2 + 6.0 * z * g3 + z2 * g4) / s4;
    w[3] = (6.0 * g1 + 18.0 * z * g2 + 9.0 * z2 * g3 + z3 * g4) / s4;
    w[4] = (6.0 + 24.0 * z * g1 + 36.0 * z2 * g2 + 12.0 * z3 * g3 + z4 * g4) / s4;
    sym_fill(2, 4, w, out);
}

} // namespace d7

#endif
