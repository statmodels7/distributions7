#ifndef D7_PT_CAUCHY_H
#define D7_PT_CAUCHY_H

// Cauchy in location and scale: one function per component and order for
// the quantities the scalar registry reads, called by the vector kernels
// (cauchy.cpp, dexpected_kernels.cpp) and by the registry
// (d7_ccallable.cpp). res = y - mu.

namespace d7 {

inline double cauchy_score_mu(double res, double s) {
    double res2 = res * res;
    double s2 = s * s;
    double den = s2 + res2;
    return (2.0 * res) / den;
}

inline double cauchy_score_sigma(double res, double s) {
    double res2 = res * res;
    double s2 = s * s;
    double den = s2 + res2;
    return (res2 - s2) / (s * den);
}

inline double cauchy_hess_mu_mu(double res, double s) {
    double res2 = res * res;
    double s2 = s * s;
    double den = s2 + res2;
    double den2 = den * den;
    return (2.0 * res2 - 2.0 * s2) / den2;
}

inline double cauchy_hess_sigma_sigma(double res, double s) {
    double res2 = res * res;
    double res4 = res2 * res2;
    double s2 = s * s;
    double s4 = s2 * s2;
    double den = s2 + res2;
    double den2 = den * den;
    return (s4 - 4.0 * s2 * res2 - res4) / (s2 * den2);
}

inline double cauchy_expected_mu_mu(double s) {
    double s2 = s * s;
    return -0.5 / s2;
}

inline double cauchy_expected_sigma_sigma(double s) {
    double s2 = s * s;
    return -0.5 / s2;
}

// E_mm and E_ss are k/sigma^2 with k = -1/2; nothing moves with mu
inline double cauchy_dexpected_mu_mu_mu() {
    return 0.0;
}

inline double cauchy_dexpected_sigma_sigma_sigma(double s) {
    const double k = -0.5;
    double u = 1.0 / s, u3 = u * u * u;
    return -2.0 * k * u3;
}

inline void cauchy_score_curv(int k, double y, const double* th,
                              double* out) {
    double res = y - th[0], s = th[1];
    if (k == 0) {
        out[0] = cauchy_score_mu(res, s);
        out[1] = cauchy_hess_mu_mu(res, s);
    } else {
        out[0] = cauchy_score_sigma(res, s);
        out[1] = cauchy_hess_sigma_sigma(res, s);
    }
}

inline void cauchy_info_dinfo(int k, double y, const double* th,
                              double* out) {
    double s = th[1];
    if (k == 0) {
        out[0] = cauchy_expected_mu_mu(s);
        out[1] = cauchy_dexpected_mu_mu_mu();
    } else {
        out[0] = cauchy_expected_sigma_sigma(s);
        out[1] = cauchy_dexpected_sigma_sigma_sigma(s);
    }
}

} // namespace d7

#endif
