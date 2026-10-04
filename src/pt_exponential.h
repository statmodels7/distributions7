#ifndef D7_PT_EXPONENTIAL_H
#define D7_PT_EXPONENTIAL_H

// Exponential in its mean, l = -log mu - y/mu: one function per component
// and order for the quantities the scalar registry reads, called by the
// vector kernels (exponential.cpp, dexpected_kernels.cpp) and by the
// registry (d7_ccallable.cpp).

namespace d7 {

inline double exponential_score_mu(double y, double m) {
    return (y - m) / (m * m);
}

inline double exponential_hess_mu_mu(double y, double m) {
    double m2 = m * m;
    return 1.0 / m2 - 2.0 * y / (m2 * m);
}

inline double exponential_expected_mu_mu(double m) {
    return -1.0 / (m * m);
}

inline double exponential_dexpected_mu_mu_mu(double m) {
    double u = 1.0 / m, u3 = u * u * u;
    return 2.0 * u3;
}

inline void exponential_score_curv(int k, double y, const double* th,
                                   double* out) {
    out[0] = exponential_score_mu(y, th[0]);
    out[1] = exponential_hess_mu_mu(y, th[0]);
}

inline void exponential_info_dinfo(int k, double y, const double* th,
                                   double* out) {
    out[0] = exponential_expected_mu_mu(th[0]);
    out[1] = exponential_dexpected_mu_mu_mu(th[0]);
}

} // namespace d7

#endif
