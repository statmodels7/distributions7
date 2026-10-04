#ifndef D7_PT_BERNOULLI_H
#define D7_PT_BERNOULLI_H

// Bernoulli in its mean: one function per component and order for the
// quantities the scalar registry reads, called by the vector kernels
// (bernoulli.cpp, dexpected_kernels.cpp) and by the registry
// (d7_ccallable.cpp).

namespace d7 {

inline double bernoulli_score_mu(double y, double m) {
    return (y - m) / (m * (1.0 - m));
}

inline double bernoulli_hess_mu_mu(double y, double m) {
    return -(y / (m * m)) - ((1.0 - y) / ((1.0 - m) * (1.0 - m)));
}

inline double bernoulli_expected_mu_mu(double m) {
    return -1.0 / (m * (1.0 - m));
}

// E = -1/q with q = m(1 - m)
inline double bernoulli_dexpected_mu_mu_mu(double m) {
    double q = m * (1.0 - m), s = 1.0 - 2.0 * m, q2 = q * q;
    return s / q2;
}

inline void bernoulli_score_curv(int k, double y, const double* th,
                                 double* out) {
    out[0] = bernoulli_score_mu(y, th[0]);
    out[1] = bernoulli_hess_mu_mu(y, th[0]);
}

inline void bernoulli_info_dinfo(int k, double y, const double* th,
                                 double* out) {
    out[0] = bernoulli_expected_mu_mu(th[0]);
    out[1] = bernoulli_dexpected_mu_mu_mu(th[0]);
}

} // namespace d7

#endif
