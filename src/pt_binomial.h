#ifndef D7_PT_BINOMIAL_H
#define D7_PT_BINOMIAL_H

// Binomial in its success probability, with the size N fixed: one function
// per component and order for the quantities the scalar registry reads,
// called by the vector kernels (binomial.cpp, dexpected_kernels.cpp).

namespace d7 {

inline double binomial_score_mu(double y, double m, double sz) {
    return (y - sz * m) / (m * (1.0 - m));
}

inline double binomial_hess_mu_mu(double y, double m, double sz) {
    double m2 = m * m;
    double one_m = 1.0 - m;
    double one_m2 = one_m * one_m;
    return -(y / m2) - ((sz - y) / one_m2);
}

inline double binomial_expected_mu_mu(double m, double sz) {
    return -sz / (m * (1.0 - m));
}

// the Bernoulli's derivative times the size
inline double binomial_dexpected_mu_mu_mu(double m, double sz) {
    double q = m * (1.0 - m), s = 1.0 - 2.0 * m, q2 = q * q;
    return sz * (s / q2);
}

// th = (mu, size): the size is a constant of the family and follows the
// parameters
inline void binomial_score_curv(int k, double y, const double* th,
                                double* out) {
    out[0] = binomial_score_mu(y, th[0], th[1]);
    out[1] = binomial_hess_mu_mu(y, th[0], th[1]);
}

inline void binomial_info_dinfo(int k, double y, const double* th,
                                double* out) {
    out[0] = binomial_expected_mu_mu(th[0], th[1]);
    out[1] = binomial_dexpected_mu_mu_mu(th[0], th[1]);
}

} // namespace d7

#endif
