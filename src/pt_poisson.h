#ifndef D7_PT_POISSON_H
#define D7_PT_POISSON_H

#include <Rcpp.h>

// Poisson in the mean, l = y log mu - mu - log y!: one function per
// component and order for the quantities the scalar registry reads, called
// by the vector kernels in poisson.cpp and by the registry in
// d7_ccallable.cpp.

namespace d7 {

inline double poisson_score_mu(double y, double m) {
    return (y - m) / m;
}

inline double poisson_hess_mu_mu(double y, double m) {
    double m2 = m * m;
    return -y / m2;
}

inline double poisson_expected_mu_mu(double m) {
    return -1.0 / m;
}

inline double poisson_dexpected_mu_mu_mu(double m) {
    double inv = 1.0 / m;
    return inv * inv;
}

inline void poisson_score_curv(int k, double y, const double* th,
                               double* out) {
    out[0] = poisson_score_mu(y, th[0]);
    out[1] = poisson_hess_mu_mu(y, th[0]);
}

inline void poisson_info_dinfo(int k, double y, const double* th,
                               double* out) {
    out[0] = poisson_expected_mu_mu(th[0]);
    out[1] = poisson_dexpected_mu_mu_mu(th[0]);
}

} // namespace d7

#endif
