#ifndef D7_PT_POISSON_H
#define D7_PT_POISSON_H

#include <Rcpp.h>

// Poisson in the mean, l = y log mu - mu - log y!: one function per
// component and order for the quantities the scalar registry reads, called
// by the vector kernels in poisson.cpp and by the registry in
// d7_ccallable.cpp.

#include "pt_discrete_pow.h"
#include "pt_parts.h"

namespace d7 {

inline double poisson_score_mu(double y, double m) {
    return (y - m) / m;
}

inline double poisson_hess_mu_mu(double y, double m) {
    return -coef_over_pow(y, m, 2);
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

// The parts route (pt_parts.h). Poisson has one parameter, the mean mu, and
// no special function, so the buffer holds the inputs only:
//   [0] level the buffer was filled to
//   [1] y                                       level 0 (input)
//   [2] mu                                      level 0 (input)
// Every order is arithmetic on y and mu. The expected buffer is the same
// with y set to NaN, since no expected order reads it.
inline void poisson_parts(double y, const double* th, int level, double* b) {
    b[0] = level; b[1] = y; b[2] = th[0];
}

inline void poisson_parts_expected(const double* th, int level, double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = th[0];
}

inline void poisson_parts_grad(const double* b, double* out) {
    out[0] = poisson_score_mu(b[1], b[2]);
}

inline void poisson_parts_hess(const double* b, double* out) {
    out[0] = poisson_hess_mu_mu(b[1], b[2]);
}

inline void poisson_parts_info(const double* b, double* out) {
    out[0] = poisson_expected_mu_mu(b[2]);
}

// poisson_deriv3_cpp and poisson_deriv4_cpp (poisson_hd.cpp)
inline void poisson_parts_d3(const double* b, double* out) {
    out[0] = 2.0 * coef_over_pow(b[1], b[2], 3);
}

inline void poisson_parts_d4(const double* b, double* out) {
    out[0] = -6.0 * coef_over_pow(b[1], b[2], 4);
}

} // namespace d7

#endif
