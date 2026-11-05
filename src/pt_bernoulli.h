#ifndef D7_PT_BERNOULLI_H
#define D7_PT_BERNOULLI_H

// Bernoulli in its mean: one function per component and order for the
// quantities the scalar registry reads, called by the vector kernels
// (bernoulli.cpp, dexpected_kernels.cpp) and by the registry
// (d7_ccallable.cpp).

#include <Rcpp.h>
#include "pt_discrete_pow.h"
#include "pt_parts.h"

namespace d7 {

inline double bernoulli_score_mu(double y, double m) {
    return (y - m) / (m * (1.0 - m));
}

inline double bernoulli_hess_mu_mu(double y, double m) {
    return -coef_over_pow(y, m, 2) - coef_over_pow(1.0 - y, 1.0 - m, 2);
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

// The parts route (pt_parts.h). Bernoulli has one parameter, the mean mu,
// and no special function, so the buffer holds the inputs only:
//   [0] level the buffer was filled to
//   [1] y                                       level 0 (input)
//   [2] mu                                      level 0 (input)
// The third and fourth orders recompute om = mu - 1 in their own
// expressions, as bernoulli_hd.cpp does. The expected buffer is the same
// with y set to NaN.
inline void bernoulli_parts(double y, const double* th, int level, double* b) {
    b[0] = level; b[1] = y; b[2] = th[0];
}

inline void bernoulli_parts_expected(const double* th, int level, double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = th[0];
}

inline void bernoulli_parts_grad(const double* b, double* out) {
    out[0] = bernoulli_score_mu(b[1], b[2]);
}

inline void bernoulli_parts_hess(const double* b, double* out) {
    out[0] = bernoulli_hess_mu_mu(b[1], b[2]);
}

inline void bernoulli_parts_info(const double* b, double* out) {
    out[0] = bernoulli_expected_mu_mu(b[2]);
}

// bernoulli_deriv3_cpp and bernoulli_deriv4_cpp (bernoulli_hd.cpp)
inline void bernoulli_parts_d3(const double* b, double* out) {
    double om = b[2] - 1.0;
    out[0] = -2.0 * coef_over_pow(b[1] - 1.0, om, 3)
           + 2.0 * coef_over_pow(b[1], b[2], 3);
}

inline void bernoulli_parts_d4(const double* b, double* out) {
    double om = b[2] - 1.0;
    out[0] = 6.0 * coef_over_pow(b[1] - 1.0, om, 4)
           - 6.0 * coef_over_pow(b[1], b[2], 4);
}

} // namespace d7

#endif
