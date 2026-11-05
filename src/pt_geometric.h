#ifndef D7_PT_GEOMETRIC_H
#define D7_PT_GEOMETRIC_H

// Geometric on {0, 1, ...} in its mean, l = y log mu - (y + 1) log(1 + mu):
// one function per component and order for the quantities the scalar
// registry reads, called by the vector kernels (geometric.cpp,
// dexpected_kernels.cpp) and by the registry (d7_ccallable.cpp).

#include <Rcpp.h>
#include "pt_discrete_pow.h"
#include "pt_parts.h"

namespace d7 {

inline double geometric_score_mu(double y, double m) {
    return (y - m) / (m * (1.0 + m));
}

inline double geometric_hess_mu_mu(double y, double m) {
    return -count_bracket(y, m, 1.0, 2);
}

inline double geometric_expected_mu_mu(double m) {
    return -1.0 / (m * (1.0 + m));
}

// E = -1/q with q = m(1 + m)
inline double geometric_dexpected_mu_mu_mu(double m) {
    double q = m * (1.0 + m), s = 1.0 + 2.0 * m, q2 = q * q;
    return s / q2;
}

inline void geometric_score_curv(int k, double y, const double* th,
                                 double* out) {
    out[0] = geometric_score_mu(y, th[0]);
    out[1] = geometric_hess_mu_mu(y, th[0]);
}

inline void geometric_info_dinfo(int k, double y, const double* th,
                                 double* out) {
    out[0] = geometric_expected_mu_mu(th[0]);
    out[1] = geometric_dexpected_mu_mu_mu(th[0]);
}

// The parts route (pt_parts.h). Geometric has one parameter, the mean mu.
// The costly quantity is L = log1p(1 / mu), which count_bracket() forms for
// the second, third and fourth orders; the buffer holds it once. The buffer:
//   [0] level the buffer was filled to
//   [1] y                                       level 0 (input)
//   [2] mu                                      level 0 (input)
//   [3] L = log1p(1 / mu)                       level 2 (read by orders 2-4)
// The expected buffer is the same with y set to NaN and L set to NaN, since
// the expected Hessian reads no log1p.

// count_bracket(y, mu, 1, k) with L = log1p(1 / mu) passed in, the same
// operations in the same order
inline double geometric_bracket_parts(double y, double m, double L, int k) {
    double one_mq = -std::expm1(-k * L);
    return coef_over_pow(y * one_mq, m, k) - coef_over_pow(1.0, m + 1.0, k);
}

inline void geometric_parts(double y, const double* th, int level, double* b) {
    const double m = th[0];
    b[0] = level; b[1] = y; b[2] = m;
    b[3] = level >= 2 ? std::log1p(1.0 / m) : R_NaN;
}

inline void geometric_parts_expected(const double* th, int level, double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = R_NaN;
}

inline void geometric_parts_grad(const double* b, double* out) {
    out[0] = geometric_score_mu(b[1], b[2]);
}

inline void geometric_parts_hess(const double* b, double* out) {
    out[0] = -geometric_bracket_parts(b[1], b[2], b[3], 2);
}

inline void geometric_parts_info(const double* b, double* out) {
    out[0] = geometric_expected_mu_mu(b[2]);
}

// geometric_deriv3_cpp and geometric_deriv4_cpp (geometric.cpp)
inline void geometric_parts_d3(const double* b, double* out) {
    out[0] = 2.0 * geometric_bracket_parts(b[1], b[2], b[3], 3);
}

inline void geometric_parts_d4(const double* b, double* out) {
    out[0] = -6.0 * geometric_bracket_parts(b[1], b[2], b[3], 4);
}

} // namespace d7

#endif
