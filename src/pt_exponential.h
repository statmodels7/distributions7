#ifndef D7_PT_EXPONENTIAL_H
#define D7_PT_EXPONENTIAL_H

// Exponential in its mean, l = -log mu - y/mu: one function per component
// and order for the quantities the scalar registry reads, called by the
// vector kernels (exponential.cpp, dexpected_kernels.cpp) and by the
// registry (d7_ccallable.cpp).

#include <Rcpp.h>
#include "pt_parts.h"

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

// The parts route (pt_parts.h). Exponential has one parameter, the mean mu,
// and no special function, so the buffer holds the inputs only:
//   [0] level the buffer was filled to
//   [1] y                                       level 0 (input)
//   [2] mu                                      level 0 (input)
// The powers of mu are cheap and are formed inside each order, as
// exponential_deriv3_cpp and exponential_deriv4_cpp form them. The expected
// buffer is the same with y set to NaN.
inline void exponential_parts(double y, const double* th, int level,
                              double* b) {
    b[0] = level; b[1] = y; b[2] = th[0];
}

inline void exponential_parts_expected(const double* th, int level,
                                       double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = th[0];
}

inline void exponential_parts_grad(const double* b, double* out) {
    out[0] = exponential_score_mu(b[1], b[2]);
}

inline void exponential_parts_hess(const double* b, double* out) {
    out[0] = exponential_hess_mu_mu(b[1], b[2]);
}

inline void exponential_parts_info(const double* b, double* out) {
    out[0] = exponential_expected_mu_mu(b[2]);
}

// exponential_deriv3_cpp and exponential_deriv4_cpp (exponential.cpp)
inline void exponential_parts_d3(const double* b, double* out) {
    const double m = b[2], y = b[1];
    double m3 = m * m * m;
    out[0] = -2.0 / m3 + 6.0 * y / (m3 * m);
}

inline void exponential_parts_d4(const double* b, double* out) {
    const double m = b[2], y = b[1];
    double m2 = m * m;
    double m4 = m2 * m2;
    out[0] = 6.0 / m4 - 24.0 * y / (m4 * m);
}

} // namespace d7

#endif
