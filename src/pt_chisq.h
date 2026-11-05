#ifndef D7_PT_CHISQ_H
#define D7_PT_CHISQ_H

#include <Rcpp.h>
#include <cmath>
#include "pt_parts.h"

// Chi-squared in its mean (the degrees of freedom): one function per
// component and order for the quantities the scalar registry reads, called
// by the vector kernels (chisq.cpp, dexpected_kernels.cpp) and by the
// registry (d7_ccallable.cpp). The caller passes log y and the polygamma
// psi^(k)(mu/2) each function reads. From the second order on no
// derivative involves y, so the expected information is the observed one.

namespace d7 {

inline double chisq_score_mu(double log_y, double dg) {
    const double log2 = std::log(2.0);
    return 0.5 * (log_y - log2 - dg);
}

inline double chisq_hess_mu_mu(double tg) {
    return -0.25 * tg;
}

inline double chisq_expected_mu_mu(double tg) {
    return chisq_hess_mu_mu(tg);
}

inline double chisq_dexpected_mu_mu_mu(double pg2) {
    return -pg2 / 8.0;
}

inline void chisq_score_curv(int k, double y, const double* th,
                             double* out) {
    double h = 0.5 * th[0];
    out[0] = chisq_score_mu(std::log(y), R::digamma(h));
    out[1] = chisq_hess_mu_mu(R::trigamma(h));
}

inline void chisq_info_dinfo(int k, double y, const double* th,
                             double* out) {
    double h = 0.5 * th[0];
    out[0] = chisq_expected_mu_mu(R::trigamma(h));
    out[1] = chisq_dexpected_mu_mu_mu(R::psigamma(h, 2));
}

// The parts route (pt_parts.h). Chi-squared has one parameter, the degrees
// of freedom mu, and h = mu / 2. The costly quantities are the logarithm of
// y and the polygammas of h, each read by the orders named below; the
// buffer holds each once:
//   [0] level the buffer was filled to
//   [1] y                                       level 0 (input)
//   [2] mu                                      level 0 (input)
//   [3] h = mu / 2                              level 0 (cheap)
//   [4] log y                                   level 1 (order 1)
//   [5] digamma(h)                              level 1 (order 1)
//   [6] trigamma(h)                             level 2 (orders 2 and info)
//   [7] psigamma(h, 2)                          level 3 (order 3)
//   [8] psigamma(h, 3)                          level 4 (order 4)
// The expected buffer is the same with y and log y set to NaN, since no
// expected order reads them, and it fills [6] at level 2 and [7], [8] at
// their levels as the observed buffer does.
inline void chisq_parts(double y, const double* th, int level, double* b) {
    const double h = 0.5 * th[0];
    b[0] = level; b[1] = y; b[2] = th[0]; b[3] = h;
    b[4] = level >= 1 ? std::log(y) : R_NaN;
    b[5] = level >= 1 ? R::digamma(h) : R_NaN;
    b[6] = level >= 2 ? R::trigamma(h) : R_NaN;
    b[7] = level >= 3 ? R::psigamma(h, 2) : R_NaN;
    b[8] = level >= 4 ? R::psigamma(h, 3) : R_NaN;
}

inline void chisq_parts_expected(const double* th, int level, double* b) {
    const double h = 0.5 * th[0];
    b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = h;
    b[4] = R_NaN; b[5] = R_NaN;
    b[6] = level >= 2 ? R::trigamma(h) : R_NaN;
    b[7] = level >= 3 ? R::psigamma(h, 2) : R_NaN;
    b[8] = level >= 4 ? R::psigamma(h, 3) : R_NaN;
}

inline void chisq_parts_grad(const double* b, double* out) {
    out[0] = chisq_score_mu(b[4], b[5]);
}

inline void chisq_parts_hess(const double* b, double* out) {
    out[0] = chisq_hess_mu_mu(b[6]);
}

inline void chisq_parts_info(const double* b, double* out) {
    out[0] = chisq_expected_mu_mu(b[6]);
}

// chisq_deriv3_cpp and chisq_deriv4_cpp (chisq.cpp)
inline void chisq_parts_d3(const double* b, double* out) {
    out[0] = -b[7] / 8.0;
}

inline void chisq_parts_d4(const double* b, double* out) {
    out[0] = -b[8] / 16.0;
}

} // namespace d7

#endif
