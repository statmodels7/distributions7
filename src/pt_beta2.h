#ifndef D7_PT_BETA2_H
#define D7_PT_BETA2_H

#include <Rcpp.h>
#include <cmath>
#include "pt_parts.h"

// The beta in its two shapes (alpha, beta), s = alpha + beta: one function
// per component and order for the quantities the scalar registry reads,
// called by the vector kernels (beta2.cpp, dexpected_kernels.cpp) and by
// the registry (d7_ccallable.cpp). The caller passes log y, log(1 - y) and
// the polygammas a_k = psi^(k)(alpha), b_k = psi^(k)(beta),
// s_k = psi^(k)(s). The second derivatives do not involve y, so the
// expected information is the observed one.

namespace d7 {

inline double beta2_score_alpha(double log_y, double a0, double s0) {
    return log_y - a0 + s0;
}

inline double beta2_score_beta(double log_1_y, double b0, double s0) {
    return log_1_y - b0 + s0;
}

inline double beta2_hess_alpha_alpha(double a1, double s1) {
    return s1 - a1;
}

inline double beta2_hess_beta_beta(double b1, double s1) {
    return s1 - b1;
}

inline double beta2_expected_alpha_alpha(double a1, double s1) {
    return beta2_hess_alpha_alpha(a1, s1);
}

inline double beta2_expected_beta_beta(double b1, double s1) {
    return beta2_hess_beta_beta(b1, s1);
}

inline double beta2_dexpected_alpha_alpha_alpha(double a2, double s2) {
    return s2 - a2;
}

inline double beta2_dexpected_beta_beta_beta(double b2, double s2) {
    return s2 - b2;
}

inline void beta2_score_curv(int k, double y, const double* th,
                             double* out) {
    double a = th[0], b = th[1];
    double s0 = R::digamma(a + b), s1 = R::trigamma(a + b);
    if (k == 0) {
        out[0] = beta2_score_alpha(std::log(y), R::digamma(a), s0);
        out[1] = beta2_hess_alpha_alpha(R::trigamma(a), s1);
    } else {
        out[0] = beta2_score_beta(std::log1p(-y), R::digamma(b), s0);
        out[1] = beta2_hess_beta_beta(R::trigamma(b), s1);
    }
}

inline void beta2_info_dinfo(int k, double y, const double* th,
                             double* out) {
    double a = th[0], b = th[1];
    double s1 = R::trigamma(a + b), s2 = R::psigamma(a + b, 2);
    if (k == 0) {
        out[0] = beta2_expected_alpha_alpha(R::trigamma(a), s1);
        out[1] = beta2_dexpected_alpha_alpha_alpha(R::psigamma(a, 2), s2);
    } else {
        out[0] = beta2_expected_beta_beta(R::trigamma(b), s1);
        out[1] = beta2_dexpected_beta_beta_beta(R::psigamma(b, 2), s2);
    }
}

// The routers of the joint registry, one per order, for (alpha, beta) in
// that order. Each fills a full P x P matrix (P = 2) column-major, out[a + P b],
// and the derivative of the expected information a P x P x P array, out[a +
// P b + P P c] being the derivative of E[l_ab] in parameter c. The mixed
// entry is psi'(alpha + beta) in every order, and so are its derivatives
// in alpha and beta, psi''(alpha + beta), as beta2_dexpected1_cpp() forms
// them.
inline void beta2_joint_grad(double y, const double* th, double* out) {
    double a = th[0], b = th[1];
    double s0 = R::digamma(a + b);
    out[0] = beta2_score_alpha(std::log(y), R::digamma(a), s0);
    out[1] = beta2_score_beta(std::log1p(-y), R::digamma(b), s0);
}

inline void beta2_joint_hess(double y, const double* th, double* out) {
    double a = th[0], b = th[1];
    double s1 = R::trigamma(a + b);
    out[0] = beta2_hess_alpha_alpha(R::trigamma(a), s1);
    out[3] = beta2_hess_beta_beta(R::trigamma(b), s1);
    out[1] = out[2] = s1;
}

inline void beta2_joint_info(double y, const double* th, double* out) {
    double a = th[0], b = th[1];
    double s1 = R::trigamma(a + b);
    out[0] = beta2_expected_alpha_alpha(R::trigamma(a), s1);
    out[3] = beta2_expected_beta_beta(R::trigamma(b), s1);
    out[1] = out[2] = s1;
}

inline void beta2_joint_dinfo(double y, const double* th, double* out) {
    double a = th[0], b = th[1];
    double ps = R::psigamma(a + b, 2), pa = R::psigamma(a, 2),
           pb = R::psigamma(b, 2);
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[0] = beta2_dexpected_alpha_alpha_alpha(pa, ps);
    out[4] = ps;
    out[2] = out[1] = ps;
    out[6] = out[5] = ps;
    out[3] = ps;
    out[7] = beta2_dexpected_beta_beta_beta(pb, ps);
}

// The parts route (pt_parts.h). The buffer, with the logarithms, the
// digammas, the trigammas and the higher polygammas of alpha, beta and
// s = alpha + beta, each computed once at the level that first reads it:
//
//   [0]  level
//   [1]  y                          NaN in the expected buffer
//   [2]  alpha
//   [3]  beta
//   [4]  log y                      level 1
//   [5]  log(1 - y)                 level 1
//   [6]  digamma(alpha)             level 1
//   [7]  digamma(beta)              level 1
//   [8]  digamma(s)                 level 1
//   [9]  trigamma(alpha)            level 2
//   [10] trigamma(beta)             level 2
//   [11] trigamma(s)                level 2
//   [12] psigamma(alpha, 2)         level 3
//   [13] psigamma(beta, 2)          level 3
//   [14] psigamma(s, 2)             level 3
//   [15] psigamma(alpha, 3)         level 4
//   [16] psigamma(beta, 3)          level 4
//   [17] psigamma(s, 3)             level 4
//
// The expected buffer fills levels 0 and 1 with NaN and level 2 with [9] to
// [11]; the expected information is the observed one.
inline void beta2_parts(double y, const double* th, int level, double* b) {
    const double a = th[0], bb = th[1];
    const double s = a + bb;
    b[0] = level; b[1] = y; b[2] = a; b[3] = bb;
    b[4] = level >= 1 ? std::log(y) : R_NaN;
    b[5] = level >= 1 ? std::log1p(-y) : R_NaN;
    b[6] = level >= 1 ? R::digamma(a) : R_NaN;
    b[7] = level >= 1 ? R::digamma(bb) : R_NaN;
    b[8] = level >= 1 ? R::digamma(s) : R_NaN;
    b[9] = level >= 2 ? R::trigamma(a) : R_NaN;
    b[10] = level >= 2 ? R::trigamma(bb) : R_NaN;
    b[11] = level >= 2 ? R::trigamma(s) : R_NaN;
    b[12] = level >= 3 ? R::psigamma(a, 2.0) : R_NaN;
    b[13] = level >= 3 ? R::psigamma(bb, 2.0) : R_NaN;
    b[14] = level >= 3 ? R::psigamma(s, 2.0) : R_NaN;
    b[15] = level >= 4 ? R::psigamma(a, 3.0) : R_NaN;
    b[16] = level >= 4 ? R::psigamma(bb, 3.0) : R_NaN;
    b[17] = level >= 4 ? R::psigamma(s, 3.0) : R_NaN;
}

inline void beta2_parts_expected(const double* th, int level, double* b) {
    const double a = th[0], bb = th[1];
    const double s = a + bb;
    b[0] = level; b[1] = R_NaN; b[2] = a; b[3] = bb;
    for (int j = 4; j <= 8; ++j) b[j] = R_NaN;
    b[9] = level >= 2 ? R::trigamma(a) : R_NaN;
    b[10] = level >= 2 ? R::trigamma(bb) : R_NaN;
    b[11] = level >= 2 ? R::trigamma(s) : R_NaN;
    for (int j = 12; j <= 17; ++j) b[j] = R_NaN;
}

inline void beta2_parts_grad(const double* b, double* out) {
    out[0] = beta2_score_alpha(b[4], b[6], b[8]);
    out[1] = beta2_score_beta(b[5], b[7], b[8]);
}

inline void beta2_parts_hess(const double* b, double* out) {
    out[0] = beta2_hess_alpha_alpha(b[9], b[11]);
    out[3] = beta2_hess_beta_beta(b[10], b[11]);
    out[1] = out[2] = b[11];
}

// from the expected parts
inline void beta2_parts_info(const double* b, double* out) {
    out[0] = beta2_expected_alpha_alpha(b[9], b[11]);
    out[3] = beta2_expected_beta_beta(b[10], b[11]);
    out[1] = out[2] = b[11];
}

// third order: no compiled kernel of its own, the R closed form of
// beta2_higher(theta, n, 3) written in C++: alpha^3 = psi''(s) - psi''(alpha),
// the two mixed components psi''(s), beta^3 = psi''(s) - psi''(beta)
inline void beta2_parts_d3(const double* b, double* out) {
    double u[4];
    u[0] = b[14] - b[12];
    u[1] = b[14];
    u[2] = b[14];
    u[3] = b[14] - b[13];
    sym_fill(2, 3, u, out);
}

// fourth order: no compiled kernel of its own, the R closed form of
// beta2_higher(theta, n, 4): alpha^4 = psi'''(s) - psi'''(alpha), the three
// mixed components psi'''(s), beta^4 = psi'''(s) - psi'''(beta)
inline void beta2_parts_d4(const double* b, double* out) {
    double u[5];
    u[0] = b[17] - b[15];
    u[1] = b[17];
    u[2] = b[17];
    u[3] = b[17];
    u[4] = b[17] - b[16];
    sym_fill(2, 4, u, out);
}

} // namespace d7

#endif
