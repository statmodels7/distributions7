#ifndef D7_PT_BETA1_H
#define D7_PT_BETA1_H

#include <Rcpp.h>
#include <cmath>
#include "pt_parts.h"

// Beta in the mean and the precision, with shapes alpha = mu phi and
// beta = (1 - mu) phi: one function per component and order for the
// quantities the scalar registry reads, called by the vector kernels in
// beta.cpp and by the registry in d7_ccallable.cpp. The functions are
// arithmetic only; the caller computes and passes the special functions,
// written a_k = psi^(k)(alpha), b_k = psi^(k)(beta), c_k = psi^(k)(phi),
// together with log y and log(1 - y).

namespace d7 {

inline double beta1_score_mu(double p, double log_y, double log_1_y,
                             double a0, double b0) {
    double log_ratio = log_y - log_1_y;
    return p * (log_ratio - a0 + b0);
}

inline double beta1_score_phi(double m, double log_y, double log_1_y,
                              double c0, double a0, double b0) {
    return c0 - m * a0 - (1.0 - m) * b0 + m * log_y + (1.0 - m) * log_1_y;
}

inline double beta1_hess_mu_mu(double p, double a1, double b1) {
    return -p * p * (a1 + b1);
}

inline double beta1_hess_phi_phi(double m, double c1, double a1, double b1) {
    return c1 - m * m * a1 - (1.0 - m) * (1.0 - m) * b1;
}

// the diagonal of the Hessian does not involve y, so the expected
// information's diagonal is the same expression
inline double beta1_expected_mu_mu(double p, double a1, double b1) {
    return beta1_hess_mu_mu(p, a1, b1);
}

inline double beta1_expected_phi_phi(double m, double c1, double a1,
                                     double b1) {
    return beta1_hess_phi_phi(m, c1, a1, b1);
}

// d alpha/d mu = phi, d beta/d mu = -phi, d alpha/d phi = mu,
// d beta/d phi = 1 - mu
inline double beta1_dexpected_mu_mu_mu(double p, double a2, double b2) {
    return -p * p * p * (a2 - b2);
}

inline double beta1_dexpected_phi_phi_phi(double m, double c2, double a2,
                                          double b2) {
    double n = 1.0 - m;
    return c2 - m * m * m * a2 - n * n * n * b2;
}

inline void beta1_score_curv(int k, double y, const double* th,
                             double* out) {
    double m = th[0], p = th[1];
    double al = m * p, be = (1.0 - m) * p;
    double log_y = std::log(y), log_1_y = std::log(1.0 - y);
    double a0 = R::digamma(al), b0 = R::digamma(be);
    double a1 = R::trigamma(al), b1 = R::trigamma(be);
    if (k == 0) {
        out[0] = beta1_score_mu(p, log_y, log_1_y, a0, b0);
        out[1] = beta1_hess_mu_mu(p, a1, b1);
    } else {
        out[0] = beta1_score_phi(m, log_y, log_1_y, R::digamma(p), a0, b0);
        out[1] = beta1_hess_phi_phi(m, R::trigamma(p), a1, b1);
    }
}

inline void beta1_info_dinfo(int k, double y, const double* th,
                             double* out) {
    double m = th[0], p = th[1];
    double al = m * p, be = (1.0 - m) * p;
    double a1 = R::trigamma(al), b1 = R::trigamma(be);
    double a2 = R::psigamma(al, 2), b2 = R::psigamma(be, 2);
    if (k == 0) {
        out[0] = beta1_expected_mu_mu(p, a1, b1);
        out[1] = beta1_dexpected_mu_mu_mu(p, a2, b2);
    } else {
        out[0] = beta1_expected_phi_phi(m, R::trigamma(p), a1, b1);
        out[1] = beta1_dexpected_phi_phi_phi(m, R::psigamma(p, 2), a2, b2);
    }
}

// The mixed observed entry, as beta_hessian_cpp() forms it: the log ratio of
// the response is log(y/(1 - y)) there, not log(y) - log(1 - y).
inline double beta1_hess_mu_phi(double y, double m, double p, double a0,
                                double b0, double a1, double b1) {
    double log_ratio = std::log(y / (1.0 - y));
    double term1 = log_ratio - a0 + b0;
    double term2 = p * (m * a1 - (1.0 - m) * b1);
    return term1 - term2;
}

// The mixed expected entry, E[l_mu,phi] = -phi (mu a1 - (1 - mu) b1).
inline double beta1_expected_mu_phi(double m, double p, double a1, double b1) {
    return -p * (m * a1 - (1.0 - m) * b1);
}

// d_mu E_mm is beta1_dexpected_mu_mu_mu(); d_phi E_mm, as beta_dexpected1_cpp()
// forms it.
inline double beta1_dexpected_mu_mu_phi(double m, double p, double a1,
                                        double b1, double a2, double b2) {
    double n = 1.0 - m;
    return -2.0 * p * (a1 + b1) - p * p * (m * a2 + n * b2);
}

// d_mu E_pp, as beta_dexpected1_cpp() forms it.
inline double beta1_dexpected_phi_phi_mu(double m, double p, double a1,
                                         double b1, double a2, double b2) {
    double n = 1.0 - m;
    return -2.0 * m * a1 - m * m * p * a2 + 2.0 * n * b1 + n * n * p * b2;
}

// d_mu E_mp, as beta_dexpected1_cpp() forms it.
inline double beta1_dexpected_mu_phi_mu(double m, double p, double a1,
                                        double b1, double a2, double b2) {
    double n = 1.0 - m;
    return -p * (a1 + b1) - p * p * (m * a2 + n * b2);
}

// d_phi E_mp, as beta_dexpected1_cpp() forms it.
inline double beta1_dexpected_mu_phi_phi(double m, double p, double a1,
                                         double b1, double a2, double b2) {
    double n = 1.0 - m;
    return -(m * a1 - n * b1) - p * (m * m * a2 - n * n * b2);
}

// The routers of the joint registry, one per order, for (mu, phi) in that
// order. Each fills a full P x P matrix (P = 2) column-major, out[a + P b],
// and the derivative of the expected information a P x P x P array, out[a +
// P b + P P c] being the derivative of E[l_ab] in parameter c. The
// polygammas are formed as the kernels form them, from alpha = mu phi and
// beta = (1 - mu) phi.
inline void beta1_joint_grad(double y, const double* th, double* out) {
    double m = th[0], p = th[1];
    double al = m * p, be = (1.0 - m) * p;
    double log_y = std::log(y), log_1_y = std::log(1.0 - y);
    double a0 = R::digamma(al), b0 = R::digamma(be);
    out[0] = beta1_score_mu(p, log_y, log_1_y, a0, b0);
    out[1] = beta1_score_phi(m, log_y, log_1_y, R::digamma(p), a0, b0);
}

inline void beta1_joint_hess(double y, const double* th, double* out) {
    double m = th[0], p = th[1];
    double al = m * p, be = (1.0 - m) * p;
    double a0 = R::digamma(al), b0 = R::digamma(be);
    double a1 = R::trigamma(al), b1 = R::trigamma(be);
    out[0] = beta1_hess_mu_mu(p, a1, b1);
    out[3] = beta1_hess_phi_phi(m, R::trigamma(p), a1, b1);
    out[1] = out[2] = beta1_hess_mu_phi(y, m, p, a0, b0, a1, b1);
}

inline void beta1_joint_info(double y, const double* th, double* out) {
    double m = th[0], p = th[1];
    double al = m * p, be = (1.0 - m) * p;
    double a1 = R::trigamma(al), b1 = R::trigamma(be);
    out[0] = beta1_expected_mu_mu(p, a1, b1);
    out[3] = beta1_expected_phi_phi(m, R::trigamma(p), a1, b1);
    out[1] = out[2] = beta1_expected_mu_phi(m, p, a1, b1);
}

inline void beta1_joint_dinfo(double y, const double* th, double* out) {
    double m = th[0], p = th[1];
    double n = 1.0 - m, al = m * p, be = n * p;
    double a1 = R::trigamma(al), b1 = R::trigamma(be);
    double a2 = R::psigamma(al, 2), b2 = R::psigamma(be, 2);
    double c2 = R::psigamma(p, 2);
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[0] = beta1_dexpected_mu_mu_mu(p, a2, b2);
    out[4] = beta1_dexpected_mu_mu_phi(m, p, a1, b1, a2, b2);
    out[2] = out[1] = beta1_dexpected_mu_phi_mu(m, p, a1, b1, a2, b2);
    out[6] = out[5] = beta1_dexpected_mu_phi_phi(m, p, a1, b1, a2, b2);
    out[3] = beta1_dexpected_phi_phi_mu(m, p, a1, b1, a2, b2);
    out[7] = beta1_dexpected_phi_phi_phi(m, c2, a2, b2);
}

// The parts route (pt_parts.h). The buffer, with mu = th[0], phi = th[1],
// al = mu phi and be = (1 - mu) phi as the joint routers form them, and
// b' = phi - mu phi as the vector kernels beta_deriv3_cpp() and
// beta_deriv4_cpp() form the second shape. The two routes give the same value
// to rounding, and each order keeps its own route, so the second shape has
// its own slots:
//
//   [0]  level
//   [1]  y                          NaN in the expected buffer
//   [2]  mu
//   [3]  phi
//   [4]  log y                      level 1
//   [5]  log(1 - y)                 level 1
//   [6]  digamma(al)                level 1
//   [7]  digamma(be)                level 1
//   [8]  digamma(phi)               level 1
//   [9]  trigamma(al)               level 2
//   [10] trigamma(be)               level 2
//   [11] trigamma(phi)              level 2
//   [12] trigamma(b')               level 3
//   [13] psigamma(al, 2)            level 3
//   [14] psigamma(b', 2)            level 3
//   [15] psigamma(phi, 2)           level 3
//   [16] psigamma(al, 3)            level 4
//   [17] psigamma(b', 3)            level 4
//   [18] psigamma(phi, 3)           level 4
//
// The expected buffer fills levels 0 and 1 with NaN and level 2 with [9] to
// [11]; the expected Hessian is the only expected order read.
inline void beta1_parts(double y, const double* th, int level, double* b) {
    const double m = th[0], p = th[1];
    const double al = m * p, be = (1.0 - m) * p, bp = p - m * p;
    b[0] = level; b[1] = y; b[2] = m; b[3] = p;
    b[4] = level >= 1 ? std::log(y) : R_NaN;
    b[5] = level >= 1 ? std::log(1.0 - y) : R_NaN;
    b[6] = level >= 1 ? R::digamma(al) : R_NaN;
    b[7] = level >= 1 ? R::digamma(be) : R_NaN;
    b[8] = level >= 1 ? R::digamma(p) : R_NaN;
    b[9] = level >= 2 ? R::trigamma(al) : R_NaN;
    b[10] = level >= 2 ? R::trigamma(be) : R_NaN;
    b[11] = level >= 2 ? R::trigamma(p) : R_NaN;
    b[12] = level >= 3 ? R::trigamma(bp) : R_NaN;
    b[13] = level >= 3 ? R::psigamma(al, 2.0) : R_NaN;
    b[14] = level >= 3 ? R::psigamma(bp, 2.0) : R_NaN;
    b[15] = level >= 3 ? R::psigamma(p, 2.0) : R_NaN;
    b[16] = level >= 4 ? R::psigamma(al, 3.0) : R_NaN;
    b[17] = level >= 4 ? R::psigamma(bp, 3.0) : R_NaN;
    b[18] = level >= 4 ? R::psigamma(p, 3.0) : R_NaN;
}

inline void beta1_parts_expected(const double* th, int level, double* b) {
    const double m = th[0], p = th[1];
    const double al = m * p, be = (1.0 - m) * p;
    b[0] = level; b[1] = R_NaN; b[2] = m; b[3] = p;
    for (int j = 4; j <= 8; ++j) b[j] = R_NaN;
    b[9] = level >= 2 ? R::trigamma(al) : R_NaN;
    b[10] = level >= 2 ? R::trigamma(be) : R_NaN;
    b[11] = level >= 2 ? R::trigamma(p) : R_NaN;
    for (int j = 12; j <= 18; ++j) b[j] = R_NaN;
}

inline void beta1_parts_grad(const double* b, double* out) {
    out[0] = beta1_score_mu(b[3], b[4], b[5], b[6], b[7]);
    out[1] = beta1_score_phi(b[2], b[4], b[5], b[8], b[6], b[7]);
}

inline void beta1_parts_hess(const double* b, double* out) {
    const double y = b[1], m = b[2], p = b[3];
    out[0] = beta1_hess_mu_mu(p, b[9], b[10]);
    out[3] = beta1_hess_phi_phi(m, b[11], b[9], b[10]);
    out[1] = out[2] = beta1_hess_mu_phi(y, m, p, b[6], b[7], b[9], b[10]);
}

// from the expected parts
inline void beta1_parts_info(const double* b, double* out) {
    const double m = b[2], p = b[3];
    out[0] = beta1_expected_mu_mu(p, b[9], b[10]);
    out[3] = beta1_expected_phi_phi(m, b[11], b[9], b[10]);
    out[1] = out[2] = beta1_expected_mu_phi(m, p, b[9], b[10]);
}

// beta_deriv3_cpp() of beta.cpp's vector kernels, the same expressions
inline void beta1_parts_d3(const double* b, double* out) {
    const double m = b[2], p = b[3];
    const double p2 = p * p, p3 = p2 * p;
    const double Ta = b[9], Tb = b[12];
    const double P2a = b[13], P2b = b[14], P2phi = b[15];
    double u[4];
    u[0] = p3 * (-P2a + P2b);
    u[1] = p * (-(m * p * P2a) + (m - 1.0) * p * P2b - 2.0 * (Ta + Tb));
    u[2] = -(m * m * p * P2a) + (m - 1.0) * (m - 1.0) * p * P2b - 2.0 * m * Ta - 2.0 * (m - 1.0) * Tb;
    u[3] = P2phi - m * m * m * P2a + (m - 1.0) * (m - 1.0) * (m - 1.0) * P2b;
    sym_fill(2, 3, u, out);
}

// beta_deriv4_cpp() of beta.cpp's vector kernels, the same expressions
inline void beta1_parts_d4(const double* b, double* out) {
    const double m = b[2], p = b[3];
    const double om = m - 1.0;
    const double p2 = p * p, p4 = p2 * p2;
    const double Ta = b[9], Tb = b[12];
    const double P2a = b[13], P2b = b[14], P2phi = b[15];
    const double P3a = b[16], P3b = b[17], P3phi = b[18];
    double u[5];
    u[0] = -(p4 * (P3a + P3b));
    u[1] = p2 * (-3.0 * P2a - m * p * P3a + 3.0 * P2b - om * p * P3b);
    u[2] = p * (-4.0 * m * P2a - m * m * p * P3a + om * (4.0 * P2b - om * p * P3b)) - 2.0 * (Ta + Tb);
    u[3] = -3.0 * m * m * P2a - m * m * m * p * P3a + om * om * (3.0 * P2b - om * p * P3b);
    u[4] = P3phi - m * m * m * m * P3a - om * om * om * om * P3b;
    sym_fill(2, 4, u, out);
}

} // namespace d7

#endif
