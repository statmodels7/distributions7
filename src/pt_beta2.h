#ifndef D7_PT_BETA2_H
#define D7_PT_BETA2_H

#include <Rcpp.h>
#include <cmath>

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

} // namespace d7

#endif
