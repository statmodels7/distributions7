#ifndef D7_PT_GAUSSIAN3_H
#define D7_PT_GAUSSIAN3_H

#include <Rcpp.h>
#include "pt_parts.h"

// Gaussian in the mean and the precision tau: one function per component
// and order for the quantities the scalar registry reads, called by the
// vector kernels (gaussian3.cpp, dexpected_kernels.cpp) and by the registry
// (d7_ccallable.cpp). r = y - mu.

namespace d7 {

inline double gaussian3_score_mu(double r, double t) {
    return t * r;
}

inline double gaussian3_score_tau(double r, double t) {
    return 0.5 / t - 0.5 * r * r;
}

inline double gaussian3_hess_mu_mu(double t) {
    return -t;
}

inline double gaussian3_hess_tau_tau(double t) {
    return -0.5 / (t * t);
}

inline double gaussian3_expected_mu_mu(double t) {
    return -t;
}

inline double gaussian3_expected_tau_tau(double t) {
    return -0.5 / (t * t);
}

// E_mm = -tau does not depend on mu
inline double gaussian3_dexpected_mu_mu_mu() {
    return 0.0;
}

inline double gaussian3_dexpected_tau_tau_tau(double t) {
    double u = 1.0 / t, u3 = u * u * u;
    return u3;
}

inline void gaussian3_score_curv(int k, double y, const double* th,
                                 double* out) {
    double r = y - th[0], t = th[1];
    if (k == 0) {
        out[0] = gaussian3_score_mu(r, t);
        out[1] = gaussian3_hess_mu_mu(t);
    } else {
        out[0] = gaussian3_score_tau(r, t);
        out[1] = gaussian3_hess_tau_tau(t);
    }
}

inline void gaussian3_info_dinfo(int k, double y, const double* th,
                                 double* out) {
    double t = th[1];
    if (k == 0) {
        out[0] = gaussian3_expected_mu_mu(t);
        out[1] = gaussian3_dexpected_mu_mu_mu();
    } else {
        out[0] = gaussian3_expected_tau_tau(t);
        out[1] = gaussian3_dexpected_tau_tau_tau(t);
    }
}

// Mixed component of the observed Hessian, d2 l / d mu d tau = r, which
// gaussian3_hessian_cpp writes as y - mu.
inline double gaussian3_hess_mu_tau(double r) {
    return r;
}

// d E[l_mm] / d tau for E_mm = -tau, the kernel's w[1].
inline double gaussian3_dexpected_mu_mu_tau() {
    return -1.0;
}

// The routers of the joint registry, one per order (see pt_gaussian2.h).
inline void gaussian3_joint_grad(double y, const double* th, double* out) {
    double r = y - th[0], t = th[1];
    out[0] = gaussian3_score_mu(r, t);
    out[1] = gaussian3_score_tau(r, t);
}

inline void gaussian3_joint_hess(double y, const double* th, double* out) {
    double r = y - th[0], t = th[1];
    out[0] = gaussian3_hess_mu_mu(t);
    out[3] = gaussian3_hess_tau_tau(t);
    out[1] = out[2] = gaussian3_hess_mu_tau(r);
}

inline void gaussian3_joint_info(double y, const double* th, double* out) {
    double t = th[1];
    out[0] = gaussian3_expected_mu_mu(t);
    out[3] = gaussian3_expected_tau_tau(t);
    out[1] = out[2] = 0.0;
}

inline void gaussian3_joint_dinfo(double y, const double* th, double* out) {
    double t = th[1];
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[0] = gaussian3_dexpected_mu_mu_mu();
    out[4] = gaussian3_dexpected_mu_mu_tau();
    out[7] = gaussian3_dexpected_tau_tau_tau(t);
}

// The parts route (pt_parts.h). The Gaussian in the precision has no special
// function, so the buffer holds only inputs that every order reads, all cheap:
//   [0] level                       always
//   [1] y                           observed parts; NaN in the expected parts
//   [2] mu                          always
//   [3] tau                         always
//   [4] r = y - mu                  observed parts; NaN in the expected parts
// The expressions are those of the joint routers above and of
// gaussian3_deriv3_cpp and gaussian3_deriv4_cpp (gaussian3.cpp).
inline void gaussian3_parts(double y, const double* th, int level, double* b) {
    const double m = th[0], t = th[1];
    b[0] = level; b[1] = y; b[2] = m; b[3] = t; b[4] = y - m;
}

inline void gaussian3_parts_expected(const double* th, int level, double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = th[1]; b[4] = R_NaN;
}

inline void gaussian3_parts_grad(const double* b, double* out) {
    out[0] = gaussian3_score_mu(b[4], b[3]);
    out[1] = gaussian3_score_tau(b[4], b[3]);
}

inline void gaussian3_parts_hess(const double* b, double* out) {
    out[0] = gaussian3_hess_mu_mu(b[3]);
    out[3] = gaussian3_hess_tau_tau(b[3]);
    out[1] = out[2] = gaussian3_hess_mu_tau(b[4]);
}

// from the expected parts
inline void gaussian3_parts_info(const double* b, double* out) {
    out[0] = gaussian3_expected_mu_mu(b[3]);
    out[3] = gaussian3_expected_tau_tau(b[3]);
    out[1] = out[2] = 0.0;
}

// third order: gaussian3_deriv3_cpp, unique components in the sorted order
// mu mu mu, mu mu tau, mu tau tau, tau tau tau (the first and the third are
// zero, and the second is the constant -1)
inline void gaussian3_parts_d3(const double* b, double* out) {
    const double t = b[3];
    double w[4];
    w[0] = 0.0;
    w[1] = -1.0;
    w[2] = 0.0;
    w[3] = 1.0 / (t * t * t);
    sym_fill(2, 3, w, out);
}

// fourth order: gaussian3_deriv4_cpp, unique components in the sorted order
// mu^4, mu^3 tau, mu^2 tau^2, mu tau^3, tau^4 (only the last is non-zero)
inline void gaussian3_parts_d4(const double* b, double* out) {
    const double t = b[3];
    double t2 = t * t;
    double w[5];
    w[0] = 0.0;
    w[1] = 0.0;
    w[2] = 0.0;
    w[3] = 0.0;
    w[4] = -3.0 / (t2 * t2);
    sym_fill(2, 4, w, out);
}

} // namespace d7

#endif
