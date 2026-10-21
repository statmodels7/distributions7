#ifndef D7_PT_GAUSSIAN3_H
#define D7_PT_GAUSSIAN3_H

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

} // namespace d7

#endif
