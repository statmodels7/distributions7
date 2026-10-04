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

} // namespace d7

#endif
