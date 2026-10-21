#ifndef D7_PT_GAUSSIAN2_H
#define D7_PT_GAUSSIAN2_H

// Gaussian in the mean and the variance v: one function per component and
// order for the quantities the scalar registry reads, called by the vector
// kernels (gaussian2.cpp, dexpected_kernels.cpp) and by the registry
// (d7_ccallable.cpp). With r = y - mu and u = 1/v. The lognormal in
// (mu, sigma2) shares the derivatives of this expected information
// (pt_lognormal1.h).

namespace d7 {

inline double gaussian2_score_mu(double r, double u) {
    return r * u;
}

inline double gaussian2_score_sigma2(double r, double v, double u) {
    double z2 = r * r / v;
    return 0.5 * (z2 - 1.0) * u;
}

inline double gaussian2_hess_mu_mu(double u) {
    return -u;
}

inline double gaussian2_hess_sigma2_sigma2(double r, double v, double u) {
    double u2 = u * u;
    double z2 = r * r / v;
    return (0.5 - z2) * u2;
}

inline double gaussian2_expected_mu_mu(double u) {
    return -u;
}

inline double gaussian2_expected_sigma2_sigma2(double u) {
    return -0.5 * u * u;
}

// E_mm = -1/v does not depend on mu
inline double gaussian2_dexpected_mu_mu_mu() {
    return 0.0;
}

inline double gaussian2_dexpected_sigma2_sigma2_sigma2(double u) {
    double u2 = u * u, u3 = u2 * u;
    return u3;
}

inline void gaussian2_score_curv(int k, double y, const double* th,
                                 double* out) {
    double v = th[1];
    double u = 1.0 / v;
    double r = y - th[0];
    if (k == 0) {
        out[0] = gaussian2_score_mu(r, u);
        out[1] = gaussian2_hess_mu_mu(u);
    } else {
        out[0] = gaussian2_score_sigma2(r, v, u);
        out[1] = gaussian2_hess_sigma2_sigma2(r, v, u);
    }
}

inline void gaussian2_info_dinfo(int k, double y, const double* th,
                                 double* out) {
    double u = 1.0 / th[1];
    if (k == 0) {
        out[0] = gaussian2_expected_mu_mu(u);
        out[1] = gaussian2_dexpected_mu_mu_mu();
    } else {
        out[0] = gaussian2_expected_sigma2_sigma2(u);
        out[1] = gaussian2_dexpected_sigma2_sigma2_sigma2(u);
    }
}

// Mixed component of the observed Hessian, d2 l / d mu d v = -r / v^2, with
// u = 1/v as gaussian2_hessian_cpp writes it.
inline double gaussian2_hess_mu_sigma2(double r, double u) {
    double u2 = u * u;
    return -r * u2;
}

// d E[l_mm] / d v for E_mm = -1/v, the kernel's w[1].
inline double gaussian2_dexpected_mu_mu_sigma2(double u) {
    double u2 = u * u;
    return u2;
}

// The routers of the joint registry, one per order: every component at one
// observation, on the parameter scale. The Hessian and the expected
// information fill a full P x P matrix and the derivative of the expected
// information a P x P x P array, column-major, out[a + P b + P^2 c] being
// the derivative of E[l_ab] in parameter c.
inline void gaussian2_joint_grad(double y, const double* th, double* out) {
    double v = th[1];
    double u = 1.0 / v;
    double r = y - th[0];
    out[0] = gaussian2_score_mu(r, u);
    out[1] = gaussian2_score_sigma2(r, v, u);
}

inline void gaussian2_joint_hess(double y, const double* th, double* out) {
    double v = th[1];
    double u = 1.0 / v;
    double r = y - th[0];
    out[0] = gaussian2_hess_mu_mu(u);
    out[3] = gaussian2_hess_sigma2_sigma2(r, v, u);
    out[1] = out[2] = gaussian2_hess_mu_sigma2(r, u);
}

inline void gaussian2_joint_info(double y, const double* th, double* out) {
    double u = 1.0 / th[1];
    out[0] = gaussian2_expected_mu_mu(u);
    out[3] = gaussian2_expected_sigma2_sigma2(u);
    out[1] = out[2] = 0.0;
}

inline void gaussian2_joint_dinfo(double y, const double* th, double* out) {
    double u = 1.0 / th[1];
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[0] = gaussian2_dexpected_mu_mu_mu();
    out[4] = gaussian2_dexpected_mu_mu_sigma2(u);
    out[7] = gaussian2_dexpected_sigma2_sigma2_sigma2(u);
}

} // namespace d7

#endif
