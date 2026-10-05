#ifndef D7_PT_INVGAUSS1_H
#define D7_PT_INVGAUSS1_H

// Inverse gaussian in the mean and the dispersion phi, Var = phi mu^3: one
// function per component and order for the quantities the scalar registry
// reads, called by the vector kernels (invgauss.cpp, dexpected_kernels.cpp)
// and by the registry (d7_ccallable.cpp).

namespace d7 {

inline double invgauss1_score_mu(double y, double m, double p) {
    double m2 = m * m;
    double m3 = m2 * m;
    double res = y - m;
    return res / (p * m3);
}

inline double invgauss1_score_phi(double y, double m, double p) {
    double m2 = m * m;
    double p2 = p * p;
    double res = y - m;
    return (res * res - y * m2 * p) / (2.0 * y * p2 * m2);
}

inline double invgauss1_hess_mu_mu(double y, double m, double p) {
    double m2 = m * m;
    double m4 = m2 * m2;
    return -(3.0 * y - 2.0 * m) / (p * m4);
}

inline double invgauss1_hess_phi_phi(double y, double m, double p) {
    double m2 = m * m;
    double p2 = p * p;
    double p3 = p2 * p;
    double res = y - m;
    double res2 = res * res;
    return (p - 2.0 * res2 / (m2 * y)) / (2.0 * p3);
}

inline double invgauss1_expected_mu_mu(double m, double p) {
    double m3 = m * m * m;
    return -1.0 / (p * m3);
}

inline double invgauss1_expected_phi_phi(double p) {
    double p2 = p * p;
    return -0.5 / p2;
}

inline double invgauss1_dexpected_mu_mu_mu(double m, double p) {
    double im = 1.0 / m, im3 = im * im * im, im4 = im3 * im;
    double ip = 1.0 / p;
    return 3.0 * ip * im4;
}

inline double invgauss1_dexpected_phi_phi_phi(double p) {
    double ip = 1.0 / p, ip2 = ip * ip, ip3 = ip2 * ip;
    return ip3;
}

inline void invgauss1_score_curv(int k, double y, const double* th,
                                 double* out) {
    double m = th[0], p = th[1];
    if (k == 0) {
        out[0] = invgauss1_score_mu(y, m, p);
        out[1] = invgauss1_hess_mu_mu(y, m, p);
    } else {
        out[0] = invgauss1_score_phi(y, m, p);
        out[1] = invgauss1_hess_phi_phi(y, m, p);
    }
}

inline void invgauss1_info_dinfo(int k, double y, const double* th,
                                 double* out) {
    double m = th[0], p = th[1];
    if (k == 0) {
        out[0] = invgauss1_expected_mu_mu(m, p);
        out[1] = invgauss1_dexpected_mu_mu_mu(m, p);
    } else {
        out[0] = invgauss1_expected_phi_phi(p);
        out[1] = invgauss1_dexpected_phi_phi_phi(p);
    }
}

} // namespace d7

#endif
