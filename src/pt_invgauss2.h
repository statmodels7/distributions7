#ifndef D7_PT_INVGAUSS2_H
#define D7_PT_INVGAUSS2_H

// Inverse gaussian in the mean and the shape lambda, Var = mu^3/lambda: one
// function per component and order for the quantities the scalar registry
// reads, called by the vector kernels (invgauss2.cpp, dexpected_kernels.cpp)
// and by the registry (d7_ccallable.cpp).

namespace d7 {

inline double invgauss2_score_mu(double y, double m, double L) {
    double r = y - m, m2 = m * m, m3 = m2 * m;
    return L * r / m3;
}

inline double invgauss2_score_lambda(double y, double m, double L) {
    double r = y - m, m2 = m * m;
    return 0.5 / L - r * r / (2.0 * m2 * y);
}

inline double invgauss2_hess_mu_mu(double y, double m, double L) {
    double m2 = m * m, m4 = m2 * m2;
    return L * (2.0 * m - 3.0 * y) / m4;
}

inline double invgauss2_hess_lambda_lambda(double L) {
    return -0.5 / (L * L);
}

inline double invgauss2_expected_mu_mu(double m, double L) {
    return -L / (m * m * m);
}

inline double invgauss2_expected_lambda_lambda(double L) {
    return -0.5 / (L * L);
}

inline double invgauss2_dexpected_mu_mu_mu(double m, double L) {
    double im = 1.0 / m, im3 = im * im * im, im4 = im3 * im;
    return 3.0 * L * im4;
}

inline double invgauss2_dexpected_lambda_lambda_lambda(double L) {
    double iL = 1.0 / L, iL3 = iL * iL * iL;
    return iL3;
}

inline void invgauss2_score_curv(int k, double y, const double* th,
                                 double* out) {
    double m = th[0], L = th[1];
    if (k == 0) {
        out[0] = invgauss2_score_mu(y, m, L);
        out[1] = invgauss2_hess_mu_mu(y, m, L);
    } else {
        out[0] = invgauss2_score_lambda(y, m, L);
        out[1] = invgauss2_hess_lambda_lambda(L);
    }
}

inline void invgauss2_info_dinfo(int k, double y, const double* th,
                                 double* out) {
    double m = th[0], L = th[1];
    if (k == 0) {
        out[0] = invgauss2_expected_mu_mu(m, L);
        out[1] = invgauss2_dexpected_mu_mu_mu(m, L);
    } else {
        out[0] = invgauss2_expected_lambda_lambda(L);
        out[1] = invgauss2_dexpected_lambda_lambda_lambda(L);
    }
}

} // namespace d7

#endif
