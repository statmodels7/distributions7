#ifndef D7_PT_INVGAUSS2_H
#define D7_PT_INVGAUSS2_H

#include <Rcpp.h>
#include "pt_parts.h"

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

// Off-diagonal entry of the observed Hessian, as invgauss2_hessian_cpp
// computes it.
inline double invgauss2_hess_mu_lambda(double y, double m) {
    double m3 = m * m * m;
    return (y - m) / m3;
}

// d E[l_mu mu] / d lambda, as invgauss2_dexpected1_cpp computes it.
inline double invgauss2_dexpected_mu_mu_lambda(double m) {
    double im = 1.0 / m, im3 = im * im * im;
    return -im3;
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

// The routers of the joint registry, one per order: every component at one
// observation, on the parameter scale (mu, lambda). The Hessian and the
// expected information fill a full P x P matrix and the derivative of the
// expected information a P x P x P array, column-major, out[a + P b + P^2 c]
// being the derivative of E[l_ab] in parameter c.
inline void invgauss2_joint_grad(double y, const double* th, double* out) {
    double m = th[0], L = th[1];
    out[0] = invgauss2_score_mu(y, m, L);
    out[1] = invgauss2_score_lambda(y, m, L);
}

inline void invgauss2_joint_hess(double y, const double* th, double* out) {
    double m = th[0], L = th[1];
    out[0] = invgauss2_hess_mu_mu(y, m, L);
    out[3] = invgauss2_hess_lambda_lambda(L);
    out[1] = out[2] = invgauss2_hess_mu_lambda(y, m);
}

inline void invgauss2_joint_info(double y, const double* th, double* out) {
    double m = th[0], L = th[1];
    out[0] = invgauss2_expected_mu_mu(m, L);
    out[3] = invgauss2_expected_lambda_lambda(L);
    out[1] = out[2] = 0.0;
}

inline void invgauss2_joint_dinfo(double y, const double* th, double* out) {
    double m = th[0], L = th[1];
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[0] = invgauss2_dexpected_mu_mu_mu(m, L);
    out[4] = invgauss2_dexpected_mu_mu_lambda(m);
    out[7] = invgauss2_dexpected_lambda_lambda_lambda(L);
}

// The parts route (pt_parts.h). Nothing in this family is costly, so the buffer
// holds the inputs every order reads, and each order is written from them:
//   [0] level                       always
//   [1] y                           observed parts; NaN in the expected parts
//   [2] mu                          always
//   [3] lambda                      always
// Orders one and two call the joint routers' component functions. Orders
// three and four are the expressions of invgauss2_deriv3_cpp and
// invgauss2_deriv4_cpp (invgauss2.cpp) with expected = FALSE, so yy = y.
inline void invgauss2_parts(double y, const double* th, int level, double* b) {
    b[0] = level; b[1] = y; b[2] = th[0]; b[3] = th[1];
}

inline void invgauss2_parts_expected(const double* th, int level, double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = th[1];
}

inline void invgauss2_parts_grad(const double* b, double* out) {
    out[0] = invgauss2_score_mu(b[1], b[2], b[3]);
    out[1] = invgauss2_score_lambda(b[1], b[2], b[3]);
}

inline void invgauss2_parts_hess(const double* b, double* out) {
    out[0] = invgauss2_hess_mu_mu(b[1], b[2], b[3]);
    out[3] = invgauss2_hess_lambda_lambda(b[3]);
    out[1] = out[2] = invgauss2_hess_mu_lambda(b[1], b[2]);
}

// from the expected parts
inline void invgauss2_parts_info(const double* b, double* out) {
    out[0] = invgauss2_expected_mu_mu(b[2], b[3]);
    out[3] = invgauss2_expected_lambda_lambda(b[3]);
    out[1] = out[2] = 0.0;
}

// third order: invgauss2_deriv3_cpp, unique components in the sorted order
// mu mu mu, mu mu lambda, mu lambda lambda, lambda lambda lambda
inline void invgauss2_parts_d3(const double* b, double* out) {
    const double yy = b[1], m = b[2], L = b[3];
    const double m2 = m * m, m4 = m2 * m2, m5 = m4 * m;
    double u[4];
    u[0] = L * (12.0 * yy - 6.0 * m) / m5;
    u[1] = (2.0 * m - 3.0 * yy) / m4;
    u[2] = 0.0;
    u[3] = 1.0 / (L * L * L);
    sym_fill(2, 3, u, out);
}

// fourth order: invgauss2_deriv4_cpp, unique components in the sorted order
// mu mu mu mu, mu mu mu lambda, mu mu lambda lambda, mu lambda lambda lambda,
// lambda lambda lambda lambda
inline void invgauss2_parts_d4(const double* b, double* out) {
    const double yy = b[1], m = b[2], L = b[3];
    const double m2 = m * m, m4 = m2 * m2, m5 = m4 * m, m6 = m5 * m;
    const double L2 = L * L;
    double u[5];
    u[0] = L * (24.0 * m - 60.0 * yy) / m6;
    u[1] = (12.0 * yy - 6.0 * m) / m5;
    u[2] = 0.0;
    u[3] = 0.0;
    u[4] = -3.0 / (L2 * L2);
    sym_fill(2, 4, u, out);
}

} // namespace d7

#endif
