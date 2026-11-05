#ifndef D7_PT_CAUCHY_H
#define D7_PT_CAUCHY_H

// Cauchy in location and scale: one function per component and order for
// the quantities the scalar registry reads, called by the vector kernels
// (cauchy.cpp, dexpected_kernels.cpp) and by the registry
// (d7_ccallable.cpp). res = y - mu.

#include <Rcpp.h>
#include <cmath>
#include <complex>
#include "pt_parts.h"

namespace d7 {

// The Cauchy through Y = res + i s, res = y - mu: log f = log s - 2 Re log Y
// - log pi and F = 1 - Im(log Y)/pi, with
//   d^n log Y / dmu^j ds^k = (-1)^j i^k (-1)^(n - 1) (n - 1)! Y^(-n),
// so every derivative of order n is one term in Y^(-n) (plus (-1)^(n-1)
// (n-1)!/s^n from log s when it is pure in s). Neither res^2 nor s^2 is
// formed, and no two terms cancel at any residual.

// 1/Y by Smith's division
inline std::complex<double> cauchy_yinv(double res, double s) {
    if (std::fabs(res) >= std::fabs(s)) {
        const double r = s / res, d = res + s * r;
        return std::complex<double>(1.0 / d, -r / d);
    }
    const double r = res / s, d = res * r + s;
    return std::complex<double>(r / d, -1.0 / d);
}

// i^k X: its real part (re = true) or imaginary part
inline double cauchy_ik_part(int k, std::complex<double> X, bool re) {
    switch (k % 4) {
    case 0: return re ? X.real() : X.imag();
    case 1: return re ? -X.imag() : X.real();
    case 2: return re ? -X.real() : -X.imag();
    default: return re ? X.imag() : -X.real();
    }
}

// the order-n derivative of log f with j derivatives in mu, n - j in s
inline double cauchy_logdens_d(int n, int j, double res, double s) {
    const std::complex<double> w = cauchy_yinv(res, s);
    std::complex<double> X(1.0, 0.0);
    double f = 1.0;
    for (int m = 0; m < n; m++) X *= w;
    for (int m = 2; m < n; m++) f *= m;
    const int k = n - j;
    const double sg = ((j + n - 1) % 2) ? -1.0 : 1.0;
    double v = -2.0 * sg * f * cauchy_ik_part(k, X, true);
    if (j == 0) {
        double is = 1.0;
        for (int m = 0; m < n; m++) is /= s;
        v += ((n - 1) % 2 ? -f : f) * is;
    }
    return v;
}

// the order-n derivative of F with j derivatives in mu, n - j in s
inline double cauchy_cdf_d(int n, int j, double res, double s) {
    const std::complex<double> w = cauchy_yinv(res, s);
    std::complex<double> X(1.0, 0.0);
    double f = 1.0;
    for (int m = 0; m < n; m++) X *= w;
    for (int m = 2; m < n; m++) f *= m;
    const double sg = ((j + n - 1) % 2) ? -1.0 : 1.0;
    return -sg * f * cauchy_ik_part(n - j, X, false) / 3.141592653589793238462643383280;
}

inline double cauchy_score_mu(double res, double s) {
    return cauchy_logdens_d(1, 1, res, s);
}

inline double cauchy_score_sigma(double res, double s) {
    return cauchy_logdens_d(1, 0, res, s);
}

inline double cauchy_hess_mu_mu(double res, double s) {
    return cauchy_logdens_d(2, 2, res, s);
}

inline double cauchy_hess_mu_sigma(double res, double s) {
    return cauchy_logdens_d(2, 1, res, s);
}

inline double cauchy_hess_sigma_sigma(double res, double s) {
    return cauchy_logdens_d(2, 0, res, s);
}

inline double cauchy_expected_mu_mu(double s) {
    double s2 = s * s;
    return -0.5 / s2;
}

inline double cauchy_expected_sigma_sigma(double s) {
    double s2 = s * s;
    return -0.5 / s2;
}

// E_mm and E_ss are k/sigma^2 with k = -1/2; nothing moves with mu
inline double cauchy_dexpected_mu_mu_mu() {
    return 0.0;
}

inline double cauchy_dexpected_sigma_sigma_sigma(double s) {
    const double k = -0.5;
    double u = 1.0 / s, u3 = u * u * u;
    return -2.0 * k * u3;
}

inline void cauchy_score_curv(int k, double y, const double* th,
                              double* out) {
    double res = y - th[0], s = th[1];
    if (k == 0) {
        out[0] = cauchy_score_mu(res, s);
        out[1] = cauchy_hess_mu_mu(res, s);
    } else {
        out[0] = cauchy_score_sigma(res, s);
        out[1] = cauchy_hess_sigma_sigma(res, s);
    }
}

inline void cauchy_info_dinfo(int k, double y, const double* th,
                              double* out) {
    double s = th[1];
    if (k == 0) {
        out[0] = cauchy_expected_mu_mu(s);
        out[1] = cauchy_dexpected_mu_mu_mu();
    } else {
        out[0] = cauchy_expected_sigma_sigma(s);
        out[1] = cauchy_dexpected_sigma_sigma_sigma(s);
    }
}

// The component the joint registry reads and the scalar one does not: the
// derivative of E_mm = -1/(2 sigma^2) in sigma, written as the kernel writes
// it, -2 k u^3 with k = -1/2 and u = 1/sigma (dexpected_kernels.cpp).
inline double cauchy_dexpected_mu_mu_sigma(double s) {
    const double k = -0.5;
    double u = 1.0 / s, u3 = u * u * u;
    return -2.0 * k * u3;
}

// The routers of the joint registry, one per order: every component at one
// observation, on the parameter scale, in (mu, sigma). The Hessian and the
// expected information fill a full 2 x 2 matrix and the derivative of the
// expected information a 2 x 2 x 2 array, column-major, out[a + 2 b + 4 c]
// being the derivative of E[l_ab] in parameter c.
inline void cauchy_joint_grad(double y, const double* th, double* out) {
    double res = y - th[0], s = th[1];
    out[0] = cauchy_score_mu(res, s);
    out[1] = cauchy_score_sigma(res, s);
}

inline void cauchy_joint_hess(double y, const double* th, double* out) {
    double res = y - th[0], s = th[1];
    out[0] = cauchy_hess_mu_mu(res, s);
    out[3] = cauchy_hess_sigma_sigma(res, s);
    out[1] = out[2] = cauchy_hess_mu_sigma(res, s);
}

inline void cauchy_joint_info(double y, const double* th, double* out) {
    double s = th[1];
    out[0] = cauchy_expected_mu_mu(s);
    out[3] = cauchy_expected_sigma_sigma(s);
    out[1] = out[2] = 0.0;
}

// E_ms = 0 for every sigma, so both of its derivatives are zero
inline void cauchy_joint_dinfo(double y, const double* th, double* out) {
    double s = th[1];
    for (int j = 0; j < 8; ++j) out[j] = 0.0;
    out[4] = cauchy_dexpected_mu_mu_sigma(s);
    out[7] = cauchy_dexpected_sigma_sigma_sigma(s);
}

// The parts route (pt_parts.h). The Cauchy has no special function. Its one
// costly quantity is the complex reciprocal 1/Y of Y = res + i s, computed by
// Smith's division, which every order reads, so the buffer holds it once:
//   [0] level                       always
//   [1] y                           observed parts; NaN in the expected parts
//   [2] mu                          observed parts; NaN in the expected parts
//   [3] sigma                       always
//   [4] res = y - mu                observed parts; NaN in the expected parts
//   [5] Re(1/Y)                     level >= 1; NaN in the expected parts
//   [6] Im(1/Y)                     level >= 1; NaN in the expected parts
// Each order is the function of the joint router with the reciprocal read
// from the buffer. cauchy_yinv() receives the same arguments that
// cauchy_logdens_d() gives it, so the values are bit-identical.

// cauchy_logdens_d() with the reciprocal w = 1/Y given: the same operations in
// the same order.
inline double cauchy_parts_logdens_d(int n, int j, std::complex<double> w,
                                     double s) {
    std::complex<double> X(1.0, 0.0);
    double f = 1.0;
    for (int m = 0; m < n; m++) X *= w;
    for (int m = 2; m < n; m++) f *= m;
    const int k = n - j;
    const double sg = ((j + n - 1) % 2) ? -1.0 : 1.0;
    double v = -2.0 * sg * f * cauchy_ik_part(k, X, true);
    if (j == 0) {
        double is = 1.0;
        for (int m = 0; m < n; m++) is /= s;
        v += ((n - 1) % 2 ? -f : f) * is;
    }
    return v;
}

inline void cauchy_parts(double y, const double* th, int level, double* b) {
    const double m = th[0], s = th[1];
    const double res = y - m;
    b[0] = level; b[1] = y; b[2] = m; b[3] = s; b[4] = res;
    if (level >= 1) {
        const std::complex<double> w = cauchy_yinv(res, s);
        b[5] = w.real(); b[6] = w.imag();
    } else {
        b[5] = R_NaN; b[6] = R_NaN;
    }
}

inline void cauchy_parts_expected(const double* th, int level, double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = R_NaN; b[3] = th[1]; b[4] = R_NaN;
    b[5] = R_NaN; b[6] = R_NaN;
}

inline std::complex<double> cauchy_parts_w(const double* b) {
    return std::complex<double>(b[5], b[6]);
}

// from the observed parts: score in mu and in sigma (the joint grad)
inline void cauchy_parts_grad(const double* b, double* out) {
    const std::complex<double> w = cauchy_parts_w(b);
    const double s = b[3];
    out[0] = cauchy_parts_logdens_d(1, 1, w, s);
    out[1] = cauchy_parts_logdens_d(1, 0, w, s);
}

inline void cauchy_parts_hess(const double* b, double* out) {
    const std::complex<double> w = cauchy_parts_w(b);
    const double s = b[3];
    out[0] = cauchy_parts_logdens_d(2, 2, w, s);
    out[3] = cauchy_parts_logdens_d(2, 0, w, s);
    out[1] = out[2] = cauchy_parts_logdens_d(2, 1, w, s);
}

// from the expected parts
inline void cauchy_parts_info(const double* b, double* out) {
    out[0] = cauchy_expected_mu_mu(b[3]);
    out[3] = cauchy_expected_sigma_sigma(b[3]);
    out[1] = out[2] = 0.0;
}

// third order: cauchy_list(3) in cauchy_hd.cpp, component with m sigma
// indices is order-3 derivative with j = 3 - m mu indices; unique components
// in the sorted order mu mu mu, mu mu sigma, mu sigma sigma, sigma sigma sigma
inline void cauchy_parts_d3(const double* b, double* out) {
    const std::complex<double> w = cauchy_parts_w(b);
    const double s = b[3];
    double u[4];
    for (int k = 0; k <= 3; k++) u[k] = cauchy_parts_logdens_d(3, 3 - k, w, s);
    sym_fill(2, 3, u, out);
}

// fourth order: cauchy_list(4) in cauchy_hd.cpp, same convention
inline void cauchy_parts_d4(const double* b, double* out) {
    const std::complex<double> w = cauchy_parts_w(b);
    const double s = b[3];
    double u[5];
    for (int k = 0; k <= 4; k++) u[k] = cauchy_parts_logdens_d(4, 4 - k, w, s);
    sym_fill(2, 4, u, out);
}

} // namespace d7

#endif
