#ifndef D7_PT_ENET_H
#define D7_PT_ENET_H

#include <Rcpp.h>
#include <cmath>
#include "pt_mills.h"
#include "pt_sqrt.h"


// The elastic-net density in (mu, lambda, alpha),
//   l = -a|z| - c z^2/2 - log Z(a, c),  z = y - mu,
//   a = lambda alpha,  c = lambda (1 - alpha),
//   log Z = log 2 + L(x) - log(c)/2,  x = a/sqrt(c),  L = log M,
// with M the Mills ratio: one function per component and order for the
// quantities the scalar registry reads, called by the vector kernels
// (enet.cpp) and by the registry (d7_ccallable.cpp).
//
// The data enter only through |z| and z^2, so every derivative beyond the
// first in (lambda, alpha) is a derivative of log Z. Its derivatives in
// (a, c) come from G = L'(x), G' = 1 + xG - G^2, G'' = G + xG' - 2GG' and
// those of x: x_a = c^-1/2, x_c = -x/(2c), x_ac = -c^-3/2/2,
// x_cc = 3x/(4c^2), x_acc = 3 c^-5/2/4, x_ccc = -15x/(8c^3), with x_aa,
// x_aaa and x_aac zero; the term -log(c)/2 adds 1/(2c^2) at second order in
// c and -1/c^3 at third. The parameter quantities are formed once per
// parameter value by enet_par() and passed as a struct.

namespace d7 {

// G(x) = d log M/dx = x - m(-x), m the Mills ratio phi/Phi (pt_mills.h);
// past |x| = 1e3 its asymptotic series
inline double enet_G(double x) {
    if (std::fabs(x) > 1e3) {
        double u = 1.0 / x, u2 = u * u, u3 = u2 * u, u5 = u3 * u2,
               u7 = u5 * u2;
        return -u + 2.0 * u3 - 10.0 * u5 + 74.0 * u7;
    }
    return x - mills_ratio(-x);
}

struct EnetPar {
    double lam, al, a, c, x, g, dg;
    double za, zc, zaa, zac, zcc;
};

inline EnetPar enet_par(double lam, double al) {
    EnetPar P;
    P.lam = lam;
    P.al = al;
    P.a = lam * al;
    P.c = lam * (1.0 - al);
    P.x = P.a / d7::sqrt_cr(P.c);
    P.g = enet_G(P.x);
    P.dg = 1.0 + P.x * P.g - P.g * P.g;
    double cc = P.c, x = P.x, g = P.g, dg = P.dg;
    P.za = g / d7::sqrt_cr(cc);
    P.zc = -(1.0 + x * g) / (2.0 * cc);
    P.zaa = dg / cc;
    P.zac = -(x * dg + g) / (2.0 * std::pow(cc, 1.5));
    P.zcc = (2.0 * (1.0 + x * g) + x * (g + x * dg)) / (4.0 * cc * cc);
    return P;
}

// the third derivatives of log Z in (a, c)
struct EnetT3 { double aaa, aac, acc, ccc; };

inline EnetT3 enet_t3(const EnetPar& P) {
    double c = P.c, x = P.x, g = P.g, dg = P.dg;
    double d2g = g + x * dg - 2.0 * g * dg;
    double rc = d7::sqrt_cr(c);
    double xa = 1.0 / rc, xc = -x / (2.0 * c);
    double xac = -0.5 / (c * rc), xcc = 3.0 * x / (4.0 * c * c);
    double xacc = 0.75 / (c * c * rc), xccc = -15.0 * x / (8.0 * c * c * c);
    EnetT3 T;
    T.aaa = d2g * xa * xa * xa;
    T.aac = d2g * xa * xa * xc + 2.0 * dg * xa * xac;
    T.acc = d2g * xa * xc * xc + dg * (2.0 * xac * xc + xcc * xa) + g * xacc;
    T.ccc = d2g * xc * xc * xc + 3.0 * dg * xcc * xc + g * xccc - 1.0 / (c * c * c);
    return T;
}

inline double enet_sign(double z) {
    return (z > 0.0) - (z < 0.0);
}

inline double enet_score_mu(const EnetPar& P, double z) {
    return P.a * enet_sign(z) + P.c * z;
}

inline double enet_score_lambda(const EnetPar& P, double z) {
    double d_a = -std::fabs(z) - P.za;
    double d_c = -z * z / 2.0 - P.zc;
    return d_a * P.al + d_c * (1.0 - P.al);
}

inline double enet_score_alpha(const EnetPar& P, double z) {
    double d_a = -std::fabs(z) - P.za;
    double d_c = -z * z / 2.0 - P.zc;
    return P.lam * (d_a - d_c);
}

inline double enet_hess_mu_mu(const EnetPar& P) {
    return -P.c;
}

inline double enet_hess_lambda_lambda(const EnetPar& P) {
    double al = P.al, bl = 1.0 - al;
    return -(P.zaa * al * al + 2.0 * P.zac * al * bl + P.zcc * bl * bl);
}

inline double enet_hess_alpha_alpha(const EnetPar& P) {
    return -P.lam * P.lam * (P.zaa - 2.0 * P.zac + P.zcc);
}

inline double enet_expected_mu_mu(const EnetPar& P) {
    return -(P.a * P.a - 2.0 * P.a * P.c * P.za - 2.0 * P.c * P.c * P.zc);
}

// the second derivatives in (lambda, alpha) are free of y
inline double enet_expected_lambda_lambda(const EnetPar& P) {
    return enet_hess_lambda_lambda(P);
}

inline double enet_expected_alpha_alpha(const EnetPar& P) {
    return enet_hess_alpha_alpha(P);
}

// E_mm depends on (lambda, alpha) alone
inline double enet_dexpected_mu_mu_mu() {
    return 0.0;
}

// with w = (alpha, 1 - alpha), E_ll = -w'Z''w and its derivative in lambda
// is -Z'''[w, w, w]
inline double enet_dexpected_lambda_lambda_lambda(const EnetPar& P,
                                                  const EnetT3& T) {
    double al = P.al, bl = 1.0 - al;
    return -(T.aaa * al * al * al + 3.0 * T.aac * al * al * bl +
             3.0 * T.acc * al * bl * bl + T.ccc * bl * bl * bl);
}

// with v = (1, -1), E_aa = -lambda^2 v'Z''v and da/dalpha = lambda,
// dc/dalpha = -lambda
inline double enet_dexpected_alpha_alpha_alpha(const EnetPar& P,
                                               const EnetT3& T) {
    double l3 = P.lam * P.lam * P.lam;
    return -l3 * (T.aaa - 3.0 * T.aac + 3.0 * T.acc - T.ccc);
}

inline void enet_score_curv(int k, double y, const double* th, double* out) {
    const EnetPar P = enet_par(th[1], th[2]);
    double z = y - th[0];
    if (k == 0) {
        out[0] = enet_score_mu(P, z);
        out[1] = enet_hess_mu_mu(P);
    } else if (k == 1) {
        out[0] = enet_score_lambda(P, z);
        out[1] = enet_hess_lambda_lambda(P);
    } else {
        out[0] = enet_score_alpha(P, z);
        out[1] = enet_hess_alpha_alpha(P);
    }
}

inline void enet_info_dinfo(int k, double y, const double* th, double* out) {
    const EnetPar P = enet_par(th[1], th[2]);
    if (k == 0) {
        out[0] = enet_expected_mu_mu(P);
        out[1] = enet_dexpected_mu_mu_mu();
    } else if (k == 1) {
        out[0] = enet_expected_lambda_lambda(P);
        out[1] = enet_dexpected_lambda_lambda_lambda(P, enet_t3(P));
    } else {
        out[0] = enet_expected_alpha_alpha(P);
        out[1] = enet_dexpected_alpha_alpha_alpha(P, enet_t3(P));
    }
}

} // namespace d7

#endif
