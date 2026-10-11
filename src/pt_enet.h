#ifndef D7_PT_ENET_H
#define D7_PT_ENET_H

#include <Rcpp.h>
#include <cmath>
#include "pt_mills.h"
#include "pt_sqrt.h"
#include "pt_enet_logz.h"


// The elastic-net density in (mu, lambda, alpha),
//   l = -a|z| - c z^2/2 - log Z(a, c),  z = y - mu,
//   a = lambda alpha,  c = lambda (1 - alpha),
//   log Z = log 2 + L(x) - log(c)/2,  x = a/sqrt(c),  L = log M,
// with M the Mills ratio: one function per component and order for the
// quantities the scalar registry reads, called by the vector kernels
// (enet.cpp) and by the registry (d7_ccallable.cpp).
//
// The data enter only through |z| and z^2, so every derivative beyond the
// first in (lambda, alpha) is a derivative of log Z. Its partial derivatives
// in (a, c) come from enet_logz_ac() (pt_enet_logz.h), which forms them
// without the cancellation that the expressions in G = L'(x) carry for large
// x. The parameter quantities are formed once per parameter value by
// enet_par() and passed as a struct.

namespace d7 {

// the third derivatives of log Z in (a, c)
struct EnetT3 { double aaa, aac, acc, ccc; };

struct EnetPar {
    double lam, al, a, c, x;
    double za, zc, zaa, zac, zcc;
    EnetT3 t3;
};

inline EnetPar enet_par(double lam, double al) {
    EnetPar P;
    P.lam = lam;
    P.al = al;
    P.a = lam * al;
    P.c = lam * (1.0 - al);
    P.x = P.a / d7::sqrt_cr(P.c);
    double d[5][5];
    enet_logz_ac(P.a, P.c, 3, d);
    P.za = d[1][0];
    P.zc = d[0][1];
    P.zaa = d[2][0];
    P.zac = d[1][1];
    P.zcc = d[0][2];
    P.t3.aaa = d[3][0];
    P.t3.aac = d[2][1];
    P.t3.acc = d[1][2];
    P.t3.ccc = d[0][3];
    return P;
}

inline EnetT3 enet_t3(const EnetPar& P) {
    return P.t3;
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
