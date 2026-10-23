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

// ---- the components only the joint registry reads ---------------------------
// Written from the kernels of enet.cpp (enet_hessian_cpp,
// enet_expected_hessian_cpp, enet_dexpected1_cpp), which keep them inline.

inline double enet_hess_mu_lambda(const EnetPar& P, double z) {
    double s = enet_sign(z), al = P.al;
    return al * s + (1.0 - al) * z;
}

inline double enet_hess_mu_alpha(const EnetPar& P, double z) {
    double s = enet_sign(z), lam = P.lam;
    return lam * (s - z);
}

inline double enet_hess_lambda_alpha(const EnetPar& P, double z) {
    double al = P.al, lam = P.lam;
    return (-std::fabs(z) + z * z / 2.0) -
        (lam * (P.zaa * al + P.zac * (1.0 - 2.0 * al) - P.zcc * (1.0 - al)) +
         P.za - P.zc);
}

inline double enet_expected_lambda_alpha(const EnetPar& P) {
    double al = P.al;
    return -P.lam * (P.zaa * al + P.zac * (1.0 - 2.0 * al) -
                     P.zcc * (1.0 - al));
}

// E_mm = -(a^2 - 2ac Z_a - 2c^2 Z_c) differentiated in a and in c, then
// carried to lambda and alpha
inline double enet_dexpected_mu_mu_lambda(const EnetPar& P) {
    double al = P.al, bl = 1.0 - al, a = P.a, c = P.c;
    double dEa = -(2.0 * a - 2.0 * c * P.za - 2.0 * a * c * P.zaa -
                   2.0 * c * c * P.zac);
    double dEc = -(-2.0 * a * P.za - 2.0 * a * c * P.zac - 4.0 * c * P.zc -
                   2.0 * c * c * P.zcc);
    return al * dEa + bl * dEc;
}

inline double enet_dexpected_mu_mu_alpha(const EnetPar& P) {
    double lam = P.lam, a = P.a, c = P.c;
    double dEa = -(2.0 * a - 2.0 * c * P.za - 2.0 * a * c * P.zaa -
                   2.0 * c * c * P.zac);
    double dEc = -(-2.0 * a * P.za - 2.0 * a * c * P.zac - 4.0 * c * P.zc -
                   2.0 * c * c * P.zcc);
    return lam * (dEa - dEc);
}

// Z'' along (w, v), and Z''' along (w, w, v)
inline double enet_dexpected_lambda_lambda_alpha(const EnetPar& P,
                                                 const EnetT3& T) {
    double lam = P.lam, al = P.al, bl = 1.0 - al;
    double Hwv = P.zaa * al + P.zac * (bl - al) - P.zcc * bl;
    double Twwv = T.aaa * al * al + T.aac * (2.0 * al * bl - al * al) +
                  T.acc * (bl * bl - 2.0 * al * bl) - T.ccc * bl * bl;
    return -(2.0 * Hwv + lam * Twwv);
}

// Z'' along (v, v), and Z''' along (v, v, w)
inline double enet_dexpected_alpha_alpha_lambda(const EnetPar& P,
                                                const EnetT3& T) {
    double lam = P.lam, al = P.al, bl = 1.0 - al;
    double Hvv = P.zaa - 2.0 * P.zac + P.zcc;
    double Tvvw = T.aaa * al + T.aac * (bl - 2.0 * al) +
                  T.acc * (al - 2.0 * bl) + T.ccc * bl;
    return -2.0 * lam * Hvv - lam * lam * Tvvw;
}

inline double enet_dexpected_lambda_alpha_lambda(const EnetPar& P,
                                                 const EnetT3& T) {
    double lam = P.lam, al = P.al, bl = 1.0 - al;
    double Hwv = P.zaa * al + P.zac * (bl - al) - P.zcc * bl;
    double Twwv = T.aaa * al * al + T.aac * (2.0 * al * bl - al * al) +
                  T.acc * (bl * bl - 2.0 * al * bl) - T.ccc * bl * bl;
    return -Hwv - lam * Twwv;
}

inline double enet_dexpected_lambda_alpha_alpha(const EnetPar& P,
                                                const EnetT3& T) {
    double lam = P.lam, al = P.al, bl = 1.0 - al;
    double Hvv = P.zaa - 2.0 * P.zac + P.zcc;
    double Tvvw = T.aaa * al + T.aac * (bl - 2.0 * al) +
                  T.acc * (al - 2.0 * bl) + T.ccc * bl;
    return -lam * (Hvv + lam * Tvvw);
}

// ---- the second derivatives of the expected information ---------------------
//
// With w = (alpha, 1 - alpha), v = (1, -1) and Z^(k)[u_1, ..., u_k] the k-th
// derivative of log Z along the directions u_i in (a, c), d/dlambda moves a
// Z^(k) to Z^(k+1)[., w] and d/dalpha to lambda Z^(k+1)[., v] while taking
// w to v. Differentiating E_ll = -Z''[w, w], E_aa = -lambda^2 Z''[v, v] and
// E_la = -lambda Z''[w, v] once more gives the nine derivatives in
// (lambda, alpha) written in enet_d2expected_block(). E_mm is differentiated
// in the form E_mm = -c - 2 a u, with u = 1/Z the density at zero and
// u_a = -Z_a u, u_c = -Z_c u, u_ab = (Z_a Z_b - Z_ab) u: it needs the
// derivatives of log Z to order two only and is a sum of terms of one sign,
// where -(a^2 - 2ac Z_a - 2c^2 Z_c) loses its second derivative in lambda to
// cancellation when alpha is small. E_ml = E_ma = 0 and E does not depend on
// mu, so every other derivative is zero.

// the density at zero, u = 1/Z = (a + sqrt(c) m_1)/2 with m_1 = M_1/M_0
inline double enet_density0(double a, double c) {
    const double sc = std::sqrt(c);
    double m[12];
    enet_s_moments(a / sc, 1, m);
    return 0.5 * (a + sc * m[1]);
}

// d[i][j] contracted along the directions u[0..order-1] (each a pair, the
// components in a and in c)
inline double enet_zk(int order, const double d[5][5], const double* const* u) {
    double tot = 0.0;
    for (int idx = 0; idx < (1 << order); ++idx) {
        double w = 1.0;
        int nc = 0;
        for (int r = 0; r < order; ++r) {
            const int bit = (idx >> r) & 1;
            w *= u[r][bit];
            nc += bit;
        }
        if (w != 0.0) tot += w * d[order - nc][nc];
    }
    return tot;
}

// The derivatives of E_mm, E_ll, E_aa and E_la (rows 0 to 3) in
// (lambda, lambda), (alpha, alpha) and (lambda, alpha) (columns 0 to 2).
inline void enet_d2expected_block(double lam, double al, double out[4][3]) {
    const double bl = 1.0 - al, a = lam * al, c = lam * bl;
    double d[5][5];
    enet_logz_ac(a, c, 4, d);
    const double w[2] = {al, bl}, v[2] = {1.0, -1.0};
    const double* const vv[2] = {v, v};
    const double* const wwv[3] = {w, w, v};
    const double* const wvv[3] = {w, v, v};
    const double* const vvv[3] = {v, v, v};
    const double* const wwww[4] = {w, w, w, w};
    const double* const wwwv[4] = {w, w, w, v};
    const double* const wwvv[4] = {w, w, v, v};
    const double* const wvvv[4] = {w, v, v, v};
    const double* const vvvv[4] = {v, v, v, v};
    const double Z_vv = enet_zk(2, d, vv);
    const double Z_wwv = enet_zk(3, d, wwv), Z_wvv = enet_zk(3, d, wvv),
                 Z_vvv = enet_zk(3, d, vvv);
    const double Z_wwww = enet_zk(4, d, wwww), Z_wwwv = enet_zk(4, d, wwwv),
                 Z_wwvv = enet_zk(4, d, wwvv), Z_wvvv = enet_zk(4, d, wvvv),
                 Z_vvvv = enet_zk(4, d, vvvv);
    // E_mm = f(a, c) = -c - 2 a u
    const double Za = d[1][0], Zc = d[0][1], Zaa = d[2][0], Zac = d[1][1],
                 Zcc = d[0][2];
    const double u = enet_density0(a, c);
    const double ua = -Za * u, uc = -Zc * u, uaa = (Za * Za - Zaa) * u,
                 uac = (Za * Zc - Zac) * u, ucc = (Zc * Zc - Zcc) * u;
    const double fa = -2.0 * u - 2.0 * a * ua, fc = -1.0 - 2.0 * a * uc;
    const double faa = -4.0 * ua - 2.0 * a * uaa,
                 fac = -2.0 * uc - 2.0 * a * uac, fcc = -2.0 * a * ucc;
    out[0][0] = al * al * faa + 2.0 * al * bl * fac + bl * bl * fcc;
    out[0][1] = lam * lam * (faa - 2.0 * fac + fcc);
    out[0][2] = lam * (al * faa + (bl - al) * fac - bl * fcc) + (fa - fc);
    // E_ll = -Z''[w, w]
    out[1][0] = -Z_wwww;
    out[1][1] = -(lam * lam * Z_wwvv + 4.0 * lam * Z_wvv + 2.0 * Z_vv);
    out[1][2] = -(lam * Z_wwwv + 3.0 * Z_wwv);
    // E_aa = -lambda^2 Z''[v, v]
    out[2][0] = -(2.0 * Z_vv + 4.0 * lam * Z_wvv + lam * lam * Z_wwvv);
    out[2][1] = -lam * lam * lam * lam * Z_vvvv;
    out[2][2] = -(3.0 * lam * lam * Z_vvv + lam * lam * lam * Z_wvvv);
    // E_la = -lambda Z''[w, v]
    out[3][0] = -(2.0 * Z_wwv + lam * Z_wwwv);
    out[3][1] = -(2.0 * lam * lam * Z_vvv + lam * lam * lam * Z_wvvv);
    out[3][2] = -(Z_vv + 3.0 * lam * Z_wvv + lam * lam * Z_wwvv);
}

// ---- the routers of the joint registry --------------------------------------
//
// One per order: every component at one observation, on the parameter scale,
// in (mu, lambda, alpha). The Hessian and the expected Hessian fill a full
// 3 x 3 matrix and the derivative of the expected Hessian a 3 x 3 x 3 array,
// column-major, out[a + 3 b + 9 c] being the derivative of E[l_ab] in
// parameter c. E_ml = E_ma = 0 and every derivative in mu is 0.

inline void enet_joint_grad(double y, const double* th, double* out) {
    const EnetPar P = enet_par(th[1], th[2]);
    double z = y - th[0];
    out[0] = enet_score_mu(P, z);
    out[1] = enet_score_lambda(P, z);
    out[2] = enet_score_alpha(P, z);
}

inline void enet_joint_hess(double y, const double* th, double* out) {
    const EnetPar P = enet_par(th[1], th[2]);
    double z = y - th[0];
    out[0] = enet_hess_mu_mu(P);
    out[4] = enet_hess_lambda_lambda(P);
    out[8] = enet_hess_alpha_alpha(P);
    out[1] = out[3] = enet_hess_mu_lambda(P, z);
    out[2] = out[6] = enet_hess_mu_alpha(P, z);
    out[5] = out[7] = enet_hess_lambda_alpha(P, z);
}

inline void enet_joint_info(double y, const double* th, double* out) {
    (void) y;
    const EnetPar P = enet_par(th[1], th[2]);
    out[0] = enet_expected_mu_mu(P);
    out[4] = enet_expected_lambda_lambda(P);
    out[8] = enet_expected_alpha_alpha(P);
    out[1] = out[3] = 0.0;
    out[2] = out[6] = 0.0;
    out[5] = out[7] = enet_expected_lambda_alpha(P);
}

inline void enet_joint_dinfo(double y, const double* th, double* out) {
    (void) y;
    const EnetPar P = enet_par(th[1], th[2]);
    const EnetT3 T = enet_t3(P);
    for (int j = 0; j < 27; ++j) out[j] = 0.0;
    // c = 1 (lambda)
    out[0 + 9] = enet_dexpected_mu_mu_lambda(P);
    out[4 + 9] = enet_dexpected_lambda_lambda_lambda(P, T);
    out[8 + 9] = enet_dexpected_alpha_alpha_lambda(P, T);
    out[5 + 9] = out[7 + 9] = enet_dexpected_lambda_alpha_lambda(P, T);
    // c = 2 (alpha)
    out[0 + 18] = enet_dexpected_mu_mu_alpha(P);
    out[4 + 18] = enet_dexpected_lambda_lambda_alpha(P, T);
    out[8 + 18] = enet_dexpected_alpha_alpha_alpha(P, T);
    out[5 + 18] = out[7 + 18] = enet_dexpected_lambda_alpha_alpha(P, T);
}

} // namespace d7

#endif
