#ifndef D7_PT_ENET_H
#define D7_PT_ENET_H

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include "pt_parts.h"
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

// ---- the parts route (piano_motore.md, D7) ----------------------------------
//
// Orders one and two and the expected Hessian call the components of the
// joint routers above. Orders three and four are transcribed from
// enet_deriv_k() and its helpers (enet.cpp), here under enet_parts_ names:
// the table of one parameter value, formed at z = 0 and s = 0, plus the one
// component that carries the data at order three, sgn(z) - z in
// (mu, lambda, alpha). The costly quantity is the table of the partial
// derivatives of log Z in (a, c), formed once by enet_logz_ac() to the order
// the level asks for (its entries of low order do not depend on that order).

// The derivatives of l = -a|z| - c z^2/2 - log Z(a, c) in (mu, a, c) along
// the index tuple ix (0 = mu, 1 = a, 2 = c) of length m, with s = sgn(z)
inline double enet_parts_l_ac(int m, const int* ix, double z, double s,
                              double a, double c, const double d[5][5]) {
    int nm = 0, na = 0, nc = 0;
    for (int r = 0; r < m; ++r) {
        if (ix[r] == 0) ++nm;
        else if (ix[r] == 1) ++na;
        else ++nc;
    }
    if (m == 1) {
        if (nm) return a * s + c * z;
        if (na) return -std::fabs(z) - d[1][0];
        return -z * z / 2.0 - d[0][1];
    }
    if (nm == 0) return -d[na][nc];
    if (m == 2) {
        if (nm == 2) return -c;
        return na ? s : z;
    }
    if (m == 3 && nm == 2 && nc == 1) return -1.0;
    return 0.0;
}

inline double enet_parts_l_form(int m, const double (*dir)[3], double z,
                                double s, double a, double c,
                                const double d[5][5]) {
    int ix[4] = {0, 0, 0, 0};
    double out = 0.0;
    int total = 1;
    for (int r = 0; r < m; ++r) total *= 3;
    for (int t = 0; t < total; ++t) {
        int u = t;
        double w = 1.0;
        for (int r = 0; r < m; ++r) {
            ix[r] = u % 3;
            u /= 3;
            w *= dir[r][ix[r]];
        }
        if (w != 0.0) out += w * enet_parts_l_ac(m, ix, z, s, a, c, d);
    }
    return out;
}

// One derivative of order k in (mu, lambda, alpha), with nm, nl and na the
// counts of each (the partition sum over the bilinear map to (a, c))
inline double enet_parts_deriv_lam_al(int nm, int nl, int na, double z,
                                      double s, double lam, double al,
                                      const double d[5][5]) {
    const double a = lam * al, c = lam * (1.0 - al);
    const double emu[3] = {1.0, 0.0, 0.0}, w[3] = {0.0, al, 1.0 - al},
                 lv[3] = {0.0, lam, -lam}, v[3] = {0.0, 1.0, -1.0};
    static const double fact[5] = {1.0, 1.0, 2.0, 6.0, 24.0};
    double out = 0.0;
    for (int j = 0; j <= std::min(nl, na); ++j) {
        const double coef = fact[nl] / (fact[j] * fact[nl - j]) *
                            fact[na] / (fact[j] * fact[na - j]) * fact[j];
        double dir[4][3];
        int m = 0;
        auto add = [&](const double* e, int times) {
            for (int t = 0; t < times; ++t, ++m)
                for (int q = 0; q < 3; ++q) dir[m][q] = e[q];
        };
        add(emu, nm);
        add(w, nl - j);
        add(lv, na - j);
        add(v, j);
        out += coef * enet_parts_l_form(m, dir, z, s, a, c, d);
    }
    return out;
}

// Buffer layout, observed and expected alike (index -> quantity -> lowest
// level; the expected buffer has NaN at 1 and reads level 2 only):
//   0  level
//   1  y            2  mu            3  lambda         4  alpha
//   5  a = lambda alpha                6  c = lambda (1 - alpha)
//   8  Z_a   9  Z_c                                           1
//  10  Z_aa  11  Z_ac  12  Z_cc                               2
//  13-16  Z_aaa, Z_aac, Z_acc, Z_ccc                          3
//  17-21  Z_aaaa, Z_aaac, Z_aacc, Z_accc, Z_cccc              4
// (Z the logarithm of the normalizing constant; slot 7 is unused.)
inline void enet_parts_fill(double y, const double* th, int level, double* b) {
    for (int i = 1; i < 22; ++i) b[i] = R_NaN;
    b[0] = (double) level;
    b[1] = y;
    b[2] = th[0];
    b[3] = th[1];
    b[4] = th[2];
    b[5] = th[1] * th[2];
    b[6] = th[1] * (1.0 - th[2]);
    if (level >= 1) {
        const int n = level > 4 ? 4 : level;
        double d[5][5];
        for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) d[i][j] = 0.0;
        enet_logz_ac(b[5], b[6], n, d);
        b[8] = d[1][0];
        b[9] = d[0][1];
        if (n >= 2) {
            b[10] = d[2][0];
            b[11] = d[1][1];
            b[12] = d[0][2];
        }
        if (n >= 3) {
            b[13] = d[3][0];
            b[14] = d[2][1];
            b[15] = d[1][2];
            b[16] = d[0][3];
        }
        if (n >= 4) {
            b[17] = d[4][0];
            b[18] = d[3][1];
            b[19] = d[2][2];
            b[20] = d[1][3];
            b[21] = d[0][4];
        }
    }
}

inline void enet_parts(double y, const double* th, int level, double* b) {
    enet_parts_fill(y, th, level, b);
}

inline void enet_parts_expected(const double* th, int level, double* b) {
    enet_parts_fill(R_NaN, th, level, b);
}

// the EnetPar the components read, from the buffer
inline EnetPar enet_parts_par(const double* b) {
    EnetPar P;
    P.lam = b[3];
    P.al = b[4];
    P.a = b[5];
    P.c = b[6];
    P.x = R_NaN;
    P.za = b[8];
    P.zc = b[9];
    P.zaa = b[10];
    P.zac = b[11];
    P.zcc = b[12];
    P.t3.aaa = b[13];
    P.t3.aac = b[14];
    P.t3.acc = b[15];
    P.t3.ccc = b[16];
    return P;
}

inline void enet_parts_grad(const double* b, double* out) {
    const EnetPar P = enet_parts_par(b);
    double z = b[1] - b[2];
    out[0] = enet_score_mu(P, z);
    out[1] = enet_score_lambda(P, z);
    out[2] = enet_score_alpha(P, z);
}

inline void enet_parts_hess(const double* b, double* out) {
    const EnetPar P = enet_parts_par(b);
    double z = b[1] - b[2];
    out[0] = enet_hess_mu_mu(P);
    out[4] = enet_hess_lambda_lambda(P);
    out[8] = enet_hess_alpha_alpha(P);
    out[1] = out[3] = enet_hess_mu_lambda(P, z);
    out[2] = out[6] = enet_hess_mu_alpha(P, z);
    out[5] = out[7] = enet_hess_lambda_alpha(P, z);
}

inline void enet_parts_info(const double* b, double* out) {
    const EnetPar P = enet_parts_par(b);
    out[0] = enet_expected_mu_mu(P);
    out[4] = enet_expected_lambda_lambda(P);
    out[8] = enet_expected_alpha_alpha(P);
    out[1] = out[3] = 0.0;
    out[2] = out[6] = 0.0;
    out[5] = out[7] = enet_expected_lambda_alpha(P);
}

// the table of partial derivatives of log Z from the buffer
inline void enet_parts_table(const double* b, double d[5][5]) {
    for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) d[i][j] = 0.0;
    d[1][0] = b[8];
    d[0][1] = b[9];
    d[2][0] = b[10];
    d[1][1] = b[11];
    d[0][2] = b[12];
    d[3][0] = b[13];
    d[2][1] = b[14];
    d[1][2] = b[15];
    d[0][3] = b[16];
    d[4][0] = b[17];
    d[3][1] = b[18];
    d[2][2] = b[19];
    d[1][3] = b[20];
    d[0][4] = b[21];
}

// the 10 unique components, tuples (i0 <= i1 <= i2) over (mu, lambda, alpha)
inline void enet_parts_d3(const double* b, double* out) {
    double d[5][5];
    enet_parts_table(b, d);
    double u[10];
    int w = 0;
    for (int i0 = 0; i0 < 3; ++i0)
        for (int i1 = i0; i1 < 3; ++i1)
            for (int i2 = i1; i2 < 3; ++i2) {
                int c3[3] = {0, 0, 0};
                ++c3[i0]; ++c3[i1]; ++c3[i2];
                double v = enet_parts_deriv_lam_al(c3[0], c3[1], c3[2], 0.0,
                                                   0.0, b[3], b[4], d);
                if (c3[0] == 1 && c3[1] == 1 && c3[2] == 1) {
                    const double z = b[1] - b[2];
                    v += enet_sign(z) - z;
                }
                u[w++] = v;
            }
    sym_fill(3, 3, u, out);
}

// the 15 unique components, tuples (i0 <= i1 <= i2 <= i3)
inline void enet_parts_d4(const double* b, double* out) {
    double d[5][5];
    enet_parts_table(b, d);
    double u[15];
    int w = 0;
    for (int i0 = 0; i0 < 3; ++i0)
        for (int i1 = i0; i1 < 3; ++i1)
            for (int i2 = i1; i2 < 3; ++i2)
                for (int i3 = i2; i3 < 3; ++i3) {
                    int c3[3] = {0, 0, 0};
                    ++c3[i0]; ++c3[i1]; ++c3[i2]; ++c3[i3];
                    u[w++] = enet_parts_deriv_lam_al(c3[0], c3[1], c3[2], 0.0,
                                                     0.0, b[3], b[4], d);
                }
    sym_fill(3, 4, u, out);
}

} // namespace d7

#endif
