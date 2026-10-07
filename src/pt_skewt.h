#ifndef D7_PT_SKEWT_H
#define D7_PT_SKEWT_H

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include "pt_loc_scale.h"
#include "pt_sqrt.h"
#include "pt_skewt_exact.h"

// Azzalini's skew t in (mu, sigma, alpha, nu): one function per component
// and order for the score and the diagonal of the Hessian, called by the
// vector kernels (skewt.cpp). With z = (y - mu)/sigma, m = nu + 1,
// s = nu + z^2, c = sqrt(m/s) and w = alpha z c,
//   l = log 2 - log sigma + log t_nu(z) + log T_{nu+1}(w).
// The components in (mu, sigma, alpha) are closed forms in the quantities
// skewt_pieces() returns (skewt_distrib.R's skewt_pieces(), written in C);
// those in nu are the exact derivatives of pt_skewt_exact.h, which the
// vector kernels (skewt.cpp) read at every order.
//
// The t distribution function is numericals7's n7_pt(), R's pt() compiled
// without its warnings, resolved once through R_GetCCallable on the calling
// thread (skewt_pt()), so that the kernels and the registry may run these
// functions on worker threads. R::dt() is silent at positive degrees of
// freedom.

namespace d7 {

typedef double (*N7Pt)(double, double, int, int);

// n7_pt, resolved on first use; the first call must be on the calling thread
N7Pt skewt_pt();

struct SkewtPieces { double z, w, c, a, da, e, b, db, q, dq; };

inline SkewtPieces skewt_pieces(double y, double mu, double sigma,
                                double alpha, double nu) {
    SkewtPieces P;
    double z = (y - mu) / sigma;
    double m = nu + 1.0;
    if (std::fabs(z) > 1e100) {
        // s = nu + z^2 overflows near |z| = 1e154, so far out the pieces are
        // written in rs = sqrt(s), a = z/rs and nu/s, which stay finite;
        // skewt_pieces() in R takes the same branch
        const double rs = std::fabs(z) * d7::sqrt_cr(1.0 + nu / (z * z));
        const double a = z / rs, irs = 1.0 / rs, b2 = nu * irs * irs;
        const double rm = d7::sqrt_cr(m);
        P.c = rm * irs;
        P.w = alpha * rm * a;
        P.a = -m * a * irs;
        P.da = -m * (b2 - a * a) * irs * irs;
        P.e = rm * b2 * irs;
        P.db = -3.0 * alpha * rm * b2 * a * irs * irs;
    } else {
        double s = nu + z * z;
        double cc = d7::sqrt_cr(m / s);
        P.c = cc;
        P.w = alpha * z * cc;
        P.a = -m * z / s;
        P.da = -m * (nu - z * z) / (s * s);
        P.e = nu * d7::sqrt_cr(m) / std::pow(s, 1.5);
        P.db = -3.0 * alpha * nu * d7::sqrt_cr(m) * z / std::pow(s, 2.5);
    }
    const double w = P.w;
    double q = std::exp(R::dt(w, m, 1) - skewt_pt()(w, m, 1, 1));
    P.z = z;
    P.b = alpha * P.e;
    P.q = q;
    P.dq = q * (-(m + 1.0) * w / (m + w * w) - q);
    return P;
}

inline double skewt_logpdf(double y, double mu, double sigma, double alpha,
                           double nu) {
    double z = (y - mu) / sigma;
    double w = alpha * z * d7::sqrt_cr((nu + 1.0) / (nu + z * z));
    if (std::fabs(z) > 1e100) {
        // as skewt_pieces(), and as distrib_pdf() in R
        const double rs = std::fabs(z) * d7::sqrt_cr(1.0 + nu / (z * z));
        w = alpha * d7::sqrt_cr(nu + 1.0) * (z / rs);
    }
    return std::log(2.0) - std::log(sigma) + R::dt(z, nu, 1) +
        skewt_pt()(w, nu + 1.0, 1, 1);
}

inline double skewt_score_mu(const SkewtPieces& P, double sigma) {
    double d = P.a + P.q * P.b;
    return -d / sigma;
}

inline double skewt_score_sigma(const SkewtPieces& P, double sigma) {
    double d = P.a + P.q * P.b;
    return -(1.0 + P.z * d) / sigma;
}

inline double skewt_score_alpha(const SkewtPieces& P) {
    return P.q * P.z * P.c;
}

inline double skewt_hess_mu_mu(const SkewtPieces& P, double sigma) {
    double dd = P.da + P.dq * P.b * P.b + P.q * P.db;
    return dd / (sigma * sigma);
}

inline double skewt_hess_sigma_sigma(const SkewtPieces& P, double sigma) {
    double d = P.a + P.q * P.b;
    double dd = P.da + P.dq * P.b * P.b + P.q * P.db;
    // z^2 overflows past |z| of 1e154 where z (z dd) does not
    const double zzdd = (std::fabs(P.z) > 1e100) ? P.z * (P.z * dd) :
        P.z * P.z * dd;
    return (1.0 + 2.0 * P.z * d + zzdd) / (sigma * sigma);
}

inline double skewt_hess_alpha_alpha(const SkewtPieces& P) {
    double zc = P.z * P.c;
    return P.dq * zc * zc;
}

// The third derivatives on the diagonal. Those in (mu, sigma, alpha) are
// skewt_msa_tower() and skewt_msa_component() (skewt_distrib.R) written in
// C: with u = z sqrt((nu+1)/(nu+z^2)), w = alpha u, Q = t_{nu+1}(w)/T_{nu+1}(w)
// and its w-derivatives q1, q2 by the Riccati recursion, L1 to L3 are the
// z-derivatives of log t_nu(z) + log T_{nu+1}(alpha u(z)). The powers are
// R_pow(), which is what R's ^ evaluates (powl() on 64-bit Windows).
struct SkewtD3 { double z, u, q2, L1, L2, L3; };

inline SkewtD3 skewt_d3_pieces(double y, double mu, double sigma,
                               double alpha, double nu) {
    SkewtD3 D;
    double z = (y - mu) / sigma;
    double m = nu + 1.0;
    double s = nu + z * z;
    double rm = d7::sqrt_cr(m);
    double u = z * rm / d7::sqrt_cr(s);
    double u1 = nu * rm * R_pow(s, -1.5);
    double u2 = -3.0 * nu * rm * z * R_pow(s, -2.5);
    double u3 = -3.0 * nu * rm * (nu - 4.0 * z * z) * R_pow(s, -3.5);
    double gz1 = -m * z / s;
    double gz2 = -m * (nu - z * z) / (s * s);
    double gz3 = 2.0 * m * z * (3.0 * nu - z * z) / R_pow(s, 3.0);
    double w = alpha * u;
    double q = std::exp(R::dt(w, m, 1) - skewt_pt()(w, m, 1, 1));
    double mw = m + w * w;
    double g0 = -(m + 1.0) * w / mw;
    double g1 = -(m + 1.0) * (m - w * w) / (mw * mw);
    double q1 = q * (g0 - q);
    double q2 = q1 * (g0 - q) + q * (g1 - q1);
    D.z = z;
    D.u = u;
    D.q2 = q2;
    D.L1 = gz1 + q * alpha * u1;
    D.L2 = gz2 + (q * alpha * u2 + q1 * (alpha * alpha) * (u1 * u1));
    D.L3 = gz3 + (q * alpha * u3 + q1 * (alpha * alpha) * (3.0 * u1 * u2) +
                  q2 * R_pow(alpha, 3.0) * R_pow(u1, 3.0));
    return D;
}

inline double skewt_d3_mu_mu_mu(const SkewtD3& D, double sigma) {
    return -D.L3 / R_pow(sigma, 3.0);
}

inline double skewt_d3_sigma_sigma_sigma(const SkewtD3& D, double sigma) {
    double z = D.z;
    double acc = 6.0 * z * D.L1 + 6.0 * (z * z) * D.L2 +
        R_pow(z, 3.0) * D.L3;
    double s3 = R_pow(sigma, 3.0);
    return -acc / s3 + -2.0 / s3;
}

inline double skewt_d3_alpha_alpha_alpha(const SkewtD3& D) {
    return R_pow(D.u, 3.0) * D.q2;
}

// the routers of the scalar registry: th = (mu, sigma, alpha, nu)
inline void skewt_score_curv(int k, double y, const double* th, double* out) {
    const double m = th[0], s = th[1], a = th[2], v = th[3];
    if (k < 3) {
        // the generated closed forms the vector kernels read, without the
        // integrals in nu, which these components do not read: the same
        // numbers to the bit
        static const double C0[7] = {0, 0, 0, 0, 0, 0, 0};
        const SkewtObs o = skewt_obs_vars(y, m, s, a, v);
        double B[6][6] = {{0}};
        skewt_btab<2>(o, C0, B, false);
        double g[4], h2[10];
        skewt_d_1(o.a, o.irs, o.isg, a, v, o.rm, o.L, C0, B, g);
        skewt_d_2(o.a, o.irs, o.isg, a, v, o.rm, o.L, C0, B, h2);
        static const int diag[3] = {0, 4, 7};
        out[0] = g[k];
        out[1] = h2[diag[k]];
        return;
    }
    // the derivatives in nu, exact (pt_skewt_exact.h), from the same
    // generated expressions as the vector kernels, so that the two agree to
    // the bit; the constants of one nu are kept per thread
    thread_local SkewtConst cc;
    cc.at(v, 2);
    const SkewtObs o = skewt_obs_vars(y, m, s, a, v);
    double B[6][6] = {{0}};
    skewt_btab<2>(o, cc.Cm, B);
    double g[4], h2[10];
    skewt_d_1(o.a, o.irs, o.isg, a, v, o.rm, o.L, cc.Cnu, B, g);
    skewt_d_2(o.a, o.irs, o.isg, a, v, o.rm, o.L, cc.Cnu, B, h2);
    out[0] = g[3];
    out[1] = h2[9];
}

// the standardized diagonal pair of parameter k at shape = (alpha, nu), by
// the rule of pt_loc_scale.h
inline void skewt_quad_diag(int k, const double* shape, double* out,
                            bool want_d) {
    const double a = shape[0], v = shape[1];
    const LocScaleRule R = loc_scale_rule();
    SkewtConst cc;
    if (k == 3) cc.at(v, want_d ? 3 : 2);
    LocScaleSum S;
    for (int j = 0; j < R.n; ++j) {
        const double z = R.x[j];
        double g = 0.0, H, T3 = 0.0;
        if (k < 3) {
            const double fw = std::exp(skewt_logpdf(z, 0.0, 1.0, a, v)) * R.w[j];
            if (fw == 0.0) continue;
            const SkewtPieces P = skewt_pieces(z, 0.0, 1.0, a, v);
            if (k == 0) H = skewt_hess_mu_mu(P, 1.0);
            else if (k == 1) H = skewt_hess_sigma_sigma(P, 1.0);
            else H = skewt_hess_alpha_alpha(P);
            if (!want_d) { S.add0(fw, H); continue; }
            const SkewtD3 D = skewt_d3_pieces(z, 0.0, 1.0, a, v);
            if (k == 0) {
                g = skewt_score_mu(P, 1.0);
                T3 = skewt_d3_mu_mu_mu(D, 1.0);
            } else if (k == 1) {
                g = skewt_score_sigma(P, 1.0);
                T3 = skewt_d3_sigma_sigma_sigma(D, 1.0);
            } else {
                g = skewt_score_alpha(P);
                T3 = skewt_d3_alpha_alpha_alpha(D);
            }
            S.add(fw, g, H, T3);
        } else {
            const double fw = std::exp(skewt_logpdf(z, 0.0, 1.0, a, v)) * R.w[j];
            if (fw == 0.0) continue;
            // the expressions the vector kernels read, to the bit
            const SkewtObs o = skewt_obs_vars(z, 0.0, 1.0, a, v);
            double B[6][6] = {{0}};
            double d1[4], d2[10], d3[20], dn[3];
            if (!want_d) {
                skewt_btab<2>(o, cc.Cm, B);
                skewt_d_2(o.a, o.irs, o.isg, a, v, o.rm, o.L, cc.Cnu, B, d2);
                S.add0(fw, d2[9]);
                continue;
            }
            skewt_btab<3>(o, cc.Cm, B);
            skewt_d_1(o.a, o.irs, o.isg, a, v, o.rm, o.L, cc.Cnu, B, d1);
            skewt_d_2(o.a, o.irs, o.isg, a, v, o.rm, o.L, cc.Cnu, B, d2);
            skewt_d_3(o.a, o.irs, o.isg, a, v, o.rm, o.L, cc.Cnu, B, d3);
            dn[0] = d1[3];
            dn[1] = d2[9];
            dn[2] = d3[19];
            S.add(fw, dn[0], dn[1], dn[2]);
        }
    }
    S.result(out);
}

inline void skewt_info_dinfo(int k, double y, const double* th, double* out) {
    (void) y;
    double std2[2];
    loc_scale_cached(1, k, th + 2, 2, [&](double* v) {
        skewt_quad_diag(k, th + 2, v, true);
    }, std2);
    loc_scale_unstandardize(k, th[1], std2, out);
}

} // namespace d7

#endif
