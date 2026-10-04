#ifndef D7_PT_SKEWT_H
#define D7_PT_SKEWT_H

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include "pt_loc_scale.h"

// Azzalini's skew t in (mu, sigma, alpha, nu): one function per component
// and order for the score and the diagonal of the Hessian, called by the
// vector kernels (skewt.cpp). With z = (y - mu)/sigma, m = nu + 1,
// s = nu + z^2, c = sqrt(m/s) and w = alpha z c,
//   l = log 2 - log sigma + log t_nu(z) + log T_{nu+1}(w).
// The components in (mu, sigma, alpha) are closed forms in the quantities
// skewt_pieces() returns (skewt_distrib.R's skewt_pieces(), written in C);
// those in nu are one five-point stencil on the log-density, at the step
// max(1e-3 |nu|, 1e-6) and with numericals7's weights at accuracy four, as
// the R methods take them. The caller evaluates the log-density at the
// stencil's nodes (skewt_logpdf) and passes the values.
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
    double s = nu + z * z;
    double cc = std::sqrt(m / s);
    double w = alpha * z * cc;
    double q = std::exp(R::dt(w, m, 1) - skewt_pt()(w, m, 1, 1));
    double e = nu * std::sqrt(m) / std::pow(s, 1.5);
    P.z = z;
    P.w = w;
    P.c = cc;
    P.a = -m * z / s;
    P.da = -m * (nu - z * z) / (s * s);
    P.e = e;
    P.b = alpha * e;
    P.db = -3.0 * alpha * nu * std::sqrt(m) * z / std::pow(s, 2.5);
    P.q = q;
    P.dq = q * (-(m + 1.0) * w / (m + w * w) - q);
    return P;
}

inline double skewt_logpdf(double y, double mu, double sigma, double alpha,
                           double nu) {
    double z = (y - mu) / sigma;
    double w = alpha * z * std::sqrt((nu + 1.0) / (nu + z * z));
    return std::log(2.0) - std::log(sigma) + R::dt(z, nu, 1) +
        skewt_pt()(w, nu + 1.0, 1, 1);
}

// the step in nu, and the stencil weights on the offsets -2, -1, 0, 1, 2,
// exactly as numericals7::fd_weights() returns them
inline double skewt_nu_step(double nu) {
    return std::max(1e-3 * std::fabs(nu), 1e-6);
}

const double kSkewtW1[5] = {0x1.5555555555555p-4, -0x1.5555555555555p-1, 0.0,
                            0x1.5555555555555p-1, -0x1.5555555555555p-4};
const double kSkewtW2[5] = {-0x1.5555555555553p-4, 0x1.5555555555553p+0,
                            -0x1.3ffffffffffffp+1, 0x1.5555555555555p+0,
                            -0x1.5555555555555p-4};

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

// lp[k] the log-density at nu + (k - 2) h, k = 0..4; lp[2] is not read
inline double skewt_score_nu(const double* lp, double h) {
    double acc = 0.0;
    for (int k = 0; k < 5; ++k) {
        if (kSkewtW1[k] == 0.0) continue;
        acc = acc + kSkewtW1[k] * lp[k];
    }
    return acc / h;
}

inline double skewt_hess_mu_mu(const SkewtPieces& P, double sigma) {
    double dd = P.da + P.dq * P.b * P.b + P.q * P.db;
    return dd / (sigma * sigma);
}

inline double skewt_hess_sigma_sigma(const SkewtPieces& P, double sigma) {
    double d = P.a + P.q * P.b;
    double dd = P.da + P.dq * P.b * P.b + P.q * P.db;
    return (1.0 + 2.0 * P.z * d + P.z * P.z * dd) / (sigma * sigma);
}

inline double skewt_hess_alpha_alpha(const SkewtPieces& P) {
    double zc = P.z * P.c;
    return P.dq * zc * zc;
}

inline double skewt_hess_nu_nu(const double* lp, double h) {
    double acc = 0.0;
    for (int k = 0; k < 5; ++k) acc = acc + kSkewtW2[k] * lp[k];
    return acc / (h * h);
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
    double rm = std::sqrt(m);
    double u = z * rm / std::sqrt(s);
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

// lp[k] the log-density at nu + (k - 2) h, k = 0..4; lp[2] is not read
inline double skewt_d3_nu_nu_nu(const double* lp, double h) {
    const double w3[5] = {-0.5, 1.0, 0.0, -1.0, 0.5};
    double acc = 0.0;
    for (int k = 0; k < 5; ++k) {
        if (w3[k] == 0.0) continue;
        acc = acc + w3[k] * lp[k];
    }
    return acc / R_pow(h, 3.0);
}

// the routers of the scalar registry: th = (mu, sigma, alpha, nu)
inline void skewt_score_curv(int k, double y, const double* th, double* out) {
    const double m = th[0], s = th[1], a = th[2], v = th[3];
    if (k < 3) {
        const SkewtPieces P = skewt_pieces(y, m, s, a, v);
        if (k == 0) {
            out[0] = skewt_score_mu(P, s);
            out[1] = skewt_hess_mu_mu(P, s);
        } else if (k == 1) {
            out[0] = skewt_score_sigma(P, s);
            out[1] = skewt_hess_sigma_sigma(P, s);
        } else {
            out[0] = skewt_score_alpha(P);
            out[1] = skewt_hess_alpha_alpha(P);
        }
        return;
    }
    const double h = skewt_nu_step(v);
    double lp[5];
    for (int j = 0; j < 5; ++j) lp[j] = skewt_logpdf(y, m, s, a, v + (j - 2) * h);
    out[0] = skewt_score_nu(lp, h);
    out[1] = skewt_hess_nu_nu(lp, h);
}

// the standardized diagonal pair of parameter k at shape = (alpha, nu), by
// the rule of pt_loc_scale.h
inline void skewt_quad_diag(int k, const double* shape, double* out,
                            bool want_d) {
    const double a = shape[0], v = shape[1];
    const LocScaleRule R = loc_scale_rule();
    const double h = skewt_nu_step(v);
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
            double lp[5];
            lp[2] = skewt_logpdf(z, 0.0, 1.0, a, v);
            const double fw = std::exp(lp[2]) * R.w[j];
            if (fw == 0.0) continue;
            for (int i = 0; i < 5; ++i) {
                if (i != 2) lp[i] = skewt_logpdf(z, 0.0, 1.0, a, v + (i - 2) * h);
            }
            H = skewt_hess_nu_nu(lp, h);
            if (!want_d) { S.add0(fw, H); continue; }
            g = skewt_score_nu(lp, h);
            T3 = skewt_d3_nu_nu_nu(lp, h);
            S.add(fw, g, H, T3);
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
