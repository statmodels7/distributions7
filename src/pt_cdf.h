#ifndef D7_PT_CDF_H
#define D7_PT_CDF_H

#include <Rcpp.h>
#include <cmath>
#include "pt_gpd.h"

// The distribution functions of the continuous families and their first
// and second derivatives in the parameters, for truncated(): one function
// per family and order, called by the scalar registry and, through
// d7_cdf_*_cpp(), by trunc_route_parts() in R, so that the retained mass
// and its derivatives are the same numbers on both routes. The derivatives
// are those of F, on the parameter scale; cdf_grad writes one value per
// parameter and cdf_hess one per pair in the order of hess_names(): the
// diagonal first, then (i, j), i < j, row by row.
//
// A family returns false when it has no compiled distribution function.

namespace d7 {

// ---- Gaussian1: (mu, sigma) ------------------------------------------------

inline double gaussian1_cdf(double q, const double* th, bool lower) {
    return R::pnorm(q, th[0], th[1], lower, 0);
}

// f and d log f / dy at q, as dnorm() and distrib_grad_y() write them
inline void gaussian1_f_ly(double q, const double* th, double* f, double* ly) {
    *f = R::dnorm(q, th[0], th[1], 0);
    *ly = -(q - th[0]) / (th[1] * th[1]);
}

// the location-scale derivatives of F from f and l_y at z = (q - mu)/s:
// F_mu = -f, F_s = -z f, F_mumu = f l_y, F_ss = f (z^2 l_y + 2 z / s),
// F_mus = f (z l_y + 1/s)
inline void loc_scale_cdf_grad(double z, double f, double* g) {
    g[0] = -f;
    g[1] = -z * f;
}
inline void loc_scale_cdf_hess(double z, double s, double f, double ly,
                               double* h) {
    h[0] = f * ly;
    h[1] = f * (z * z * ly + 2 * z / s);
    h[2] = f * (z * ly + 1 / s);
}

inline void gaussian1_cdf_grad(double q, const double* th, double* g) {
    double f, ly;
    gaussian1_f_ly(q, th, &f, &ly);
    loc_scale_cdf_grad((q - th[0]) / th[1], f, g);
}

inline void gaussian1_cdf_hess(double q, const double* th, double* h) {
    double f, ly;
    gaussian1_f_ly(q, th, &f, &ly);
    loc_scale_cdf_hess((q - th[0]) / th[1], th[1], f, ly, h);
}

// ---- location-scale families in (mu, sigma) --------------------------------

// Cauchy: f = 1/(pi s (1 + z^2)), l_y = -2 z / (s (1 + z^2))
inline double cauchy_cdf(double q, const double* th, bool lower) {
    return R::pcauchy(q, th[0], th[1], lower, 0);
}
inline void cauchy_f_ly(double q, const double* th, double* z, double* f,
                        double* ly) {
    *z = (q - th[0]) / th[1];
    const double u = 1 + *z * *z;
    *f = 1 / (M_PI * th[1] * u);
    *ly = -2 * *z / (th[1] * u);
}

// logistic: f = dlogis, l_y = -tanh(z/2)/s
inline double logistic_cdf(double q, const double* th, bool lower) {
    return R::plogis(q, th[0], th[1], lower, 0);
}
inline void logistic_f_ly(double q, const double* th, double* z, double* f,
                          double* ly) {
    *z = (q - th[0]) / th[1];
    *f = R::dlogis(q, th[0], th[1], 0);
    *ly = -std::tanh(*z / 2) / th[1];
}

// Gumbel (maxima): F = exp(-w), w = exp(-z), f = w F / s, l_y = (w - 1)/s
inline double gumbel_cdf(double q, const double* th, bool lower) {
    const double w = std::exp(-(q - th[0]) / th[1]);
    return lower ? std::exp(-w) : -std::expm1(-w);
}
inline void gumbel_f_ly(double q, const double* th, double* z, double* f,
                        double* ly) {
    *z = (q - th[0]) / th[1];
    const double w = std::exp(-*z);
    *f = w * std::exp(-w) / th[1];
    *ly = (w - 1) / th[1];
}

// Laplace: f = exp(-|z|)/(2 s), l_y = -sign(z)/s; both tails directly
inline double laplace_cdf(double q, const double* th, bool lower) {
    const double z = (q - th[0]) / th[1];
    const double t = 0.5 * std::exp(-std::fabs(z));
    if (lower) return (z < 0) ? t : 1 - t;
    return (z < 0) ? 1 - t : t;
}
inline void laplace_f_ly(double q, const double* th, double* z, double* f,
                         double* ly) {
    *z = (q - th[0]) / th[1];
    *f = std::exp(-std::fabs(*z)) / (2 * th[1]);
    *ly = -((*z > 0) - (*z < 0)) / th[1];
}

#define D7_LOC_SCALE_CDF(fam)                                              \
inline void fam##_cdf_grad(double q, const double* th, double* g) {        \
    double z, f, ly;                                                       \
    fam##_f_ly(q, th, &z, &f, &ly);                                        \
    loc_scale_cdf_grad(z, f, g);                                           \
}                                                                          \
inline void fam##_cdf_hess(double q, const double* th, double* h) {        \
    double z, f, ly;                                                       \
    fam##_f_ly(q, th, &z, &f, &ly);                                        \
    loc_scale_cdf_hess(z, th[1], f, ly, h);                                \
}

D7_LOC_SCALE_CDF(cauchy)
D7_LOC_SCALE_CDF(logistic)
D7_LOC_SCALE_CDF(gumbel)
D7_LOC_SCALE_CDF(laplace)

// ---- Gaussian2 (mu, sigma2) and Gaussian3 (mu, tau): Gaussian1 at
// sigma = sqrt(sigma2) or 1/sqrt(tau), by the chain rule in the scale, with
// d1 and d2 the first and second derivatives of sigma in the parameter

inline void gaussian_scale_chain(double q, double mu, double sg, double d1,
                                 double d2, double* g, double* h) {
    const double th1[2] = {mu, sg};
    double g1[2], h1[3];
    gaussian1_cdf_grad(q, th1, g1);
    if (g) { g[0] = g1[0]; g[1] = g1[1] * d1; }
    if (h) {
        gaussian1_cdf_hess(q, th1, h1);
        h[0] = h1[0];
        h[1] = h1[1] * d1 * d1 + g1[1] * d2;
        h[2] = h1[2] * d1;
    }
}

inline double gaussian2_cdf(double q, const double* th, bool lower) {
    return R::pnorm(q, th[0], std::sqrt(th[1]), lower, 0);
}
inline void gaussian2_cdf_grad(double q, const double* th, double* g) {
    const double sg = std::sqrt(th[1]);
    gaussian_scale_chain(q, th[0], sg, 0.5 / sg, 0, g, nullptr);
}
inline void gaussian2_cdf_hess(double q, const double* th, double* h) {
    const double sg = std::sqrt(th[1]);
    gaussian_scale_chain(q, th[0], sg, 0.5 / sg, -0.25 / (sg * th[1]),
                         nullptr, h);
}

inline double gaussian3_cdf(double q, const double* th, bool lower) {
    return R::pnorm(q, th[0], 1 / std::sqrt(th[1]), lower, 0);
}
inline void gaussian3_cdf_grad(double q, const double* th, double* g) {
    const double r = std::sqrt(th[1]);
    gaussian_scale_chain(q, th[0], 1 / r, -0.5 / (th[1] * r), 0, g, nullptr);
}
inline void gaussian3_cdf_hess(double q, const double* th, double* h) {
    const double r = std::sqrt(th[1]);
    gaussian_scale_chain(q, th[0], 1 / r, -0.5 / (th[1] * r),
                         0.75 / (th[1] * th[1] * r), nullptr, h);
}

// ---- Laplace2 (mu, lambda): r = q - mu, f = lambda/2 exp(-lambda |r|);
// F_mu = -f, F_lambda = r f / lambda, F_mumu = -lambda sign(r) f,
// F_lambdalambda = -sign(r) r^2 f / lambda, F_mulambda = -(1/lambda - |r|) f

inline double laplace2_cdf(double q, const double* th, bool lower) {
    const double r = q - th[0];
    const double t = 0.5 * std::exp(-th[1] * std::fabs(r));
    if (lower) return (r < 0) ? t : 1 - t;
    return (r < 0) ? 1 - t : t;
}
inline void laplace2_cdf_grad(double q, const double* th, double* g) {
    const double r = q - th[0], lam = th[1];
    const double f = 0.5 * lam * std::exp(-lam * std::fabs(r));
    g[0] = -f;
    g[1] = r * f / lam;
}
inline void laplace2_cdf_hess(double q, const double* th, double* h) {
    const double r = q - th[0], lam = th[1];
    const double f = 0.5 * lam * std::exp(-lam * std::fabs(r));
    const double sg = (r > 0) - (r < 0);
    h[0] = -lam * sg * f;
    h[1] = -sg * r * r * f / lam;
    h[2] = -(1 / lam - std::fabs(r)) * f;
}

// ---- Lognormal1 (mu, sigma2 = v): z = (log q - mu)/sigma, phi = dnorm(z);
// F_mu = -phi/sigma, F_v = -z phi/(2 v), F_mumu = -z phi/v,
// F_vv = z phi (3 - z^2)/(4 v^2), F_muv = phi (1 - z^2)/(2 sigma^3)

inline double lognormal1_cdf(double q, const double* th, bool lower) {
    return R::plnorm(q, th[0], std::sqrt(th[1]), lower, 0);
}
inline void lognormal1_cdf_grad(double q, const double* th, double* g) {
    if (!(q > 0)) { g[0] = g[1] = 0.0; return; }
    const double v = th[1], sg = std::sqrt(v);
    const double z = (std::log(q) - th[0]) / sg, ph = R::dnorm(z, 0.0, 1.0, 0);
    g[0] = -ph / sg;
    g[1] = -z * ph / (2 * v);
}
inline void lognormal1_cdf_hess(double q, const double* th, double* h) {
    if (!(q > 0)) { h[0] = h[1] = h[2] = 0.0; return; }
    const double v = th[1], sg = std::sqrt(v);
    const double z = (std::log(q) - th[0]) / sg, ph = R::dnorm(z, 0.0, 1.0, 0);
    h[0] = -z * ph / v;
    h[1] = z * ph * (3 - z * z) / (4 * v * v);
    h[2] = ph * (1 - z * z) / (2 * v * sg);
}

// ---- survival families: S = exp(L), F_i = -S L_i, F_ij = -S (L_ij + L_i L_j)

// exponential (mu): L = -q/mu
inline double exponential_cdf(double q, const double* th, bool lower) {
    if (!(q > 0)) return lower ? 0.0 : 1.0;
    const double L = -q / th[0];
    return lower ? -std::expm1(L) : std::exp(L);
}
inline void exponential_cdf_grad(double q, const double* th, double* g) {
    if (!(q > 0)) { g[0] = 0.0; return; }
    const double m = th[0], S = std::exp(-q / m);
    g[0] = -S * q / (m * m);
}
inline void exponential_cdf_hess(double q, const double* th, double* h) {
    if (!(q > 0)) { h[0] = 0.0; return; }
    const double m = th[0], S = std::exp(-q / m);
    const double Li = q / (m * m), Lii = -2 * q / (m * m * m);
    h[0] = -S * (Lii + Li * Li);
}

// Weibull1 (mu scale, sigma shape): L = -e^h, h = sigma log(q/mu)
inline double weibull1_cdf(double q, const double* th, bool lower) {
    return R::pweibull(q, th[1], th[0], lower, 0);
}
inline void weibull1_L(double q, const double* th, double* L1, double* L2) {
    const double mu = th[0], sg = th[1], lq = std::log(q / mu);
    const double eh = std::exp(sg * lq);
    const double hm = -sg / mu, hs = lq;
    L1[0] = -eh * hm;
    L1[1] = -eh * hs;
    if (L2) {
        L2[0] = -eh * (sg / (mu * mu) + hm * hm);
        L2[1] = -eh * (hs * hs);
        L2[2] = -eh * (-1 / mu + hm * hs);
    }
}
inline void weibull1_cdf_grad(double q, const double* th, double* g) {
    if (!(q > 0)) { g[0] = g[1] = 0.0; return; }
    double L1[2];
    weibull1_L(q, th, L1, nullptr);
    const double S = R::pweibull(q, th[1], th[0], 0, 0);
    g[0] = -S * L1[0];
    g[1] = -S * L1[1];
}
inline void weibull1_cdf_hess(double q, const double* th, double* h) {
    if (!(q > 0)) { h[0] = h[1] = h[2] = 0.0; return; }
    double L1[2], L2[3];
    weibull1_L(q, th, L1, L2);
    const double S = R::pweibull(q, th[1], th[0], 0, 0);
    h[0] = -S * (L2[0] + L1[0] * L1[0]);
    h[1] = -S * (L2[1] + L1[1] * L1[1]);
    h[2] = -S * (L2[2] + L1[0] * L1[1]);
}

// ---- GPD (sigma, xi): S = exp(-W), W = z phi_0(xi z), z = q/sigma, with
// pt_gpd.h's phi_j (a series near xi z = 0); F = 1 - S, so the derivatives
// of F are minus those of S, written as in gpd_*_surv_cpp(). Below the
// support F = 0 and past its end F = 1, both constant.

inline bool gpd_inside(double q, const double* th, double* t) {
    *t = gpd_t(q, th[0], th[1]);
    return q > 0.0 && *t > 0.0;
}
inline double gpd_cdf(double q, const double* th, bool lower) {
    double t;
    if (!gpd_inside(q, th, &t)) {
        const bool above = q > 0.0;
        return (above == lower) ? 1.0 : 0.0;
    }
    const double z = q / th[0], u = th[1] * z;
    const double LT = gpd_logt(u, t);
    const double W = z * gpd_phi(0, u, u / t, t, LT);
    return lower ? -std::expm1(-W) : std::exp(-W);
}
inline void gpd_cdf_grad(double q, const double* th, double* g) {
    double t;
    if (!gpd_inside(q, th, &t)) { g[0] = g[1] = 0.0; return; }
    const double sv = th[0], z = q / sv, u = th[1] * z;
    const double T1 = 1.0 / t, LT = gpd_logt(u, t), v = u / t;
    const double F1 = gpd_phi(1, u, v, t, LT);
    const double Sv = std::exp(-z * gpd_phi(0, u, v, t, LT));
    g[0] = -(Sv * T1 * z) / sv;
    g[1] = F1 * Sv * z * z;
}
inline void gpd_cdf_hess(double q, const double* th, double* h) {
    double t;
    if (!gpd_inside(q, th, &t)) { h[0] = h[1] = h[2] = 0.0; return; }
    const double sv = th[0], x = th[1], z = q / sv, u = x * z;
    const double T1 = 1.0 / t, LT = gpd_logt(u, t), v = u / t;
    const double F1 = gpd_phi(1, u, v, t, LT), F2 = gpd_phi(2, u, v, t, LT);
    const double Sv = std::exp(-z * gpd_phi(0, u, v, t, LT));
    const double w0 = T1 * z;
    h[0] = -(Sv * w0 * (w0 * x + w0 - 2)) / (sv * sv);
    h[1] = -(Sv * z * z * z * (F1 * F1 * z - F2));
    h[2] = (Sv * T1 * z * z * (F1 * z + T1)) / sv;
}

// ---- Weibull3 (mean m, shape s): the scale is b = m / Gamma(1 + 1/s), and
// S = exp(-e^h), h = s log(q/b). With x = 1 + 1/s,
// h_m = -s/m, h_s = log(q/b) - psi(x)/s, h_mm = s/m^2, h_ms = -1/m,
// h_ss = psi'(x)/s^3; L = -e^h, L_i = -e^h h_i, L_ij = -e^h (h_ij + h_i h_j)

inline double weibull3_scale(const double* th) {
    return th[0] / std::exp(R::lgammafn(1.0 + 1.0 / th[1]));
}
inline double weibull3_cdf(double q, const double* th, bool lower) {
    return R::pweibull(q, th[1], weibull3_scale(th), lower, 0);
}
inline void weibull3_L(double q, const double* th, double* L1, double* L2,
                       double* S) {
    const double m = th[0], s = th[1], x = 1.0 + 1.0 / s;
    const double b = m / std::exp(R::lgammafn(x));
    const double lq = std::log(q / b), eh = std::exp(s * lq);
    const double p0 = R::digamma(x);
    const double hm = -s / m, hs = lq - p0 / s;
    *S = std::exp(-eh);
    L1[0] = -eh * hm;
    L1[1] = -eh * hs;
    if (L2) {
        const double p1 = R::trigamma(x);
        L2[0] = -eh * (s / (m * m) + hm * hm);
        L2[1] = -eh * (p1 / (s * s * s) + hs * hs);
        L2[2] = -eh * (-1 / m + hm * hs);
    }
}
inline void weibull3_cdf_grad(double q, const double* th, double* g) {
    if (!(q > 0)) { g[0] = g[1] = 0.0; return; }
    double L1[2], S;
    weibull3_L(q, th, L1, nullptr, &S);
    g[0] = -S * L1[0];
    g[1] = -S * L1[1];
}
inline void weibull3_cdf_hess(double q, const double* th, double* h) {
    if (!(q > 0)) { h[0] = h[1] = h[2] = 0.0; return; }
    double L1[2], L2[3], S;
    weibull3_L(q, th, L1, L2, &S);
    h[0] = -S * (L2[0] + L1[0] * L1[0]);
    h[1] = -S * (L2[1] + L1[1] * L1[1]);
    h[2] = -S * (L2[2] + L1[0] * L1[1]);
}

// ---- Lognormal2 (mean m, var v): Lognormal1 at meanlog mu = log m - S/2
// and sdlog^2 = S = log1p(v/m^2), by the chain rule in (mu, S):
// S_m = -2v/(m(m^2+v)), S_v = 1/(m^2+v), S_mm = 2v(3m^2+v)/(m(m^2+v))^2,
// S_vv = -1/(m^2+v)^2, S_mv = -2m/(m^2+v)^2; mu_m = 1/m - S_m/2,
// mu_v = -S_v/2, mu_mm = -1/m^2 - S_mm/2, mu_mv = -S_mv/2, mu_vv = -S_vv/2

inline void lognormal2_map(const double* th, double* l1, double* u1,
                           double* u2) {
    const double m = th[0], v = th[1], w = m * m + v;
    const double S = std::log1p(v / (m * m));
    l1[0] = std::log(m) - S / 2;
    l1[1] = S;
    const double Sm = -2 * v / (m * w), Sv = 1 / w;
    const double Smm = 2 * v * (3 * m * m + v) / (m * w * m * w);
    const double Svv = -1 / (w * w), Smv = -2 * m / (w * w);
    // u1[k + 2 a]: d u_k / d theta_a, u = (mu, S)
    u1[0] = 1 / m - Sm / 2; u1[1] = Sm;
    u1[2] = -Sv / 2;        u1[3] = Sv;
    if (u2) {
        // u2[k + 2 r], r over (mm, vv, mv)
        u2[0] = -1 / (m * m) - Smm / 2; u2[1] = Smm;
        u2[2] = -Svv / 2;               u2[3] = Svv;
        u2[4] = -Smv / 2;               u2[5] = Smv;
    }
}
inline double lognormal2_cdf(double q, const double* th, bool lower) {
    double l1[2], u1[4];
    lognormal2_map(th, l1, u1, nullptr);
    return R::plnorm(q, l1[0], std::sqrt(l1[1]), lower, 0);
}
inline void lognormal2_cdf_grad(double q, const double* th, double* g) {
    double l1[2], u1[4], G[2];
    lognormal2_map(th, l1, u1, nullptr);
    lognormal1_cdf_grad(q, l1, G);
    g[0] = G[0] * u1[0] + G[1] * u1[1];
    g[1] = G[0] * u1[2] + G[1] * u1[3];
}
inline void lognormal2_cdf_hess(double q, const double* th, double* h) {
    double l1[2], u1[4], u2[6], G[2], H[3];
    lognormal2_map(th, l1, u1, u2);
    lognormal1_cdf_grad(q, l1, G);
    lognormal1_cdf_hess(q, l1, H);
    // H in (mu, S): [mumu, SS, muS]
    auto quad = [&](const double* ua, const double* ub) {
        return H[0] * ua[0] * ub[0] + H[1] * ua[1] * ub[1] +
            H[2] * (ua[0] * ub[1] + ua[1] * ub[0]);
    };
    const double* um = u1;
    const double* uv = u1 + 2;
    h[0] = quad(um, um) + G[0] * u2[0] + G[1] * u2[1];
    h[1] = quad(uv, uv) + G[0] * u2[2] + G[1] * u2[3];
    h[2] = quad(um, uv) + G[0] * u2[4] + G[1] * u2[5];
}

// ---- InvGauss1 (mu, phi): F = Phi(a) + B, B = exp(2/(phi mu)) Phi(b),
// a = (q/mu - 1)/sqrt(phi q), b = -(q/mu + 1)/sqrt(phi q); with
// exp(2/(phi mu)) phi(b) = phi(a) the derivatives are
// F_mu = -2B/(phi mu^2), F_phi = phi(a)/(phi r) - 2B/(phi^2 mu), r = sqrt(phi q),
// B_mu = -2B/(phi mu^2) + phi(a) b_mu, B_phi = -2B/(phi^2 mu) - phi(a) b/(2 phi),
// b_mu = q/(mu^2 r); F_mumu = -2/phi (B_mu/mu^2 - 2B/mu^3),
// F_phiphi = T (a^2 - 3)/(2 phi) - 2/mu (B_phi/phi^2 - 2B/phi^3),
// T = phi(a)/(phi r), F_muphi = -2/mu^2 (B_phi/phi - B/phi^2)

struct IG1 { double a, b, B, pa, r; bool in; };
inline IG1 invgauss1_parts(double q, const double* th) {
    IG1 p;
    p.in = q > 0;
    if (!p.in) return p;
    const double mu = th[0], ph = th[1];
    p.r = std::sqrt(ph * q);
    p.a = (q / mu - 1) / p.r;
    p.b = -(q / mu + 1) / p.r;
    p.B = std::exp(2 / (ph * mu) + R::pnorm(p.b, 0.0, 1.0, 1, 1));
    p.pa = R::dnorm(p.a, 0.0, 1.0, 0);
    return p;
}
inline double invgauss1_cdf(double q, const double* th, bool lower) {
    const IG1 p = invgauss1_parts(q, th);
    if (!p.in) return lower ? 0.0 : 1.0;
    return lower ? R::pnorm(p.a, 0.0, 1.0, 1, 0) + p.B
                 : R::pnorm(p.a, 0.0, 1.0, 0, 0) - p.B;
}
inline void invgauss1_cdf_grad(double q, const double* th, double* g) {
    const IG1 p = invgauss1_parts(q, th);
    if (!p.in) { g[0] = g[1] = 0.0; return; }
    const double mu = th[0], ph = th[1];
    g[0] = -2 * p.B / (ph * mu * mu);
    g[1] = p.pa / (ph * p.r) - 2 * p.B / (ph * ph * mu);
}
inline void invgauss1_cdf_hess(double q, const double* th, double* h) {
    const IG1 p = invgauss1_parts(q, th);
    if (!p.in) { h[0] = h[1] = h[2] = 0.0; return; }
    const double mu = th[0], ph = th[1];
    const double bmu = q / (mu * mu * p.r);
    const double Bmu = -2 * p.B / (ph * mu * mu) + p.pa * bmu;
    const double Bph = -2 * p.B / (ph * ph * mu) - p.pa * p.b / (2 * ph);
    const double T = p.pa / (ph * p.r);
    h[0] = -2 / ph * (Bmu / (mu * mu) - 2 * p.B / (mu * mu * mu));
    h[1] = T * (p.a * p.a - 3) / (2 * ph) -
        2 / mu * (Bph / (ph * ph) - 2 * p.B / (ph * ph * ph));
    h[2] = -2 / (mu * mu) * (Bph / ph - p.B / (ph * ph));
}

// ---- InvGauss2 (mu, lambda): InvGauss1 at phi = 1/lambda
inline double invgauss2_cdf(double q, const double* th, bool lower) {
    const double t1[2] = {th[0], 1 / th[1]};
    return invgauss1_cdf(q, t1, lower);
}
inline void invgauss2_cdf_grad(double q, const double* th, double* g) {
    const double lam = th[1], t1[2] = {th[0], 1 / lam};
    double G[2];
    invgauss1_cdf_grad(q, t1, G);
    g[0] = G[0];
    g[1] = -G[1] / (lam * lam);
}
inline void invgauss2_cdf_hess(double q, const double* th, double* h) {
    const double lam = th[1], t1[2] = {th[0], 1 / lam};
    double G[2], H[3];
    invgauss1_cdf_grad(q, t1, G);
    invgauss1_cdf_hess(q, t1, H);
    const double l2 = lam * lam;
    h[0] = H[0];
    h[1] = H[1] / (l2 * l2) + 2 * G[1] / (l2 * lam);
    h[2] = -H[2] / l2;
}

// ---- Enet (mu, lambda, alpha): with a = lambda alpha, c = lambda (1 - alpha),
// s = sqrt(c), x = a/s and z = q - mu, the distribution function is
// Phi(u)/(2 Q) for z <= 0 and 1 - Phi(v)/(2 Q) above, u = s z - x,
// v = -s z - x, Q = Phi(-x), each term formed on the log scale. A term
// P = Phi(w)/(2Q) has d_i log P = D_i = R(w) w_i + G(x) x_i and
// d_ij P = P (D_i D_j + D_ij), D_ij = R'(w) w_i w_j + R(w) w_ij
// + G'(x) x_i x_j + G(x) x_ij, R = phi/Phi, R' = -R (w + R),
// G = phi(x)/Phi(-x), G' = G (G - x); the derivatives in (mu, a, c) are
// then carried to (mu, lambda, alpha).

struct EnetTerm { double P, D[3], DD[6]; };  // DD: 00 11 22 01 02 12

inline EnetTerm enet_term(double q, const double* th, bool low) {
    const double mu = th[0], lam = th[1], al = th[2];
    const double a = lam * al, c = lam * (1 - al), s = std::sqrt(c);
    const double x = a / s, z = q - mu;
    const double sg = low ? 1.0 : -1.0;
    const double w = sg * s * z - x;
    const double lQ = R::pnorm(-x, 0.0, 1.0, 1, 1);
    EnetTerm t;
    t.P = std::exp(R::pnorm(w, 0.0, 1.0, 1, 1) - lQ) / 2;
    const double Rw = std::exp(R::dnorm(w, 0.0, 1.0, 1) - R::pnorm(w, 0.0, 1.0, 1, 1));
    const double Gx = std::exp(R::dnorm(x, 0.0, 1.0, 1) - lQ);
    const double R1 = -Rw * (w + Rw), G1 = Gx * (Gx - x);
    // partials of w and x in (mu, a, c)
    const double w1[3] = {-sg * s, -1 / s, sg * z / (2 * s) + a / (2 * c * s)};
    const double x1[3] = {0.0, 1 / s, -a / (2 * c * s)};
    const double w2[6] = {0.0, 0.0, -sg * z / (4 * c * s) - 3 * a / (4 * c * c * s),
                          0.0, -sg / (2 * s), 1 / (2 * c * s)};
    const double x2[6] = {0.0, 0.0, 3 * a / (4 * c * c * s), 0.0, 0.0,
                          -1 / (2 * c * s)};
    for (int i = 0; i < 3; ++i) t.D[i] = Rw * w1[i] + Gx * x1[i];
    const int I[6] = {0, 1, 2, 0, 0, 1}, J[6] = {0, 1, 2, 1, 2, 2};
    for (int r = 0; r < 6; ++r)
        t.DD[r] = R1 * w1[I[r]] * w1[J[r]] + Rw * w2[r] +
            G1 * x1[I[r]] * x1[J[r]] + Gx * x2[r];
    return t;
}

inline double enet_cdf(double q, const double* th, bool lower) {
    const bool low = q - th[0] <= 0;
    const double P = enet_term(q, th, low).P;
    return (low == lower) ? P : 1 - P;
}

// the derivatives of F in (mu, a, c): g3[3] and h3[6] (00 11 22 01 02 12)
inline void enet_cdf_mac(double q, const double* th, double* g3, double* h3) {
    const bool low = q - th[0] <= 0;
    const EnetTerm t = enet_term(q, th, low);
    const double sg = low ? 1.0 : -1.0;
    const int I[6] = {0, 1, 2, 0, 0, 1}, J[6] = {0, 1, 2, 1, 2, 2};
    for (int i = 0; i < 3; ++i) g3[i] = sg * t.P * t.D[i];
    if (h3)
        for (int r = 0; r < 6; ++r)
            h3[r] = sg * t.P * (t.D[I[r]] * t.D[J[r]] + t.DD[r]);
}

inline void enet_cdf_grad(double q, const double* th, double* g) {
    const double lam = th[1], al = th[2];
    double g3[3];
    enet_cdf_mac(q, th, g3, nullptr);
    g[0] = g3[0];
    g[1] = g3[1] * al + g3[2] * (1 - al);
    g[2] = lam * (g3[1] - g3[2]);
}

inline void enet_cdf_hess(double q, const double* th, double* h) {
    const double lam = th[1], al = th[2], be = 1 - al;
    double g3[3], h3[6];
    enet_cdf_mac(q, th, g3, h3);
    const double Faa = h3[1], Fcc = h3[2], Fac = h3[5];
    h[0] = h3[0];
    h[1] = al * al * Faa + 2 * al * be * Fac + be * be * Fcc;
    h[2] = lam * lam * (Faa - 2 * Fac + Fcc);
    h[3] = al * h3[3] + be * h3[4];
    h[4] = lam * (h3[3] - h3[4]);
    h[5] = lam * al * Faa + lam * (1 - 2 * al) * Fac - lam * be * Fcc +
        g3[1] - g3[2];
}

} // namespace d7

#endif
