// The generator, the density and the quantile of the pseudo-Huber family in
// (mu, sigma, nu), value for value as the R methods compute them. The cores
// (pseudohuber_*) serve pt_expect_pseudohuber2.h as well, which includes this
// file.
//
// Sources transcribed:
//   R/pseudohuber_distrib.R:380 (rng -> pseudohuber_rng_cpp,
//     src/pseudohuber_rng.cpp), :161 (pdf), :283 (quantile),
//     :309 (symmetric_log_quantile()), R/moments.R:973 (variance(), the
//     bracket's spread);
//   the log distribution function of the method distrib_cdf =
//     compiled_cdf (R/cdf_compiled.R:67, compiled_log_cdf() :105,
//     compiled_far_tail() :84), d7_cdf_cpp -> cont_cdf() -> quad_cdf() ->
//     quad_cdf_sums() of src/d7_ccallable.cpp (family ids 36 and 37), and
//     quad_log_tail_prob()/quad_tail_parts() of R/log_tail.R:403, 413,
//     log_tail_nodes() :460;
//   R/utility_functions.R:620 (quantile_log_tails), :536 (log1mexp_neg).

#ifndef D7_PT_EXPECT_PSEUDOHUBER_H
#define D7_PT_EXPECT_PSEUDOHUBER_H

#include <algorithm>
#include "pt_expect_common.h"
#include "pt_sqrt.h"
#include "pt_pseudohuber.h"
#include "pt_pseudohuber2.h"
#include "pt_trunc_rule.h"

namespace d7x {

// ---------------------------------------------------------------------------
// The generator: the per-draw body of pseudohuber_rng_cpp() (src/
// pseudohuber_rng.cpp), copied unchanged. The box cache keyed on omega does
// not change a value.

struct pseudohuber_RouBox {
    double um, up;
};

inline double pseudohuber_cubic(double x, double a) {
    return ((x + a) * x - 1.0) * x + 1.0;
}
inline double pseudohuber_cubic_d(double x, double a) {
    return (3.0 * x + 2.0 * a) * x - 1.0;
}

inline pseudohuber_RouBox pseudohuber_rou_box(double omega) {
    const double a = -(1.0 + 4.0 / omega), b = -1.0, c = 1.0;
    // depressed cubic t^3 + p t + q, x = t - a/3
    const double p = b - a * a / 3.0;
    const double q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
    double arg = -q / 2.0 * d7::sqrt_cr(-27.0 / (p * p * p));
    arg = std::max(-1.0, std::min(1.0, arg));
    const double phi = std::acos(arg);
    const double r = d7::sqrt_cr(-4.0 * p / 3.0);
    double xs[3];
    for (int k = 0; k < 3; k++)
        xs[k] = r * std::cos((phi + 2.0 * M_PI * k) / 3.0) - a / 3.0;
    double xm = NA_REAL, xp = NA_REAL;
    for (int k = 0; k < 3; k++) {
        if (xs[k] > 0.0 && xs[k] < 1.0) xm = xs[k];
        if (xs[k] > 1.0) xp = xs[k];
    }
    // two Newton steps on the cubic sharpen each root
    for (int it = 0; it < 2; it++) {
        xm -= pseudohuber_cubic(xm, a) / pseudohuber_cubic_d(xm, a);
        xp -= pseudohuber_cubic(xp, a) / pseudohuber_cubic_d(xp, a);
    }
    pseudohuber_RouBox box;
    box.um = (xm - 1.0) * std::exp(-omega * (xm - 1.0) * (xm - 1.0) / (4.0 * xm));
    box.up = (xp - 1.0) * std::exp(-omega * (xp - 1.0) * (xp - 1.0) / (4.0 * xp));
    return box;
}

inline double pseudohuber_gig1_draw(double omega, const pseudohuber_RouBox& box) {
    while (true) {
        const double u = box.um + (box.up - box.um) * unif_rand();
        const double v = unif_rand();
        const double x = u / v + 1.0;
        if (x <= 0.0) continue;
        // accept when v^2 <= h(x), on the log scale
        if (2.0 * std::log(v) <= -omega * (x - 1.0) * (x - 1.0) / (2.0 * x)) return x;
    }
}

// pseudohuber_rng_cpp(n, mu, sigma1, nu) at scalars
inline void pseudohuber_rng_core(double m, double s1, double v, int n,
                                 double* out) {
    double last_omega = NA_REAL;
    pseudohuber_RouBox box = {0.0, 0.0};
    for (int i = 0; i < n; i++) {
        const double omega = d7::sqrt_cr(v);
        double w;
        if (omega < 1e-8) {
            w = 2.0 * s1 * s1 * exp_rand();
        } else {
            if (!(omega == last_omega)) {
                box = pseudohuber_rou_box(omega);
                last_omega = omega;
            }
            w = s1 * s1 * omega * pseudohuber_gig1_draw(omega, box);
        }
        out[i] = m + d7::sqrt_cr(w) * norm_rand();
    }
}

// th = (mu, sigma, nu): sigma and nu finite and positive, mu finite (R signals
// an error or returns NaN otherwise)
inline int pseudohuber_rng_ok(const double* th) {
    return (R_FINITE(th[0]) && th[1] > 0 && R_FINITE(th[1]) && th[2] > 0 &&
            R_FINITE(th[2])) ? 0 : 1;
}

// R/pseudohuber_distrib.R:380: pseudohuber_rng_cpp(n, mu, sigma, nu)
inline void pseudohuber_rng(const double* th, int n, double* out) {
    pseudohuber_rng_core(th[0], th[1], th[2], n, out);
}

// ---------------------------------------------------------------------------
// The density

// R's pmin() of two doubles: NaN if either is NaN
inline double pseudohuber_r_pmin(double a, double b) {
    if (ISNAN(a)) return a;
    if (ISNAN(b)) return b;
    return a <= b ? a : b;
}

// R/pseudohuber_distrib.R:161, log = TRUE; besselK(x, 1, expon.scaled = TRUE)
// is R::bessel_k(x, 1, 2)
inline double pseudohuber_logpdf_r(double y, double mu, double sigma, double nu) {
    const double sq_nu = d7::sqrt_cr(nu);
    const double az = std::fabs((y - mu) / sigma);
    const double big = n7::r_pmax(az, sq_nu);
    const double rr = pseudohuber_r_pmin(az, sq_nu) / big;
    const double D = big * d7::sqrt_cr(1 + rr * rr);
    double log_norm = std::log(2.0) + std::log(sigma);
    log_norm = log_norm + 0.5 * std::log(nu);
    log_norm = log_norm + std::log(R::bessel_k(sq_nu, 1.0, 2.0));
    log_norm = log_norm - sq_nu;
    return -D - log_norm;
}

// R/pseudohuber_distrib.R:161, log = FALSE
inline void pseudohuber_pdf(const double* y, int n, const double* th, int vb,
                            double* out) {
    (void) vb;
    for (int i = 0; i < n; ++i)
        out[i] = r_exp(pseudohuber_logpdf_r(y[i], th[0], th[1], th[2]));
}

// ---------------------------------------------------------------------------
// The distribution function on the log scale, as compiled_cdf(log.p = TRUE)
// computes it for one quantile at the lower tail. id 36 is pseudohuber and
// 37 pseudohuber2; th = (mu, sigma, nu). The family's bounds are (-Inf, Inf).

// d7_logpdf(id, y, th) (pt_logpdf.h)
inline double pseudohuber_cdf_logpdf(int id, double y, const double* th) {
    if (id == 36)
        return d7::pseudohuber_logpdf(y, th[0], th[1], th[2], d7::ph_nu(th[2]));
    return d7::pseudohuber2_logpdf(y, th[0], th[1], th[2],
                                   d7::pseudohuber2_terms(th[2]));
}

// quad_cdf_sums() with mask 0 for a family on the real line, location mu and
// scale sigma (center_scale(), src/d7_ccallable.cpp:335), no kinks, no head:
// the integral of the density over the side of q, `upper` true for [q, Inf)
inline double pseudohuber_cdf_sums(int id, const double* th, double q,
                                   bool upper) {
    const double slo = R_NegInf, shi = R_PosInf, c = th[0], s = th[1];
    std::vector<double> ys, ws;
    const double a = upper ? std::max(q, slo) : slo;
    const double b = upper ? shi : std::min(q, shi);
    double kk[1] = {0.0};
    if (b > a)
        d7::trunc_rule(a, b, c, s, kk, 0, ys, ws, d7::kTsStep, d7::kDoublings,
                       R_NaN);
    long double S0 = 0.0L;
    for (std::size_t j = 0; j < ys.size(); ++j) {
        const double y = ys[j];
        const double fw = std::exp(pseudohuber_cdf_logpdf(id, y, th)) * ws[j];
        if (fw == 0.0) continue;
        S0 += fw;
    }
    return (double) S0;
}

// quad_cdf(): the side is q > center, flipped when the integral exceeds 1/2
inline double pseudohuber_cdf_quad(int id, const double* th, double q,
                                   bool lower) {
    bool up = q > th[0];
    double I0 = pseudohuber_cdf_sums(id, th, q, up);
    if (I0 > 0.5) {
        up = !up;
        I0 = pseudohuber_cdf_sums(id, th, q, up);
    }
    if (up) return lower ? 1 - I0 : I0;
    return lower ? I0 : 1 - I0;
}

// distrib_grad_y() of the two families (R/pseudohuber_distrib.R:802,
// R/pseudohuber2_distrib.R:387)
inline double pseudohuber_cdf_grad_y(int id, double y, const double* th) {
    const double mu = th[0], sigma = th[1], nu = th[2];
    if (id == 36) {
        const double r = y - mu;
        const double sigma2 = sigma * sigma;
        const double D = d7::sqrt_cr(nu + (r * r) / sigma2);
        return (-r) / (sigma2 * D);
    }
    const double R = d7::pseudohuber2_terms(nu).R;
    const double r = y - mu;
    const double D = d7::sqrt_cr(nu + (R * (r * r)) / (sigma * sigma));
    return ((-R) * r) / ((sigma * sigma) * D);
}

// distrib_pdf(log = TRUE) of the two families (R/pseudohuber_distrib.R:161,
// R/pseudohuber2_distrib.R:125); R2 and h2 are the family-2 terms of nu
inline double pseudohuber_cdf_logpdf_r(int id, double y, const double* th,
                                       double R2, double h2) {
    if (id == 36) return pseudohuber_logpdf_r(y, th[0], th[1], th[2]);
    const double sigma = th[1], nu = th[2];
    const double az = d7::sqrt_cr(R2) * std::fabs((y - th[0]) / sigma);
    const double sq_nu = d7::sqrt_cr(nu);
    const double big = n7::r_pmax(az, sq_nu);
    const double rr = pseudohuber_r_pmin(az, sq_nu) / big;
    const double D = big * d7::sqrt_cr(1 + rr * rr);
    return (((-D) - std::log(2.0)) - std::log(sigma)) + h2;
}

// quad_log_tail_prob(distrib, q, th, lower.tail = TRUE) for one q
// (R/log_tail.R:403, quad_tail_parts() :413), bounds (-Inf, Inf)
inline double pseudohuber_cdf_log_tail(int id, const double* th, double q) {
    // log_tail_nodes() (R/log_tail.R:460): GK15 on 16 intervals
    static std::vector<double> U, W;
    if (U.empty()) {
        const double br[17] = {0, 0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32,
                               48, 64, 96, 128};
        const n7::GK15& g = n7::gk15();
        U.resize(240);
        W.resize(240);
        for (int j = 0; j < 16; ++j) {
            const double a = br[j];
            const double h = (br[j + 1] - br[j]) / 2;
            for (int i = 0; i < 15; ++i) {
                U[j * 15 + i] = g.node[i] * h + (a + h);
                W[j * 15 + i] = g.kw[i] * h;
            }
        }
    }
    const int K = 240;
    const double sg = -1.0;
    const double gy = sg * pseudohuber_cdf_grad_y(id, q, th);
    double s = (R_FINITE(gy) && gy < 0) ? (-1.0) / gy : R_PosInf;
    if (!R_FINITE(s)) s = 1.0;
    double R2 = 0.0, h2 = 0.0;
    if (id == 37) {
        const d7::Ph2Terms T = d7::pseudohuber2_terms(th[2]);
        R2 = T.R;
        h2 = T.h;
    }
    const double lfq = pseudohuber_cdf_logpdf_r(id, q, th, R2, h2);
    long double tot = 0.0L;
    const double sgs = sg * s;
    for (int k = 0; k < K; ++k) {
        double Y = q + sgs * r_expm1(U[k]);
        const double jac = s * r_exp(U[k]);
        const bool ok = R_FINITE(Y) && (Y < q) && Y > R_NegInf && Y < R_PosInf;
        if (!ok) Y = q;
        const double lf = pseudohuber_cdf_logpdf_r(id, Y, th, R2, h2);
        double Wt = r_exp(lf - lfq) * jac * W[k];
        if (!ok || !R_FINITE(Wt)) Wt = 0.0;
        tot += Wt;
    }
    return lfq + std::log((double) tot);
}

// compiled_cdf(distrib, q, th, lower.tail = TRUE, log.p = TRUE) for one q
// (R/cdf_compiled.R:67, compiled_log_cdf() :105)
inline double pseudohuber_cdf_log(int id, const double* th, double q) {
    const double out = pseudohuber_cdf_quad(id, th, q, true);
    double lo = std::log(out);
    const bool big = !ISNAN(out) && out > 0.5;
    if (big) {
        const double other = pseudohuber_cdf_quad(id, th, q, false);
        lo = r_log1p(-other);
    }
    // compiled_far_tail(): the tail below exp(log_tail_cut = -5), inside
    // the support
    const bool far = !ISNAN(out) && out < r_exp(-5.0) && R_FINITE(q) &&
        q > R_NegInf && q < R_PosInf;
    if (far) lo = pseudohuber_cdf_log_tail(id, th, q);
    return lo;
}

// ---------------------------------------------------------------------------
// The quantile

// log1mexp_neg(a) (R/utility_functions.R:536)
inline double pseudohuber_log1mexp_neg(double a) {
    double out = r_log1p(-r_exp(a));
    const bool big = !ISNAN(a) && a > -std::log(2.0);
    if (big) out = std::log(-r_expm1(a));
    return out;
}

// symmetric_log_quantile(distrib, p, theta, lower.tail = TRUE, log.p = FALSE,
// spread) for one p (R/pseudohuber_distrib.R:309), the cdf of family `id`.
// Returns nonzero where uniroot() or the bracketing loop of R would signal an
// error or reach a NaN test.
inline int pseudohuber_sym_quantile(int id, double p, const double* th,
                                    double spread, double* out) {
    // quantile_log_tails(p, TRUE, FALSE)
    double lo_t = std::log(p);
    if (!ISNAN(lo_t) && lo_t > 0) lo_t = R_NaN;
    const double up_t = pseudohuber_log1mexp_neg(lo_t);
    if (ISNAN(lo_t) || ISNAN(up_t)) { *out = R_NaN; return 0; }
    const bool reflect = up_t < lo_t;
    const double ls = reflect ? up_t : lo_t;
    const double m = th[0];
    if (ISNAN(ls) || ls > 0) { *out = R_NaN; return 0; }
    if (ls == R_NegInf) { *out = reflect ? R_PosInf : R_NegInf; return 0; }
    if (ls >= std::log(0.5)) { *out = m; return 0; }
    auto f = [&](double q) -> double {
        return pseudohuber_cdf_log(id, th, q) - ls;
    };
    double lo = m - 10 * spread;
    int it = 0;
    for (;;) {
        const double fl = f(lo);
        if (ISNAN(fl)) return 3;                 // `NA && ...` in R's while()
        if (!(fl > 0 && it < 2000)) break;
        lo = m - 2 * (m - lo);
        it = it + 1;
    }
    double q_low;
    const double flo = f(lo), fm = f(m);
    const int rc = n7::uniroot(f, lo, m, flo, fm, R_POW(DBL_EPSILON, 0.5),
                               1000, &q_low);
    if (rc) return 4;
    *out = reflect ? 2 * m - q_low : q_low;
    return 0;
}

// the spread of the bracket: sqrt(variance(distrib, th)), R/moments.R:973
inline double pseudohuber_spread(const double* th) {
    const double sq_nu = d7::sqrt_cr(th[2]);
    const double ratio = R::bessel_k(sq_nu, 2.0, 2.0) / R::bessel_k(sq_nu, 1.0, 2.0);
    const double v = (th[1] * th[1]) * sq_nu * ratio + 0.0;
    return d7::sqrt_cr(v);
}

// R/pseudohuber_distrib.R:283. Parameters must be valid (sigma, nu > 0).
inline int pseudohuber_quantile(const double* p, int np, const double* th,
                                int vb, double* out) {
    (void) vb;
    if (pseudohuber_rng_ok(th)) return 1;
    const double spread = pseudohuber_spread(th);
    for (int k = 0; k < np; ++k)
        if (pseudohuber_sym_quantile(36, p[k], th, spread, &out[k])) return 1;
    return 0;
}

}  // namespace d7x

#endif
