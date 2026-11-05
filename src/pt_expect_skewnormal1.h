// skewnormal1: generator, density and quantile for the compiled expected
// information by simulation and by integration (see pt_expect_common.h).
// th = (mu, sigma, alpha).
//
// The quantile is the numerical fallback of pt_expect_qfallback.h over this
// family's density and distribution function. The distribution function is
// transcribed from R/skewnormal1_distrib.R: distrib_cdf() (:238-253),
// sn_log_lower() (:255-261), the rescue log_tail_integral()
// (R/cdf_compiled.R:124-157, stats::integrate() through Rdqagi) with
// distrib_grad_y() (R/skewnormal1_distrib.R:807-812), and numericals7's
// owen_t() (R/special.R:134-181) whose quad_vec() row is n7::quad_row().

#ifndef D7_PT_EXPECT_SKEWNORMAL1_H
#define D7_PT_EXPECT_SKEWNORMAL1_H

#include "pt_expect_common.h"
#include "pt_expect_qfallback.h"

namespace d7x {

// numericals7::owen_t(h, a) for one pair (R/special.R:134-181); a failed
// quad_vec() row is NA, as quad_vec() returns it
inline double sn_owen_t(double h0, double a0) {
    const double h = std::fabs(h0);
    const double a = a0;
    double out = NA_REAL;
    const bool ok = !ISNAN(h) && !ISNAN(a);
    const double sgn = r_sign(a);
    const double aa = std::fabs(a);
    if (ok && aa == 0) out = 0;
    if (ok && std::isinf(h) && aa > 0) out = 0;
    if (ok && R_FINITE(h) && std::isinf(aa))
        out = sgn * 0.5 * R::pnorm(h, 0.0, 1.0, 0, 0);
    const bool fin = ok && R_FINITE(h) && R_FINITE(aa) && aa > 0;
    if (fin) {
        const bool small = aa <= 1;
        const double hh = small ? h : aa * h;
        const double sl = small ? aa : 1 / aa;
        auto f = [hh](const double* x, int nx, double* o) {
            for (int i = 0; i < nx; ++i) {
                const double xv = x[i];
                o[i] = r_exp(-(hh * hh) * (xv * xv) / 2) / (1 + xv * xv);
            }
        };
        double Q = 0;
        const int bad = n7::quad_row(f, 0.0, sl, 0.0, 1e-13, 4096, &Q);
        const double v = bad ? NA_REAL
                             : r_exp(-(hh * hh) / 2) * Q / (2 * M_PI);
        if (small) {
            out = sgn * v;
        } else {
            const double q1 = R::pnorm(h, 0.0, 1.0, 0, 0);
            const double q2 = R::pnorm(aa * h, 0.0, 1.0, 0, 0);
            out = sgn * ((q1 + q2) / 2 - q1 * q2 - v);
        }
    }
    return out;
}

// sn_log_lower(z, alpha) (R/skewnormal1_distrib.R:255-261)
inline double sn_log_lower(double z, double alpha) {
    const double a = r_pmin(r_pmax(alpha, -1e50), 1e50);
    double p = R::pnorm(z, 0.0, 1.0, 1, 0) - 2 * sn_owen_t(z, a);
    p = r_pmin(r_pmax(p, 0.0), 1.0);
    const double up = R::pnorm(-z, 0.0, 1.0, 1, 0) + 2 * sn_owen_t(z, a);
    if (ISNAN(p)) return NA_REAL;
    if (p > 0.5) return r_log1p(-r_pmin(r_pmax(up, 0.0), 1.0));
    return std::log(p);
}

// the family at one setting, in (mu, sigma, alpha)
struct Sn1Fam {
    double lb, ub;
    bool fail;
    double mu, sg, al;

    Sn1Fam(double mu_, double sg_, double al_)
        : lb(R_NegInf), ub(R_PosInf), fail(false), mu(mu_), sg(sg_), al(al_) {}

    // distrib_pdf(log = TRUE) (:163-172)
    double lpdf(double y) {
        const double z = (y - mu) / sg;
        return std::log(2.0) - std::log(sg) + R::dnorm(z, 0.0, 1.0, 1) +
               R::pnorm(al * z, 0.0, 1.0, 1, 1);
    }
    double pdf(double y) { return r_exp(lpdf(y)); }

    // distrib_grad_y() (:807-812), with mills_ratio()$r
    double grad_y(double y) {
        const double z = (y - mu) / sg;
        const double t = al * z;
        const double r = r_exp(R::dnorm(t, 0.0, 1.0, 1) - R::pnorm(t, 0.0, 1.0, 1, 1));
        return (al * r - z) / sg;
    }

    // log_tail_integral(distrib, q, th, lower.tail) (R/cdf_compiled.R:124)
    double log_tail_integral(double q, bool lower) {
        const double lq = lpdf(q);
        if (!R_FINITE(lq)) return R_NegInf;
        const double dir = lower ? -1.0 : 1.0;
        const double room = lower ? q - lb : ub - q;
        const double gy = std::fabs(grad_y(q));
        double h = (R_FINITE(gy) && gy > 0) ? 1 / gy : 1.0;
        for (int it = 1; it <= 200; ++it) {
            if (R_FINITE(room) && h > room) h = room;
            const double drop = lq - lpdf(q + dir * h);
            if (ISNAN(drop)) break;
            if (drop > 50) h = h * 50 / drop;
            else if (drop < 0.5 && !(R_FINITE(room) && h >= room)) h = h * 2;
            else break;
        }
        auto g = [&](double t) {
            double v = r_exp(lpdf(q + dir * h * t) - lq);
            if (!R_FINITE(v)) v = 0;
            return v;
        };
        double val = NA_REAL;
        int bad;
        if (R_FINITE(room)) {
            bad = r_integrate_fin(g, 0.0, room / h, 1e-12, 1e-12, 1000, &val);
        } else {
            bad = r_integrate_inf1(g, 0.0, 1e-12, 1e-12, 1000, &val);
        }
        if (bad) val = NA_REAL;
        if (!R_FINITE(val) || val <= 0) {
            return (R_FINITE(gy) && gy > 0) ? lq - std::log(gy) : R_NegInf;
        }
        return lq + std::log(h) + std::log(val);
    }

    // distrib_cdf() (:238-253)
    double cdf_any(double q, bool lower, bool logp) {
        const double z = (q - mu) / sg;
        const double sgn = lower ? 1.0 : -1.0;
        double lp = sn_log_lower(sgn * z, sgn * al);
        const bool redo = ISNAN(lp) ||
            lp < std::log(1e-3) + R::pnorm(sgn * z, 0.0, 1.0, 1, 1);
        if (redo && R_FINITE(q)) lp = log_tail_integral(q, lower);
        return logp ? lp : r_exp(lp);
    }
    double lcdf(double x, bool lower) { return cdf_any(x, lower, true); }
    double cdf(double x) { return cdf_any(x, true, false); }
};

// the generator of distrib_rng(SkewNormal1Distrib) (:334-344) at (mu, sigma,
// alpha): u = rnorm(n), then v = rnorm(n)
inline void sn1_rng_dp(double mu, double sigma, double alpha, int n, double* out) {
    const bool big = std::fabs(alpha) > 1;
    const double ia = 1 / std::fabs(alpha);
    const double delta = ifelse(big, r_sign(alpha) / std::sqrt(1 + ia * ia),
                                alpha / std::sqrt(1 + alpha * alpha));
    const double cdelta = ifelse(big, ia / std::sqrt(1 + ia * ia),
                                 1 / std::sqrt(1 + alpha * alpha));
    std::vector<double> u(n), v(n);
    for (int i = 0; i < n; ++i) u[i] = R::rnorm(0.0, 1.0);
    for (int i = 0; i < n; ++i) v[i] = R::rnorm(0.0, 1.0);
    for (int i = 0; i < n; ++i) {
        const double z = delta * std::fabs(u[i]) + cdelta * v[i];
        out[i] = mu + sigma * z;
    }
}

// density at (mu, sigma, alpha), log = FALSE (:163-172)
inline void sn1_pdf_dp(const double* y, int n, double mu, double sigma,
                       double alpha, double* out) {
    for (int i = 0; i < n; ++i) {
        const double z = (y[i] - mu) / sigma;
        const double log_d = std::log(2.0) - std::log(sigma) +
                             R::dnorm(z, 0.0, 1.0, 1) +
                             R::pnorm(alpha * z, 0.0, 1.0, 1, 1);
        out[i] = r_exp(log_d);
    }
}

// rng_ok: the generator has no case in which R signals an error.
inline int skewnormal1_rng_ok(const double* th) {
    (void) th;
    return 0;
}

inline void skewnormal1_rng(const double* th, int n, double* out) {
    sn1_rng_dp(th[0], th[1], th[2], n, out);
}

inline void skewnormal1_pdf(const double* y, int n, const double* th, int vb,
                            double* out) {
    (void) vb;
    sn1_pdf_dp(y, n, th[0], th[1], th[2], out);
}

inline int skewnormal1_quantile(const double* p, int np, const double* th,
                                int vb, double* out) {
    (void) vb;
    Sn1Fam F(th[0], th[1], th[2]);
    return qf_quantile(p, np, F, out);
}

}  // namespace d7x

#endif
