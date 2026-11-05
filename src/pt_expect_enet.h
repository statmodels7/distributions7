// The elastic-net family's generator, density and quantile for the compiled
// expected information (see pt_expect_common.h). th = (mu, lambda, alpha).
#ifndef D7_PT_EXPECT_ENET_H
#define D7_PT_EXPECT_ENET_H

#include "pt_expect_common.h"

namespace d7x {

// R/utility_functions.R:536-542, log1mexp_neg() at a scalar.
inline double enet_log1mexp_neg(double a) {
    double out = r_log1p(-r_exp(a));
    bool big = !ISNAN(a) && a > -std::log(2.0);
    if (big) out = std::log(-r_expm1(a));
    return out;
}

// R/enet_distrib.R:32-48, .enet_logM() at a scalar (the branch is on x).
inline double enet_logM(double x) {
    if (!(x > 30)) {
        double xb = x;
        return R::pnorm(-xb, 0.0, 1.0, 1, 1) + (xb * xb) / 2 +
               0.5 * std::log(2 * M_PI);
    }
    double u = 1 / x;
    double u2 = u * u;
    return -std::log(x) +
           r_log1p(u2 * (-1 + u2 * (3 + u2 * (-15 + u2 * 105))));
}

// Generator: see enet_rng() below. It is always reproducible.
inline int enet_rng_ok(const double* th) { (void) th; return 0; }

// R/enet_distrib.R:372-388, distrib_quantile (lower.tail = TRUE, log.p = FALSE),
// with quantile_log_tails() (R/utility_functions.R:620-626), .enet_parts()
// (:53-62, only x and c are used) and .enet_logQ() (:50).
inline double enet_quantile_one(double p, const double* th) {
    // quantile_log_tails(p, TRUE, FALSE)
    double lp = std::log(p);
    if (!ISNAN(lp) && lp > 0) lp = R_NaN;
    double lq0 = enet_log1mexp_neg(lp);
    double lower = lp, upper = lq0;
    // .enet_parts
    double lam = th[1], al = th[2];
    double a = lam * al;
    double cc = lam * (1 - al);
    double x = a / std::sqrt(cc);
    double mu = th[0];
    double s = std::sqrt(cc);
    double lq = R::pnorm(-x, 0.0, 1.0, 1, 1);
    double lo = (x + R::qnorm(lower + std::log(2.0) + lq, 0.0, 1.0, 1, 1)) / s;
    double hi = -(x + R::qnorm(upper + std::log(2.0) + lq, 0.0, 1.0, 1, 1)) / s;
    bool pick = !ISNAN(lower) && lower <= upper;
    if (pick) hi = lo;
    return mu + hi;
}

// R/enet_distrib.R:372-388. Never signals an error in R; returns 0.
inline int enet_quantile(const double* p, int np, const double* th, int vb,
                         double* out) {
    (void) vb;
    for (int i = 0; i < np; i++) out[i] = enet_quantile_one(p[i], th);
    return 0;
}

// R/enet_distrib.R:425-427: distrib_quantile(distrib, stats::runif(n), theta).
// All n uniforms are drawn first, then inverted.
inline void enet_rng(const double* th, int n, double* out) {
    for (int i = 0; i < n; i++) out[i] = unif_rand();
    for (int i = 0; i < n; i++) out[i] = enet_quantile_one(out[i], th);
}

// R/enet_distrib.R:233-240, distrib_pdf (log = FALSE).
inline void enet_pdf(const double* y, int n, const double* th, int vb,
                     double* out) {
    (void) vb;
    double lam = th[1], al = th[2];
    double a = lam * al;
    double cc = lam * (1 - al);
    double x = a / std::sqrt(cc);
    double mu = th[0];
    double lM = enet_logM(x);
    for (int i = 0; i < n; i++) {
        double z = y[i] - mu;
        double log_z = std::log(2.0) - 0.5 * std::log(cc) + lM;
        double log_d = -a * std::fabs(z) - cc * (z * z) / 2 - log_z;
        out[i] = r_exp(log_d);
    }
}

}  // namespace d7x

#endif
