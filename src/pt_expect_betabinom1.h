// The beta-binomial in (mu, sigma): generator and density for the compiled
// expected information (see pt_expect_common.h). th = (mu, sigma, size).
#ifndef D7_PT_EXPECT_BETABINOM1_H
#define D7_PT_EXPECT_BETABINOM1_H

#include "pt_expect_common.h"
#include "pt_betabinom.h"

namespace d7x {

inline int betabinom1_rng_ok(const double* th) { (void) th; return 0; }

// R/betabinom1_distrib.R:380-384, distrib_rng: rbeta(n, a, b) for all n,
// then rbinom(n, size, prob = p).
inline void betabinom1_rng(const double* th, int n, double* out) {
    double a = th[0] / th[1];
    double b = (1 - th[0]) / th[1];
    std::vector<double> p(n);
    for (int i = 0; i < n; i++) p[i] = R::rbeta(a, b);
    for (int i = 0; i < n; i++) out[i] = R::rbinom(th[2], p[i]);
}

// R/betabinom1_distrib.R:143-146, distrib_pdf (log = FALSE), through
// betabinom_logpmf_cpp (src/betabinom.cpp:136-155, per-element body).
inline void betabinom1_pdf(const double* y, int n, const double* th, int vb,
                           double* out) {
    (void) vb;
    double m = th[0], s = th[1], size = th[2];
    double A = m / s, B = (1.0 - m) / s;
    for (int i = 0; i < n; i++) {
        double o;
        if (ISNAN(y[i])) o = y[i];
        else if (y[i] < 0 || y[i] > size || y[i] != std::floor(y[i])) o = R_NegInf;
        else o = d7::bb_log_mass(y[i], A, B, size);
        out[i] = r_exp(o);
    }
}

// Discrete family: the quantile is never called.
inline int betabinom1_quantile(const double* p, int np, const double* th,
                               int vb, double* out) {
    (void) p; (void) np; (void) th; (void) vb; (void) out;
    return 1;
}

}  // namespace d7x

#endif
