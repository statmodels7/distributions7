// The beta-binomial in the shapes: generator and density for the compiled
// expected information (see pt_expect_common.h). th = (alpha, beta, size).
#ifndef D7_PT_EXPECT_BETABINOM2_H
#define D7_PT_EXPECT_BETABINOM2_H

#include "pt_expect_common.h"

namespace d7x {

// R/betabinom2_distrib.R:352-358: the branch test of betabinom_log_mass() at
// one setting. It depends on the shapes and the size only, not on the counts.
inline int betabinom2_vb(const double* th) {
    double s = th[0] + th[1];
    return (R_FINITE(s) &&
            R::lgammafn(s + th[2]) * DBL_EPSILON < 1e-8) ? 1 : 0;
}

inline int betabinom2_rng_ok(const double* th) { (void) th; return 0; }

// R/betabinom2_distrib.R:412-415, distrib_rng: rbeta(n, a, b) for all n,
// then rbinom(n, size, prob = p).
inline void betabinom2_rng(const double* th, int n, double* out) {
    std::vector<double> p(n);
    for (int i = 0; i < n; i++) p[i] = R::rbeta(th[0], th[1]);
    for (int i = 0; i < n; i++) out[i] = R::rbinom(th[2], p[i]);
}

// R/betabinom2_distrib.R:352-374, betabinom_log_mass() at a count on the
// support; vb is its whole-vector test.
inline double betabinom2_log_mass(double y, double a, double b, double n,
                                  int vb) {
    if (vb) {
        return R::lchoose(n, y) + (R::lbeta(y + a, n - y + b) - R::lbeta(a, b));
    }
    double s1 = 0.0, s2 = 0.0, s3 = 0.0;
    int N = (int) n;
    for (int j = 0; j < N; j++) {
        s1 = s1 + (double) (j < y) * std::log(a + j);
        s2 = s2 + (double) (j < n - y) * std::log(b + j);
        s3 = s3 + std::log(a + b + j);
    }
    return R::lchoose(n, y) + s1 + s2 - s3;
}

// R/betabinom2_distrib.R:291-308, distrib_pdf (log = FALSE).
inline void betabinom2_pdf(const double* y, int n, const double* th, int vb,
                           double* out) {
    double a = th[0], b = th[1], size = th[2];
    for (int i = 0; i < n; i++) {
        double o;
        if (ISNAN(y[i])) o = y[i];
        else if (y[i] >= 0 && y[i] <= size && y[i] == std::nearbyint(y[i]))
            o = betabinom2_log_mass(y[i], a, b, size, vb);
        else o = R_NegInf;
        out[i] = r_exp(o);
    }
}

// Discrete family: the quantile is never called.
inline int betabinom2_quantile(const double* p, int np, const double* th,
                               int vb, double* out) {
    (void) p; (void) np; (void) th; (void) vb; (void) out;
    return 1;
}

}  // namespace d7x

#endif
