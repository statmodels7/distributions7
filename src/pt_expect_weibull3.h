// The Weibull in the mean and shape: generator, density and quantile for the
// compiled expected information (see pt_expect_common.h). th = (mean, sigma).
// R/weibull3_distrib.R tests weibull3_scale_ok() over the whole vector; here
// the test is weibull3_vb() at one setting.
#ifndef D7_PT_EXPECT_WEIBULL3_H
#define D7_PT_EXPECT_WEIBULL3_H

#include "pt_expect_common.h"

namespace d7x {

// R/weibull3_distrib.R:91-93, weibull3_log_scale()
inline double weibull3_log_scale(const double* th) {
    return std::log(th[0]) - R::lgammafn(1 + 1 / th[1]);
}

// R/weibull3_distrib.R:97-100, weibull3_scale_ok() at one setting: 1 where
// the scale exp(lb) is finite and positive.
inline int weibull3_vb(const double* th) {
    double b = r_exp(weibull3_log_scale(th));
    return (R_FINITE(b) && b > 0) ? 1 : 0;
}

// R/weibull3_distrib.R:194-202; the branch is per setting.
inline int weibull3_rng_ok(const double* th) { (void) th; return 0; }

// R/weibull3_distrib.R:194-202, distrib_rng.
inline void weibull3_rng(const double* th, int n, double* out) {
    double lb = weibull3_log_scale(th);
    double b = r_exp(lb);
    if (R_FINITE(b) && b > 0) {
        for (int i = 0; i < n; i++) out[i] = R::rweibull(th[1], b);
        return;
    }
    for (int i = 0; i < n; i++) out[i] = unif_rand();
    for (int i = 0; i < n; i++)
        out[i] = r_exp(lb + std::log(-std::log(out[i])) / th[1]);
}

// R/weibull3_distrib.R:136-150, distrib_pdf (log = FALSE). vb = 1: the
// dweibull route. A NaN y is not handled as R's logical subscripting would
// (the closed-form route is reproduced for y not NaN).
inline void weibull3_pdf(const double* y, int n, const double* th, int vb,
                         double* out) {
    double lb = weibull3_log_scale(th);
    double s = th[1];
    if (vb) {
        double sc = r_exp(lb);
        for (int i = 0; i < n; i++) out[i] = R::dweibull(y[i], s, sc, 0);
        return;
    }
    for (int i = 0; i < n; i++) {
        double lz = std::log(y[i]) - lb;
        double o = std::log(s) - lb + (s - 1) * lz - r_exp(s * lz);
        if (y[i] < 0 || y[i] == R_PosInf) o = R_NegInf;
        if (y[i] == 0) o = ifelse(s < 1, R_PosInf,
                                  ifelse(s > 1, R_NegInf, -lb));
        out[i] = r_exp(o);
    }
}

// R/weibull3_distrib.R:174-192, distrib_quantile (lower, not log).
inline int weibull3_quantile(const double* p, int np, const double* th, int vb,
                             double* out) {
    double lb = weibull3_log_scale(th);
    double s = th[1];
    if (vb) {
        double sc = r_exp(lb);
        for (int i = 0; i < np; i++) out[i] = R::qweibull(p[i], s, sc, 1, 0);
        return 0;
    }
    for (int i = 0; i < np; i++) {
        double H = -r_log1p(-p[i]);
        out[i] = r_exp(lb + std::log(H) / s);
    }
    return 0;
}

}  // namespace d7x

#endif
