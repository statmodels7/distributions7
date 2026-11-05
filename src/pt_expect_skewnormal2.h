// skewnormal2: generator, density and quantile for the compiled expected
// information by simulation and by integration (see pt_expect_common.h).
// th = (mu, sigma, gamma1), the centered parametrization. Every method is
// skewnormal1's at the direct parameters sn2_theta(theta)
// (R/skewnormal2_distrib.R:337-343 with sn_cp_to_dp() :192-205,
// sn_one_minus_delta2() :133-135, sn_max_skew() :79-82): distrib_pdf :475,
// distrib_rng :611-613, distrib_quantile :572 (the numerical fallback of
// skewnormal1 at dp).
//
// sn2_dp() returns nonzero where sn2_reject_unmappable() (:387-414) signals
// an error. The generator and the quantile then return nonzero; the density
// writes NaN (R errors, and the engine reaches the density only after the
// quantile or generator check).

#ifndef D7_PT_EXPECT_SKEWNORMAL2_H
#define D7_PT_EXPECT_SKEWNORMAL2_H

#include "pt_expect_common.h"
#include "pt_expect_skewnormal1.h"

namespace d7x {

inline int sn2_dp(const double* th, double* dp) {
    const double g = th[2];
    const double s = g >= 0 ? 1.0 : -1.0;
    const double b = std::sqrt(2 / M_PI);
    const double smax = (4 - M_PI) / 2 * R_POW(b / std::sqrt(1 - b * b), 3);
    const double cc = s * R_POW(2 * (s * g) / (4 - M_PI), 1.0 / 3);
    const double omd2 = -r_expm1((2.0 / 3) * std::log(s * g / smax));
    dp[0] = th[0] - th[1] * cc;
    dp[1] = th[1] * std::sqrt(1 + cc * cc);
    dp[2] = cc / (b * std::sqrt(omd2));
    // sn2_reject_unmappable(): all finite, sigma > 0
    if (R_FINITE(dp[0]) && R_FINITE(dp[1]) && dp[1] > 0 && R_FINITE(dp[2])) return 0;
    return 1;
}

// nonzero at gamma1 == 0 (the observed Hessian kernel errors there,
// sn2_reject_zero(), :353-366) and where the map is not finite
inline int skewnormal2_rng_ok(const double* th) {
    if (th[2] == 0) return 1;
    double dp[3];
    return sn2_dp(th, dp);
}

inline void skewnormal2_rng(const double* th, int n, double* out) {
    double dp[3];
    sn2_dp(th, dp);
    sn1_rng_dp(dp[0], dp[1], dp[2], n, out);
}

inline void skewnormal2_pdf(const double* y, int n, const double* th, int vb,
                            double* out) {
    (void) vb;
    double dp[3];
    if (sn2_dp(th, dp)) {
        for (int i = 0; i < n; ++i) out[i] = R_NaN;
        return;
    }
    sn1_pdf_dp(y, n, dp[0], dp[1], dp[2], out);
}

inline int skewnormal2_quantile(const double* p, int np, const double* th,
                                int vb, double* out) {
    (void) vb;
    double dp[3];
    if (sn2_dp(th, dp)) return 1;
    Sn1Fam F(dp[0], dp[1], dp[2]);
    return qf_quantile(p, np, F, out);
}

}  // namespace d7x

#endif
