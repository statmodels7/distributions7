// The generator, the density and the quantile of the pseudo-Huber family in
// (mu, sigma, nu) with sigma the standard deviation. The generator is that of
// pt_expect_pseudohuber.h at the scale sigma / sqrt(R(nu)); the quantile is
// symmetric_log_quantile() with the spread sigma on the cdf of family id 37.
//
// Sources: R/pseudohuber2_distrib.R:229 (rng), :125 (pdf), :198 (quantile),
// :165 (distrib_cdf = compiled_cdf); pseudohuber2_nu_terms_cpp()
// (src/pseudohuber2.cpp:219 -> d7::pseudohuber2_terms()).

#ifndef D7_PT_EXPECT_PSEUDOHUBER2_H
#define D7_PT_EXPECT_PSEUDOHUBER2_H

#include "pt_expect_pseudohuber.h"

namespace d7x {

inline int pseudohuber2_rng_ok(const double* th) {
    return pseudohuber_rng_ok(th);
}

// s1 <- sigma / sqrt(pseudohuber2_nu_terms_cpp(nu)$R);
// pseudohuber_rng_cpp(n, mu, s1, nu)
inline void pseudohuber2_rng(const double* th, int n, double* out) {
    const double s1 = th[1] / d7::sqrt_cr(d7::pseudohuber2_terms(th[2]).R);
    pseudohuber_rng_core(th[0], s1, th[2], n, out);
}

// R/pseudohuber2_distrib.R:125, log = FALSE
inline void pseudohuber2_pdf(const double* y, int n, const double* th, int vb,
                             double* out) {
    (void) vb;
    const d7::Ph2Terms T = d7::pseudohuber2_terms(th[2]);
    for (int i = 0; i < n; ++i)
        out[i] = r_exp(pseudohuber_cdf_logpdf_r(37, y[i], th, T.R, T.h));
}

// R/pseudohuber2_distrib.R:198, spread = th[[2]]
inline int pseudohuber2_quantile(const double* p, int np, const double* th,
                                 int vb, double* out) {
    (void) vb;
    if (pseudohuber_rng_ok(th)) return 1;
    for (int k = 0; k < np; ++k)
        if (pseudohuber_sym_quantile(37, p[k], th, th[1], &out[k])) return 1;
    return 0;
}

}  // namespace d7x

#endif
