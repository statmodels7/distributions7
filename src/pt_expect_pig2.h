// The generator and the density of the orthogonal Poisson-inverse Gaussian
// family in (mu, alpha), a discrete family: its quantile is never called and
// returns nonzero.
//
// Sources: R/pig2_distrib.R:472 (rng), :118 (pig2_sigma), :176 (pdf,
// exp(pig2_pdf_cpp(...))); src/pig_hd.cpp:383 (pig2_pdf_cpp, which calls
// pig2_row() of pt_pig_rows.h with order 0, only = true).

#ifndef D7_PT_EXPECT_PIG2_H
#define D7_PT_EXPECT_PIG2_H

#include "pt_expect_pig1.h"

namespace d7x {

// R/pig2_distrib.R:118: (mu + sqrt(mu^2 + alpha^2)) / alpha^2
inline double pig2_sigma_x(double mu, double alpha) {
    return (mu + d7::sqrt_cr(mu * mu + alpha * alpha)) / (alpha * alpha);
}

// th = (mu, alpha); R: sg = pig2_sigma(mu, alpha), shape = mu / sg
inline int pig2_rng_ok(const double* th) {
    const double sg = pig2_sigma_x(th[0], th[1]);
    return invgauss_params_ok(th[0], 1 / (th[0] / sg)) ? 0 : 1;
}
inline void pig2_rng(const double* th, int n, double* out) {
    const double sg = pig2_sigma_x(th[0], th[1]);
    pig_rng_core(th[0], th[0] / sg, n, out);
}

// R/pig2_distrib.R:176: exp(pig2_pdf_cpp(y, mu, alpha, threads))
inline void pig2_pdf(const double* y, int n, const double* th, int vb,
                     double* out) {
    (void) vb;
    for (int i = 0; i < n; ++i) {
        double v[1];
        d7pigrow::pig2_row(y[i], th[0], th[1], 0, true, v);
        out[i] = r_exp(v[0]);
    }
}

// discrete: never called
inline int pig2_quantile(const double* p, int np, const double* th, int vb,
                         double* out) {
    (void) p; (void) np; (void) th; (void) vb; (void) out;
    return 1;
}

}  // namespace d7x

#endif
