// The generator, the density and the quantile of the inverse Gaussian family
// in (mean, shape): the cores of pt_expect_invgauss1.h (statmod's rinvgauss,
// dinvgauss and qinvgauss) at the dispersion 1 / lambda.
//
// Sources: R/invgauss2_distrib.R:144 (pdf), :259 (quantile), :303 (rng).

#ifndef D7_PT_EXPECT_INVGAUSS2_H
#define D7_PT_EXPECT_INVGAUSS2_H

#include "pt_expect_invgauss1.h"

namespace d7x {

// th = (mu, lambda)

// R/invgauss2_distrib.R:303: statmod::rinvgauss(n, mean = mu,
// dispersion = 1 / lambda); the dispersion must be finite and positive
inline int invgauss2_rng_ok(const double* th) {
    return invgauss_params_ok(th[0], 1 / th[1]) ? 0 : 1;
}
inline void invgauss2_rng(const double* th, int n, double* out) {
    invgauss_rng_core(th[0], 1 / th[1], n, out);
}

// R/invgauss2_distrib.R:144: statmod::dinvgauss(y, mean = mu,
// dispersion = 1 / lambda, log = FALSE)
inline void invgauss2_pdf(const double* y, int n, const double* th, int vb,
                          double* out) {
    (void) vb;
    const double phi = 1 / th[1];
    for (int i = 0; i < n; ++i) out[i] = invgauss_pdf_one(y[i], th[0], phi);
}

// R/invgauss2_distrib.R:259 -> invgauss_quantile() (R/invgauss1_distrib.R:297)
inline int invgauss2_quantile(const double* p, int np, const double* th,
                              int vb, double* out) {
    (void) vb;
    return invgauss_quantile_core(p, np, th[0], 1 / th[1], out);
}

}  // namespace d7x

#endif
