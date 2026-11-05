// The generator and the density of the Poisson-inverse Gaussian family in
// (mu, sigma), a discrete family: its quantile is never called and returns
// nonzero. The pdf takes the counts as doubles.
//
// Sources: R/pig1_distrib.R:443 (rng), :142 (pdf, exp(pig1_pdf_cpp(...)));
// src/pig_hd.cpp:226 (pig1_pdf_cpp, which calls pig1_row() of
// pt_pig_rows.h with order 0, only = true); statmod::rinvgauss (see
// pt_expect_invgauss1.h) and stats::rpois (Rmath rpois per element).

#ifndef D7_PT_EXPECT_PIG1_H
#define D7_PT_EXPECT_PIG1_H

#include "pt_expect_invgauss1.h"
#include "pt_pig_rows.h"

namespace d7x {

// shared by pig1 and pig2: lam <- statmod::rinvgauss(n, mean = mu, shape =
// shape) (dispersion = 1/shape, as rinvgauss computes it), then
// stats::rpois(n, lam): n normals, n uniforms, then n x rpois(lam_i)
inline void pig_rng_core(double mu, double shape, int n, double* out) {
    const double phi = 1 / shape;
    std::vector<double> lam(n);
    invgauss_rng_core(mu, phi, n, lam.data());
    for (int i = 0; i < n; ++i) out[i] = R::rpois(lam[i]);
}

// th = (mu, sigma); R: mean = rep_len(mu), shape = mu / sigma
inline int pig1_rng_ok(const double* th) {
    return invgauss_params_ok(th[0], 1 / (th[0] / th[1])) ? 0 : 1;
}
inline void pig1_rng(const double* th, int n, double* out) {
    pig_rng_core(th[0], th[0] / th[1], n, out);
}

// R/pig1_distrib.R:142: exp(pig1_pdf_cpp(y, mu, sigma, threads))
inline void pig1_pdf(const double* y, int n, const double* th, int vb,
                     double* out) {
    (void) vb;
    for (int i = 0; i < n; ++i) {
        double v[1];
        d7pigrow::pig1_row(y[i], th[0], th[1], 0, true, v);
        out[i] = r_exp(v[0]);
    }
}

// discrete: never called
inline int pig1_quantile(const double* p, int np, const double* th, int vb,
                         double* out) {
    (void) p; (void) np; (void) th; (void) vb; (void) out;
    return 1;
}

}  // namespace d7x

#endif
