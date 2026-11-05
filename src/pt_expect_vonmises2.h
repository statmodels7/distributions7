// vonmises2: generator, density and quantile for the compiled expected
// information by simulation and by integration (see pt_expect_common.h).
// th = (mu, rho). Every method of R/vonmises2_distrib.R is vonmises1's at
// kappa = numericals7::bessel_i_ratio_inverse(rho), which is
// bessel_ratio_inverse_cpp() -> n7_bessel_ratio_inverse() per element
// (numericals7 src/bessel_ratio.cpp), resolved here through R_GetCCallable.
// distrib_pdf: :146-152; distrib_rng: :189-193; distrib_cdf: :800-805; the
// quantile is the numerical fallback over those.

#ifndef D7_PT_EXPECT_VONMISES2_H
#define D7_PT_EXPECT_VONMISES2_H

#include "pt_expect_common.h"
#include "pt_expect_vonmises1.h"

namespace d7x {

inline double vm2_kappa_of_rho(double rho) {
    typedef double (*fn_t)(double);
    static fn_t p = nullptr;
    if (p == nullptr) p = (fn_t) R_GetCCallable("numericals7", "n7_bessel_ratio_inverse");
    return p(rho);
}

inline int vonmises2_rng_ok(const double* th) {
    (void) th;
    return 0;
}

inline void vonmises2_rng(const double* th, int n, double* out) {
    vm1_rng_kappa(th[0], vm2_kappa_of_rho(th[1]), n, out);
}

inline void vonmises2_pdf(const double* y, int n, const double* th, int vb,
                          double* out) {
    (void) vb;
    vm1_pdf_kappa(y, n, th[0], vm2_kappa_of_rho(th[1]), out);
}

inline int vonmises2_quantile(const double* p, int np, const double* th, int vb,
                              double* out) {
    (void) vb;
    Vm1Fam F(th[0], vm2_kappa_of_rho(th[1]));
    return qf_quantile(p, np, F, out);
}

}  // namespace d7x

#endif
