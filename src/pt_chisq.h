#ifndef D7_PT_CHISQ_H
#define D7_PT_CHISQ_H

#include <Rcpp.h>
#include <cmath>

// Chi-squared in its mean (the degrees of freedom): one function per
// component and order for the quantities the scalar registry reads, called
// by the vector kernels (chisq.cpp, dexpected_kernels.cpp) and by the
// registry (d7_ccallable.cpp). The caller passes log y and the polygamma
// psi^(k)(mu/2) each function reads. From the second order on no
// derivative involves y, so the expected information is the observed one.

namespace d7 {

inline double chisq_score_mu(double log_y, double dg) {
    const double log2 = std::log(2.0);
    return 0.5 * (log_y - log2 - dg);
}

inline double chisq_hess_mu_mu(double tg) {
    return -0.25 * tg;
}

inline double chisq_expected_mu_mu(double tg) {
    return chisq_hess_mu_mu(tg);
}

inline double chisq_dexpected_mu_mu_mu(double pg2) {
    return -pg2 / 8.0;
}

inline void chisq_score_curv(int k, double y, const double* th,
                             double* out) {
    double h = 0.5 * th[0];
    out[0] = chisq_score_mu(std::log(y), R::digamma(h));
    out[1] = chisq_hess_mu_mu(R::trigamma(h));
}

inline void chisq_info_dinfo(int k, double y, const double* th,
                             double* out) {
    double h = 0.5 * th[0];
    out[0] = chisq_expected_mu_mu(R::trigamma(h));
    out[1] = chisq_dexpected_mu_mu_mu(R::psigamma(h, 2));
}

} // namespace d7

#endif
