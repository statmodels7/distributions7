#ifndef D7_PT_WEIBULL1_H
#define D7_PT_WEIBULL1_H

#include <Rcpp.h>
#include <cmath>
#include "pt_constants.h"

// The Weibull in its scale mu and shape sigma, with z = y/mu, lz = log z
// and u = z^sigma = exp(sigma lz): one function per component and order
// for the quantities the scalar registry reads, called by the vector
// kernels (weibull1.cpp, dexpected_kernels.cpp) and by the registry
// (d7_ccallable.cpp). The caller passes lz and u.

namespace d7 {

inline double weibull1_score_mu(double u, double m, double s) {
    return s * (u - 1.0) / m;
}

inline double weibull1_score_sigma(double lz, double u, double s) {
    return 1.0 / s + (1.0 - u) * lz;
}

inline double weibull1_hess_mu_mu(double u, double m, double s) {
    return s * (1.0 - (1.0 + s) * u) / (m * m);
}

inline double weibull1_hess_sigma_sigma(double lz, double u, double s) {
    return -1.0 / (s * s) - u * lz * lz;
}

// E_mm = -s^2/m^2, E_ss = -((1 - gamma)^2 + pi^2/6)/s^2, E_ms = (1 - gamma)/m
inline double weibull1_expected_mu_mu(double m, double s) {
    return -(s * s) / (m * m);
}

inline double weibull1_expected_sigma_sigma(double s) {
    return -kGumbelInfo / (s * s);
}

inline double weibull1_dexpected_mu_mu_mu(double m, double s) {
    double im = 1.0 / m, im2 = im * im, im3 = im2 * im;
    return 2.0 * s * s * im3;
}

inline double weibull1_dexpected_sigma_sigma_sigma(double s) {
    double is = 1.0 / s, is3 = is * is * is;
    return 2.0 * kGumbelInfo * is3;
}

inline void weibull1_score_curv(int k, double y, const double* th,
                                double* out) {
    double m = th[0], s = th[1];
    double lz = std::log(y / m), u = std::exp(s * lz);
    if (k == 0) {
        out[0] = weibull1_score_mu(u, m, s);
        out[1] = weibull1_hess_mu_mu(u, m, s);
    } else {
        out[0] = weibull1_score_sigma(lz, u, s);
        out[1] = weibull1_hess_sigma_sigma(lz, u, s);
    }
}

inline void weibull1_info_dinfo(int k, double y, const double* th,
                                double* out) {
    double m = th[0], s = th[1];
    if (k == 0) {
        out[0] = weibull1_expected_mu_mu(m, s);
        out[1] = weibull1_dexpected_mu_mu_mu(m, s);
    } else {
        out[0] = weibull1_expected_sigma_sigma(s);
        out[1] = weibull1_dexpected_sigma_sigma_sigma(s);
    }
}

} // namespace d7

#endif
