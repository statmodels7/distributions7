#ifndef D7_PT_SKEWNORMAL1_H
#define D7_PT_SKEWNORMAL1_H

#include <Rcpp.h>
#include <cmath>
#include "pt_mills.h"

// The skew normal in its direct parametrization (mu, sigma, alpha),
// l = log 2 - log sigma + log phi(z) + log Phi(alpha z), z = (y - mu)/sigma:
// one function per component and order for the score and the diagonal of
// the Hessian, called by the vector kernels (skewnormal1.cpp). The caller
// passes z, the Mills ratio r at t = alpha z (mills_ratio) and its
// derivative dr = -r (t + r).

namespace d7 {

inline double skewnormal1_score_mu(double z, double s, double a, double r) {
    return (z - a * r) / s;
}

inline double skewnormal1_score_sigma(double z, double s, double a, double r) {
    return (z * z - 1.0 - a * z * r) / s;
}

inline double skewnormal1_score_alpha(double z, double r) {
    return z * r;
}

inline double skewnormal1_hess_mu_mu(double s, double a, double dr) {
    return (a * a * dr - 1.0) / (s * s);
}

inline double skewnormal1_hess_sigma_sigma(double z, double s, double a,
                                           double r, double dr) {
    return (1.0 - 3.0 * z * z + 2.0 * a * z * r + a * a * z * z * dr) / (s * s);
}

inline double skewnormal1_hess_alpha_alpha(double z, double dr) {
    return z * z * dr;
}

} // namespace d7

#endif
