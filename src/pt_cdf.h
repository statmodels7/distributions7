#ifndef D7_PT_CDF_H
#define D7_PT_CDF_H

#include <Rcpp.h>
#include <cmath>

// The distribution functions of the continuous families and their first
// and second derivatives in the parameters, for truncated(): one function
// per family and order, called by the scalar registry and, through
// d7_cdf_*_cpp(), by trunc_route_parts() in R, so that the retained mass
// and its derivatives are the same numbers on both routes. The derivatives
// are those of F, on the parameter scale; cdf_grad writes one value per
// parameter and cdf_hess one per pair in the order of hess_names(): the
// diagonal first, then (i, j), i < j, row by row.
//
// A family returns false when it has no compiled distribution function.

namespace d7 {

// ---- Gaussian1: (mu, sigma) ------------------------------------------------

inline double gaussian1_cdf(double q, const double* th, bool lower) {
    return R::pnorm(q, th[0], th[1], lower, 0);
}

// f and d log f / dy at q, as dnorm() and distrib_grad_y() write them
inline void gaussian1_f_ly(double q, const double* th, double* f, double* ly) {
    *f = R::dnorm(q, th[0], th[1], 0);
    *ly = -(q - th[0]) / (th[1] * th[1]);
}

// the location-scale derivatives of F from f and l_y at z = (q - mu)/s:
// F_mu = -f, F_s = -z f, F_mumu = f l_y, F_ss = f (z^2 l_y + 2 z / s),
// F_mus = f (z l_y + 1/s)
inline void loc_scale_cdf_grad(double z, double f, double* g) {
    g[0] = -f;
    g[1] = -z * f;
}
inline void loc_scale_cdf_hess(double z, double s, double f, double ly,
                               double* h) {
    h[0] = f * ly;
    h[1] = f * (z * z * ly + 2 * z / s);
    h[2] = f * (z * ly + 1 / s);
}

inline void gaussian1_cdf_grad(double q, const double* th, double* g) {
    double f, ly;
    gaussian1_f_ly(q, th, &f, &ly);
    loc_scale_cdf_grad((q - th[0]) / th[1], f, g);
}

inline void gaussian1_cdf_hess(double q, const double* th, double* h) {
    double f, ly;
    gaussian1_f_ly(q, th, &f, &ly);
    loc_scale_cdf_hess((q - th[0]) / th[1], th[1], f, ly, h);
}

} // namespace d7

#endif
