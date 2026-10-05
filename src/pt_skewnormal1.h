#ifndef D7_PT_SKEWNORMAL1_H
#define D7_PT_SKEWNORMAL1_H

#include <Rcpp.h>
#include <cmath>
#include "pt_mills.h"
#include "pt_loc_scale.h"

// The skew normal in its direct parametrization (mu, sigma, alpha),
// l = log 2 - log sigma + log phi(z) + log Phi(alpha z), z = (y - mu)/sigma:
// one function per component and order for the score, the diagonal of the
// Hessian and the diagonal of the third derivatives, called by the vector
// kernels (skewnormal1.cpp, skewnormal_hd.cpp), by the quadrature of the
// expected information and by the scalar registry. The caller
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

// the third derivatives on the diagonal; with t = alpha z,
// g3 = -dr (t + r) - r (1 + dr) is the third derivative of log Phi at t, and
// h1, h2, h3 those of h(z) = -z^2/2 + log Phi(alpha z) in z
inline double skewnormal1_g3(double t, double r, double dr) {
    return -dr * (t + r) - r * (1.0 + dr);
}

inline double skewnormal1_d3_mu_mu_mu(double s, double a, double g3) {
    double h3 = a * a * a * g3;
    double s3 = s * s * s;
    return -h3 / s3;
}

inline double skewnormal1_d3_sigma_sigma_sigma(double z, double s, double a,
                                               double r, double dr,
                                               double g3) {
    double h1 = -z + a * r;
    double h2 = -1.0 + a * a * dr;
    double h3 = a * a * a * g3;
    double z2 = z * z, z3 = z2 * z, s3 = s * s * s;
    return -(2.0 + 6.0 * z * h1 + 6.0 * z2 * h2 + z3 * h3) / s3;
}

inline double skewnormal1_d3_alpha_alpha_alpha(double z, double g3) {
    double z2 = z * z, z3 = z2 * z;
    return z3 * g3;
}

// log f, as distrib_pdf() forms it
inline double skewnormal1_logpdf(double y, double mu, double sigma,
                                 double alpha) {
    double z = (y - mu) / sigma;
    return std::log(2.0) - std::log(sigma) + R::dnorm4(z, 0.0, 1.0, 1) +
        R::pnorm5(alpha * z, 0.0, 1.0, 1, 1);
}

// the routers of the scalar registry: th = (mu, sigma, alpha)
inline void skewnormal1_score_curv(int k, double y, const double* th,
                                   double* out) {
    const double m = th[0], s = th[1], a = th[2];
    const double z = (y - m) / s, t = a * z;
    const double r = mills_ratio(t), dr = -r * (t + r);
    if (k == 0) {
        out[0] = skewnormal1_score_mu(z, s, a, r);
        out[1] = skewnormal1_hess_mu_mu(s, a, dr);
    } else if (k == 1) {
        out[0] = skewnormal1_score_sigma(z, s, a, r);
        out[1] = skewnormal1_hess_sigma_sigma(z, s, a, r, dr);
    } else {
        out[0] = skewnormal1_score_alpha(z, r);
        out[1] = skewnormal1_hess_alpha_alpha(z, dr);
    }
}

// the standardized diagonal pair of parameter k at shape = (alpha), by the
// rule of pt_loc_scale.h
inline void skewnormal1_quad_diag(int k, const double* shape, double* out,
                                  bool want_d) {
    const double a = shape[0];
    const LocScaleRule R = loc_scale_rule();
    LocScaleSum S;
    for (int j = 0; j < R.n; ++j) {
        const double z = R.x[j], t = a * z;
        const double fw = std::exp(skewnormal1_logpdf(z, 0.0, 1.0, a)) * R.w[j];
        const double r = mills_ratio(t), dr = -r * (t + r);
        double H;
        if (k == 0) H = skewnormal1_hess_mu_mu(1.0, a, dr);
        else if (k == 1) H = skewnormal1_hess_sigma_sigma(z, 1.0, a, r, dr);
        else H = skewnormal1_hess_alpha_alpha(z, dr);
        if (!want_d) { S.add0(fw, H); continue; }
        const double g3 = skewnormal1_g3(t, r, dr);
        double g, T3;
        if (k == 0) {
            g = skewnormal1_score_mu(z, 1.0, a, r);
            T3 = skewnormal1_d3_mu_mu_mu(1.0, a, g3);
        } else if (k == 1) {
            g = skewnormal1_score_sigma(z, 1.0, a, r);
            T3 = skewnormal1_d3_sigma_sigma_sigma(z, 1.0, a, r, dr, g3);
        } else {
            g = skewnormal1_score_alpha(z, r);
            T3 = skewnormal1_d3_alpha_alpha_alpha(z, g3);
        }
        S.add(fw, g, H, T3);
    }
    S.result(out);
}

inline void skewnormal1_info_dinfo(int k, double y, const double* th,
                                   double* out) {
    (void) y;
    double std2[2];
    loc_scale_cached(0, k, th + 2, 1, [&](double* v) {
        skewnormal1_quad_diag(k, th + 2, v, true);
    }, std2);
    loc_scale_unstandardize(k, th[1], std2, out);
}

} // namespace d7

#endif
