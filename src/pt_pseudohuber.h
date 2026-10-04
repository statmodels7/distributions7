#ifndef D7_PT_PSEUDOHUBER_H
#define D7_PT_PSEUDOHUBER_H

#include <Rcpp.h>
#include <cmath>
#include "pt_bessel_k.h"
#include "pt_loc_scale.h"

// The pseudo-Huber family (mu, sigma, nu): with res = y - mu and
// D = sqrt(nu + res^2/sigma^2),
//   l = -D - log 2 - log sigma - log(nu)/2 - log K_1(sqrt(nu)).
// One function per component and order for the score, the diagonal of the
// Hessian and the diagonal of the third derivatives, called by the vector
// kernels (pseudohuber.cpp, pseudohuber_hd.cpp), by the quadrature of the
// expected information and by the scalar registry. nu enters the
// normalizing constant only through the ratios of the scaled Bessel
// functions K_m(sqrt(nu)), which PhNu carries.

namespace d7 {

// the nu-only quantities: K_0 to K_3 scaled at sqrt(nu),
// r1 = K_1'/K_1 = -(K_0 + K_2)/(2 K_1), r2 = K_1''/K_1 = (3 K_1 + K_3)/(4 K_1),
// and the nu-only part of l_nu,nu
struct PhNu { double sv, k0, k1, k2, k3, r1, r2, nu_const; };

inline PhNu ph_nu(double v) {
    PhNu P;
    P.sv = std::sqrt(v);
    P.k0 = bessel_k_scaled(P.sv, 0.0);
    P.k1 = bessel_k_scaled(P.sv, 1.0);
    P.k2 = bessel_k_scaled(P.sv, 2.0);
    P.k3 = bessel_k_scaled(P.sv, 3.0);
    P.r1 = -(P.k0 + P.k2) / (2.0 * P.k1);
    P.r2 = (3.0 * P.k1 + P.k3) / (4.0 * P.k1);
    P.nu_const = 0.5 / (v * v) +
        0.25 * (P.r1 / (v * P.sv) + P.r1 * P.r1 / v - P.r2 / v);
    return P;
}

inline double pseudohuber_score_mu(double res, double s2, double D) {
    return res / (s2 * D);
}

inline double pseudohuber_score_sigma(double res2, double s, double s2,
                                      double D) {
    return (res2 / (s2 * D) - 1.0) / s;
}

inline double pseudohuber_score_nu(double v, double D, const PhNu& P) {
    return -0.5 * (1.0 / v + 1.0 / D + P.r1 / P.sv);
}

inline double pseudohuber_hess_mu_mu(double v, double s2, double D3) {
    return -v / (s2 * D3);
}

inline double pseudohuber_hess_sigma_sigma(double res2, double s2, double s4,
                                           double D, double D3) {
    return (s4 - 3.0 * s2 * res2 / D + res2 * res2 / D3) / (s4 * s2);
}

inline double pseudohuber_hess_nu_nu(double D3, const PhNu& P) {
    return 0.25 / D3 + P.nu_const;
}

// the third derivatives on the diagonal, in r = mu - y and
// S = nu sigma^2 + r^2, so that D^k = S^(k/2) / sigma^k
inline double pseudohuber_d3_mu_mu_mu(double v, double s, double r,
                                      double S) {
    double S12 = std::sqrt(S), S52 = S * S * S12;
    return 3.0 * v * s * r / S52;
}

inline double pseudohuber_d3_sigma_sigma_sigma(double s, double r2,
                                               double S) {
    double s2 = s * s, s3 = s2 * s, s4 = s2 * s2;
    double r4 = r2 * r2, r6 = r4 * r2;
    double S12 = std::sqrt(S), S32 = S * S12, S52 = S * S * S12;
    return -2.0 / s3 + 12.0 * r2 / (s4 * S12)
        - 9.0 * r4 / (s4 * S32) + 3.0 * r6 / (s4 * S52);
}

// k4 the scaled K_4(sqrt(nu))
inline double pseudohuber_d3_nu_nu_nu(double v, double s, double S,
                                      const PhNu& P, double k4) {
    double s2 = s * s, s4 = s2 * s2, s5 = s4 * s;
    double S12 = std::sqrt(S), S52 = S * S * S12;
    double sv = P.sv, v2 = v * v, v32 = v * sv;
    double k0 = P.k0, k1 = P.k1, k2 = P.k2, k3 = P.k3;
    double A = k0 + k2;
    return (
        -32.0 / (v2 * v)
        - 12.0 * s5 / S52
        + 6.0 * A / (v2 * sv * k1)
        - 3.0 * A * A / (v2 * k1 * k1)
        + A * A * A / (v32 * k1 * k1 * k1)
        - 3.0 * A * (3.0 * k1 + k3) / (2.0 * v32 * k1 * k1)
        + 2.0 * (3.0 + k3 / k1) / v2
        + (2.0 * (3.0 * k1 + k3) + sv * (3.0 * k0 + 4.0 * k2 + k4)) / (2.0 * v2 * k1)
    ) / 32.0;
}

// log f, as distrib_pdf() forms it
inline double pseudohuber_logpdf(double y, double mu, double sigma, double v,
                                 const PhNu& P) {
    double z = (y - mu) / sigma;
    double D = std::sqrt(v + z * z);
    double log_norm = std::log(2.0) + std::log(sigma) + 0.5 * std::log(v) +
        std::log(P.k1) - P.sv;
    return -D - log_norm;
}

// the routers of the scalar registry: th = (mu, sigma, nu)
inline void pseudohuber_score_curv(int k, double y, const double* th,
                                   double* out) {
    const double m = th[0], s = th[1], v = th[2];
    const double s2 = s * s;
    const double res = y - m, res2 = res * res;
    const double D = std::sqrt(v + res2 / s2), D3 = D * D * D;
    if (k == 0) {
        out[0] = pseudohuber_score_mu(res, s2, D);
        out[1] = pseudohuber_hess_mu_mu(v, s2, D3);
    } else if (k == 1) {
        out[0] = pseudohuber_score_sigma(res2, s, s2, D);
        out[1] = pseudohuber_hess_sigma_sigma(res2, s2, s2 * s2, D, D3);
    } else {
        const PhNu P = ph_nu(v);
        out[0] = pseudohuber_score_nu(v, D, P);
        out[1] = pseudohuber_hess_nu_nu(D3, P);
    }
}

// the standardized diagonal pair of parameter k at shape = (nu), by the rule
// of pt_loc_scale.h
inline void pseudohuber_quad_diag(int k, const double* shape, double* out,
                                  bool want_d) {
    const double v = shape[0];
    const PhNu P = ph_nu(v);
    const double k4 = (want_d && k == 2) ? bessel_k_scaled(P.sv, 4.0) : 0.0;
    const LocScaleRule R = loc_scale_rule();
    LocScaleSum S;
    for (int j = 0; j < R.n; ++j) {
        const double z = R.x[j];
        const double fw = std::exp(pseudohuber_logpdf(z, 0.0, 1.0, v, P)) * R.w[j];
        const double res = z, res2 = res * res;
        const double D = std::sqrt(v + res2 / 1.0), D3 = D * D * D;
        double H;
        if (k == 0) H = pseudohuber_hess_mu_mu(v, 1.0, D3);
        else if (k == 1) H = pseudohuber_hess_sigma_sigma(res2, 1.0, 1.0, D, D3);
        else H = pseudohuber_hess_nu_nu(D3, P);
        if (!want_d) { S.add0(fw, H); continue; }
        const double r = -z, r2 = r * r, Sq = v + r2;
        double g, T3;
        if (k == 0) {
            g = pseudohuber_score_mu(res, 1.0, D);
            T3 = pseudohuber_d3_mu_mu_mu(v, 1.0, r, Sq);
        } else if (k == 1) {
            g = pseudohuber_score_sigma(res2, 1.0, 1.0, D);
            T3 = pseudohuber_d3_sigma_sigma_sigma(1.0, r2, Sq);
        } else {
            g = pseudohuber_score_nu(v, D, P);
            T3 = pseudohuber_d3_nu_nu_nu(v, 1.0, Sq, P, k4);
        }
        S.add(fw, g, H, T3);
    }
    S.result(out);
}

inline void pseudohuber_info_dinfo(int k, double y, const double* th,
                                   double* out) {
    (void) y;
    double std2[2];
    loc_scale_cached(2, k, th + 2, 1, [&](double* v) {
        pseudohuber_quad_diag(k, th + 2, v, true);
    }, std2);
    loc_scale_unstandardize(k, th[1], std2, out);
}

} // namespace d7

#endif
