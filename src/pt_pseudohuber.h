#ifndef D7_PT_PSEUDOHUBER_H
#define D7_PT_PSEUDOHUBER_H

#include <Rcpp.h>
#include <cmath>
#include "pt_bessel_k.h"
#include "pt_loc_scale.h"
#include "pt_sqrt.h"

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
    P.sv = d7::sqrt_cr(v);
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

// With r = mu - y, A = sqrt(nu sigma^2 + r^2) = sigma D is formed by hypot()
// so that r^2 never overflows, and every component below is written in
// rho = r/A and kap2 = nu sigma^2/A^2, both bounded (rho^2 + kap2 = 1),
// times a power of 1/A or a single factor r: written in r^k / S^(m/2) the
// fourth derivatives were NaN from |y| of about 1e40.
struct PhA { double A, rho, kap2; };

inline PhA ph_a(double r, double s, double v) {
    PhA P;
    const double c = s * d7::sqrt_cr(v);
    P.A = std::hypot(r, c);
    P.rho = r / P.A;
    const double k = c / P.A;
    P.kap2 = k * k;
    return P;
}

inline double pseudohuber_score_mu(double s, const PhA& P) {
    return -P.rho / s;
}

inline double pseudohuber_score_sigma(double r, double s, const PhA& P) {
    return (r * P.rho / s - 1.0) / s;
}

inline double pseudohuber_score_nu(double v, double s, const PhA& P,
                                   const PhNu& N) {
    return -0.5 * (1.0 / v + s / P.A + N.r1 / N.sv);
}

inline double pseudohuber_hess_mu_mu(double s, const PhA& P) {
    return -P.kap2 / (s * P.A);
}

inline double pseudohuber_hess_sigma_sigma(double r, double s, const PhA& P) {
    return (1.0 + (r * P.rho / s) * (P.rho * P.rho - 3.0)) / (s * s);
}

inline double pseudohuber_hess_nu_nu(double s, const PhA& P, const PhNu& N) {
    const double q = s / P.A;
    return 0.25 * q * q * q + N.nu_const;
}

// the third derivatives on the diagonal
inline double pseudohuber_d3_mu_mu_mu(double s, const PhA& P) {
    return 3.0 * P.kap2 * P.rho / (s * P.A * P.A);
}

inline double pseudohuber_d3_sigma_sigma_sigma(double r, double s,
                                               const PhA& P) {
    const double p = P.rho, p2 = p * p;
    return -2.0 / (s * s * s) +
        (r / (s * s * s * s)) * p * (12.0 - 9.0 * p2 + 3.0 * p2 * p2);
}

// k4 the scaled K_4(sqrt(nu))
inline double pseudohuber_d3_nu_nu_nu(double v, double s, const PhA& PA,
                                      const PhNu& P, double k4) {
    const double q = s / PA.A, q2 = q * q;
    double sv = P.sv, v2 = v * v, v32 = v * sv;
    double k0 = P.k0, k1 = P.k1, k2 = P.k2, k3 = P.k3;
    double A = k0 + k2;
    return (
        -32.0 / (v2 * v)
        - 12.0 * q2 * q2 * q
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
    double D = d7::sqrt_cr(v + z * z);
    double log_norm = std::log(2.0) + std::log(sigma) + 0.5 * std::log(v) +
        std::log(P.k1) - P.sv;
    return -D - log_norm;
}

// the routers of the scalar registry: th = (mu, sigma, nu)
inline void pseudohuber_score_curv(int k, double y, const double* th,
                                   double* out) {
    const double m = th[0], s = th[1], v = th[2];
    const double r = m - y;
    const PhA A = ph_a(r, s, v);
    if (k == 0) {
        out[0] = pseudohuber_score_mu(s, A);
        out[1] = pseudohuber_hess_mu_mu(s, A);
    } else if (k == 1) {
        out[0] = pseudohuber_score_sigma(r, s, A);
        out[1] = pseudohuber_hess_sigma_sigma(r, s, A);
    } else {
        const PhNu P = ph_nu(v);
        out[0] = pseudohuber_score_nu(v, s, A, P);
        out[1] = pseudohuber_hess_nu_nu(s, A, P);
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
        const double r = -z;
        const PhA A = ph_a(r, 1.0, v);
        double H;
        if (k == 0) H = pseudohuber_hess_mu_mu(1.0, A);
        else if (k == 1) H = pseudohuber_hess_sigma_sigma(r, 1.0, A);
        else H = pseudohuber_hess_nu_nu(1.0, A, P);
        if (!want_d) { S.add0(fw, H); continue; }
        double g, T3;
        if (k == 0) {
            g = pseudohuber_score_mu(1.0, A);
            T3 = pseudohuber_d3_mu_mu_mu(1.0, A);
        } else if (k == 1) {
            g = pseudohuber_score_sigma(r, 1.0, A);
            T3 = pseudohuber_d3_sigma_sigma_sigma(r, 1.0, A);
        } else {
            g = pseudohuber_score_nu(v, 1.0, A, P);
            T3 = pseudohuber_d3_nu_nu_nu(v, 1.0, A, P, k4);
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
