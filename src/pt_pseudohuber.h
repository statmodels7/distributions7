#ifndef D7_PT_PSEUDOHUBER_H
#define D7_PT_PSEUDOHUBER_H

#include <Rcpp.h>
#include <cmath>
#include "pt_bessel_k.h"
#include "pt_loc_scale.h"
#include "pt_sqrt.h"
#include "pt_ph_nu.h"

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
// the nu-only part of l_nu,nu, and dB[k] = B^(k)(nu), k = 1..4, with
// B = log(sqrt(nu) K_1(sqrt(nu))) from ph_nu_derivs() (pt_ph_nu.h): the
// log-density carries -B, and its pure nu derivatives are
// -c_k D^(1 - 2k) - B^(k), c_k the k-th derivative of nu^(1/2) times
// nu^(k - 1/2)
struct PhNu { double sv, k0, k1, k2, k3, r1, r2, nu_const; double dB[6]; };

inline PhNu ph_nu(double v) {
    PhNu P;
    P.sv = d7::sqrt_cr(v);
    P.k0 = bessel_k_scaled(P.sv, 0.0);
    P.k1 = bessel_k_scaled(P.sv, 1.0);
    P.k2 = bessel_k_scaled(P.sv, 2.0);
    P.k3 = bessel_k_scaled(P.sv, 3.0);
    P.r1 = -(P.k0 + P.k2) / (2.0 * P.k1);
    P.r2 = (3.0 * P.k1 + P.k3) / (4.0 * P.k1);
    double dlR[6], R0;
    ph_nu_derivs(v, 4, P.dB, dlR, R0);
    P.nu_const = -P.dB[2];
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
    (void) v;
    return -0.5 * s / P.A - N.dB[1];
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

// d^3 l/dnu^3 = -(3/8) D^-5 - B'''(nu), with q = s/A = 1/D; k4 is no
// longer read
inline double pseudohuber_d3_nu_nu_nu(double v, double s, const PhA& PA,
                                      const PhNu& P, double k4) {
    (void) v; (void) k4;
    const double q = s / PA.A, q2 = q * q;
    return -0.375 * q2 * q2 * q - P.dB[3];
}

// the off-diagonal components of the Hessian, as pseudohuber_hessian_cpp
// forms them (s2 = sigma^2)
inline double pseudohuber_hess_mu_sigma(double s, const PhA& P) {
    const double s2 = s * s;
    return P.rho * (2.0 * P.kap2 + P.rho * P.rho) / s2;
}

inline double pseudohuber_hess_mu_nu(double s, const PhA& P) {
    return P.rho * (s / P.A) / (2.0 * P.A);
}

inline double pseudohuber_hess_sigma_nu(const PhA& P) {
    return -P.rho * P.rho / (2.0 * P.A);
}

// the off-diagonal third derivatives, as pseudohuber_deriv3_cpp forms them
// (s2 = sigma^2, s3 = sigma^3, q = sigma/A)
inline double pseudohuber_d3_mu_mu_sigma(double s, const PhA& P) {
    const double s2 = s * s;
    const double p2 = P.rho * P.rho, k = P.kap2;
    return k * (2.0 * k - p2) / (s2 * P.A);
}

inline double pseudohuber_d3_mu_mu_nu(double s, const PhA& P) {
    const double p2 = P.rho * P.rho, k = P.kap2, q = s / P.A;
    return q * (k - 2.0 * p2) / (2.0 * P.A * P.A);
}

inline double pseudohuber_d3_mu_sigma_sigma(double s, const PhA& P) {
    const double s2 = s * s, s3 = s2 * s;
    const double p = P.rho, p2 = p * p, k = P.kap2;
    return -p * (6.0 - 7.0 * k * p2 - 4.0 * p2 * p2) / s3;
}

inline double pseudohuber_d3_mu_sigma_nu(const PhA& P) {
    const double p = P.rho, p2 = p * p, k = P.kap2;
    return p * (p2 - 2.0 * k) / (2.0 * P.A * P.A);
}

inline double pseudohuber_d3_mu_nu_nu(double s, const PhA& P) {
    const double p = P.rho, q = s / P.A;
    return -3.0 * q * q * q * p / (4.0 * P.A);
}

inline double pseudohuber_d3_sigma_sigma_nu(double s, const PhA& P) {
    const double p2 = P.rho * P.rho, k = P.kap2;
    return 3.0 * k * p2 / (2.0 * s * P.A);
}

inline double pseudohuber_d3_sigma_nu_nu(double s, const PhA& P) {
    const double p2 = P.rho * P.rho, q = s / P.A;
    return 3.0 * q * q * p2 / (4.0 * P.A);
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

// The routers of the joint registry for the score and the Hessian: every
// component at one observation, on the parameter scale, th = (mu, sigma,
// nu), P = 3. The Hessian fills a full 3 x 3 matrix, column-major,
// out[a + 3 b]. The expected information is an exp-sinh quadrature; its
// routers, joint_info and joint_dinfo, follow pseudohuber_quad_full() below.
inline void pseudohuber_joint_grad(double y, const double* th, double* out) {
    const double m = th[0], s = th[1], v = th[2];
    const double r = m - y;
    const PhA A = ph_a(r, s, v);
    const PhNu P = ph_nu(v);
    out[0] = pseudohuber_score_mu(s, A);
    out[1] = pseudohuber_score_sigma(r, s, A);
    out[2] = pseudohuber_score_nu(v, s, A, P);
}

inline void pseudohuber_joint_hess(double y, const double* th, double* out) {
    const double m = th[0], s = th[1], v = th[2];
    const double r = m - y;
    const PhA A = ph_a(r, s, v);
    const PhNu P = ph_nu(v);
    out[0] = pseudohuber_hess_mu_mu(s, A);
    out[4] = pseudohuber_hess_sigma_sigma(r, s, A);
    out[8] = pseudohuber_hess_nu_nu(s, A, P);
    out[1] = out[3] = pseudohuber_hess_mu_sigma(s, A);
    out[2] = out[6] = pseudohuber_hess_mu_nu(s, A);
    out[5] = out[7] = pseudohuber_hess_sigma_nu(A);
}

// sets T3 (3 x 3 x 3, index a + 3 b + 9 c) at the six permutations of (a, b, c)
inline void pseudohuber_t3_set(double* T3, int a, int b, int c, double x) {
    T3[a + 3 * b + 9 * c] = T3[a + 3 * c + 9 * b] = T3[b + 3 * a + 9 * c] =
        T3[b + 3 * c + 9 * a] = T3[c + 3 * a + 9 * b] =
        T3[c + 3 * b + 9 * a] = x;
}

// The full standardized pair of the expected information at shape = (nu), by
// the rule of pt_loc_scale.h: at each node the score g, the Hessian H and the
// third derivative T3 of all three parameters (index 0 mu, 1 sigma, 2 nu), at
// location 0 and scale 1, accumulated by LocScaleFull. The pure entries are
// the functions pseudohuber_quad_diag() calls, the others those of
// pseudohuber_hessian_cpp and pseudohuber_deriv3_cpp. I is 3 x 3 and D is
// 3 x 3 x 3, column-major, both triangles filled.
inline void pseudohuber_quad_full(const double* shape, double* I, double* D) {
    const double v = shape[0];
    const PhNu P = ph_nu(v);
    const double k4 = bessel_k_scaled(P.sv, 4.0);
    const LocScaleRule R = loc_scale_rule();
    LocScaleFull S(3);
    for (int j = 0; j < R.n; ++j) {
        const double z = R.x[j];
        const double fw = std::exp(pseudohuber_logpdf(z, 0.0, 1.0, v, P)) * R.w[j];
        if (fw == 0.0) continue;
        const double r = -z;
        const PhA A = ph_a(r, 1.0, v);
        double g[3], H[9], T3[27];
        g[0] = pseudohuber_score_mu(1.0, A);
        g[1] = pseudohuber_score_sigma(r, 1.0, A);
        g[2] = pseudohuber_score_nu(v, 1.0, A, P);
        H[0] = pseudohuber_hess_mu_mu(1.0, A);
        H[4] = pseudohuber_hess_sigma_sigma(r, 1.0, A);
        H[8] = pseudohuber_hess_nu_nu(1.0, A, P);
        H[1] = H[3] = pseudohuber_hess_mu_sigma(1.0, A);
        H[2] = H[6] = pseudohuber_hess_mu_nu(1.0, A);
        H[5] = H[7] = pseudohuber_hess_sigma_nu(A);
        pseudohuber_t3_set(T3, 0, 0, 0, pseudohuber_d3_mu_mu_mu(1.0, A));
        pseudohuber_t3_set(T3, 0, 0, 1, pseudohuber_d3_mu_mu_sigma(1.0, A));
        pseudohuber_t3_set(T3, 0, 0, 2, pseudohuber_d3_mu_mu_nu(1.0, A));
        pseudohuber_t3_set(T3, 0, 1, 1, pseudohuber_d3_mu_sigma_sigma(1.0, A));
        pseudohuber_t3_set(T3, 0, 1, 2, pseudohuber_d3_mu_sigma_nu(A));
        pseudohuber_t3_set(T3, 0, 2, 2, pseudohuber_d3_mu_nu_nu(1.0, A));
        pseudohuber_t3_set(T3, 1, 1, 1,
                           pseudohuber_d3_sigma_sigma_sigma(r, 1.0, A));
        pseudohuber_t3_set(T3, 1, 1, 2, pseudohuber_d3_sigma_sigma_nu(1.0, A));
        pseudohuber_t3_set(T3, 1, 2, 2, pseudohuber_d3_sigma_nu_nu(1.0, A));
        pseudohuber_t3_set(T3, 2, 2, 2,
                           pseudohuber_d3_nu_nu_nu(v, 1.0, A, P, k4));
        S.add(fw, g, H, T3);
    }
    S.result(I, D);
}

// The routers of the joint registry for the expected Hessian E[l_ab] (3 x 3,
// column-major, out[a + 3 b]) and its derivative in parameter c (3 x 3 x 3,
// out[a + 3 b + 9 c]), both triangles filled, th = (mu, sigma, nu). The
// standardized pair is computed once per distinct nu on each thread and
// scaled by sigma.
inline void pseudohuber_joint_info(double y, const double* th, double* out) {
    (void) y;
    const double* sh = th + 2;
    const LocScaleFullSlot& s = loc_scale_full_cached(2, sh, 1,
        [&](double* I, double* D) { pseudohuber_quad_full(sh, I, D); });
    loc_scale_full_info(s, 3, th[1], out);
}

inline void pseudohuber_joint_dinfo(double y, const double* th, double* out) {
    (void) y;
    const double* sh = th + 2;
    const LocScaleFullSlot& s = loc_scale_full_cached(2, sh, 1,
        [&](double* I, double* D) { pseudohuber_quad_full(sh, I, D); });
    loc_scale_full_dinfo(s, 3, th[1], out);
}

} // namespace d7

#endif
