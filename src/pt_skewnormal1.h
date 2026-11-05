#ifndef D7_PT_SKEWNORMAL1_H
#define D7_PT_SKEWNORMAL1_H

#include <Rcpp.h>
#include <cmath>
#include "pt_mills.h"
#include "pt_loc_scale.h"
#include "pt_logphi.h"
#include "pt_parts.h"

// The skew normal in its direct parametrization (mu, sigma, alpha),
// l = log 2 - log sigma + log phi(z) + log Phi(alpha z), z = (y - mu)/sigma:
// one function per component and order for the score, the diagonal of the
// Hessian and the diagonal of the third derivatives, called by the vector
// kernels (skewnormal1.cpp, skewnormal_hd.cpp), by the quadrature of the
// expected information and by the scalar registry. The caller
// passes z, the Mills ratio r = phi/Phi at t = alpha z and its derivative
// dr, both from logphi_derivs() (pt_logphi.h), which forms them without the
// cancellation of dr = -r (t + r) far below t = 0.

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
    double gg[6];
    logphi_derivs(t, 2, gg);
    const double r = gg[1], dr = gg[2];
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
        double gg[6];
        logphi_derivs(t, want_d ? 3 : 2, gg);
        const double r = gg[1], dr = gg[2];
        double H;
        if (k == 0) H = skewnormal1_hess_mu_mu(1.0, a, dr);
        else if (k == 1) H = skewnormal1_hess_sigma_sigma(z, 1.0, a, r, dr);
        else H = skewnormal1_hess_alpha_alpha(z, dr);
        if (!want_d) { S.add0(fw, H); continue; }
        const double g3 = gg[3];
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

// The off-diagonal components of the Hessian, from skewnormal1_hessian_cpp
// (skewnormal1.cpp): the caller passes z, the Mills ratio r and dr.
inline double skewnormal1_hess_mu_sigma(double z, double s, double a,
                                        double r, double dr) {
    const double s2 = s * s;
    return (a * a * z * dr - 2.0 * z + a * r) / s2;
}

inline double skewnormal1_hess_mu_alpha(double z, double s, double a,
                                        double r, double dr) {
    return -(r + a * z * dr) / s;
}

inline double skewnormal1_hess_sigma_alpha(double z, double s, double a,
                                           double r, double dr) {
    return -(z * r + a * z * z * dr) / s;
}

// The routers of the joint registry, in (mu, sigma, alpha). The skew normal's
// expected information has no closed form (it is a quadrature in the shape,
// see skewnormal1_info_dinfo()), so only the gradient and the Hessian are
// written here.
inline void skewnormal1_joint_grad(double y, const double* th, double* out) {
    const double m = th[0], s = th[1], a = th[2];
    const double z = (y - m) / s;
    double gg[6];
    logphi_derivs(a * z, 1, gg);
    const double r = gg[1];
    out[0] = skewnormal1_score_mu(z, s, a, r);
    out[1] = skewnormal1_score_sigma(z, s, a, r);
    out[2] = skewnormal1_score_alpha(z, r);
}

inline void skewnormal1_joint_hess(double y, const double* th, double* out) {
    const double m = th[0], s = th[1], a = th[2];
    const double z = (y - m) / s;
    double gg[6];
    logphi_derivs(a * z, 2, gg);
    const double r = gg[1], dr = gg[2];
    out[0] = skewnormal1_hess_mu_mu(s, a, dr);
    out[4] = skewnormal1_hess_sigma_sigma(z, s, a, r, dr);
    out[8] = skewnormal1_hess_alpha_alpha(z, dr);
    out[1] = out[3] = skewnormal1_hess_mu_sigma(z, s, a, r, dr);
    out[2] = out[6] = skewnormal1_hess_mu_alpha(z, s, a, r, dr);
    out[5] = out[7] = skewnormal1_hess_sigma_alpha(z, s, a, r, dr);
}

// The third derivatives with a mixed index, from skewnormal_deriv3_cpp
// (skewnormal_hd.cpp). With t = alpha z, h1, h2, h3 as above and g1, g2, g3
// the derivatives of log Phi at t (g1 = r, g2 = dr).
inline double skewnormal1_d3_mu_mu_sigma(double z, double s, double h2,
                                         double h3) {
    return -(z * h3 + 2.0 * h2) / (s * s * s);
}

inline double skewnormal1_d3_mu_mu_alpha(double t, double s, double a,
                                         double g2, double g3) {
    return a * (2.0 * g2 + t * g3) / (s * s);
}

inline double skewnormal1_d3_mu_sigma_sigma(double z, double s, double h1,
                                            double h2, double h3) {
    return -(z * z * h3 + 4.0 * z * h2 + 2.0 * h1) / (s * s * s);
}

inline double skewnormal1_d3_mu_sigma_alpha(double t, double s, double g1,
                                            double g2, double g3) {
    return (g1 + 3.0 * t * g2 + t * t * g3) / (s * s);
}

inline double skewnormal1_d3_mu_alpha_alpha(double z, double t, double s,
                                            double g2, double g3) {
    return -z * (2.0 * g2 + t * g3) / s;
}

inline double skewnormal1_d3_sigma_sigma_alpha(double z, double t, double s,
                                               double g1, double g2,
                                               double g3) {
    return z * (2.0 * g1 + 4.0 * t * g2 + t * t * g3) / (s * s);
}

inline double skewnormal1_d3_sigma_alpha_alpha(double z, double t, double s,
                                               double g2, double g3) {
    return -z * z * (2.0 * g2 + t * g3) / s;
}

// the ten distinct third derivatives, in the order (mmm, mms, mma, mss, msa,
// maa, sss, ssa, saa, aaa), into the P = 3 array T[a + 3 b + 9 c] with every
// permutation of the indices filled
inline void skewnormal1_fill_t3(const double* v, double* T) {
    static const int idx[10][3] = {{0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {0, 1, 1},
                                   {0, 1, 2}, {0, 2, 2}, {1, 1, 1}, {1, 1, 2},
                                   {1, 2, 2}, {2, 2, 2}};
    static const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                                   {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (int q = 0; q < 10; ++q)
        for (int p = 0; p < 6; ++p)
            T[idx[q][perm[p][0]] + 3 * idx[q][perm[p][1]] +
              9 * idx[q][perm[p][2]]] = v[q];
}

// The full standardized pair at shape = (alpha), by the rule of
// pt_loc_scale.h: at each node the score g, the Hessian H and the third
// derivatives T3 in (mu, sigma, alpha) at location 0 and scale 1, added with
// the density weight. I is 3 x 3 and D is 3 x 3 x 3, column-major.
inline void skewnormal1_quad_full(const double* shape, double* I, double* D) {
    const double a = shape[0];
    const LocScaleRule R = loc_scale_rule();
    LocScaleFull S(3);
    for (int j = 0; j < R.n; ++j) {
        const double z = R.x[j], t = a * z;
        const double fw = std::exp(skewnormal1_logpdf(z, 0.0, 1.0, a)) * R.w[j];
        if (fw == 0.0) continue;
        double gg[6];
        logphi_derivs(t, 3, gg);
        const double r = gg[1], dr = gg[2], g3 = gg[3];
        double g[3], H[9], T3[27];
        g[0] = skewnormal1_score_mu(z, 1.0, a, r);
        g[1] = skewnormal1_score_sigma(z, 1.0, a, r);
        g[2] = skewnormal1_score_alpha(z, r);
        H[0] = skewnormal1_hess_mu_mu(1.0, a, dr);
        H[4] = skewnormal1_hess_sigma_sigma(z, 1.0, a, r, dr);
        H[8] = skewnormal1_hess_alpha_alpha(z, dr);
        H[1] = H[3] = skewnormal1_hess_mu_sigma(z, 1.0, a, r, dr);
        H[2] = H[6] = skewnormal1_hess_mu_alpha(z, 1.0, a, r, dr);
        H[5] = H[7] = skewnormal1_hess_sigma_alpha(z, 1.0, a, r, dr);
        const double h1 = -z + a * r;
        const double h2 = -1.0 + a * a * dr;
        const double h3 = a * a * a * g3;
        double v[10];
        v[0] = skewnormal1_d3_mu_mu_mu(1.0, a, g3);
        v[1] = skewnormal1_d3_mu_mu_sigma(z, 1.0, h2, h3);
        v[2] = skewnormal1_d3_mu_mu_alpha(t, 1.0, a, dr, g3);
        v[3] = skewnormal1_d3_mu_sigma_sigma(z, 1.0, h1, h2, h3);
        v[4] = skewnormal1_d3_mu_sigma_alpha(t, 1.0, r, dr, g3);
        v[5] = skewnormal1_d3_mu_alpha_alpha(z, t, 1.0, dr, g3);
        v[6] = skewnormal1_d3_sigma_sigma_sigma(z, 1.0, a, r, dr, g3);
        v[7] = skewnormal1_d3_sigma_sigma_alpha(z, t, 1.0, r, dr, g3);
        v[8] = skewnormal1_d3_sigma_alpha_alpha(z, t, 1.0, dr, g3);
        v[9] = skewnormal1_d3_alpha_alpha_alpha(z, g3);
        skewnormal1_fill_t3(v, T3);
        S.add(fw, g, H, T3);
    }
    S.result(I, D);
}

// The routers of the joint registry for the expected Hessian and its
// derivative, in (mu, sigma, alpha): P x P and P x P x P, column-major, both
// triangles filled. The standardized pair is the full quadrature, kept per
// thread for the last shape (tag 0 of loc_scale_diag_cpp()).
inline void skewnormal1_joint_info(double y, const double* th, double* out) {
    (void) y;
    const LocScaleFullSlot& s = loc_scale_full_cached(
        0, th + 2, 1, [&](double* I, double* D) {
            skewnormal1_quad_full(th + 2, I, D);
        });
    loc_scale_full_info(s, 3, th[1], out);
}

inline void skewnormal1_joint_dinfo(double y, const double* th, double* out) {
    (void) y;
    const LocScaleFullSlot& s = loc_scale_full_cached(
        0, th + 2, 1, [&](double* I, double* D) {
            skewnormal1_quad_full(th + 2, I, D);
        });
    loc_scale_full_dinfo(s, 3, th[1], out);
}

// The parts route (pt_parts.h), P = 3, (mu, sigma, alpha). The buffer:
// [0] level, [1] y, [2] mu, [3] sigma, [4] alpha, then the derivatives of
// log Phi at t = alpha z, from one call of logphi_derivs():
// [5] g1 = r (level 1), [6] g2 = dr (level 2), [7] g3 (level 3), [8] g4
// (level 4). The expected parts hold y as NaN, the parameters and no g.
inline void skewnormal1_parts(double y, const double* th, int level,
                              double* b) {
    const double m = th[0], s = th[1], a = th[2];
    b[0] = level; b[1] = y; b[2] = m; b[3] = s; b[4] = a;
    for (int k = 5; k <= 8; ++k) b[k] = R_NaN;
    if (level >= 1) {
        const double z = (y - m) / s;
        const int n = level > 4 ? 4 : level;
        double gg[6];
        logphi_derivs(a * z, n, gg);
        for (int k = 1; k <= n; ++k) b[4 + k] = gg[k];
    }
}

inline void skewnormal1_parts_expected(const double* th, int level,
                                       double* b) {
    b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = th[1]; b[4] = th[2];
    for (int k = 5; k <= 8; ++k) b[k] = R_NaN;
}

inline void skewnormal1_parts_grad(const double* b, double* out) {
    const double s = b[3], a = b[4];
    const double z = (b[1] - b[2]) / s;
    const double r = b[5];
    out[0] = skewnormal1_score_mu(z, s, a, r);
    out[1] = skewnormal1_score_sigma(z, s, a, r);
    out[2] = skewnormal1_score_alpha(z, r);
}

inline void skewnormal1_parts_hess(const double* b, double* out) {
    const double s = b[3], a = b[4];
    const double z = (b[1] - b[2]) / s;
    const double r = b[5], dr = b[6];
    out[0] = skewnormal1_hess_mu_mu(s, a, dr);
    out[4] = skewnormal1_hess_sigma_sigma(z, s, a, r, dr);
    out[8] = skewnormal1_hess_alpha_alpha(z, dr);
    out[1] = out[3] = skewnormal1_hess_mu_sigma(z, s, a, r, dr);
    out[2] = out[6] = skewnormal1_hess_mu_alpha(z, s, a, r, dr);
    out[5] = out[7] = skewnormal1_hess_sigma_alpha(z, s, a, r, dr);
}

// the expected Hessian, from the expected parts: the full quadrature in the
// shape, as skewnormal1_joint_info()
inline void skewnormal1_parts_info(const double* b, double* out) {
    const double shape[1] = {b[4]};
    const LocScaleFullSlot& s = loc_scale_full_cached(
        0, shape, 1, [&](double* I, double* D) {
            skewnormal1_quad_full(shape, I, D);
        });
    loc_scale_full_info(s, 3, b[3], out);
}

// the ten unique components, as skewnormal_deriv3_cpp
inline void skewnormal1_parts_d3(const double* b, double* out) {
    const double s = b[3], a = b[4];
    const double s2 = s * s, s3 = s2 * s;
    const double z = (b[1] - b[2]) / s;
    const double t = a * z;
    const double z2 = z * z, t2 = t * t;
    const double g1 = b[5], g2 = b[6], g3 = b[7];
    const double h1 = -z + a * g1;
    const double h2 = -1.0 + a * a * g2;
    const double h3 = a * a * a * g3;
    double u[10];
    u[0] = skewnormal1_d3_mu_mu_mu(s, a, g3);
    u[1] = -(z * h3 + 2.0 * h2) / s3;
    u[2] = a * (2.0 * g2 + t * g3) / s2;
    u[3] = -(z2 * h3 + 4.0 * z * h2 + 2.0 * h1) / s3;
    u[4] = (g1 + 3.0 * t * g2 + t2 * g3) / s2;
    u[5] = -z * (2.0 * g2 + t * g3) / s;
    u[6] = skewnormal1_d3_sigma_sigma_sigma(z, s, a, g1, g2, g3);
    u[7] = z * (2.0 * g1 + 4.0 * t * g2 + t2 * g3) / s2;
    u[8] = -z2 * (2.0 * g2 + t * g3) / s;
    u[9] = skewnormal1_d3_alpha_alpha_alpha(z, g3);
    sym_fill(3, 3, u, out);
}

// the fifteen unique components, as skewnormal_deriv4_cpp
inline void skewnormal1_parts_d4(const double* b, double* out) {
    const double s = b[3], a = b[4];
    const double s2 = s * s, s3 = s2 * s, s4 = s2 * s2;
    const double z = (b[1] - b[2]) / s;
    const double t = a * z;
    const double z2 = z * z, z3 = z2 * z, z4 = z2 * z2;
    const double t2 = t * t, t3 = t2 * t;
    const double g1 = b[5], g2 = b[6], g3 = b[7], g4 = b[8];
    const double h1 = -z + a * g1;
    const double h2 = -1.0 + a * a * g2;
    const double h3 = a * a * a * g3;
    const double h4 = a * a * a * a * g4;
    double u[15];
    u[0] = h4 / s4;
    u[1] = (z * h4 + 3.0 * h3) / s4;
    u[2] = -a * a * (3.0 * g3 + t * g4) / s3;
    u[3] = (z2 * h4 + 6.0 * z * h3 + 6.0 * h2) / s4;
    u[4] = -a * (4.0 * g2 + 5.0 * t * g3 + t2 * g4) / s3;
    u[5] = (2.0 * g2 + 4.0 * t * g3 + t2 * g4) / s2;
    u[6] = (z3 * h4 + 9.0 * z2 * h3 + 18.0 * z * h2 + 6.0 * h1) / s4;
    u[7] = -(2.0 * g1 + 10.0 * t * g2 + 7.0 * t2 * g3 + t3 * g4) / s3;
    u[8] = z * (4.0 * g2 + 5.0 * t * g3 + t2 * g4) / s2;
    u[9] = -z2 * (3.0 * g3 + t * g4) / s;
    u[10] = (6.0 + 24.0 * z * h1 + 36.0 * z2 * h2 + 12.0 * z3 * h3
             + z4 * h4) / s4;
    u[11] = -z * (6.0 * g1 + 18.0 * t * g2 + 9.0 * t2 * g3 + t3 * g4) / s3;
    u[12] = z2 * (6.0 * g2 + 6.0 * t * g3 + t2 * g4) / s2;
    u[13] = -z3 * (3.0 * g3 + t * g4) / s;
    u[14] = z4 * g4;
    sym_fill(3, 4, u, out);
}

} // namespace d7

#endif
