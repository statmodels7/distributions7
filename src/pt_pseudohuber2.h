#ifndef D7_PT_PSEUDOHUBER2_H
#define D7_PT_PSEUDOHUBER2_H

#include <Rcpp.h>
#include <cmath>
#include <cstdlib>
#include "pt_bessel_k.h"
#include "pt_loc_scale.h"
#include "pt_sqrt.h"
#include "pt_ph_nu.h"

// The pseudo-Huber family in its standard-deviation parametrization
// (mu, sigma, nu): with R(nu) = sqrt(nu) K_2(sqrt(nu)) / K_1(sqrt(nu)),
// r = y - mu, w = r^2/sigma^2, Q = nu + R(nu) w and D = sqrt(Q),
//   l = -D - log 2 - log sigma + h(nu)
// (pseudohuber2.cpp gives h). One function per component and order for the
// score, the diagonal of the Hessian and the diagonal of the third
// derivatives, called by the vector kernels (pseudohuber2.cpp), by the
// quadrature of the expected information and by the scalar registry. The
// nu-only quantities R^(k) and h^(k) come from nu_part(), which the
// higher-order kernels of pseudohuber2.cpp read as well, and which takes them
// from ph_nu_derivs() (pt_ph_nu.h).

namespace d7 {

// nu-only quantities to order K: R^(k) and h^(k), k = 0..K (h[0] unused).
struct NuPart {
  double R[6];
  double h[6];
};

inline void nu_part(double v, int K, NuPart &o) {
  // B^(k) and (log R)^(k) from ph_nu_derivs() (pt_ph_nu.h); R^(k) is R times
  // the complete Bell polynomial of the (log R)^(k), and h = log(R)/2 - B
  double dB[6] = {0}, lam[6] = {0}, R0;
  ph_nu_derivs(v, K > 0 ? K : 1, dB, lam, R0);
  o.R[0] = R0;
  if (K == 0) return;
  o.R[1] = o.R[0] * lam[1];
  if (K >= 2) o.R[2] = o.R[0] * (lam[1] * lam[1] + lam[2]);
  if (K >= 3) o.R[3] = o.R[0] * (std::pow(lam[1], 3) + 3.0 * lam[1] * lam[2] +
                                 lam[3]);
  if (K >= 4) o.R[4] = o.R[0] * (std::pow(lam[1], 4) +
                                 6.0 * lam[1] * lam[1] * lam[2] +
                                 4.0 * lam[1] * lam[3] + 3.0 * lam[2] * lam[2] +
                                 lam[4]);
  if (K >= 5) o.R[5] = o.R[0] * (std::pow(lam[1], 5) +
                                 10.0 * std::pow(lam[1], 3) * lam[2] +
                                 15.0 * lam[1] * lam[2] * lam[2] +
                                 10.0 * lam[1] * lam[1] * lam[3] +
                                 10.0 * lam[2] * lam[3] + 5.0 * lam[1] * lam[4] +
                                 lam[5]);
  for (int k = 1; k <= K; k++) o.h[k] = 0.5 * lam[k] - dB[k];
}

// R(nu) and h(nu) themselves, as the density reads them
struct Ph2Terms { double R, h; };

inline Ph2Terms pseudohuber2_terms(double v) {
    Ph2Terms T;
    const double t = d7::sqrt_cr(v);
    const double k1 = bessel_k_scaled(t, 1.0), k2 = bessel_k_scaled(t, 2.0);
    T.R = t * k2 / k1;
    // log K_n(t) = log(scaled K_n(t)) - t
    T.h = -0.25 * std::log(v) + 0.5 * std::log(k2) - 1.5 * std::log(k1) + t;
    return T;
}

// D = sqrt(nu + R z^2) and zD = z/D, with z = (y - mu)/sigma, formed as
// |z| sqrt(R + nu/z^2) above |z| = 1 so that z^2 never overflows. The
// components of orders one and two are written in z, D and zD: each
// Q_x/D is then of the order of |z| at most and each Q_x/D^2 bounded,
// where Q_x Q_y/(4 D^3) taken apart overflowed from |y| of about 1e80.
struct Ph2Z { double D, zD; };

inline Ph2Z ph2_z(double z, double v, double R0) {
    Ph2Z P;
    if (std::fabs(z) > 1.0) {
        const double a = d7::sqrt_cr(R0 + v / (z * z));
        P.D = std::fabs(z) * a;
        P.zD = (z > 0 ? 1.0 : -1.0) / a;
    } else {
        P.D = d7::sqrt_cr(v + R0 * (z * z));
        P.zD = z / P.D;
    }
    return P;
}

inline double pseudohuber2_score_mu(double s, const Ph2Z& P,
                                    const NuPart& np) {
    return np.R[0] * P.zD / s;
}

inline double pseudohuber2_score_sigma(double z, double s, const Ph2Z& P,
                                       const NuPart& np) {
    return (np.R[0] * z * P.zD - 1.0) / s;
}

inline double pseudohuber2_score_nu(double z, const Ph2Z& P,
                                    const NuPart& np) {
    return -0.5 * (1.0 / P.D + np.R[1] * z * P.zD) + np.h[1];
}

// l_xx = -(Q_xx / (2 D) - (Q_x / D)^2 / (4 D)) + the terms free of Q
// 1 - R zD^2 is nu/D^2 exactly (D^2 = nu + R z^2), and is written so: as a
// difference it cancels as |z| grows, and l_mu,mu read 2.7e-26 for -5.7e-31
// at |z| of 1e10
inline double pseudohuber2_hess_mu_mu(double s, const Ph2Z& P,
                                      const NuPart& np, double v) {
    const double R0 = np.R[0], iD = 1.0 / P.D;
    return -(R0 / (s * s * P.D)) * (v * iD * iD);
}

inline double pseudohuber2_hess_sigma_sigma(double z, double s,
                                            const Ph2Z& P, const NuPart& np,
                                            double v) {
    const double R0 = np.R[0], x = R0 * z * P.zD, iD = 1.0 / P.D;
    return (1.0 - 2.0 * x - x * (v * iD * iD)) / (s * s);
}

// Q_nu / D and Q_nu / D^2
inline void ph2_qn(double z, const Ph2Z& P, const NuPart& np, double* qn,
                   double* qnD) {
    const double iD = 1.0 / P.D;
    *qn = iD + np.R[1] * z * P.zD;
    *qnD = iD * iD + np.R[1] * P.zD * P.zD;
}

inline double pseudohuber2_hess_nu_nu(double z, const Ph2Z& P,
                                      const NuPart& np) {
    double qn, qnD;
    ph2_qn(z, P, np, &qn, &qnD);
    return -(0.5 * np.R[2] * z * P.zD - 0.25 * qn * qnD) + np.h[2];
}

// the off-diagonal components of the Hessian, as pseudohuber2_hessian_cpp
// forms them; np carries R and h to order two. With vD2 = nu/D^2,
// 2 - R zD^2 = 1 + vD2 and 1 - R zD^2 / 2 = (1 + vD2)/2 exactly
inline double pseudohuber2_hess_mu_sigma(double s, const Ph2Z& P,
                                         const NuPart& np, double v) {
    const double R0 = np.R[0], s2 = s * s;
    const double vD2 = v / (P.D * P.D);
    return -R0 * P.zD * (1.0 + vD2) / s2;
}

inline double pseudohuber2_hess_mu_nu(double s, const Ph2Z& P,
                                      const NuPart& np, double v) {
    const double R0 = np.R[0], R1 = np.R[1];
    const double vD2 = v / (P.D * P.D);
    return P.zD * (0.5 * R1 * (1.0 + vD2) - 0.5 * R0 / (P.D * P.D)) / s;
}

inline double pseudohuber2_hess_sigma_nu(double z, double s, const Ph2Z& P,
                                         const NuPart& np, double v) {
    const double R0 = np.R[0], R1 = np.R[1];
    const double zz = z * P.zD;
    const double vD2 = v / (P.D * P.D);
    return zz * (0.5 * R1 * (1.0 + vD2) - 0.5 * R0 / (P.D * P.D)) / s;
}

// the third derivatives on the diagonal: with (d/dQ)^m sqrt(Q) = c_m D/Q^m,
// c = (1/2, -1/4, 3/8),
//   (d/dx)^3 D = c_1 Q_xxx D/Q + 3 c_2 Q_x Q_xx D/Q^2 + c_3 Q_x^3 D/Q^3;
// Q_mumumu = 0, and np must carry R and h to order 3 for nu
inline double ph2_d3_sqrt(double Q, double D, double Qx, double Qxx,
                          double Qxxx) {
    const double dq1 = D / Q, dq2 = dq1 / Q, dq3 = dq2 / Q;
    return 0.5 * Qxxx * dq1 - 0.75 * Qx * Qxx * dq2 + 0.375 * Qx * Qx * Qx * dq3;
}

inline double pseudohuber2_d3_mu_mu_mu(double r, double s2, double Q,
                                       double D, const NuPart& np) {
    const double R0 = np.R[0];
    const double Qm = -2.0 * R0 * r / s2, Qmm = 2.0 * R0 / s2;
    return -ph2_d3_sqrt(Q, D, Qm, Qmm, 0.0);
}

inline double pseudohuber2_d3_sigma_sigma_sigma(double s, double s2,
                                                double w, double Q, double D,
                                                const NuPart& np) {
    const double R0 = np.R[0];
    const double Qs = -2.0 * R0 * w / s, Qss = 6.0 * R0 * w / s2,
                 Qsss = -24.0 * R0 * w / (s2 * s);
    return -ph2_d3_sqrt(Q, D, Qs, Qss, Qsss) - 2.0 / (s2 * s);
}

inline double pseudohuber2_d3_nu_nu_nu(double w, double Q, double D,
                                       const NuPart& np) {
    const double Qn = 1.0 + np.R[1] * w, Qnn = np.R[2] * w,
                 Qnnn = np.R[3] * w;
    return -ph2_d3_sqrt(Q, D, Qn, Qnn, Qnnn) + np.h[3];
}

// log f, as distrib_pdf() forms it
inline double pseudohuber2_logpdf(double y, double mu, double sigma, double v,
                                  const Ph2Terms& T) {
    const double z = (y - mu) / sigma;
    const double D = d7::sqrt_cr(v + T.R * (z * z));
    return -D - std::log(2.0) - std::log(sigma) + T.h;
}

// the routers of the scalar registry: th = (mu, sigma, nu)
inline void pseudohuber2_score_curv(int k, double y, const double* th,
                                    double* out) {
    const double m = th[0], s = th[1], v = th[2];
    NuPart np = {{0}, {0}};
    nu_part(v, k == 2 ? 2 : 0, np);
    const double z = (y - m) / s;
    const Ph2Z P = ph2_z(z, v, np.R[0]);
    if (k == 0) {
        out[0] = pseudohuber2_score_mu(s, P, np);
        out[1] = pseudohuber2_hess_mu_mu(s, P, np, v);
    } else if (k == 1) {
        out[0] = pseudohuber2_score_sigma(z, s, P, np);
        out[1] = pseudohuber2_hess_sigma_sigma(z, s, P, np, v);
    } else {
        out[0] = pseudohuber2_score_nu(z, P, np);
        out[1] = pseudohuber2_hess_nu_nu(z, P, np);
    }
}

// the standardized diagonal pair of parameter k at shape = (nu), by the rule
// of pt_loc_scale.h
inline void pseudohuber2_quad_diag(int k, const double* shape, double* out,
                                   bool want_d) {
    const double v = shape[0];
    NuPart np = {{0}, {0}};
    nu_part(v, want_d && k == 2 ? 3 : 2, np);
    const Ph2Terms T = pseudohuber2_terms(v);
    const LocScaleRule R = loc_scale_rule();
    LocScaleSum S;
    for (int j = 0; j < R.n; ++j) {
        const double z = R.x[j];
        const double fw = std::exp(pseudohuber2_logpdf(z, 0.0, 1.0, v, T)) * R.w[j];
        const double r = z, w = r * r / 1.0;
        const double Q = v + np.R[0] * w, D = d7::sqrt_cr(Q);
        const Ph2Z P = ph2_z(z, v, np.R[0]);
        double H;
        if (k == 0) H = pseudohuber2_hess_mu_mu(1.0, P, np, v);
        else if (k == 1) H = pseudohuber2_hess_sigma_sigma(z, 1.0, P, np, v);
        else H = pseudohuber2_hess_nu_nu(z, P, np);
        if (!want_d) { S.add0(fw, H); continue; }
        double g, T3;
        if (k == 0) {
            g = pseudohuber2_score_mu(1.0, P, np);
            T3 = pseudohuber2_d3_mu_mu_mu(r, 1.0, Q, D, np);
        } else if (k == 1) {
            g = pseudohuber2_score_sigma(z, 1.0, P, np);
            T3 = pseudohuber2_d3_sigma_sigma_sigma(1.0, 1.0, w, Q, D, np);
        } else {
            g = pseudohuber2_score_nu(z, P, np);
            T3 = pseudohuber2_d3_nu_nu_nu(w, Q, D, np);
        }
        S.add(fw, g, H, T3);
    }
    S.result(out);
}

inline void pseudohuber2_info_dinfo(int k, double y, const double* th,
                                    double* out) {
    (void) y;
    double std2[2];
    loc_scale_cached(3, k, th + 2, 1, [&](double* v) {
        pseudohuber2_quad_diag(k, th + 2, v, true);
    }, std2);
    loc_scale_unstandardize(k, th[1], std2, out);
}

// The routers of the joint registry for the score and the Hessian: every
// component at one observation, on the parameter scale, th = (mu, sigma,
// nu), P = 3. The Hessian fills a full 3 x 3 matrix, column-major,
// out[a + 3 b]. The expected information of this family is an exp-sinh
// quadrature (the diagonal compiled in pseudohuber2_quad_diag(), the
// sigma-nu entry only in R), so there is no joint_info and no joint_dinfo.
inline void pseudohuber2_joint_grad(double y, const double* th, double* out) {
    const double m = th[0], s = th[1], v = th[2];
    NuPart np = {{0}, {0}};
    nu_part(v, 1, np);
    const double z = (y - m) / s;
    const Ph2Z P = ph2_z(z, v, np.R[0]);
    out[0] = pseudohuber2_score_mu(s, P, np);
    out[1] = pseudohuber2_score_sigma(z, s, P, np);
    out[2] = pseudohuber2_score_nu(z, P, np);
}

inline void pseudohuber2_joint_hess(double y, const double* th, double* out) {
    const double m = th[0], s = th[1], v = th[2];
    NuPart np = {{0}, {0}};
    nu_part(v, 2, np);
    const double z = (y - m) / s;
    const Ph2Z P = ph2_z(z, v, np.R[0]);
    out[0] = pseudohuber2_hess_mu_mu(s, P, np, v);
    out[4] = pseudohuber2_hess_sigma_sigma(z, s, P, np, v);
    out[8] = pseudohuber2_hess_nu_nu(z, P, np);
    out[1] = out[3] = pseudohuber2_hess_mu_sigma(s, P, np, v);
    out[2] = out[6] = pseudohuber2_hess_mu_nu(s, P, np, v);
    out[5] = out[7] = pseudohuber2_hess_sigma_nu(z, s, P, np, v);
}

} // namespace d7

#endif
