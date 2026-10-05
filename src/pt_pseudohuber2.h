#ifndef D7_PT_PSEUDOHUBER2_H
#define D7_PT_PSEUDOHUBER2_H

#include <Rcpp.h>
#include <cmath>
#include <cstdlib>
#include "pt_bessel_k.h"
#include "pt_loc_scale.h"
#include "pt_sqrt.h"

// The pseudo-Huber family in its standard-deviation parametrization
// (mu, sigma, nu): with R(nu) = sqrt(nu) K_2(sqrt(nu)) / K_1(sqrt(nu)),
// r = y - mu, w = r^2/sigma^2, Q = nu + R(nu) w and D = sqrt(Q),
//   l = -D - log 2 - log sigma + h(nu)
// (pseudohuber2.cpp gives h). One function per component and order for the
// score, the diagonal of the Hessian and the diagonal of the third
// derivatives, called by the vector kernels (pseudohuber2.cpp), by the
// quadrature of the expected information and by the scalar registry. The
// nu-only quantities R^(k) and h^(k) come from nu_part(), which the
// higher-order kernels of pseudohuber2.cpp read as well.

namespace d7 {

// nu-only quantities to order K: R^(k) and h^(k), k = 0..K (h[0] unused).
struct NuPart {
  double R[6];
  double h[6];
};

inline void nu_part(double v, int K, NuPart &o) {
  const double t = d7::sqrt_cr(v);
  double km[8] = {0};
  for (int m = 0; m <= K + 2; m++) km[m] = bessel_k_scaled(t, (double) m);
  o.R[0] = t * km[2] / km[1];
  if (K == 0) return;
  // derivatives of t(nu) = nu^(1/2)
  double tk[6] = {0};
  {
    double c = 1.0;
    for (int k = 1; k <= K; k++) {
      c *= (0.5 - (k - 1));
      tk[k] = c * std::pow(v, 0.5 - k);
    }
  }
  double L[3][6] = {{0}};
  for (int n = 1; n <= 2; n++) {
    // rho_j = K_n^(j)(t) / K_n(t)
    double rho[6] = {0};
    for (int j = 1; j <= K; j++) {
      double s = 0.0, binom = 1.0;
      for (int i = 0; i <= j; i++) {
        if (i > 0) binom = binom * (j - i + 1) / i;
        s += binom * km[std::abs(n - j + 2 * i)];
      }
      rho[j] = std::pow(-0.5, j) * s / km[n];
    }
    // g_k = (d/dt)^k log K_n(t)
    double g[6] = {0};
    const double r1 = rho[1];
    g[1] = r1;
    if (K >= 2) g[2] = rho[2] - r1 * r1;
    if (K >= 3) g[3] = rho[3] - 3.0 * r1 * rho[2] + 2.0 * r1 * r1 * r1;
    if (K >= 4) g[4] = rho[4] - 4.0 * r1 * rho[3] - 3.0 * rho[2] * rho[2] +
                       12.0 * r1 * r1 * rho[2] - 6.0 * std::pow(r1, 4);
    if (K >= 5) g[5] = rho[5] - 5.0 * r1 * rho[4] - 10.0 * rho[2] * rho[3] +
                       20.0 * r1 * r1 * rho[3] + 30.0 * r1 * rho[2] * rho[2] -
                       60.0 * std::pow(r1, 3) * rho[2] + 24.0 * std::pow(r1, 5);
    L[n][1] = g[1] * tk[1];
    if (K >= 2) L[n][2] = g[2] * tk[1] * tk[1] + g[1] * tk[2];
    if (K >= 3) L[n][3] = g[3] * std::pow(tk[1], 3) + 3.0 * g[2] * tk[1] * tk[2] +
                          g[1] * tk[3];
    if (K >= 4) L[n][4] = g[4] * std::pow(tk[1], 4) +
                          6.0 * g[3] * tk[1] * tk[1] * tk[2] +
                          g[2] * (3.0 * tk[2] * tk[2] + 4.0 * tk[1] * tk[3]) +
                          g[1] * tk[4];
    if (K >= 5) L[n][5] = g[5] * std::pow(tk[1], 5) +
                          10.0 * g[4] * std::pow(tk[1], 3) * tk[2] +
                          g[3] * (15.0 * tk[1] * tk[2] * tk[2] +
                                  10.0 * tk[1] * tk[1] * tk[3]) +
                          g[2] * (10.0 * tk[2] * tk[3] + 5.0 * tk[1] * tk[4]) +
                          g[1] * tk[5];
  }
  // (d/dnu)^k log nu = (-1)^(k-1) (k-1)! / nu^k
  double lg[6] = {0};
  {
    double f = 1.0;
    for (int k = 1; k <= K; k++) {
      if (k > 1) f *= -(k - 1);
      lg[k] = f / std::pow(v, k);
    }
  }
  double lam[6] = {0};
  for (int k = 1; k <= K; k++) lam[k] = 0.5 * lg[k] + L[2][k] - L[1][k];
  // R^(k) = R Y_k(lambda), Y_k the complete Bell polynomial
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
  for (int k = 1; k <= K; k++)
    o.h[k] = -0.25 * lg[k] + 0.5 * L[2][k] - 1.5 * L[1][k];
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

} // namespace d7

#endif
