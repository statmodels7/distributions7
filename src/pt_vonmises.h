#ifndef D7_PT_VONMISES_H
#define D7_PT_VONMISES_H

#include <Rcpp.h>
#include <cmath>

// The von Mises families, vonmises1 in (mu, kappa) and vonmises2 in
// (mu, rho), rho = A(kappa) = I1(kappa)/I0(kappa): one function per
// component and order for the quantities the scalar registry reads, called
// by the vector kernels (vonmises.cpp) and by the registry
// (d7_ccallable.cpp). vonmises.cpp records the derivation. The caller
// passes sin(y - mu), cos(y - mu), the ratio A and its derivatives a_j for
// vonmises1, and kappa(rho) with its derivatives k_j for vonmises2.

namespace d7 {

// A, its derivatives and the inverse come from numericals7's compiled
// functions (src/bessel_ratio.cpp there), resolved once through
// R_GetCCallable on the calling thread, before any worker starts.
struct N7Bessel {
  double (*A)(double);
  double (*d1)(double);
  double (*d2)(double);
  double (*d3)(double);
  void (*upto)(double, int, double*);
  double (*inv)(double);
};

// defined in vonmises.cpp; resolves on first use
const N7Bessel& vm_bessel();

// the derivatives of the inverse kappa(rho) to order n from A' ... A^(n)
// at kappa, into k[1..n]
inline void vm_inverse_derivs(const double* a, int n, double* k) {
  const double p1 = a[1], ip1 = 1.0 / p1;
  k[1] = ip1;
  if (n < 2) return;
  const double p2 = a[2], ip3 = ip1 * ip1 * ip1;
  k[2] = -p2 * ip3;
  if (n < 3) return;
  const double p3 = a[3], ip5 = ip3 * ip1 * ip1;
  k[3] = (3.0 * p2 * p2 - p1 * p3) * ip5;
  if (n < 4) return;
  const double p4 = a[4];
  k[4] = (-15.0 * p2 * p2 * p2 + 10.0 * p1 * p2 * p3 - p1 * p1 * p4) *
    ip5 * ip1 * ip1;
}

// kappa(rho) and its derivatives to order n (1 <= n <= 4), k[0] = kappa
inline void vm2_kappa(const N7Bessel& B, double rho, int n, double* k) {
  double a[5];
  k[0] = B.inv(rho);
  B.upto(k[0], n, a);
  vm_inverse_derivs(a, n, k);
}

// --- vonmises1 ---------------------------------------------------------------

inline double vonmises1_score_mu(double kappa, double sin_d) {
  return kappa * sin_d;
}

inline double vonmises1_score_kappa(double cos_d, double A) {
  return cos_d - A;
}

inline double vonmises1_hess_mu_mu(double kappa, double cos_d) {
  return -kappa * cos_d;
}

inline double vonmises1_hess_kappa_kappa(double A1) {
  return -A1;
}

inline double vonmises1_expected_mu_mu(double kappa, double A) {
  return -kappa * A;
}

inline double vonmises1_expected_kappa_kappa(double A1) {
  return -A1;
}

// E_mm = -kappa A does not depend on mu
inline double vonmises1_dexpected_mu_mu_mu() {
  return 0.0;
}

inline double vonmises1_dexpected_kappa_kappa_kappa(double A2) {
  return -A2;
}

inline void vonmises1_score_curv(int k, double y, const double* th,
                                 double* out) {
  const N7Bessel& B = vm_bessel();
  const double kappa = th[1], d = y - th[0];
  if (k == 0) {
    out[0] = vonmises1_score_mu(kappa, std::sin(d));
    out[1] = vonmises1_hess_mu_mu(kappa, std::cos(d));
  } else {
    out[0] = vonmises1_score_kappa(std::cos(d), B.A(kappa));
    out[1] = vonmises1_hess_kappa_kappa(B.d1(kappa));
  }
}

// the expected information reads A and A' from one evaluation of the pair,
// as its vector kernel does
inline void vonmises1_info_dinfo(int k, double y, const double* th,
                                 double* out) {
  const N7Bessel& B = vm_bessel();
  const double kappa = th[1];
  if (k == 0) {
    double a[2];
    B.upto(kappa, 1, a);
    out[0] = vonmises1_expected_mu_mu(kappa, a[0]);
    out[1] = vonmises1_dexpected_mu_mu_mu();
  } else {
    double a[3];
    B.upto(kappa, 2, a);
    out[0] = vonmises1_expected_kappa_kappa(a[1]);
    out[1] = vonmises1_dexpected_kappa_kappa_kappa(a[2]);
  }
}

// --- vonmises2 ---------------------------------------------------------------

inline double vonmises2_score_mu(double k0, double sin_d) {
  return k0 * sin_d;
}

inline double vonmises2_score_rho(double cos_d, double r, double k1) {
  return (cos_d - r) * k1;
}

inline double vonmises2_hess_mu_mu(double k0, double cos_d) {
  return -k0 * cos_d;
}

inline double vonmises2_hess_rho_rho(double cos_d, double r, double k1,
                                     double k2) {
  return (cos_d - r) * k2 - k1;
}

inline double vonmises2_expected_mu_mu(double k0, double r) {
  return -k0 * r;
}

inline double vonmises2_expected_rho_rho(double k1) {
  return -k1;
}

// E_mm = -rho kappa(rho) does not depend on mu
inline double vonmises2_dexpected_mu_mu_mu() {
  return 0.0;
}

inline double vonmises2_dexpected_rho_rho_rho(double k2) {
  return -k2;
}

inline void vonmises2_score_curv(int k, double y, const double* th,
                                 double* out) {
  const N7Bessel& B = vm_bessel();
  const double r = th[1], d = y - th[0];
  if (k == 0) {
    double kk[2];
    vm2_kappa(B, r, 1, kk);
    out[0] = vonmises2_score_mu(kk[0], std::sin(d));
    out[1] = vonmises2_hess_mu_mu(kk[0], std::cos(d));
  } else {
    double kk[3];
    vm2_kappa(B, r, 2, kk);
    const double c = std::cos(d);
    out[0] = vonmises2_score_rho(c, r, kk[1]);
    out[1] = vonmises2_hess_rho_rho(c, r, kk[1], kk[2]);
  }
}

inline void vonmises2_info_dinfo(int k, double y, const double* th,
                                 double* out) {
  const N7Bessel& B = vm_bessel();
  const double r = th[1];
  if (k == 0) {
    double kk[2];
    vm2_kappa(B, r, 1, kk);
    out[0] = vonmises2_expected_mu_mu(kk[0], r);
    out[1] = vonmises2_dexpected_mu_mu_mu();
  } else {
    double kk[3];
    vm2_kappa(B, r, 2, kk);
    out[0] = vonmises2_expected_rho_rho(kk[1]);
    out[1] = vonmises2_dexpected_rho_rho_rho(kk[2]);
  }
}

} // namespace d7

#endif
