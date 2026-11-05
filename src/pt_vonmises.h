#ifndef D7_PT_VONMISES_H
#define D7_PT_VONMISES_H

#include <Rcpp.h>
#include <cmath>
#include "pt_parts.h"

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

// l_{mu kappa} = sin(y - mu)
inline double vonmises1_hess_mu_kappa(double sin_d) {
  return sin_d;
}

// E_mm = -kappa A, so its derivative in kappa is -(A + kappa A')
inline double vonmises1_dexpected_mu_mu_kappa(double kappa, double A,
                                              double A1) {
  return -(A + kappa * A1);
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

// l_{mu rho} = sin(y - mu) kappa'(rho)
inline double vonmises2_hess_mu_rho(double sin_d, double k1) {
  return sin_d * k1;
}

// E_mm = -rho kappa(rho), so its derivative in rho is -(kappa + rho kappa')
inline double vonmises2_dexpected_mu_mu_rho(double r, double k0, double k1) {
  return -(k0 + r * k1);
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

// The routers of the joint registry, one per order: every component at one
// observation, on the parameter scale, th = (mu, kappa) for vonmises1 and
// (mu, rho) for vonmises2. The Hessian and the expected information fill a
// full 2 x 2 matrix and the derivative of the expected information a
// 2 x 2 x 2 array, column-major, out[a + 2 b + 4 c] being the derivative of
// E[l_ab] in parameter c. The expected information does not depend on mu and
// E[l_12] = 0, so only the entries in the second parameter are non-zero.
inline void vonmises1_joint_grad(double y, const double* th, double* out) {
  const N7Bessel& B = vm_bessel();
  const double kappa = th[1], d = y - th[0];
  out[0] = vonmises1_score_mu(kappa, std::sin(d));
  out[1] = vonmises1_score_kappa(std::cos(d), B.A(kappa));
}

inline void vonmises1_joint_hess(double y, const double* th, double* out) {
  const N7Bessel& B = vm_bessel();
  const double kappa = th[1], d = y - th[0];
  out[0] = vonmises1_hess_mu_mu(kappa, std::cos(d));
  out[3] = vonmises1_hess_kappa_kappa(B.d1(kappa));
  out[1] = out[2] = vonmises1_hess_mu_kappa(std::sin(d));
}

inline void vonmises1_joint_info(double y, const double* th, double* out) {
  const N7Bessel& B = vm_bessel();
  const double kappa = th[1];
  double a[2];
  B.upto(kappa, 1, a);
  out[0] = vonmises1_expected_mu_mu(kappa, a[0]);
  out[3] = vonmises1_expected_kappa_kappa(a[1]);
  out[1] = out[2] = 0.0;
}

inline void vonmises1_joint_dinfo(double y, const double* th, double* out) {
  const N7Bessel& B = vm_bessel();
  const double kappa = th[1];
  double a[3];
  B.upto(kappa, 2, a);
  for (int j = 0; j < 8; ++j) out[j] = 0.0;
  out[4] = vonmises1_dexpected_mu_mu_kappa(kappa, a[0], a[1]);
  out[7] = vonmises1_dexpected_kappa_kappa_kappa(a[2]);
}

inline void vonmises2_joint_grad(double y, const double* th, double* out) {
  const N7Bessel& B = vm_bessel();
  const double r = th[1], d = y - th[0];
  double kk[2];
  vm2_kappa(B, r, 1, kk);
  out[0] = vonmises2_score_mu(kk[0], std::sin(d));
  out[1] = vonmises2_score_rho(std::cos(d), r, kk[1]);
}

inline void vonmises2_joint_hess(double y, const double* th, double* out) {
  const N7Bessel& B = vm_bessel();
  const double r = th[1], d = y - th[0];
  double kk[3];
  vm2_kappa(B, r, 2, kk);
  const double c = std::cos(d);
  out[0] = vonmises2_hess_mu_mu(kk[0], c);
  out[3] = vonmises2_hess_rho_rho(c, r, kk[1], kk[2]);
  out[1] = out[2] = vonmises2_hess_mu_rho(std::sin(d), kk[1]);
}

inline void vonmises2_joint_info(double y, const double* th, double* out) {
  const N7Bessel& B = vm_bessel();
  const double r = th[1];
  double kk[2];
  vm2_kappa(B, r, 1, kk);
  out[0] = vonmises2_expected_mu_mu(kk[0], r);
  out[3] = vonmises2_expected_rho_rho(kk[1]);
  out[1] = out[2] = 0.0;
}

inline void vonmises2_joint_dinfo(double y, const double* th, double* out) {
  const N7Bessel& B = vm_bessel();
  const double r = th[1];
  double kk[3];
  vm2_kappa(B, r, 2, kk);
  for (int j = 0; j < 8; ++j) out[j] = 0.0;
  out[4] = vonmises2_dexpected_mu_mu_rho(r, kk[0], kk[1]);
  out[7] = vonmises2_dexpected_rho_rho_rho(kk[2]);
}

// The parts route (pt_parts.h), vonmises1. The buffer, with mu = th[0] and
// kappa = th[1]; the Bessel ratio and its derivatives are read once, each
// from the numericals7 function the vector kernels use for that order:
//
//   [0]  level
//   [1]  y                          NaN in the expected buffer
//   [2]  mu
//   [3]  kappa
//   [4]  A(kappa)                   level 1 (expected buffer: level 2)
//   [5]  A'(kappa)                  level 2 (expected buffer: level 2)
//   [6]  A''(kappa)                 level 3
//   [7]  A'''(kappa)                level 4
//
// The expected buffer reads A and A' from one evaluation of the pair
// (upto), as vonmises1_joint_info() does, and fills [6] and [7] with NaN;
// the expected Hessian is the only expected order read.
inline void vonmises1_parts(double y, const double* th, int level, double* b) {
  const double kappa = th[1];
  b[0] = level; b[1] = y; b[2] = th[0]; b[3] = kappa;
  if (level >= 1) {
    const N7Bessel& B = vm_bessel();
    b[4] = level >= 1 ? B.A(kappa) : R_NaN;
    b[5] = level >= 2 ? B.d1(kappa) : R_NaN;
    b[6] = level >= 3 ? B.d2(kappa) : R_NaN;
    b[7] = level >= 4 ? B.d3(kappa) : R_NaN;
  } else {
    for (int j = 4; j <= 7; ++j) b[j] = R_NaN;
  }
}

inline void vonmises1_parts_expected(const double* th, int level, double* b) {
  const double kappa = th[1];
  b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = kappa;
  if (level >= 2) {
    const N7Bessel& B = vm_bessel();
    double a[2];
    B.upto(kappa, 1, a);
    b[4] = a[0];
    b[5] = a[1];
  } else {
    b[4] = R_NaN;
    b[5] = R_NaN;
  }
  b[6] = R_NaN;
  b[7] = R_NaN;
}

inline void vonmises1_parts_grad(const double* b, double* out) {
  const double kappa = b[3], d = b[1] - b[2];
  out[0] = vonmises1_score_mu(kappa, std::sin(d));
  out[1] = vonmises1_score_kappa(std::cos(d), b[4]);
}

inline void vonmises1_parts_hess(const double* b, double* out) {
  const double kappa = b[3], d = b[1] - b[2];
  out[0] = vonmises1_hess_mu_mu(kappa, std::cos(d));
  out[3] = vonmises1_hess_kappa_kappa(b[5]);
  out[1] = out[2] = vonmises1_hess_mu_kappa(std::sin(d));
}

// from the expected parts
inline void vonmises1_parts_info(const double* b, double* out) {
  const double kappa = b[3];
  out[0] = vonmises1_expected_mu_mu(kappa, b[4]);
  out[3] = vonmises1_expected_kappa_kappa(b[5]);
  out[1] = out[2] = 0.0;
}

// vonmises1_deriv3_cpp() of vonmises.cpp, the same expressions
inline void vonmises1_parts_d3(const double* b, double* out) {
  const double k = b[3], d = b[1] - b[2];
  double u[4];
  u[0] = -k * std::sin(d);
  u[1] = -std::cos(d);
  u[2] = 0.0;
  u[3] = -b[6];
  sym_fill(2, 3, u, out);
}

// vonmises1_deriv4_cpp() of vonmises.cpp, the same expressions
inline void vonmises1_parts_d4(const double* b, double* out) {
  const double k = b[3], d = b[1] - b[2];
  double u[5];
  u[0] = k * std::cos(d);
  u[1] = -std::sin(d);
  u[2] = 0.0;
  u[3] = 0.0;
  u[4] = -b[7];
  sym_fill(2, 4, u, out);
}

// The parts route, vonmises2, with mu = th[0] and rho = th[1]. kappa(rho)
// and its derivatives come from vm2_kappa() once, to the order the level
// asks for (1 to 4):
//
//   [0]  level
//   [1]  y                          NaN in the expected buffer
//   [2]  mu
//   [3]  rho
//   [4]  kappa(rho)                 level 1 (expected buffer: level 2)
//   [5]  kappa'(rho)                level 1 (expected buffer: level 2)
//   [6]  kappa''(rho)               level 2
//   [7]  kappa'''(rho)              level 3
//   [8]  kappa''''(rho)             level 4
//
// The expected buffer evaluates kappa to the first derivative, as
// vonmises2_joint_info() does, and fills [6] to [8] with NaN. Orders 1 to 4
// are all closed form here; each calls vm2_kappa() with the order of the
// buffer's level, whereas each vector kernel calls it with its own order.
inline void vonmises2_parts(double y, const double* th, int level, double* b) {
  const double r = th[1];
  b[0] = level; b[1] = y; b[2] = th[0]; b[3] = r;
  if (level >= 1) {
    const N7Bessel& B = vm_bessel();
    double kk[5];
    vm2_kappa(B, r, level < 4 ? level : 4, kk);
    b[4] = kk[0];
    b[5] = kk[1];
    b[6] = level >= 2 ? kk[2] : R_NaN;
    b[7] = level >= 3 ? kk[3] : R_NaN;
    b[8] = level >= 4 ? kk[4] : R_NaN;
  } else {
    for (int j = 4; j <= 8; ++j) b[j] = R_NaN;
  }
}

inline void vonmises2_parts_expected(const double* th, int level, double* b) {
  const double r = th[1];
  b[0] = level; b[1] = R_NaN; b[2] = th[0]; b[3] = r;
  if (level >= 2) {
    const N7Bessel& B = vm_bessel();
    double kk[2];
    vm2_kappa(B, r, 1, kk);
    b[4] = kk[0];
    b[5] = kk[1];
  } else {
    b[4] = R_NaN;
    b[5] = R_NaN;
  }
  for (int j = 6; j <= 8; ++j) b[j] = R_NaN;
}

inline void vonmises2_parts_grad(const double* b, double* out) {
  const double r = b[3], d = b[1] - b[2];
  out[0] = vonmises2_score_mu(b[4], std::sin(d));
  out[1] = vonmises2_score_rho(std::cos(d), r, b[5]);
}

inline void vonmises2_parts_hess(const double* b, double* out) {
  const double r = b[3], d = b[1] - b[2];
  const double c = std::cos(d);
  out[0] = vonmises2_hess_mu_mu(b[4], c);
  out[3] = vonmises2_hess_rho_rho(c, r, b[5], b[6]);
  out[1] = out[2] = vonmises2_hess_mu_rho(std::sin(d), b[5]);
}

// from the expected parts
inline void vonmises2_parts_info(const double* b, double* out) {
  out[0] = vonmises2_expected_mu_mu(b[4], b[3]);
  out[3] = vonmises2_expected_rho_rho(b[5]);
  out[1] = out[2] = 0.0;
}

// vonmises2_deriv3_cpp() of vonmises.cpp, the same expressions
inline void vonmises2_parts_d3(const double* b, double* out) {
  const double r = b[3], d = b[1] - b[2];
  const double c = std::cos(d), s = std::sin(d);
  double u[4];
  u[0] = -b[4] * s;
  u[1] = -c * b[5];
  u[2] = s * b[6];
  u[3] = (c - r) * b[7] - 2.0 * b[6];
  sym_fill(2, 3, u, out);
}

// vonmises2_deriv4_cpp() of vonmises.cpp, the same expressions
inline void vonmises2_parts_d4(const double* b, double* out) {
  const double r = b[3], d = b[1] - b[2];
  const double c = std::cos(d), s = std::sin(d);
  double u[5];
  u[0] = b[4] * c;
  u[1] = -s * b[5];
  u[2] = -c * b[6];
  u[3] = s * b[7];
  u[4] = (c - r) * b[8] - 3.0 * b[7];
  sym_fill(2, 4, u, out);
}

} // namespace d7

#endif
