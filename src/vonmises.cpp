#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <algorithm>
#include <cmath>
#include "d7_par.h"
#include "pt_vonmises.h"
using namespace Rcpp;

// The von Mises families, one kernel per derivative order.
//
// vonmises1 (mu, kappa): log f = kappa cos(y - mu) - log(2 pi) - log I0(kappa),
// so every derivative in kappa alone is -A^(n-1)(kappa), with
// A = I1/I0, and kappa enters the mixed ones linearly.
//
// vonmises2 (mu, rho), rho = A(kappa): with kappa = kappa(rho) the inverse,
// d/drho log I0(kappa(rho)) = A(kappa) kappa' = rho kappa', so
//   l_rho^n = (cos(y - mu) - rho) kappa^(n) - (n - 1) kappa^(n-1),
// and only the derivatives of the inverse are needed, each from A' ... A^(n)
// by the inverse function rule. (cos(y - mu) - rho) is formed from rho
// itself, not from A(kappa(rho)).
//
// A, its derivatives and the inverse come from numericals7's compiled
// functions (src/bessel_ratio.cpp there), resolved once through
// R_GetCCallable on the calling thread, before any worker starts.

namespace d7 {

// Resolved on first use, on the calling thread. Not a guarded static: a
// lookup that fails raises an R error, and a longjmp out of a static
// initializer leaves its guard held, so every later call would block.
// d7_scalar_id() resolves it for the registry, on the consumer's thread.
static N7Bessel n7b_ptrs = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

const N7Bessel& vm_bessel() {
  if (n7b_ptrs.inv == nullptr) {
    N7Bessel b;
    b.A = (double (*)(double)) R_GetCCallable("numericals7", "n7_bessel_ratio");
    b.d1 = (double (*)(double)) R_GetCCallable("numericals7", "n7_bessel_ratio_d1");
    b.d2 = (double (*)(double)) R_GetCCallable("numericals7", "n7_bessel_ratio_d2");
    b.d3 = (double (*)(double)) R_GetCCallable("numericals7", "n7_bessel_ratio_d3");
    b.upto = (void (*)(double, int, double*)) R_GetCCallable("numericals7",
                                                             "n7_bessel_ratio_upto");
    b.inv = (double (*)(double)) R_GetCCallable("numericals7",
                                                "n7_bessel_ratio_inverse");
    n7b_ptrs = b;
  }
  return n7b_ptrs;
}

}  // namespace d7

namespace {
using d7::N7Bessel;
using d7::vm2_kappa;
inline const N7Bessel& n7b() { return d7::vm_bessel(); }
}  // namespace

// --- vonmises1 -------------------------------------------------------------

// [[Rcpp::export]]
List vonmises1_gradient_cpp(NumericVector y, NumericVector mu,
                            NumericVector kappa, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n_y = y.size(), n_mu = mu.size(), n_k = kappa.size();
  const int n = std::max(n_y, std::max(n_mu, n_k));
  NumericVector o_mu(n), o_kappa(n);
  const bool ks = n_k == 1;
  const double A0 = ks ? B.A(kappa[0]) : 0.0;
  d7::par_for(n, threads, ks ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const double k = kappa[i % n_k], d = y[i % n_y] - mu[i % n_mu];
    const double A = ks ? A0 : B.A(k);
    o_mu[i] = d7::vonmises1_score_mu(k, std::sin(d));
    o_kappa[i] = d7::vonmises1_score_kappa(std::cos(d), A);
  });
  return List::create(Named("mu") = o_mu, Named("kappa") = o_kappa);
}

// [[Rcpp::export]]
List vonmises1_hessian_cpp(NumericVector y, NumericVector mu,
                           NumericVector kappa, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n_y = y.size(), n_mu = mu.size(), n_k = kappa.size();
  const int n = std::max(n_y, std::max(n_mu, n_k));
  NumericVector o_mm(n), o_mk(n), o_kk(n);
  const bool ks = n_k == 1;
  const double D0 = ks ? B.d1(kappa[0]) : 0.0;
  d7::par_for(n, threads, ks ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const double k = kappa[i % n_k], d = y[i % n_y] - mu[i % n_mu];
    o_mm[i] = d7::vonmises1_hess_mu_mu(k, std::cos(d));
    o_mk[i] = std::sin(d);
    o_kk[i] = d7::vonmises1_hess_kappa_kappa(ks ? D0 : B.d1(k));
  });
  return List::create(Named("mu_mu") = o_mm, Named("mu_kappa") = o_mk,
                      Named("kappa_kappa") = o_kk);
}

// The expected information needs A and A', both used, so one evaluation
// of the pair serves them.
// [[Rcpp::export]]
List vonmises1_expected_hessian_cpp(int n, NumericVector kappa, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n_k = kappa.size();
  NumericVector o_mm(n), o_mk(n), o_kk(n);
  const bool ks = n_k == 1;
  double a0[2] = {0.0, 0.0};
  if (ks) B.upto(kappa[0], 1, a0);
  d7::par_for(n, threads, ks ? d7::kMinTiny : d7::kMinCostly, [&](std::size_t i) {
    const double k = kappa[i % n_k];
    double a[2];
    if (ks) { a[0] = a0[0]; a[1] = a0[1]; } else B.upto(k, 1, a);
    o_mm[i] = d7::vonmises1_expected_mu_mu(k, a[0]);
    o_mk[i] = 0.0;
    o_kk[i] = d7::vonmises1_expected_kappa_kappa(a[1]);
  });
  return List::create(Named("mu_mu") = o_mm, Named("mu_kappa") = o_mk,
                      Named("kappa_kappa") = o_kk);
}

// [[Rcpp::export]]
List vonmises1_deriv3_cpp(NumericVector y, NumericVector mu,
                          NumericVector kappa, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n_y = y.size(), n_mu = mu.size(), n_k = kappa.size();
  const int n = std::max(n_y, std::max(n_mu, n_k));
  NumericVector o_mmm(n), o_mmk(n), o_mkk(n), o_kkk(n);
  const bool ks = n_k == 1;
  const double D0 = ks ? B.d2(kappa[0]) : 0.0;
  d7::par_for(n, threads, ks ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const double k = kappa[i % n_k], d = y[i % n_y] - mu[i % n_mu];
    o_mmm[i] = -k * std::sin(d);
    o_mmk[i] = -std::cos(d);
    o_mkk[i] = 0.0;
    o_kkk[i] = -(ks ? D0 : B.d2(k));
  });
  return List::create(Named("mu_mu_mu") = o_mmm, Named("mu_mu_kappa") = o_mmk,
                      Named("mu_kappa_kappa") = o_mkk,
                      Named("kappa_kappa_kappa") = o_kkk);
}

// [[Rcpp::export]]
List vonmises1_deriv4_cpp(NumericVector y, NumericVector mu,
                          NumericVector kappa, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n_y = y.size(), n_mu = mu.size(), n_k = kappa.size();
  const int n = std::max(n_y, std::max(n_mu, n_k));
  NumericVector o_mmmm(n), o_mmmk(n), o_mmkk(n), o_mkkk(n), o_kkkk(n);
  const bool ks = n_k == 1;
  const double D0 = ks ? B.d3(kappa[0]) : 0.0;
  d7::par_for(n, threads, ks ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const double k = kappa[i % n_k], d = y[i % n_y] - mu[i % n_mu];
    o_mmmm[i] = k * std::cos(d);
    o_mmmk[i] = -std::sin(d);
    o_mmkk[i] = 0.0;
    o_mkkk[i] = 0.0;
    o_kkkk[i] = -(ks ? D0 : B.d3(k));
  });
  return List::create(Named("mu_mu_mu_mu") = o_mmmm,
                      Named("mu_mu_mu_kappa") = o_mmmk,
                      Named("mu_mu_kappa_kappa") = o_mmkk,
                      Named("mu_kappa_kappa_kappa") = o_mkkk,
                      Named("kappa_kappa_kappa_kappa") = o_kkkk);
}

// First derivatives of the expected information in kappa:
// E_mm = -kappa A gives -(A + kappa A'), E_kk = -A' gives -A''.
// [[Rcpp::export]]
List vonmises1_dexpected1_cpp(NumericVector kappa, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n = kappa.size();
  NumericVector f(n), g(n);
  d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
    const double k = kappa[i];
    double a[3];
    B.upto(k, 2, a);
    f[i] = -(a[0] + k * a[1]);
    g[i] = d7::vonmises1_dexpected_kappa_kappa_kappa(a[2]);
  });
  return List::create(Named("f") = f, Named("g") = g);
}

// Second derivatives: -(2 A' + kappa A'') and -A'''.
// [[Rcpp::export]]
List vonmises1_dexpected2_cpp(NumericVector kappa, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n = kappa.size();
  NumericVector f(n), g(n);
  d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
    const double k = kappa[i];
    double a[4];
    B.upto(k, 3, a);
    f[i] = -(2.0 * a[1] + k * a[2]);
    g[i] = -a[3];
  });
  return List::create(Named("f") = f, Named("g") = g);
}

// --- vonmises2 -------------------------------------------------------------

// [[Rcpp::export]]
List vonmises2_gradient_cpp(NumericVector y, NumericVector mu,
                            NumericVector rho, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n_y = y.size(), n_mu = mu.size(), n_r = rho.size();
  const int n = std::max(n_y, std::max(n_mu, n_r));
  NumericVector o_mu(n), o_rho(n);
  const bool rs = n_r == 1;
  double k0[2] = {0.0, 0.0};
  if (rs) vm2_kappa(B, rho[0], 1, k0);
  d7::par_for(n, threads, rs ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const double r = rho[i % n_r], d = y[i % n_y] - mu[i % n_mu];
    double k[2];
    if (rs) { k[0] = k0[0]; k[1] = k0[1]; } else vm2_kappa(B, r, 1, k);
    o_mu[i] = d7::vonmises2_score_mu(k[0], std::sin(d));
    o_rho[i] = d7::vonmises2_score_rho(std::cos(d), r, k[1]);
  });
  return List::create(Named("mu") = o_mu, Named("rho") = o_rho);
}

// [[Rcpp::export]]
List vonmises2_hessian_cpp(NumericVector y, NumericVector mu,
                           NumericVector rho, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n_y = y.size(), n_mu = mu.size(), n_r = rho.size();
  const int n = std::max(n_y, std::max(n_mu, n_r));
  NumericVector o_mm(n), o_rr(n), o_mr(n);
  const bool rs = n_r == 1;
  double k0[3] = {0.0, 0.0, 0.0};
  if (rs) vm2_kappa(B, rho[0], 2, k0);
  d7::par_for(n, threads, rs ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const double r = rho[i % n_r], d = y[i % n_y] - mu[i % n_mu];
    double k[3];
    if (rs) { k[0] = k0[0]; k[1] = k0[1]; k[2] = k0[2]; } else vm2_kappa(B, r, 2, k);
    const double c = std::cos(d);
    o_mm[i] = d7::vonmises2_hess_mu_mu(k[0], c);
    o_rr[i] = d7::vonmises2_hess_rho_rho(c, r, k[1], k[2]);
    o_mr[i] = std::sin(d) * k[1];
  });
  return List::create(Named("mu_mu") = o_mm, Named("rho_rho") = o_rr,
                      Named("mu_rho") = o_mr);
}

// [[Rcpp::export]]
List vonmises2_expected_hessian_cpp(int n, NumericVector rho, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n_r = rho.size();
  NumericVector o_mm(n), o_rr(n), o_mr(n);
  const bool rs = n_r == 1;
  double k0[2] = {0.0, 0.0};
  if (rs) vm2_kappa(B, rho[0], 1, k0);
  d7::par_for(n, threads, rs ? d7::kMinTiny : d7::kMinCostly, [&](std::size_t i) {
    const double r = rho[i % n_r];
    double k[2];
    if (rs) { k[0] = k0[0]; k[1] = k0[1]; } else vm2_kappa(B, r, 1, k);
    o_mm[i] = d7::vonmises2_expected_mu_mu(k[0], r);
    o_rr[i] = d7::vonmises2_expected_rho_rho(k[1]);
    o_mr[i] = 0.0;
  });
  return List::create(Named("mu_mu") = o_mm, Named("rho_rho") = o_rr,
                      Named("mu_rho") = o_mr);
}

// [[Rcpp::export]]
List vonmises2_deriv3_cpp(NumericVector y, NumericVector mu,
                          NumericVector rho, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n_y = y.size(), n_mu = mu.size(), n_r = rho.size();
  const int n = std::max(n_y, std::max(n_mu, n_r));
  NumericVector o_mmm(n), o_mmr(n), o_mrr(n), o_rrr(n);
  const bool rs = n_r == 1;
  double k0[4] = {0.0, 0.0, 0.0, 0.0};
  if (rs) vm2_kappa(B, rho[0], 3, k0);
  d7::par_for(n, threads, rs ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const double r = rho[i % n_r], d = y[i % n_y] - mu[i % n_mu];
    double k[4];
    if (rs) { for (int j = 0; j < 4; ++j) k[j] = k0[j]; } else vm2_kappa(B, r, 3, k);
    const double c = std::cos(d), s = std::sin(d);
    o_mmm[i] = -k[0] * s;
    o_mmr[i] = -c * k[1];
    o_mrr[i] = s * k[2];
    o_rrr[i] = (c - r) * k[3] - 2.0 * k[2];
  });
  return List::create(Named("mu_mu_mu") = o_mmm, Named("mu_mu_rho") = o_mmr,
                      Named("mu_rho_rho") = o_mrr, Named("rho_rho_rho") = o_rrr);
}

// [[Rcpp::export]]
List vonmises2_deriv4_cpp(NumericVector y, NumericVector mu,
                          NumericVector rho, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n_y = y.size(), n_mu = mu.size(), n_r = rho.size();
  const int n = std::max(n_y, std::max(n_mu, n_r));
  NumericVector o_mmmm(n), o_mmmr(n), o_mmrr(n), o_mrrr(n), o_rrrr(n);
  const bool rs = n_r == 1;
  double k0[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  if (rs) vm2_kappa(B, rho[0], 4, k0);
  d7::par_for(n, threads, rs ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    const double r = rho[i % n_r], d = y[i % n_y] - mu[i % n_mu];
    double k[5];
    if (rs) { for (int j = 0; j < 5; ++j) k[j] = k0[j]; } else vm2_kappa(B, r, 4, k);
    const double c = std::cos(d), s = std::sin(d);
    o_mmmm[i] = k[0] * c;
    o_mmmr[i] = -s * k[1];
    o_mmrr[i] = -c * k[2];
    o_mrrr[i] = s * k[3];
    o_rrrr[i] = (c - r) * k[4] - 3.0 * k[3];
  });
  return List::create(Named("mu_mu_mu_mu") = o_mmmm,
                      Named("mu_mu_mu_rho") = o_mmmr,
                      Named("mu_mu_rho_rho") = o_mmrr,
                      Named("mu_rho_rho_rho") = o_mrrr,
                      Named("rho_rho_rho_rho") = o_rrrr);
}

// First derivatives of the expected information in rho:
// E_mm = -rho kappa gives -(kappa + rho kappa'), E_rr = -kappa' gives -kappa''.
// [[Rcpp::export]]
List vonmises2_dexpected1_cpp(NumericVector rho, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n = rho.size();
  NumericVector f(n), g(n);
  d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
    const double r = rho[i];
    double k[3];
    vm2_kappa(B, r, 2, k);
    f[i] = -(k[0] + r * k[1]);
    g[i] = d7::vonmises2_dexpected_rho_rho_rho(k[2]);
  });
  return List::create(Named("f") = f, Named("g") = g);
}

// Second derivatives: -(2 kappa' + rho kappa'') and -kappa'''.
// [[Rcpp::export]]
List vonmises2_dexpected2_cpp(NumericVector rho, int threads = 1) {
  const N7Bessel& B = n7b();
  const int n = rho.size();
  NumericVector f(n), g(n);
  d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
    const double r = rho[i];
    double k[4];
    vm2_kappa(B, r, 3, k);
    f[i] = -(2.0 * k[1] + r * k[2]);
    g[i] = -k[3];
  });
  return List::create(Named("f") = f, Named("g") = g);
}
