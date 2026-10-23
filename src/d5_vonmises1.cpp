#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "d7_par.h"
#include "pt_vonmises.h"
using namespace Rcpp;

// The fifth derivatives of the vonmises1 log-density in (mu, kappa):
//   l = kappa cos(y - mu) - log(2 pi) - log I0(kappa).
// The mu-derivatives are those of kappa cos(y - mu); kappa enters the
// log-density linearly and log I0 has derivative A = I1/I0, so every
// derivative with two or more kappa indices from the kappa cosine is zero,
// and the pure kappa derivatives are -A^(4)(kappa). A^(j) come from
// numericals7's compiled ratio, to order four (N7Bessel::upto), as in
// vonmises1_deriv4_cpp(). The loop body is written out per component; it is
// derived in misure_gas_congiunto/c/fase4/D5/vonmises1/derive.py.

// [[Rcpp::export]]
List vonmises1_deriv5_cpp(NumericVector y, NumericVector mu,
                          NumericVector kappa, int threads = 1) {
  const d7::N7Bessel& B = d7::vm_bessel();
  const int n_y = y.size(), n_mu = mu.size(), n_k = kappa.size();
  const int n = std::max(n_y, std::max(n_mu, n_k));
  NumericVector o_mu_mu_mu_mu_mu(n), o_mu_mu_mu_mu_kappa(n),
    o_mu_mu_mu_kappa_kappa(n), o_mu_mu_kappa_kappa_kappa(n),
    o_mu_kappa_kappa_kappa_kappa(n), o_kappa_kappa_kappa_kappa_kappa(n);
  const bool ks = n_k == 1;
  double Av0[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  if (ks) B.upto(kappa[0], 4, Av0);
  d7::par_for(n, threads, ks ? d7::kMinMid : d7::kMinCostly, [&](std::size_t i) {
    double Av[5];
    if (ks) {
      for (int j = 0; j < 5; ++j) Av[j] = Av0[j];
    } else {
      B.upto(kappa[i % n_k], 4, Av);
    }
    const double yy = y[i % n_y];
    const double mm = mu[i % n_mu];
    const double kk = kappa[i % n_k];
    const double A4 = Av[4];
    const double t0 = mm - yy;
    o_mu_mu_mu_mu_mu[i] = -kk*std::sin(t0);
    o_mu_mu_mu_mu_kappa[i] = std::cos(t0);
    o_mu_mu_mu_kappa_kappa[i] = 0;
    o_mu_mu_kappa_kappa_kappa[i] = 0;
    o_mu_kappa_kappa_kappa_kappa[i] = 0;
    o_kappa_kappa_kappa_kappa_kappa[i] = -A4;
  });
  return List::create(Named("mu_mu_mu_mu_mu") = o_mu_mu_mu_mu_mu,
                      Named("mu_mu_mu_mu_kappa") = o_mu_mu_mu_mu_kappa,
                      Named("mu_mu_mu_kappa_kappa") = o_mu_mu_mu_kappa_kappa,
                      Named("mu_mu_kappa_kappa_kappa") = o_mu_mu_kappa_kappa_kappa,
                      Named("mu_kappa_kappa_kappa_kappa") = o_mu_kappa_kappa_kappa_kappa,
                      Named("kappa_kappa_kappa_kappa_kappa") = o_kappa_kappa_kappa_kappa_kappa);
}
