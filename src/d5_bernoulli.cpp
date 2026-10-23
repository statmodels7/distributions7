#include <Rcpp.h>
#include "d7_par.h"
#include "pt_discrete_pow.h"
using namespace Rcpp;

// Fifth derivative of the Bernoulli log-mass in its mean mu, the closed form
// of the fourth-order kernels of bernoulli_hd.cpp taken one order further.
// With om = mu - 1,
//   d^5 l / d mu^5 = 24 y / mu^5 - 24 (y - 1) / om^5,
// each term written as coef_over_pow (c / d^k by repeated division).

// [[Rcpp::export]]
List bernoulli_deriv5_cpp(NumericVector y, NumericVector mu,
                          int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double om = m - 1.0;
        mu_mu_mu_mu_mu[i] = 24.0 * d7::coef_over_pow(y[i], m, 5)
                          - 24.0 * d7::coef_over_pow(y[i] - 1.0, om, 5);
    });
    return List::create(Named("mu_mu_mu_mu_mu") = mu_mu_mu_mu_mu);
}
