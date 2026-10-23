#include <Rcpp.h>
#include "d7_par.h"
#include "pt_discrete_pow.h"
using namespace Rcpp;

// Fifth derivative of the Binomial log-mass in its success probability mu,
// the size N fixed. With om = mu - 1 and the log-mass
// l = y log mu + (N - y) log(1 - mu) + const,
//   d^5 l / d mu^5 = 24 y / mu^5 + 24 (N - y) / om^5,
// each term written as coef_over_pow (c / d^k by repeated division).

// [[Rcpp::export]]
List binomial_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector size,
                         int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool size_is_scalar = (size.size() == 1);
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double sz = size_is_scalar ? size[0] : size[i];
        double om = m - 1.0;
        mu_mu_mu_mu_mu[i] = 24.0 * d7::coef_over_pow(y[i], m, 5)
                          + 24.0 * d7::coef_over_pow(sz - y[i], om, 5);
    });
    return List::create(Named("mu_mu_mu_mu_mu") = mu_mu_mu_mu_mu);
}
