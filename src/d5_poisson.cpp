#include <Rcpp.h>
#include "d7_par.h"
#include "pt_discrete_pow.h"
using namespace Rcpp;

// Fifth derivative of the Poisson log-mass in its mean mu. The log-mass
// y log mu - mu - log y! has no mu-dependence beyond the first two terms, so
// for k >= 2 its k-th derivative is (-1)^k (k-1)! y / mu^k, and at k = 5
//   d^5 l / d mu^5 = 24 y / mu^5,
// written as coef_over_pow (c / d^k by repeated division).

// [[Rcpp::export]]
List poisson_deriv5_cpp(NumericVector y, NumericVector mu, int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        mu_mu_mu_mu_mu[i] = 24.0 * d7::coef_over_pow(y[i], m, 5);
    });
    return List::create(Named("mu_mu_mu_mu_mu") = mu_mu_mu_mu_mu);
}
