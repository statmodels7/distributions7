#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// Chi-squared in the MEAN parametrization: the mean IS the degrees of
// freedom, so mu = nu and
//   l = (mu/2 - 1) log y - y/2 - (mu/2) log 2 - lgamma(mu/2).
// Only the score carries the data:
//   l' = (log y - log 2 - psi(mu/2)) / 2,
// and from the second order on the derivatives do not involve y at all,
//   l^(k) = -psi^(k-2)(mu/2) / 2^k,
// because the family is a one-parameter exponential family in log y. The
// observed information is therefore EXACTLY the expected information, and
// the same holds at third and fourth order -- there is nothing to average.
// E[log y] = psi(mu/2) + log 2 is what makes the score mean zero.

// [[Rcpp::export]]
List chisq_gradient_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector g(n);
    bool mu_is_scalar = (mu.size() == 1);
    const double log2 = std::log(2.0);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and
        // written inside it is shared once the iterations are split
        const double m = mu_is_scalar ? mu[0] : mu[i];
        g[i] = 0.5 * (std::log(y[i]) - log2 - R::digamma(0.5 * m));
    });
    return List::create(Named("mu") = g);
}

// [[Rcpp::export]]
List chisq_hessian_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector h(n);
    bool mu_is_scalar = (mu.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and
        // written inside it is shared once the iterations are split
        const double m = mu_is_scalar ? mu[0] : mu[i];
        h[i] = -0.25 * R::trigamma(0.5 * m);
    });
    return List::create(Named("mu_mu") = h);
}

// [[Rcpp::export]]
List chisq_deriv3_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector d(n);
    bool mu_is_scalar = (mu.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and
        // written inside it is shared once the iterations are split
        const double m = mu_is_scalar ? mu[0] : mu[i];
        d[i] = -R::psigamma(0.5 * m, 2) / 8.0;
    });
    return List::create(Named("mu_mu_mu") = d);
}

// [[Rcpp::export]]
List chisq_deriv4_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector d(n);
    bool mu_is_scalar = (mu.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and
        // written inside it is shared once the iterations are split
        const double m = mu_is_scalar ? mu[0] : mu[i];
        d[i] = -R::psigamma(0.5 * m, 3) / 16.0;
    });
    return List::create(Named("mu_mu_mu_mu") = d);
}
