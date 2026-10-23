#include <Rcpp.h>
#include "d7_par.h"
#include "pt_discrete_pow.h"
using namespace Rcpp;

// Fifth derivative of the geometric log-mass in its mean mu, on the support
// 0, 1, ... The log-mass is y log mu - (y + 1) log(1 + mu), and its k-th
// derivative for k >= 2 is (-1)^(k-1) (k-1)! [y/mu^k - (y + 1)/(1 + mu)^k].
// At k = 5 the sign is positive, and the bracket is count_bracket with c = 1,
// which carries 1 - q^k through expm1 and log1p (pt_discrete_pow.h):
//   d^5 l / d mu^5 = 24 count_bracket(y, mu, 1, 5).

// [[Rcpp::export]]
List geometric_deriv5_cpp(NumericVector y, NumericVector mu, int threads = 1) {
    int n = y.size();
    NumericVector d(n);
    bool mu_is_scalar = (mu.size() == 1);
    d7::par_for(n, threads, d7::kMinTiny, [&](std::size_t i) {
        // LOCAL to the region: a scalar hoisted out of the loop and
        // written inside it is shared once the iterations are split
        const double m = mu_is_scalar ? mu[0] : mu[i];
        d[i] = 24.0 * d7::count_bracket(y[i], m, 1.0, 5);
    });
    return List::create(Named("mu_mu_mu_mu_mu") = d);
}
