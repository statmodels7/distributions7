// The third and fourth derivatives of the expected information of the
// Bernoulli in its probability mu, on the parameter scale. The expected
// Hessian is E = -1/q with q = mu (1 - mu). With s = 2 mu - 1 the third
// derivative is -6 s (1 - 2 q) / q^4 and the fourth is -24 (1 - 5 q + 5 q^2) /
// q^5, both obtained by differentiating the partial fractions of 1/q.

#include <Rcpp.h>
#include "d7_par.h"
#include "d34e_run.h"
using namespace Rcpp;

// [[Rcpp::export]]
List bernoulli_dexpected3_cpp(NumericVector y, NumericVector mu,
                              int threads = 1) {
    d34e::Par M(mu);
    return d34e::run(y.size(), threads, d7::kMinCheap, 1,
                     [&](std::size_t i, d34e::Row o) {
        const double m = M[i];
        const double q = m * (1.0 - m);
        const double s = 2.0 * m - 1.0;
        const double q2 = q * q;
        o[0] = -6.0 * s * (1.0 - 2.0 * q) / (q2 * q2);
    });
}

// [[Rcpp::export]]
List bernoulli_dexpected4_cpp(NumericVector y, NumericVector mu,
                              int threads = 1) {
    d34e::Par M(mu);
    return d34e::run(y.size(), threads, d7::kMinCheap, 1,
                     [&](std::size_t i, d34e::Row o) {
        const double m = M[i];
        const double q = m * (1.0 - m);
        const double q2 = q * q;
        o[0] = -24.0 * (1.0 - 5.0 * q + 5.0 * q2) / (q2 * q2 * q);
    });
}
