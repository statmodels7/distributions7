// The third and fourth derivatives of the expected information of the
// exponential in its mean mu, on the parameter scale. The expected information
// is E = -1/mu^2, so the n-th derivative is (-1)^(n+1) (n+1)! mu^(-2-n):
// 24/mu^5 at order 3 and -120/mu^6 at order 4.

#include <Rcpp.h>
#include "d7_par.h"
#include "d34e_run.h"
using namespace Rcpp;

// [[Rcpp::export]]
List exponential_dexpected3_cpp(NumericVector y, NumericVector mu,
                                int threads = 1) {
    d34e::Par M(mu);
    return d34e::run(y.size(), threads, d7::kMinCheap, 1,
                     [&](std::size_t i, d34e::Row o) {
        const double u = 1.0 / M[i];
        o[0] = 24.0 * u * u * u * u * u;
    });
}

// [[Rcpp::export]]
List exponential_dexpected4_cpp(NumericVector y, NumericVector mu,
                                int threads = 1) {
    d34e::Par M(mu);
    return d34e::run(y.size(), threads, d7::kMinCheap, 1,
                     [&](std::size_t i, d34e::Row o) {
        const double u = 1.0 / M[i];
        o[0] = -120.0 * u * u * u * u * u * u;
    });
}
