// The third and fourth derivatives of the expected information of the
// Poisson in its mean mu, on the parameter scale. The expected Hessian is
// E = -1/mu, so the third derivative is 6/mu^4 and the fourth is -24/mu^5.

// R evaluates these expressions one operation at a time; clang (arm64
// macOS) would otherwise fuse a*b + c into one rounding
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
#include <Rcpp.h>
#include "d7_par.h"
#include "d34e_run.h"
using namespace Rcpp;

// [[Rcpp::export]]
List poisson_dexpected3_cpp(NumericVector y, NumericVector mu,
                            int threads = 1) {
    d34e::Par M(mu);
    return d34e::run(y.size(), threads, d7::kMinCheap, 1,
                     [&](std::size_t i, d34e::Row o) {
        const double u = 1.0 / M[i];
        o[0] = 6.0 * u * u * u * u;
    });
}

// [[Rcpp::export]]
List poisson_dexpected4_cpp(NumericVector y, NumericVector mu,
                            int threads = 1) {
    d34e::Par M(mu);
    return d34e::run(y.size(), threads, d7::kMinCheap, 1,
                     [&](std::size_t i, d34e::Row o) {
        const double u = 1.0 / M[i];
        o[0] = -24.0 * u * u * u * u * u;
    });
}
