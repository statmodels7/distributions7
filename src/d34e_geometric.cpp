// The third and fourth derivatives of the expected information of the
// geometric in its mean mu, on the parameter scale. The expected Hessian is
// E = -1/r with r = mu (1 + mu). With t = 2 mu + 1 the third derivative is
// 6 t (1 + 2 r) / r^4 and the fourth is -24 (1 + 5 r + 5 r^2) / r^5.

#include <Rcpp.h>
#include "d7_par.h"
#include "d34e_run.h"
using namespace Rcpp;

// [[Rcpp::export]]
List geometric_dexpected3_cpp(NumericVector y, NumericVector mu,
                              int threads = 1) {
    d34e::Par M(mu);
    return d34e::run(y.size(), threads, d7::kMinCheap, 1,
                     [&](std::size_t i, d34e::Row o) {
        const double m = M[i];
        const double r = m * (1.0 + m);
        const double t = 2.0 * m + 1.0;
        const double r2 = r * r;
        o[0] = 6.0 * t * (1.0 + 2.0 * r) / (r2 * r2);
    });
}

// [[Rcpp::export]]
List geometric_dexpected4_cpp(NumericVector y, NumericVector mu,
                              int threads = 1) {
    d34e::Par M(mu);
    return d34e::run(y.size(), threads, d7::kMinCheap, 1,
                     [&](std::size_t i, d34e::Row o) {
        const double m = M[i];
        const double r = m * (1.0 + m);
        const double r2 = r * r;
        o[0] = -24.0 * (1.0 + 5.0 * r + 5.0 * r2) / (r2 * r2 * r);
    });
}
