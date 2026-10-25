// The third and fourth derivatives of the expected information of the
// chi-squared in its degrees of freedom mu, on the parameter scale. The
// expected Hessian is E = -psi'(mu/2)/4, where psi' is the trigamma function,
// so with x = mu/2 the third derivative is -psi^(4)(x)/32 and the fourth is
// -psi^(5)(x)/64. Each is one polygamma call of R, the same function the
// order-1 and order-2 kernels use.

#include <Rcpp.h>
#include "d7_par.h"
#include "d34e_run.h"
using namespace Rcpp;

// [[Rcpp::export]]
List chisq_dexpected3_cpp(NumericVector y, NumericVector mu,
                          int threads = 1) {
    d34e::Par M(mu);
    return d34e::run(y.size(), threads, d7::kMinCostly, 1,
                     [&](std::size_t i, d34e::Row o) {
        o[0] = -R::psigamma(0.5 * M[i], 4) / 32.0;
    });
}

// [[Rcpp::export]]
List chisq_dexpected4_cpp(NumericVector y, NumericVector mu,
                          int threads = 1) {
    d34e::Par M(mu);
    return d34e::run(y.size(), threads, d7::kMinCostly, 1,
                     [&](std::size_t i, d34e::Row o) {
        o[0] = -R::psigamma(0.5 * M[i], 5) / 64.0;
    });
}
