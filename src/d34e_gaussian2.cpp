// The third and fourth derivatives of the expected information of gaussian2 on
// the parameter scale, in closed form, one kernel per order. The closed
// form of E_ab is differentiated symbolically (misure_d34e/agents/gaussian2/
// d34e_elem.py) and emitted with common subexpressions. Components are in
// the order of d3expected_names() and d4expected_names().

#include <Rcpp.h>
#include "d34e_run.h"
using namespace Rcpp;

// [[Rcpp::export]]
List gaussian2_dexpected3_cpp(NumericVector y, NumericVector mu, NumericVector sigma2, int threads = 1) {
    d34e::Par P_mu(mu), P_sigma2(sigma2);
    return d34e::run(y.size(), threads, d7::kMinCheap, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double mu = P_mu[i];
        const double sigma2 = P_sigma2[i];
        o[0] = 0;
        o[1] = 0;
        o[2] = 0;
        o[3] = 6/std::pow(sigma2, 4);
        o[4] = 0;
        o[5] = 0;
        o[6] = 0;
        o[7] = 12/std::pow(sigma2, 5);
        o[8] = 0;
        o[9] = 0;
        o[10] = 0;
        o[11] = 0;
    });
}

// [[Rcpp::export]]
List gaussian2_dexpected4_cpp(NumericVector y, NumericVector mu, NumericVector sigma2, int threads = 1) {
    d34e::Par P_mu(mu), P_sigma2(sigma2);
    return d34e::run(y.size(), threads, d7::kMinCheap, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double mu = P_mu[i];
        const double sigma2 = P_sigma2[i];
        o[0] = 0;
        o[1] = 0;
        o[2] = 0;
        o[3] = 0;
        o[4] = -24/std::pow(sigma2, 5);
        o[5] = 0;
        o[6] = 0;
        o[7] = 0;
        o[8] = 0;
        o[9] = -60/std::pow(sigma2, 6);
        o[10] = 0;
        o[11] = 0;
        o[12] = 0;
        o[13] = 0;
        o[14] = 0;
    });
}
