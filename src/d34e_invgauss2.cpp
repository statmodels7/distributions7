// The third and fourth derivatives of the expected information of invgauss2 on
// the parameter scale, in closed form, one kernel per order. The closed
// form of E_ab is differentiated symbolically (misure_d34e/agents/invgauss2/
// d34e_elem.py) and emitted with common subexpressions. Components are in
// the order of d3expected_names() and d4expected_names().

#include <Rcpp.h>
#include "d34e_run.h"
using namespace Rcpp;

// [[Rcpp::export]]
List invgauss2_dexpected3_cpp(NumericVector y, NumericVector mu, NumericVector lambda, int threads = 1) {
    d34e::Par P_mu(mu), P_lambda(lambda);
    return d34e::run(y.size(), threads, d7::kMinCheap, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double mu = P_mu[i];
        const double lambda = P_lambda[i];
        o[0] = 60*lambda/std::pow(mu, 6);
        o[1] = -12/std::pow(mu, 5);
        o[2] = 0;
        o[3] = 0;
        o[4] = 0;
        o[5] = 0;
        o[6] = 0;
        o[7] = 12/std::pow(lambda, 5);
        o[8] = 0;
        o[9] = 0;
        o[10] = 0;
        o[11] = 0;
    });
}

// [[Rcpp::export]]
List invgauss2_dexpected4_cpp(NumericVector y, NumericVector mu, NumericVector lambda, int threads = 1) {
    d34e::Par P_mu(mu), P_lambda(lambda);
    return d34e::run(y.size(), threads, d7::kMinCheap, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double mu = P_mu[i];
        const double lambda = P_lambda[i];
        o[0] = -360*lambda/std::pow(mu, 7);
        o[1] = 60/std::pow(mu, 6);
        o[2] = 0;
        o[3] = 0;
        o[4] = 0;
        o[5] = 0;
        o[6] = 0;
        o[7] = 0;
        o[8] = 0;
        o[9] = -60/std::pow(lambda, 6);
        o[10] = 0;
        o[11] = 0;
        o[12] = 0;
        o[13] = 0;
        o[14] = 0;
    });
}
