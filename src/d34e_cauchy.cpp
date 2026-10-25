// The third and fourth derivatives of the expected information of the cauchy
// family on the parameter scale, in closed form. The order-2 kernel of
// distributions7 states E_ab in closed form; the kernels below differentiate
// that closed form, one order per kernel, and write each component in the
// order of d3expected_names() and d4expected_names() for cauchy.

#include <Rcpp.h>
#include "d7_par.h"
#include "pt_constants.h"
#include "d34e_run.h"
using namespace Rcpp;

// [[Rcpp::export]]
List cauchy_dexpected3_cpp(NumericVector y, NumericVector mu, NumericVector sigma, int threads = 1) {
    d34e::Par S(sigma);
    return d34e::run(y.size(), threads, d7::kMinCheap, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double s = S[i];
        o[0] = 0;   // mu_mu_mu_mu_mu
        o[1] = 0;   // mu_mu_mu_mu_sigma
        o[2] = 0;   // mu_mu_mu_sigma_sigma
        o[3] = 12/pow(s, 5);   // mu_mu_sigma_sigma_sigma
        o[4] = 0;   // sigma_sigma_mu_mu_mu
        o[5] = 0;   // sigma_sigma_mu_mu_sigma
        o[6] = 0;   // sigma_sigma_mu_sigma_sigma
        o[7] = 12/pow(s, 5);   // sigma_sigma_sigma_sigma_sigma
        o[8] = 0;   // mu_sigma_mu_mu_mu
        o[9] = 0;   // mu_sigma_mu_mu_sigma
        o[10] = 0;   // mu_sigma_mu_sigma_sigma
        o[11] = 0;   // mu_sigma_sigma_sigma_sigma
    });
}

// [[Rcpp::export]]
List cauchy_dexpected4_cpp(NumericVector y, NumericVector mu, NumericVector sigma, int threads = 1) {
    d34e::Par S(sigma);
    return d34e::run(y.size(), threads, d7::kMinCheap, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double s = S[i];
        o[0] = 0;   // mu_mu_mu_mu_mu_mu
        o[1] = 0;   // mu_mu_mu_mu_mu_sigma
        o[2] = 0;   // mu_mu_mu_mu_sigma_sigma
        o[3] = 0;   // mu_mu_mu_sigma_sigma_sigma
        o[4] = -60/pow(s, 6);   // mu_mu_sigma_sigma_sigma_sigma
        o[5] = 0;   // sigma_sigma_mu_mu_mu_mu
        o[6] = 0;   // sigma_sigma_mu_mu_mu_sigma
        o[7] = 0;   // sigma_sigma_mu_mu_sigma_sigma
        o[8] = 0;   // sigma_sigma_mu_sigma_sigma_sigma
        o[9] = -60/pow(s, 6);   // sigma_sigma_sigma_sigma_sigma_sigma
        o[10] = 0;   // mu_sigma_mu_mu_mu_mu
        o[11] = 0;   // mu_sigma_mu_mu_mu_sigma
        o[12] = 0;   // mu_sigma_mu_mu_sigma_sigma
        o[13] = 0;   // mu_sigma_mu_sigma_sigma_sigma
        o[14] = 0;   // mu_sigma_sigma_sigma_sigma_sigma
    });
}
