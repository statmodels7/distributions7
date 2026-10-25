// The third and fourth derivatives of the expected information of the laplace2
// family on the parameter scale, in closed form. The order-2 kernel of
// distributions7 states E_ab in closed form; the kernels below differentiate
// that closed form, one order per kernel, and write each component in the
// order of d3expected_names() and d4expected_names() for laplace2.

#include <Rcpp.h>
#include "d7_par.h"
#include "d34e_run.h"
using namespace Rcpp;

// [[Rcpp::export]]
List laplace2_dexpected3_cpp(NumericVector y, NumericVector mu, NumericVector lambda, int threads = 1) {
    d34e::Par LAM(lambda);
    return d34e::run(y.size(), threads, d7::kMinCheap, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double lam = LAM[i];
        o[0] = 0;   // mu_mu_mu_mu_mu
        o[1] = 0;   // mu_mu_mu_mu_lambda
        o[2] = 0;   // mu_mu_mu_lambda_lambda
        o[3] = 0;   // mu_mu_lambda_lambda_lambda
        o[4] = 0;   // lambda_lambda_mu_mu_mu
        o[5] = 0;   // lambda_lambda_mu_mu_lambda
        o[6] = 0;   // lambda_lambda_mu_lambda_lambda
        o[7] = 24/pow(lam, 5);   // lambda_lambda_lambda_lambda_lambda
        o[8] = 0;   // mu_lambda_mu_mu_mu
        o[9] = 0;   // mu_lambda_mu_mu_lambda
        o[10] = 0;   // mu_lambda_mu_lambda_lambda
        o[11] = 0;   // mu_lambda_lambda_lambda_lambda
    });
}

// [[Rcpp::export]]
List laplace2_dexpected4_cpp(NumericVector y, NumericVector mu, NumericVector lambda, int threads = 1) {
    d34e::Par LAM(lambda);
    return d34e::run(y.size(), threads, d7::kMinCheap, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double lam = LAM[i];
        o[0] = 0;   // mu_mu_mu_mu_mu_mu
        o[1] = 0;   // mu_mu_mu_mu_mu_lambda
        o[2] = 0;   // mu_mu_mu_mu_lambda_lambda
        o[3] = 0;   // mu_mu_mu_lambda_lambda_lambda
        o[4] = 0;   // mu_mu_lambda_lambda_lambda_lambda
        o[5] = 0;   // lambda_lambda_mu_mu_mu_mu
        o[6] = 0;   // lambda_lambda_mu_mu_mu_lambda
        o[7] = 0;   // lambda_lambda_mu_mu_lambda_lambda
        o[8] = 0;   // lambda_lambda_mu_lambda_lambda_lambda
        o[9] = -120/pow(lam, 6);   // lambda_lambda_lambda_lambda_lambda_lambda
        o[10] = 0;   // mu_lambda_mu_mu_mu_mu
        o[11] = 0;   // mu_lambda_mu_mu_mu_lambda
        o[12] = 0;   // mu_lambda_mu_mu_lambda_lambda
        o[13] = 0;   // mu_lambda_mu_lambda_lambda_lambda
        o[14] = 0;   // mu_lambda_lambda_lambda_lambda_lambda
    });
}
