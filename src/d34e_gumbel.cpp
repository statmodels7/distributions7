// The third and fourth derivatives of the expected information of the gumbel
// family on the parameter scale, in closed form. The order-2 kernel of
// distributions7 states E_ab in closed form; the kernels below differentiate
// that closed form, one order per kernel, and write each component in the
// order of d3expected_names() and d4expected_names() for gumbel.

#include <Rcpp.h>
#include "d7_par.h"
#include "pt_constants.h"
#include "d34e_run.h"
using namespace Rcpp;

namespace {
// Euler's constant, the only constant that is not elementary
const double kEG = d7::kEulerGamma;
}  // namespace

// [[Rcpp::export]]
List gumbel_dexpected3_cpp(NumericVector y, NumericVector mu, NumericVector sigma, int threads = 1) {
    d34e::Par S(sigma);
    return d34e::run(y.size(), threads, d7::kMinCheap, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double s = S[i];
        o[0] = 0;   // mu_mu_mu_mu_mu
        o[1] = 0;   // mu_mu_mu_mu_sigma
        o[2] = 0;   // mu_mu_mu_sigma_sigma
        o[3] = 24/pow(s, 5);   // mu_mu_sigma_sigma_sigma
        o[4] = 0;   // sigma_sigma_mu_mu_mu
        o[5] = 0;   // sigma_sigma_mu_mu_sigma
        o[6] = 0;   // sigma_sigma_mu_sigma_sigma
        o[7] = 4*(6*pow(kEG, 2) - 12*kEG + 6 + pow(M_PI, 2))/pow(s, 5);   // sigma_sigma_sigma_sigma_sigma
        o[8] = 0;   // mu_sigma_mu_mu_mu
        o[9] = 0;   // mu_sigma_mu_mu_sigma
        o[10] = 0;   // mu_sigma_mu_sigma_sigma
        o[11] = 24*(kEG - 1)/pow(s, 5);   // mu_sigma_sigma_sigma_sigma
    });
}

// [[Rcpp::export]]
List gumbel_dexpected4_cpp(NumericVector y, NumericVector mu, NumericVector sigma, int threads = 1) {
    d34e::Par S(sigma);
    return d34e::run(y.size(), threads, d7::kMinCheap, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double s = S[i];
        o[0] = 0;   // mu_mu_mu_mu_mu_mu
        o[1] = 0;   // mu_mu_mu_mu_mu_sigma
        o[2] = 0;   // mu_mu_mu_mu_sigma_sigma
        o[3] = 0;   // mu_mu_mu_sigma_sigma_sigma
        o[4] = -120/pow(s, 6);   // mu_mu_sigma_sigma_sigma_sigma
        o[5] = 0;   // sigma_sigma_mu_mu_mu_mu
        o[6] = 0;   // sigma_sigma_mu_mu_mu_sigma
        o[7] = 0;   // sigma_sigma_mu_mu_sigma_sigma
        o[8] = 0;   // sigma_sigma_mu_sigma_sigma_sigma
        o[9] = -20*(6*pow(kEG, 2) - 12*kEG + 6 + pow(M_PI, 2))/pow(s, 6);   // sigma_sigma_sigma_sigma_sigma_sigma
        o[10] = 0;   // mu_sigma_mu_mu_mu_mu
        o[11] = 0;   // mu_sigma_mu_mu_mu_sigma
        o[12] = 0;   // mu_sigma_mu_mu_sigma_sigma
        o[13] = 0;   // mu_sigma_mu_sigma_sigma_sigma
        o[14] = -120*(kEG - 1)/pow(s, 6);   // mu_sigma_sigma_sigma_sigma_sigma
    });
}
