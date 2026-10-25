// The third and fourth derivatives of the expected information of invgauss1 on
// the parameter scale, in closed form, one kernel per order. The closed
// form of E_ab is differentiated symbolically (misure_d34e/agents/invgauss1/
// d34e_elem.py) and emitted with common subexpressions. Components are in
// the order of d3expected_names() and d4expected_names().

#include <Rcpp.h>
#include "d34e_run.h"
using namespace Rcpp;

// [[Rcpp::export]]
List invgauss1_dexpected3_cpp(NumericVector y, NumericVector mu, NumericVector phi, int threads = 1) {
    d34e::Par P_mu(mu), P_phi(phi);
    return d34e::run(y.size(), threads, d7::kMinCheap, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double mu = P_mu[i];
        const double phi = P_phi[i];
        o[0] = 60/(std::pow(mu, 6)*phi);
        o[1] = 12/(std::pow(mu, 5)*std::pow(phi, 2));
        o[2] = 6/(std::pow(mu, 4)*std::pow(phi, 3));
        o[3] = 6/(std::pow(mu, 3)*std::pow(phi, 4));
        o[4] = 0;
        o[5] = 0;
        o[6] = 0;
        o[7] = 12/std::pow(phi, 5);
        o[8] = 0;
        o[9] = 0;
        o[10] = 0;
        o[11] = 0;
    });
}

// [[Rcpp::export]]
List invgauss1_dexpected4_cpp(NumericVector y, NumericVector mu, NumericVector phi, int threads = 1) {
    d34e::Par P_mu(mu), P_phi(phi);
    return d34e::run(y.size(), threads, d7::kMinCheap, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double mu = P_mu[i];
        const double phi = P_phi[i];
        o[0] = -360/(std::pow(mu, 7)*phi);
        o[1] = -60/(std::pow(mu, 6)*std::pow(phi, 2));
        o[2] = -24/(std::pow(mu, 5)*std::pow(phi, 3));
        o[3] = -18/(std::pow(mu, 4)*std::pow(phi, 4));
        o[4] = -24/(std::pow(mu, 3)*std::pow(phi, 5));
        o[5] = 0;
        o[6] = 0;
        o[7] = 0;
        o[8] = 0;
        o[9] = -60/std::pow(phi, 6);
        o[10] = 0;
        o[11] = 0;
        o[12] = 0;
        o[13] = 0;
        o[14] = 0;
    });
}
