#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// Third/fourth-order derivatives of the Bernoulli log-mass, transcribed from the
// Wolfram output. Single parameter mu; om = mu - 1.

// [[Rcpp::export]]
List bernoulli_deriv3_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double om = m - 1.0;
        mu_mu_mu[i] = -2.0 * (y[i] - 1.0) / (om * om * om) + 2.0 * y[i] / (m * m * m);
    });
    return List::create(Named("mu_mu_mu") = mu_mu_mu);
}

// [[Rcpp::export]]
List bernoulli_deriv3_expected_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double om = m - 1.0;
        mu_mu_mu[i] = -2.0 / (om * om) + 2.0 / (m * m);
    });
    return List::create(Named("mu_mu_mu") = mu_mu_mu);
}

// [[Rcpp::export]]
List bernoulli_deriv4_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double om = m - 1.0;
        mu_mu_mu_mu[i] = 6.0 * (y[i] - 1.0) / (om * om * om * om) - 6.0 * y[i] / (m * m * m * m);
    });
    return List::create(Named("mu_mu_mu_mu") = mu_mu_mu_mu);
}

// [[Rcpp::export]]
List bernoulli_deriv4_expected_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double om = m - 1.0;
        mu_mu_mu_mu[i] = 6.0 / (om * om * om) - 6.0 / (m * m * m);
    });
    return List::create(Named("mu_mu_mu_mu") = mu_mu_mu_mu);
}
