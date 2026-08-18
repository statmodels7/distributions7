#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// Third- and fourth-order derivatives of the Gaussian log-density (sd
// parameterization). Written in z = (y - mu)/sigma and inv = 1/sigma for the
// reason given in gaussian.cpp: a form carrying sigma^5 or sigma^6 in a
// denominator loses the whole component once that power overflows, at
// sigma = 4.4e61 for the fourth order, while the component itself is
// representable far beyond.

// [[Rcpp::export]]
List gaussian_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                          int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double inv = 1.0 / s, inv3 = inv * inv * inv;
        double z = (y[i] - m) / s;

        mu_mu_mu[i] = 0.0;
        mu_mu_sigma[i] = 2.0 * inv3;
        mu_sigma_sigma[i] = 6.0 * z * inv3;
        sigma_sigma_sigma[i] = -2.0 * (1.0 - 6.0 * z * z) * inv3;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List gaussian_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                          int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);

    bool sigma_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double inv = 1.0 / s, inv3 = inv * inv * inv;

        mu_mu_mu[i] = 0.0;
        mu_mu_sigma[i] = 2.0 * inv3;
        mu_sigma_sigma[i] = 0.0;
        sigma_sigma_sigma[i] = 10.0 * inv3;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List gaussian_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                          int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n),
                  mu_sigma_sigma_sigma(n), sigma_sigma_sigma_sigma(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double inv = 1.0 / s, inv2 = inv * inv, inv4 = inv2 * inv2;
        double z = (y[i] - m) / s;

        mu_mu_mu_mu[i] = 0.0;
        mu_mu_mu_sigma[i] = 0.0;
        mu_mu_sigma_sigma[i] = -6.0 * inv4;
        mu_sigma_sigma_sigma[i] = -24.0 * z * inv4;
        sigma_sigma_sigma_sigma[i] = 6.0 * (1.0 - 10.0 * z * z) * inv4;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List gaussian_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                          int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n),
                  mu_sigma_sigma_sigma(n), sigma_sigma_sigma_sigma(n);

    bool sigma_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double inv = 1.0 / s, inv2 = inv * inv, inv4 = inv2 * inv2;

        mu_mu_mu_mu[i] = 0.0;
        mu_mu_mu_sigma[i] = 0.0;
        mu_mu_sigma_sigma[i] = -6.0 * inv4;
        mu_sigma_sigma_sigma[i] = 0.0;
        sigma_sigma_sigma_sigma[i] = -54.0 * inv4;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma
    );
}
