#include <Rcpp.h>
#include "d7_par.h"
#include <cmath>
using namespace Rcpp;

// Third/fourth-order derivatives of the Lognormal log-density (mean/variance
// on the log scale). The parameter derivatives coincide with those of a Gaussian
// with variance sigma2 evaluated at log(y): r = mu - log(y). Transcribed from the
// Wolfram "gaussian2" output with y -> log(y).

// [[Rcpp::export]]
List lognormal_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma2(n), mu_sigma2_sigma2(n), sigma2_sigma2_sigma2(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool s_is_scalar = (sigma2.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s2 = s_is_scalar ? sigma2[0] : sigma2[i];
        double s2_2 = s2 * s2, s2_3 = s2_2 * s2, s2_4 = s2_2 * s2_2;
        double r = m - std::log(y[i]);

        mu_mu_mu[i] = 0.0;
        mu_mu_sigma2[i] = 1.0 / s2_2;
        mu_sigma2_sigma2[i] = -2.0 * r / s2_3;
        sigma2_sigma2_sigma2[i] = (-s2 + 3.0 * r * r) / s2_4;
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma2") = mu_mu_sigma2,
        Named("mu_sigma2_sigma2") = mu_sigma2_sigma2,
        Named("sigma2_sigma2_sigma2") = sigma2_sigma2_sigma2
    );
}

// [[Rcpp::export]]
List lognormal_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma2(n), mu_sigma2_sigma2(n), sigma2_sigma2_sigma2(n);
    bool s_is_scalar = (sigma2.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double s2 = s_is_scalar ? sigma2[0] : sigma2[i];
        mu_mu_mu[i] = 0.0;
        mu_mu_sigma2[i] = 1.0 / (s2 * s2);
        mu_sigma2_sigma2[i] = 0.0;
        sigma2_sigma2_sigma2[i] = 2.0 / (s2 * s2 * s2);
    });

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma2") = mu_mu_sigma2,
        Named("mu_sigma2_sigma2") = mu_sigma2_sigma2,
        Named("sigma2_sigma2_sigma2") = sigma2_sigma2_sigma2
    );
}

// [[Rcpp::export]]
List lognormal_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma2(n), mu_mu_sigma2_sigma2(n),
                  mu_sigma2_sigma2_sigma2(n), sigma2_sigma2_sigma2_sigma2(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool s_is_scalar = (sigma2.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s2 = s_is_scalar ? sigma2[0] : sigma2[i];
        double s2_3 = s2 * s2 * s2, s2_4 = s2_3 * s2, s2_5 = s2_4 * s2;
        double r = m - std::log(y[i]);

        mu_mu_mu_mu[i] = 0.0;
        mu_mu_mu_sigma2[i] = 0.0;
        mu_mu_sigma2_sigma2[i] = -2.0 / s2_3;
        mu_sigma2_sigma2_sigma2[i] = 6.0 * r / s2_4;
        sigma2_sigma2_sigma2_sigma2[i] = 3.0 * (s2 - 4.0 * r * r) / s2_5;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma2") = mu_mu_mu_sigma2,
        Named("mu_mu_sigma2_sigma2") = mu_mu_sigma2_sigma2,
        Named("mu_sigma2_sigma2_sigma2") = mu_sigma2_sigma2_sigma2,
        Named("sigma2_sigma2_sigma2_sigma2") = sigma2_sigma2_sigma2_sigma2
    );
}

// [[Rcpp::export]]
List lognormal_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma2(n), mu_mu_sigma2_sigma2(n),
                  mu_sigma2_sigma2_sigma2(n), sigma2_sigma2_sigma2_sigma2(n);
    bool s_is_scalar = (sigma2.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double s2 = s_is_scalar ? sigma2[0] : sigma2[i];
        double s2_3 = s2 * s2 * s2, s2_4 = s2_3 * s2;
        mu_mu_mu_mu[i] = 0.0;
        mu_mu_mu_sigma2[i] = 0.0;
        mu_mu_sigma2_sigma2[i] = -2.0 / s2_3;
        mu_sigma2_sigma2_sigma2[i] = 0.0;
        sigma2_sigma2_sigma2_sigma2[i] = -9.0 / s2_4;
    });

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma2") = mu_mu_mu_sigma2,
        Named("mu_mu_sigma2_sigma2") = mu_mu_sigma2_sigma2,
        Named("mu_sigma2_sigma2_sigma2") = mu_sigma2_sigma2_sigma2,
        Named("sigma2_sigma2_sigma2_sigma2") = sigma2_sigma2_sigma2_sigma2
    );
}
