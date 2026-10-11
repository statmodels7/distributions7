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

// order 5 of log f in (mu, sigma): -log sigma - r^2 / (2 sigma^2), r = y - mu
// [[Rcpp::export]]
List gaussian_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                          int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n), o3(n), o4(n), o5(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool p_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = p_is_scalar ? sigma[0] : sigma[i];
        double r = y[i] - m;
        double inv = 1.0 / s, inv2 = inv * inv, inv5 = inv2 * inv2 * inv;
        double z = r * inv;
        o0[i] = 0.0;
        o1[i] = 0.0;
        o2[i] = 0.0;
        o3[i] = 24.0 * inv5;
        o4[i] = 120.0 * z * inv5;
        o5[i] = (360.0 * z * z - 24.0) * inv5;
    });

    return List::create(
        Named("mu_mu_mu_mu_mu") = o0,
        Named("mu_mu_mu_mu_sigma") = o1,
        Named("mu_mu_mu_sigma_sigma") = o2,
        Named("mu_mu_sigma_sigma_sigma") = o3,
        Named("mu_sigma_sigma_sigma_sigma") = o4,
        Named("sigma_sigma_sigma_sigma_sigma") = o5
    );
}

// order 5 of log f in (mu, sigma2): -log(sigma2) / 2 - r^2 / (2 sigma2)
// [[Rcpp::export]]
List gaussian2_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                          int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n), o3(n), o4(n), o5(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool p_is_scalar = (sigma2.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = p_is_scalar ? sigma2[0] : sigma2[i];
        double r = y[i] - m;
        double inv = 1.0 / s, inv2 = inv * inv, inv4 = inv2 * inv2, inv5 = inv4 * inv;
        o0[i] = 0.0;
        o1[i] = 0.0;
        o2[i] = 0.0;
        o3[i] = 6.0 * inv4;
        o4[i] = 24.0 * r * inv5;
        o5[i] = (60.0 * r * r * inv - 12.0) * inv5;
    });

    return List::create(
        Named("mu_mu_mu_mu_mu") = o0,
        Named("mu_mu_mu_mu_sigma2") = o1,
        Named("mu_mu_mu_sigma2_sigma2") = o2,
        Named("mu_mu_sigma2_sigma2_sigma2") = o3,
        Named("mu_sigma2_sigma2_sigma2_sigma2") = o4,
        Named("sigma2_sigma2_sigma2_sigma2_sigma2") = o5
    );
}

// order 5 of log f in (mu, tau): log(tau) / 2 - tau r^2 / 2
// [[Rcpp::export]]
List gaussian3_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector tau,
                          int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n), o3(n), o4(n), o5(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool p_is_scalar = (tau.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = p_is_scalar ? tau[0] : tau[i];
        double r = y[i] - m;
        double inv = 1.0 / s, inv2 = inv * inv, inv5 = inv2 * inv2 * inv;
        (void) r;
        o0[i] = 0.0;
        o1[i] = 0.0;
        o2[i] = 0.0;
        o3[i] = 0.0;
        o4[i] = 0.0;
        o5[i] = 12.0 * inv5;
    });

    return List::create(
        Named("mu_mu_mu_mu_mu") = o0,
        Named("mu_mu_mu_mu_tau") = o1,
        Named("mu_mu_mu_tau_tau") = o2,
        Named("mu_mu_tau_tau_tau") = o3,
        Named("mu_tau_tau_tau_tau") = o4,
        Named("tau_tau_tau_tau_tau") = o5
    );
}
