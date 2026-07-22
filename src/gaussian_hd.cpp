#include <Rcpp.h>
using namespace Rcpp;

// Third- and fourth-order derivatives of the Gaussian log-density (sd
// parameterization), transcribed from the Wolfram output. r = mu - y.

// [[Rcpp::export]]
List gaussian_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s2 = s * s, s3 = s2 * s, s4 = s2 * s2, s5 = s4 * s;
        double r = m - y[i];

        mu_mu_mu[i] = 0.0;
        mu_mu_sigma[i] = 2.0 / s3;
        mu_sigma_sigma[i] = 6.0 * (-m + y[i]) / s4;
        sigma_sigma_sigma[i] = -2.0 * (s2 - 6.0 * r * r) / s5;
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List gaussian_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);

    bool sigma_is_scalar = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s3 = s * s * s;

        mu_mu_mu[i] = 0.0;
        mu_mu_sigma[i] = 2.0 / s3;
        mu_sigma_sigma[i] = 0.0;
        sigma_sigma_sigma[i] = 10.0 / s3;
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List gaussian_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n),
                  mu_sigma_sigma_sigma(n), sigma_sigma_sigma_sigma(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s2 = s * s, s4 = s2 * s2, s5 = s4 * s, s6 = s4 * s2;
        double r = m - y[i];

        mu_mu_mu_mu[i] = 0.0;
        mu_mu_mu_sigma[i] = 0.0;
        mu_mu_sigma_sigma[i] = -6.0 / s4;
        mu_sigma_sigma_sigma[i] = 24.0 * (m - y[i]) / s5;
        sigma_sigma_sigma_sigma[i] = 6.0 * (s2 - 10.0 * r * r) / s6;
    }

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List gaussian_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n),
                  mu_sigma_sigma_sigma(n), sigma_sigma_sigma_sigma(n);

    bool sigma_is_scalar = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s4 = s * s * s * s;

        mu_mu_mu_mu[i] = 0.0;
        mu_mu_mu_sigma[i] = 0.0;
        mu_mu_sigma_sigma[i] = -6.0 / s4;
        mu_sigma_sigma_sigma[i] = 0.0;
        sigma_sigma_sigma_sigma[i] = -54.0 / s4;
    }

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma
    );
}
