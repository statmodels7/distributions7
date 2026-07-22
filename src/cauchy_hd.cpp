#include <Rcpp.h>
using namespace Rcpp;

// Third/fourth-order derivatives of the Cauchy log-density, transcribed from the
// Wolfram output. r = mu - y, D = sigma^2 + r^2.

// [[Rcpp::export]]
List cauchy_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s2 = s * s, s3 = s2 * s;
        double r = m - y[i], r2 = r * r;
        double D = s2 + r2, D2 = D * D, D3 = D2 * D;

        mu_mu_mu[i] = -4.0 * r * (r2 - 3.0 * s2) / D3;
        mu_mu_sigma[i] = 4.0 * s * (s2 - 3.0 * r2) / D3;
        mu_sigma_sigma[i] = 4.0 * r * (r2 - 3.0 * s2) / D3;
        sigma_sigma_sigma[i] = 2.0 / s3 - 16.0 * s3 / D3 + 12.0 * s / D2;
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List cauchy_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu(n), mu_mu_sigma(n), mu_sigma_sigma(n), sigma_sigma_sigma(n);
    bool sigma_is_scalar = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s3 = s * s * s;
        mu_mu_mu[i] = 0.0;
        mu_mu_sigma[i] = 1.0 / (2.0 * s3);
        mu_sigma_sigma[i] = 0.0;
        sigma_sigma_sigma[i] = 3.0 / (2.0 * s3);
    }

    return List::create(
        Named("mu_mu_mu") = mu_mu_mu,
        Named("mu_mu_sigma") = mu_mu_sigma,
        Named("mu_sigma_sigma") = mu_sigma_sigma,
        Named("sigma_sigma_sigma") = sigma_sigma_sigma
    );
}

// [[Rcpp::export]]
List cauchy_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n),
                  mu_sigma_sigma_sigma(n), sigma_sigma_sigma_sigma(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s2 = s * s, s4 = s2 * s2;
        double r = m - y[i], r2 = r * r;
        double D = s2 + r2, D2 = D * D, D3 = D2 * D, D4 = D2 * D2;

        // (mu^2+2 mu sigma-sigma^2-2(mu+sigma)y+y^2)(mu^2-sigma^2+2 sigma y+y^2-2 mu(sigma+y))
        double f1 = r2 + 2.0 * s * r - s2;
        double f2 = r2 - 2.0 * s * r - s2;
        double cross = s * r * (s2 - r2);   // sigma*(mu-y)*(mu+sigma-y)*(-mu+sigma+y)

        mu_mu_mu_mu[i] = 12.0 * f1 * f2 / D4;
        mu_mu_mu_sigma[i] = -48.0 * cross / D4;
        mu_mu_sigma_sigma[i] = -12.0 * f1 * f2 / D4;
        mu_sigma_sigma_sigma[i] = 48.0 * cross / D4;
        sigma_sigma_sigma_sigma[i] = -6.0 / s4 + 96.0 * s4 / D4 - 96.0 * s2 / D3 + 12.0 / D2;
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
List cauchy_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n), mu_mu_mu_sigma(n), mu_mu_sigma_sigma(n),
                  mu_sigma_sigma_sigma(n), sigma_sigma_sigma_sigma(n);
    bool sigma_is_scalar = (sigma.size() == 1);

    for (int i = 0; i < n; i++) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s4 = s * s * s * s;
        mu_mu_mu_mu[i] = 3.0 / (4.0 * s4);
        mu_mu_mu_sigma[i] = 0.0;
        mu_mu_sigma_sigma[i] = -3.0 / (4.0 * s4);
        mu_sigma_sigma_sigma[i] = 0.0;
        sigma_sigma_sigma_sigma[i] = -21.0 / (4.0 * s4);
    }

    return List::create(
        Named("mu_mu_mu_mu") = mu_mu_mu_mu,
        Named("mu_mu_mu_sigma") = mu_mu_mu_sigma,
        Named("mu_mu_sigma_sigma") = mu_mu_sigma_sigma,
        Named("mu_sigma_sigma_sigma") = mu_sigma_sigma_sigma,
        Named("sigma_sigma_sigma_sigma") = sigma_sigma_sigma_sigma
    );
}
