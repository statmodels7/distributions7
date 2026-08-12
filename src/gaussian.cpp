#include <Rcpp.h>
using namespace Rcpp;

// Gaussian in the mean and the STANDARD DEVIATION. Every component is written
// in z = (y - mu)/sigma and inv = 1/sigma, and never in a positive power of
// the scale. The algebraically equivalent forms -- (res^2 - sigma^2)/sigma^3
// for the score, (sigma^2 - 3 res^2)/sigma^4 for the second derivative --
// overflow in the DENOMINATOR before the ratio does, so the score returns
// exactly 0 above sigma = 5.6e102 and NaN above 1.3e154 where its true value,
// -1/sigma on this scale and -1 on the link scale, stays representable. A
// score of zero is what a stopping rule reads as stationarity, so an
// optimizer that wanders out there is told it has arrived.

// [[Rcpp::export]]
List gaussian_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();

    NumericVector grad_mu(n);
    NumericVector grad_sigma(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    for(int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];

        double inv = 1.0 / s;
        double z = (y[i] - m) / s;

        grad_mu[i] = z * inv;
        grad_sigma[i] = (z * z - 1.0) * inv;
    }

    return List::create(
        Named("mu") = grad_mu,
        Named("sigma") = grad_sigma
    );
}

// [[Rcpp::export]]
List gaussian_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();

    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma_sigma(n);
    NumericVector hess_mu_sigma(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    for(int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];

        double inv = 1.0 / s;
        double inv2 = inv * inv;
        double z = (y[i] - m) / s;

        hess_mu_mu[i] = -inv2;
        hess_sigma_sigma[i] = (1.0 - 3.0 * z * z) * inv2;
        hess_mu_sigma[i] = -2.0 * z * inv2;
    }
    return List::create(
        Named("mu_mu") = hess_mu_mu,
        Named("sigma_sigma") = hess_sigma_sigma,
        Named("mu_sigma") = hess_mu_sigma
    );
}

// [[Rcpp::export]]
List gaussian_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();

    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma_sigma(n);
    NumericVector hess_mu_sigma(n);

    bool sigma_is_scalar = (sigma.size() == 1);

    for(int i = 0; i < n; i++) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double inv = 1.0 / s;
        double inv2 = inv * inv;

        hess_mu_mu[i] = -inv2;
        hess_sigma_sigma[i] = -2.0 * inv2;
        hess_mu_sigma[i] = 0.0;
    }
    return List::create(
        Named("mu_mu") = hess_mu_mu,
        Named("sigma_sigma") = hess_sigma_sigma,
        Named("mu_sigma") = hess_mu_sigma
    );
}
