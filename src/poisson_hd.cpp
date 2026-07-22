#include <Rcpp.h>
using namespace Rcpp;

// Third/fourth-order derivatives of the Poisson log-mass. ell = y log(mu) - mu
// - log(y!), so d^k ell / d mu^k = (-1)^k (k-1)! y / mu^k for k >= 2, with
// expectation obtained via E[y] = mu.

// [[Rcpp::export]]
List poisson_deriv3_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        mu_mu_mu[i] = 2.0 * y[i] / (m * m * m);
    }
    return List::create(Named("mu_mu_mu") = mu_mu_mu);
}

// [[Rcpp::export]]
List poisson_deriv3_expected_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        mu_mu_mu[i] = 2.0 / (m * m);
    }
    return List::create(Named("mu_mu_mu") = mu_mu_mu);
}

// [[Rcpp::export]]
List poisson_deriv4_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        mu_mu_mu_mu[i] = -6.0 * y[i] / (m * m * m * m);
    }
    return List::create(Named("mu_mu_mu_mu") = mu_mu_mu_mu);
}

// [[Rcpp::export]]
List poisson_deriv4_expected_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        mu_mu_mu_mu[i] = -6.0 / (m * m * m);
    }
    return List::create(Named("mu_mu_mu_mu") = mu_mu_mu_mu);
}
