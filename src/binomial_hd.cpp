#include <Rcpp.h>
using namespace Rcpp;

// Third/fourth-order derivatives of the Binomial log-mass, transcribed from the
// Wolfram output. Single free parameter mu (probability); om = mu - 1.

// [[Rcpp::export]]
List binomial_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector size) {
    int n = y.size();
    NumericVector mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool size_is_scalar = (size.size() == 1);
    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double sz = size_is_scalar ? size[0] : size[i];
        double om = m - 1.0;
        mu_mu_mu[i] = 2.0 * (sz - y[i]) / (om * om * om) + 2.0 * y[i] / (m * m * m);
    }
    return List::create(Named("mu_mu_mu") = mu_mu_mu);
}

// [[Rcpp::export]]
List binomial_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector size) {
    int n = y.size();
    NumericVector mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool size_is_scalar = (size.size() == 1);
    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double sz = size_is_scalar ? size[0] : size[i];
        double om = m - 1.0;
        mu_mu_mu[i] = (2.0 - 4.0 * m) * sz / (om * om * m * m);
    }
    return List::create(Named("mu_mu_mu") = mu_mu_mu);
}

// [[Rcpp::export]]
List binomial_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector size) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool size_is_scalar = (size.size() == 1);
    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double sz = size_is_scalar ? size[0] : size[i];
        double om = m - 1.0;
        mu_mu_mu_mu[i] = -6.0 * (sz - y[i]) / (om * om * om * om) - 6.0 * y[i] / (m * m * m * m);
    }
    return List::create(Named("mu_mu_mu_mu") = mu_mu_mu_mu);
}

// [[Rcpp::export]]
List binomial_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector size) {
    int n = y.size();
    NumericVector mu_mu_mu_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    bool size_is_scalar = (size.size() == 1);
    for (int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double sz = size_is_scalar ? size[0] : size[i];
        double om = m - 1.0;
        mu_mu_mu_mu[i] = 6.0 * (1.0 + 3.0 * om * m) * sz / (om * om * om * m * m * m);
    }
    return List::create(Named("mu_mu_mu_mu") = mu_mu_mu_mu);
}
