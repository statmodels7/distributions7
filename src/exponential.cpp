#include <Rcpp.h>
using namespace Rcpp;

// Exponential in the MEAN parametrisation: f(y) = exp(-y/mu)/mu, so
//   l = -log(mu) - y/mu
// and each part contributes one term whose every order is elementary:
//   d^k/dmu^k [-log mu]  = (-1)^k (k-1)! / mu^k
//   d^k/dmu^k [-y / mu]  = (-1)^(k+1) k! y / mu^(k+1)
// Under E[y] = mu the two collapse into one term,
//   E[l^(k)] = (-1)^k (k-1)! (1 - k) / mu^k,
// which vanishes at k = 1 (the first Bartlett identity) and gives the
// information 1/mu^2 at k = 2.

// [[Rcpp::export]]
List exponential_gradient_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector grad_mu(n);
    bool mu_is_scalar = (mu.size() == 1);
    double m = mu_is_scalar ? mu[0] : 0.0;

    for (int i = 0; i < n; i++) {
        if (!mu_is_scalar) m = mu[i];
        grad_mu[i] = (y[i] - m) / (m * m);
    }
    return List::create(Named("mu") = grad_mu);
}

// [[Rcpp::export]]
List exponential_hessian_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector h(n);
    bool mu_is_scalar = (mu.size() == 1);
    double m = mu_is_scalar ? mu[0] : 0.0;

    for (int i = 0; i < n; i++) {
        if (!mu_is_scalar) m = mu[i];
        double m2 = m * m;
        h[i] = 1.0 / m2 - 2.0 * y[i] / (m2 * m);
    }
    return List::create(Named("mu_mu") = h);
}

// [[Rcpp::export]]
List exponential_expected_hessian_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector h(n);
    bool mu_is_scalar = (mu.size() == 1);
    double m = mu_is_scalar ? mu[0] : 0.0;

    for (int i = 0; i < n; i++) {
        if (!mu_is_scalar) m = mu[i];
        h[i] = -1.0 / (m * m);
    }
    return List::create(Named("mu_mu") = h);
}

// [[Rcpp::export]]
List exponential_deriv3_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector d(n);
    bool mu_is_scalar = (mu.size() == 1);
    double m = mu_is_scalar ? mu[0] : 0.0;

    for (int i = 0; i < n; i++) {
        if (!mu_is_scalar) m = mu[i];
        double m3 = m * m * m;
        d[i] = -2.0 / m3 + 6.0 * y[i] / (m3 * m);
    }
    return List::create(Named("mu_mu_mu") = d);
}

// [[Rcpp::export]]
List exponential_deriv3_expected_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector d(n);
    bool mu_is_scalar = (mu.size() == 1);
    double m = mu_is_scalar ? mu[0] : 0.0;

    for (int i = 0; i < n; i++) {
        if (!mu_is_scalar) m = mu[i];
        d[i] = 4.0 / (m * m * m);
    }
    return List::create(Named("mu_mu_mu") = d);
}

// [[Rcpp::export]]
List exponential_deriv4_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector d(n);
    bool mu_is_scalar = (mu.size() == 1);
    double m = mu_is_scalar ? mu[0] : 0.0;

    for (int i = 0; i < n; i++) {
        if (!mu_is_scalar) m = mu[i];
        double m2 = m * m;
        double m4 = m2 * m2;
        d[i] = 6.0 / m4 - 24.0 * y[i] / (m4 * m);
    }
    return List::create(Named("mu_mu_mu_mu") = d);
}

// [[Rcpp::export]]
List exponential_deriv4_expected_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector d(n);
    bool mu_is_scalar = (mu.size() == 1);
    double m = mu_is_scalar ? mu[0] : 0.0;

    for (int i = 0; i < n; i++) {
        if (!mu_is_scalar) m = mu[i];
        double m2 = m * m;
        d[i] = -18.0 / (m2 * m2);
    }
    return List::create(Named("mu_mu_mu_mu") = d);
}
