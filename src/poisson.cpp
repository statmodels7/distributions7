#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List poisson_gradient_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector grad_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    double m = 0;
    
    if (mu_is_scalar) m = mu[0];

    for(int i = 0; i < n; i++) {
        if (!mu_is_scalar) {
            m = mu[i];
        }
        grad_mu[i] = (y[i] - m) / m;
    }
    
    return List::create(Named("mu") = grad_mu);
}

// [[Rcpp::export]]
List poisson_hessian_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    double m = 0, m2 = 0;
    
    if (mu_is_scalar) {
        m = mu[0];
        m2 = m * m;
    }

    for(int i = 0; i < n; i++) {
        if (!mu_is_scalar) {
            m2 = mu[i] * mu[i];
        }
        hess_mu_mu[i] = -y[i] / m2;
    }
    
    return List::create(Named("mu_mu") = hess_mu_mu);
}

// [[Rcpp::export]]
List poisson_expected_hessian_cpp(NumericVector y, NumericVector mu) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    double val = 0;
    
    if (mu_is_scalar) {
        val = -1.0 / mu[0];
    }

    for(int i = 0; i < n; i++) {
        if (!mu_is_scalar) {
            val = -1.0 / mu[i];
        }
        hess_mu_mu[i] = val;
    }
    
    return List::create(Named("mu_mu") = hess_mu_mu);
}