#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List cauchy_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector grad_mu(n);
    NumericVector grad_sigma(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    
    for(int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        
        double res = y[i] - m;
        double res2 = res * res;
        double s2 = s * s;
        double den = s2 + res2;
        
        grad_mu[i] = (2.0 * res) / den;
        grad_sigma[i] = (res2 - s2) / (s * den);
    }
    
    return List::create(Named("mu") = grad_mu, Named("sigma") = grad_sigma);
}

// [[Rcpp::export]]
List cauchy_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma_sigma(n);
    NumericVector hess_mu_sigma(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    for(int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        
        double res = y[i] - m;
        double res2 = res * res;
        double res4 = res2 * res2;
        double s2 = s * s;
        double s4 = s2 * s2;
        double den = s2 + res2;
        double den2 = den * den;
        
        hess_mu_mu[i] = (2.0 * res2 - 2.0 * s2) / den2;
        hess_sigma_sigma[i] = (s4 - 4.0 * s2 * res2 - res4) / (s2 * den2);
        hess_mu_sigma[i] = -4.0 * s * res / den2;
    }
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("mu_sigma") = hess_mu_sigma);
}

// [[Rcpp::export]]
List cauchy_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma_sigma(n);
    NumericVector hess_mu_sigma(n);
    
    bool sigma_is_scalar = (sigma.size() == 1);

    for(int i = 0; i < n; i++) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s2 = s * s;
        double val = -0.5 / s2;
        
        hess_mu_mu[i] = val;
        hess_sigma_sigma[i] = val;
        hess_mu_sigma[i] = 0.0;
    }
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("mu_sigma") = hess_mu_sigma);
}