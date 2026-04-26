#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List binomial_gradient_cpp(NumericVector y, NumericVector mu, NumericVector size) {
    int n = y.size();
    NumericVector grad_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool size_is_scalar = (size.size() == 1);
    
    for(int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double sz = size_is_scalar ? size[0] : size[i];
        
        grad_mu[i] = (y[i] - sz * m) / (m * (1.0 - m));
    }
    
    return List::create(Named("mu") = grad_mu);
}

// [[Rcpp::export]]
List binomial_hessian_cpp(NumericVector y, NumericVector mu, NumericVector size) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool size_is_scalar = (size.size() == 1);

    for(int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double sz = size_is_scalar ? size[0] : size[i];
        
        double m2 = m * m;
        double one_m = 1.0 - m;
        double one_m2 = one_m * one_m;
        
        hess_mu_mu[i] = -(y[i] / m2) - ((sz - y[i]) / one_m2);
    }
    
    return List::create(Named("mu_mu") = hess_mu_mu);
}

// [[Rcpp::export]]
List binomial_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector size) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool size_is_scalar = (size.size() == 1);

    for(int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double sz = size_is_scalar ? size[0] : size[i];
        
        hess_mu_mu[i] = -sz / (m * (1.0 - m));
    }
    return List::create(Named("mu_mu") = hess_mu_mu);
}