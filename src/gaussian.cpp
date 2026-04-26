#include <Rcpp.h>
using namespace Rcpp;

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
        
        double s2 = s * s;
        double s3 = s2 * s;
        double res = y[i] - m;
        
        grad_mu[i] = res / s2;
        grad_sigma[i] = (res * res - s2) / s3;
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
        
        double s2 = s * s;
        double s3 = s2 * s;
        double s4 = s2 * s2;
        double res = y[i] - m;
        
        hess_mu_mu[i] = -1.0 / s2;
        hess_sigma_sigma[i] = (s2 - 3.0 * res * res) / s4;
        hess_mu_sigma[i] = -2.0 * res / s3;
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
        double s2 = s * s;
        
        hess_mu_mu[i] = -1.0 / s2;
        hess_sigma_sigma[i] = -2.0 / s2;
        hess_mu_sigma[i] = 0.0;
    }
    return List::create(
        Named("mu_mu") = hess_mu_mu,
        Named("sigma_sigma") = hess_sigma_sigma,
        Named("mu_sigma") = hess_mu_sigma
    );
}
