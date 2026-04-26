#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List logistic_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector grad_mu(n);
    NumericVector grad_sigma(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    
    for(int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        
        double res = y[i] - m;
        double tanh_z = std::tanh(0.5 * res / s);
        
        grad_mu[i] = tanh_z / s;
        grad_sigma[i] = (res * tanh_z - s) / (s * s);
    }
    
    return List::create(Named("mu") = grad_mu, Named("sigma") = grad_sigma);
}

// [[Rcpp::export]]
List logistic_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
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
        double res = y[i] - m;
        double z_half = 0.5 * res / s;
        
        double tanh_z = std::tanh(z_half);
        double sech2_z = 1.0 - tanh_z * tanh_z;
        
        hess_mu_mu[i] = -sech2_z / (2.0 * s2);
        hess_sigma_sigma[i] = (1.0 - 4.0 * z_half * tanh_z - 2.0 * z_half * z_half * sech2_z) / s2;
        hess_mu_sigma[i] = -(tanh_z + z_half * sech2_z) / s2;
    }
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("mu_sigma") = hess_mu_sigma);
}

// [[Rcpp::export]]
List logistic_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma_sigma(n);
    NumericVector hess_mu_sigma(n);
    
    bool sigma_is_scalar = (sigma.size() == 1);

    for(int i = 0; i < n; i++) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        double s2 = s * s;
        
        hess_mu_mu[i] = -1.0 / (3.0 * s2);
        // M_PI è predefinito in Rmath.h (incluso tramite Rcpp.h)
        hess_sigma_sigma[i] = -(3.0 + M_PI * M_PI) / (9.0 * s2);
        hess_mu_sigma[i] = 0.0;
    }
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("mu_sigma") = hess_mu_sigma);
}