#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List lognormal_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector grad_mu(n);
    NumericVector grad_sigma2(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma2_is_scalar = (sigma2.size() == 1);
    
    for(int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s2 = sigma2_is_scalar ? sigma2[0] : sigma2[i];
        
        double log_y = std::log(y[i]);
        double res = log_y - m;
        
        grad_mu[i] = res / s2;
        grad_sigma2[i] = (res * res - s2) / (2.0 * s2 * s2);
    }
    
    return List::create(Named("mu") = grad_mu, Named("sigma2") = grad_sigma2);
}

// [[Rcpp::export]]
List lognormal_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma2_sigma2(n);
    NumericVector hess_mu_sigma2(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma2_is_scalar = (sigma2.size() == 1);

    for(int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s2 = sigma2_is_scalar ? sigma2[0] : sigma2[i];
        
        double s2_sq = s2 * s2;
        double s2_cub = s2_sq * s2;
        double log_y = std::log(y[i]);
        double res = log_y - m;
        double res2 = res * res;
        
        hess_mu_mu[i] = -1.0 / s2;
        hess_sigma2_sigma2[i] = -res2 / s2_cub + 0.5 / s2_sq;
        hess_mu_sigma2[i] = -res / s2_sq;
    }
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma2_sigma2") = hess_sigma2_sigma2, Named("mu_sigma2") = hess_mu_sigma2);
}

// [[Rcpp::export]]
List lognormal_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma2_sigma2(n);
    NumericVector hess_mu_sigma2(n);
    
    bool sigma2_is_scalar = (sigma2.size() == 1);

    for(int i = 0; i < n; i++) {
        double s2 = sigma2_is_scalar ? sigma2[0] : sigma2[i];
        
        double s2_sq = s2 * s2;
        
        hess_mu_mu[i] = -1.0 / s2;
        hess_sigma2_sigma2[i] = -0.5 / s2_sq;
        hess_mu_sigma2[i] = 0.0;
    }
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma2_sigma2") = hess_sigma2_sigma2, Named("mu_sigma2") = hess_mu_sigma2);
}
