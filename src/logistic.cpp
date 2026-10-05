#include <Rcpp.h>
#include "d7_par.h"
#include "pt_logistic.h"
using namespace Rcpp;

// [[Rcpp::export]]
List logistic_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    int n = y.size();
    NumericVector grad_mu(n);
    NumericVector grad_sigma(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        
        double res = y[i] - m;
        double tanh_z = std::tanh(0.5 * res / s);
        
        grad_mu[i] = d7::logistic_score_mu(tanh_z, s);
        grad_sigma[i] = d7::logistic_score_sigma(res, tanh_z, s);
    });
    
    return List::create(Named("mu") = grad_mu, Named("sigma") = grad_sigma);
}

// [[Rcpp::export]]
List logistic_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma_sigma(n);
    NumericVector hess_mu_sigma(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        
        double s2 = s * s;
        double res = y[i] - m;
        double z_half = 0.5 * res / s;
        
        double tanh_z = std::tanh(z_half);
        double sech2_z = 1.0 - tanh_z * tanh_z;
        
        hess_mu_mu[i] = d7::logistic_hess_mu_mu(tanh_z, s);
        hess_sigma_sigma[i] = d7::logistic_hess_sigma_sigma(res, tanh_z, s);
        hess_mu_sigma[i] = -(tanh_z + z_half * sech2_z) / s2;
    });
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("mu_sigma") = hess_mu_sigma);
}

// [[Rcpp::export]]
List logistic_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma_sigma(n);
    NumericVector hess_mu_sigma(n);
    
    bool sigma_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        hess_mu_mu[i] = d7::logistic_expected_mu_mu(s);
        hess_sigma_sigma[i] = d7::logistic_expected_sigma_sigma(s);
        hess_mu_sigma[i] = 0.0;
    });
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("mu_sigma") = hess_mu_sigma);
}