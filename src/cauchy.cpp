#include <Rcpp.h>
#include "d7_par.h"
#include "pt_cauchy.h"
using namespace Rcpp;

// [[Rcpp::export]]
List cauchy_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    int n = y.size();
    NumericVector grad_mu(n);
    NumericVector grad_sigma(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        
        double res = y[i] - m;
        grad_mu[i] = d7::cauchy_score_mu(res, s);
        grad_sigma[i] = d7::cauchy_score_sigma(res, s);
    });
    
    return List::create(Named("mu") = grad_mu, Named("sigma") = grad_sigma);
}

// [[Rcpp::export]]
List cauchy_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma_sigma(n);
    NumericVector hess_mu_sigma(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        
        double res = y[i] - m;
        double res2 = res * res;
        double s2 = s * s;
        double den = s2 + res2;
        double den2 = den * den;
        
        hess_mu_mu[i] = d7::cauchy_hess_mu_mu(res, s);
        hess_sigma_sigma[i] = d7::cauchy_hess_sigma_sigma(res, s);
        hess_mu_sigma[i] = -4.0 * s * res / den2;
    });
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("mu_sigma") = hess_mu_sigma);
}

// [[Rcpp::export]]
List cauchy_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma_sigma(n);
    NumericVector hess_mu_sigma(n);
    
    bool sigma_is_scalar = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double s = sigma_is_scalar ? sigma[0] : sigma[i];
        hess_mu_mu[i] = d7::cauchy_expected_mu_mu(s);
        hess_sigma_sigma[i] = d7::cauchy_expected_sigma_sigma(s);
        hess_mu_sigma[i] = 0.0;
    });
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("mu_sigma") = hess_mu_sigma);
}