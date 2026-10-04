#include <Rcpp.h>
#include "d7_par.h"
#include "pt_binomial.h"
using namespace Rcpp;

// [[Rcpp::export]]
List binomial_gradient_cpp(NumericVector y, NumericVector mu, NumericVector size,
                        int threads = 1) {
    int n = y.size();
    NumericVector grad_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool size_is_scalar = (size.size() == 1);
    
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double sz = size_is_scalar ? size[0] : size[i];
        
        grad_mu[i] = d7::binomial_score_mu(y[i], m, sz);
    });
    
    return List::create(Named("mu") = grad_mu);
}

// [[Rcpp::export]]
List binomial_hessian_cpp(NumericVector y, NumericVector mu, NumericVector size,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool size_is_scalar = (size.size() == 1);

    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double sz = size_is_scalar ? size[0] : size[i];
        
        hess_mu_mu[i] = d7::binomial_hess_mu_mu(y[i], m, sz);
    });
    
    return List::create(Named("mu_mu") = hess_mu_mu);
}

// [[Rcpp::export]]
List binomial_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector size,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool size_is_scalar = (size.size() == 1);

    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double sz = size_is_scalar ? size[0] : size[i];
        
        hess_mu_mu[i] = d7::binomial_expected_mu_mu(m, sz);
    });
    return List::create(Named("mu_mu") = hess_mu_mu);
}