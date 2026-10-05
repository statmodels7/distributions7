#include <Rcpp.h>
#include "d7_par.h"
#include "pt_bernoulli.h"
using namespace Rcpp;

// [[Rcpp::export]]
List bernoulli_gradient_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector grad_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        
        grad_mu[i] = d7::bernoulli_score_mu(y[i], m);
    });
    
    return List::create(
        Named("mu") = grad_mu
    );
}

// [[Rcpp::export]]
List bernoulli_hessian_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);

    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        
        hess_mu_mu[i] = d7::bernoulli_hess_mu_mu(y[i], m);
    });
    return List::create(
        Named("mu_mu") = hess_mu_mu
    );
}

// [[Rcpp::export]]
List bernoulli_expected_hessian_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);

    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        hess_mu_mu[i] = d7::bernoulli_expected_mu_mu(m);
    });
    return List::create(
        Named("mu_mu") = hess_mu_mu
    );
}