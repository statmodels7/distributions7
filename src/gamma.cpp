#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// [[Rcpp::export]]
List gamma_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    int n = y.size();
    NumericVector grad_mu(n);
    NumericVector grad_sigma2(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma2_is_scalar = (sigma2.size() == 1);
    bool both_scalar = mu_is_scalar && sigma2_is_scalar;
    
    double m0 = 0, s20 = 0, alpha0 = 0, lambda0 = 0;
    double digamma_alpha0 = 0, log_lambda0 = 0;
    
    if (both_scalar) {
        m0 = mu[0];
        s20 = sigma2[0];
        alpha0 = m0 * m0 / s20;
        lambda0 = m0 / s20;
        digamma_alpha0 = R::digamma(alpha0);
        log_lambda0 = std::log(lambda0);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m0, s2 = s20, alpha = alpha0, lambda = lambda0, digamma_alpha = digamma_alpha0, log_lambda = log_lambda0;
        if (!both_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            s2 = sigma2_is_scalar ? sigma2[0] : sigma2[i];
            alpha = m * m / s2;
            lambda = m / s2;
            digamma_alpha = R::digamma(alpha);
            log_lambda = std::log(lambda);
        }
        
        double log_y = std::log(y[i]);
        
        grad_mu[i] = (-2.0 * m * digamma_alpha + 2.0 * m * log_lambda + m + 2.0 * m * log_y - y[i]) / s2;
        grad_sigma2[i] = -(m * (-m * digamma_alpha + m + m * (log_lambda + log_y) - y[i])) / (s2 * s2);
    });
    
    return List::create(Named("mu") = grad_mu, Named("sigma2") = grad_sigma2);
}

// [[Rcpp::export]]
List gamma_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma2_sigma2(n);
    NumericVector hess_mu_sigma2(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma2_is_scalar = (sigma2.size() == 1);
    bool both_scalar = mu_is_scalar && sigma2_is_scalar;
    
    double m0 = 0, s20 = 0, alpha0 = 0, lambda0 = 0;
    double digamma_alpha0 = 0, trigamma_alpha0 = 0, log_lambda0 = 0;

    if (both_scalar) {
        m0 = mu[0];
        s20 = sigma2[0];
        alpha0 = m0 * m0 / s20;
        lambda0 = m0 / s20;
        digamma_alpha0 = R::digamma(alpha0);
        trigamma_alpha0 = R::trigamma(alpha0);
        log_lambda0 = std::log(lambda0);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m0, s2 = s20, alpha = alpha0, lambda = lambda0, digamma_alpha = digamma_alpha0, trigamma_alpha = trigamma_alpha0, log_lambda = log_lambda0;
        if (!both_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            s2 = sigma2_is_scalar ? sigma2[0] : sigma2[i];
            alpha = m * m / s2;
            lambda = m / s2;
            digamma_alpha = R::digamma(alpha);
            trigamma_alpha = R::trigamma(alpha);
            log_lambda = std::log(lambda);
        }
        
        double log_y = std::log(y[i]);
        
        hess_mu_mu[i] = (-(4.0 * m * m * trigamma_alpha) / s2 - 2.0 * digamma_alpha + 2.0 * log_lambda + 2.0 * log_y + 3.0) / s2;
        hess_sigma2_sigma2[i] = -(m * (2.0 * m * s2 * digamma_alpha + m * m * m * trigamma_alpha + s2 * (-2.0 * m * log_lambda - 3.0 * m - 2.0 * m * log_y + 2.0 * y[i]))) / (s2 * s2 * s2 * s2);
        hess_mu_sigma2[i] = (2.0 * m * s2 * digamma_alpha + 2.0 * m * m * m * trigamma_alpha + s2 * (-2.0 * m * log_lambda - 3.0 * m - 2.0 * m * log_y + y[i])) / (s2 * s2 * s2);
    });
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma2_sigma2") = hess_sigma2_sigma2, Named("mu_sigma2") = hess_mu_sigma2);
}

// [[Rcpp::export]]
List gamma_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma2_sigma2(n);
    NumericVector hess_mu_sigma2(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma2_is_scalar = (sigma2.size() == 1);
    bool both_scalar = mu_is_scalar && sigma2_is_scalar;

    double m0 = 0, s20 = 0, alpha0 = 0, trigamma_alpha0 = 0;
    
    if (both_scalar) {
        m0 = mu[0];
        s20 = sigma2[0];
        alpha0 = m0 * m0 / s20;
        trigamma_alpha0 = R::trigamma(alpha0);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m0, s2 = s20, alpha = alpha0, trigamma_alpha = trigamma_alpha0;
        if (!both_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            s2 = sigma2_is_scalar ? sigma2[0] : sigma2[i];
            alpha = m * m / s2;
            trigamma_alpha = R::trigamma(alpha);
        }
        
        hess_mu_mu[i] = (3.0 * s2 - 4.0 * m * m * trigamma_alpha) / (s2 * s2);
        hess_sigma2_sigma2[i] = -(m * m * (m * m * trigamma_alpha - s2)) / (s2 * s2 * s2 * s2);
        hess_mu_sigma2[i] = 2.0 * m * (m * m * trigamma_alpha - s2) / (s2 * s2 * s2);
    });
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma2_sigma2") = hess_sigma2_sigma2, Named("mu_sigma2") = hess_mu_sigma2);
}