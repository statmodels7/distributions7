#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List gamma_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector grad_mu(n);
    NumericVector grad_sigma2(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma2_is_scalar = (sigma2.size() == 1);
    bool both_scalar = mu_is_scalar && sigma2_is_scalar;
    
    double m = 0, s2 = 0, alpha = 0, lambda = 0;
    double digamma_alpha = 0, log_lambda = 0;
    
    if (both_scalar) {
        m = mu[0];
        s2 = sigma2[0];
        alpha = m * m / s2;
        lambda = m / s2;
        digamma_alpha = R::digamma(alpha);
        log_lambda = std::log(lambda);
    }

    for(int i = 0; i < n; i++) {
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
    }
    
    return List::create(Named("mu") = grad_mu, Named("sigma2") = grad_sigma2);
}

// [[Rcpp::export]]
List gamma_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma2_sigma2(n);
    NumericVector hess_mu_sigma2(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma2_is_scalar = (sigma2.size() == 1);
    bool both_scalar = mu_is_scalar && sigma2_is_scalar;
    
    double m = 0, s2 = 0, alpha = 0, lambda = 0;
    double digamma_alpha = 0, trigamma_alpha = 0, log_lambda = 0;

    if (both_scalar) {
        m = mu[0];
        s2 = sigma2[0];
        alpha = m * m / s2;
        lambda = m / s2;
        digamma_alpha = R::digamma(alpha);
        trigamma_alpha = R::trigamma(alpha);
        log_lambda = std::log(lambda);
    }

    for(int i = 0; i < n; i++) {
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
    }
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma2_sigma2") = hess_sigma2_sigma2, Named("mu_sigma2") = hess_mu_sigma2);
}

// [[Rcpp::export]]
List gamma_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma2) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma2_sigma2(n);
    NumericVector hess_mu_sigma2(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma2_is_scalar = (sigma2.size() == 1);
    bool both_scalar = mu_is_scalar && sigma2_is_scalar;

    double m = 0, s2 = 0, alpha = 0, trigamma_alpha = 0;
    
    if (both_scalar) {
        m = mu[0];
        s2 = sigma2[0];
        alpha = m * m / s2;
        trigamma_alpha = R::trigamma(alpha);
    }

    for(int i = 0; i < n; i++) {
        if (!both_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            s2 = sigma2_is_scalar ? sigma2[0] : sigma2[i];
            alpha = m * m / s2;
            trigamma_alpha = R::trigamma(alpha);
        }
        
        hess_mu_mu[i] = (3.0 * s2 - 4.0 * m * m * trigamma_alpha) / (s2 * s2);
        hess_sigma2_sigma2[i] = -(m * m * (m * m * trigamma_alpha - s2)) / (s2 * s2 * s2 * s2);
        hess_mu_sigma2[i] = 2.0 * m * (m * m * trigamma_alpha - s2) / (s2 * s2 * s2);
    }
    return List::create(Named("mu_mu") = hess_mu_mu, Named("sigma2_sigma2") = hess_sigma2_sigma2, Named("mu_sigma2") = hess_mu_sigma2);
}