#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// [[Rcpp::export]]
List beta_gradient_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int threads = 1) {
    int n = y.size();
    NumericVector grad_mu(n);
    NumericVector grad_phi(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);
    bool both_scalar = mu_is_scalar && phi_is_scalar;
    
    double m0 = 0, p0 = 0, alpha0 = 0, beta_shape0 = 0;
    double digamma_alpha0 = 0, digamma_beta0 = 0, digamma_p0 = 0;
    
    if (both_scalar) {
        m0 = mu[0];
        p0 = phi[0];
        alpha0 = m0 * p0;
        beta_shape0 = (1.0 - m0) * p0;
        digamma_alpha0 = R::digamma(alpha0);
        digamma_beta0 = R::digamma(beta_shape0);
        digamma_p0 = R::digamma(p0);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m0, p = p0, alpha = alpha0, beta_shape = beta_shape0, digamma_alpha = digamma_alpha0, digamma_beta = digamma_beta0, digamma_p = digamma_p0;
        if (!both_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            p = phi_is_scalar ? phi[0] : phi[i];
            alpha = m * p;
            beta_shape = (1.0 - m) * p;
            digamma_alpha = R::digamma(alpha);
            digamma_beta = R::digamma(beta_shape);
            digamma_p = R::digamma(p);
        }
        
        double log_y = std::log(y[i]);
        double log_1_y = std::log(1.0 - y[i]);
        double log_ratio = log_y - log_1_y;
        
        grad_mu[i] = p * (log_ratio - digamma_alpha + digamma_beta);
        grad_phi[i] = digamma_p - m * digamma_alpha - (1.0 - m) * digamma_beta + m * log_y + (1.0 - m) * log_1_y;
    });
    
    return List::create(Named("mu") = grad_mu, Named("phi") = grad_phi);
}

// [[Rcpp::export]]
List beta_hessian_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_phi_phi(n);
    NumericVector hess_mu_phi(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);
    bool both_scalar = mu_is_scalar && phi_is_scalar;
    
    double m0 = 0, p0 = 0, alpha0 = 0, beta_shape0 = 0;
    double digamma_alpha0 = 0, digamma_beta0 = 0;
    double trigamma_alpha0 = 0, trigamma_beta0 = 0, trigamma_p0 = 0;

    if (both_scalar) {
        m0 = mu[0];
        p0 = phi[0];
        alpha0 = m0 * p0;
        beta_shape0 = (1.0 - m0) * p0;
        digamma_alpha0 = R::digamma(alpha0);
        digamma_beta0 = R::digamma(beta_shape0);
        trigamma_alpha0 = R::trigamma(alpha0);
        trigamma_beta0 = R::trigamma(beta_shape0);
        trigamma_p0 = R::trigamma(p0);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m0, p = p0, alpha = alpha0, beta_shape = beta_shape0, digamma_alpha = digamma_alpha0, digamma_beta = digamma_beta0, trigamma_alpha = trigamma_alpha0, trigamma_beta = trigamma_beta0, trigamma_p = trigamma_p0;
        if (!both_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            p = phi_is_scalar ? phi[0] : phi[i];
            alpha = m * p;
            beta_shape = (1.0 - m) * p;
            digamma_alpha = R::digamma(alpha);
            digamma_beta = R::digamma(beta_shape);
            trigamma_alpha = R::trigamma(alpha);
            trigamma_beta = R::trigamma(beta_shape);
            trigamma_p = R::trigamma(p);
        }
        
        double log_ratio = std::log(y[i] / (1.0 - y[i]));
        
        hess_mu_mu[i] = -p * p * (trigamma_alpha + trigamma_beta);
        hess_phi_phi[i] = trigamma_p - m * m * trigamma_alpha - (1.0 - m) * (1.0 - m) * trigamma_beta;
        
        double term1 = log_ratio - digamma_alpha + digamma_beta;
        double term2 = p * (m * trigamma_alpha - (1.0 - m) * trigamma_beta);
        hess_mu_phi[i] = term1 - term2;
    });
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("phi_phi") = hess_phi_phi, Named("mu_phi") = hess_mu_phi);
}

// [[Rcpp::export]]
List beta_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_phi_phi(n);
    NumericVector hess_mu_phi(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);
    bool both_scalar = mu_is_scalar && phi_is_scalar;

    double m0 = 0, p0 = 0, alpha0 = 0, beta_shape0 = 0;
    double trigamma_alpha0 = 0, trigamma_beta0 = 0, trigamma_p0 = 0;
    
    if (both_scalar) {
        m0 = mu[0];
        p0 = phi[0];
        alpha0 = m0 * p0;
        beta_shape0 = (1.0 - m0) * p0;
        trigamma_alpha0 = R::trigamma(alpha0);
        trigamma_beta0 = R::trigamma(beta_shape0);
        trigamma_p0 = R::trigamma(p0);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m0, p = p0, alpha = alpha0, beta_shape = beta_shape0, trigamma_alpha = trigamma_alpha0, trigamma_beta = trigamma_beta0, trigamma_p = trigamma_p0;
        if (!both_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            p = phi_is_scalar ? phi[0] : phi[i];
            alpha = m * p;
            beta_shape = (1.0 - m) * p;
            trigamma_alpha = R::trigamma(alpha);
            trigamma_beta = R::trigamma(beta_shape);
            trigamma_p = R::trigamma(p);
        }
        
        hess_mu_mu[i] = -p * p * (trigamma_alpha + trigamma_beta);
        hess_phi_phi[i] = trigamma_p - m * m * trigamma_alpha - (1.0 - m) * (1.0 - m) * trigamma_beta;
        hess_mu_phi[i] = -p * (m * trigamma_alpha - (1.0 - m) * trigamma_beta);
    });
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("phi_phi") = hess_phi_phi, Named("mu_phi") = hess_mu_phi);
}