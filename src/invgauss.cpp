#include <Rcpp.h>
#include "d7_par.h"
#include "pt_invgauss1.h"
using namespace Rcpp;

// [[Rcpp::export]]
List invgauss_gradient_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int threads = 1) {
    int n = y.size();
    NumericVector grad_mu(n);
    NumericVector grad_phi(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);
    
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double p = phi_is_scalar ? phi[0] : phi[i];
        
        grad_mu[i] = d7::invgauss1_score_mu(y[i], m, p);
        grad_phi[i] = d7::invgauss1_score_phi(y[i], m, p);
    });
    
    return List::create(Named("mu") = grad_mu, Named("phi") = grad_phi);
}

// [[Rcpp::export]]
List invgauss_hessian_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_phi_phi(n);
    NumericVector hess_mu_phi(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double p = phi_is_scalar ? phi[0] : phi[i];
        
        double m3 = m * m * m;
        double p2 = p * p;
        double res = y[i] - m;
        hess_mu_mu[i] = d7::invgauss1_hess_mu_mu(y[i], m, p);
        hess_phi_phi[i] = d7::invgauss1_hess_phi_phi(y[i], m, p);
        hess_mu_phi[i] = -res / (p2 * m3);
    });
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("phi_phi") = hess_phi_phi, Named("mu_phi") = hess_mu_phi);
}

// [[Rcpp::export]]
List invgauss_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_phi_phi(n);
    NumericVector hess_mu_phi(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);
    bool both_scalar = mu_is_scalar && phi_is_scalar;
    
    double hmm0 = 0, hpp0 = 0;

    if (both_scalar) {
        double m = mu[0];
        double p = phi[0];
        hmm0 = d7::invgauss1_expected_mu_mu(m, p);
        hpp0 = d7::invgauss1_expected_phi_phi(p);
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double hmm = hmm0, hpp = hpp0;
        if (!both_scalar) {
            double m = mu_is_scalar ? mu[0] : mu[i];
            double p = phi_is_scalar ? phi[0] : phi[i];
            hmm = d7::invgauss1_expected_mu_mu(m, p);
            hpp = d7::invgauss1_expected_phi_phi(p);
        }
        
        hess_mu_mu[i] = hmm;
        hess_phi_phi[i] = hpp;
        hess_mu_phi[i] = 0.0;
    });
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("phi_phi") = hess_phi_phi, Named("mu_phi") = hess_mu_phi);
}