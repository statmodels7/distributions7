#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List invgauss_gradient_cpp(NumericVector y, NumericVector mu, NumericVector phi) {
    int n = y.size();
    NumericVector grad_mu(n);
    NumericVector grad_phi(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);
    
    for(int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double p = phi_is_scalar ? phi[0] : phi[i];
        
        double m2 = m * m;
        double m3 = m2 * m;
        double p2 = p * p;
        
        double res = y[i] - m;
        
        grad_mu[i] = res / (p * m3);
        grad_phi[i] = (res * res - y[i] * m2 * p) / (2.0 * y[i] * p2 * m2);
    }
    
    return List::create(Named("mu") = grad_mu, Named("phi") = grad_phi);
}

// [[Rcpp::export]]
List invgauss_hessian_cpp(NumericVector y, NumericVector mu, NumericVector phi) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_phi_phi(n);
    NumericVector hess_mu_phi(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);

    for(int i = 0; i < n; i++) {
        double m = mu_is_scalar ? mu[0] : mu[i];
        double p = phi_is_scalar ? phi[0] : phi[i];
        
        double m2 = m * m;
        double m3 = m2 * m;
        double m4 = m2 * m2;
        double p2 = p * p;
        double p3 = p2 * p;
        
        double res = y[i] - m;
        double res2 = res * res;
        
        hess_mu_mu[i] = -(3.0 * y[i] - 2.0 * m) / (p * m4);
        hess_phi_phi[i] = (p - 2.0 * res2 / (m2 * y[i])) / (2.0 * p3);
        hess_mu_phi[i] = -res / (p2 * m3);
    }
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("phi_phi") = hess_phi_phi, Named("mu_phi") = hess_mu_phi);
}

// [[Rcpp::export]]
List invgauss_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector phi) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    NumericVector hess_phi_phi(n);
    NumericVector hess_mu_phi(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool phi_is_scalar = (phi.size() == 1);
    bool both_scalar = mu_is_scalar && phi_is_scalar;
    
    double hmm = 0, hpp = 0;

    if (both_scalar) {
        double m = mu[0];
        double p = phi[0];
        double m3 = m * m * m;
        double p2 = p * p;
        
        hmm = -1.0 / (p * m3);
        hpp = -0.5 / p2;
    }

    for(int i = 0; i < n; i++) {
        if (!both_scalar) {
            double m = mu_is_scalar ? mu[0] : mu[i];
            double p = phi_is_scalar ? phi[0] : phi[i];
            double m3 = m * m * m;
            double p2 = p * p;
            
            hmm = -1.0 / (p * m3);
            hpp = -0.5 / p2;
        }
        
        hess_mu_mu[i] = hmm;
        hess_phi_phi[i] = hpp;
        hess_mu_phi[i] = 0.0;
    }
    
    return List::create(Named("mu_mu") = hess_mu_mu, Named("phi_phi") = hess_phi_phi, Named("mu_phi") = hess_mu_phi);
}