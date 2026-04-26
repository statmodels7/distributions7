#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List student_t_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu) {
    int n = y.size();
    NumericVector grad_mu(n), grad_sigma(n), grad_nu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    bool nu_is_scalar = (nu.size() == 1);
    bool all_scalar = mu_is_scalar && sigma_is_scalar && nu_is_scalar;

    double m = 0, s = 0, v = 0, s2 = 0, v_s2 = 0;
    double digamma_v_half = 0, digamma_v1_half = 0;

    if (all_scalar) {
        m = mu[0]; s = sigma[0]; v = nu[0];
        s2 = s * s; v_s2 = v * s2;
        digamma_v_half = R::digamma(0.5 * v);
        digamma_v1_half = R::digamma(0.5 * (v + 1.0));
    }

    for(int i = 0; i < n; i++) {
        if (!all_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            s = sigma_is_scalar ? sigma[0] : sigma[i];
            v = nu_is_scalar ? nu[0] : nu[i];
            s2 = s * s; v_s2 = v * s2;
            digamma_v_half = R::digamma(0.5 * v);
            digamma_v1_half = R::digamma(0.5 * (v + 1.0));
        }
        
        double res = y[i] - m;
        double res2 = res * res;
        double den = v_s2 + res2;
        
        grad_mu[i] = ((v + 1.0) * res) / den;
        grad_sigma[i] = (v * (res2 - s2)) / (s * den);
        grad_nu[i] = 0.5 * (-1.0 / v - digamma_v_half + digamma_v1_half + 
                            ((v + 1.0) * res2) / (v * den) - 
                            std::log(res2 / v_s2 + 1.0));
    }
    
    return List::create(Named("mu") = grad_mu, Named("sigma") = grad_sigma, Named("nu") = grad_nu);
}

// [[Rcpp::export]]
List student_t_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu) {
    int n = y.size();
    NumericVector hess_mu_mu(n), hess_sigma_sigma(n), hess_nu_nu(n);
    NumericVector hess_mu_sigma(n), hess_mu_nu(n), hess_sigma_nu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    bool nu_is_scalar = (nu.size() == 1);
    bool all_scalar = mu_is_scalar && sigma_is_scalar && nu_is_scalar;

    double m = 0, s = 0, v = 0, s2 = 0, s4 = 0, v_s2 = 0, v_s4 = 0;
    double trigamma_v_half = 0, trigamma_v1_half = 0;

    if (all_scalar) {
        m = mu[0]; s = sigma[0]; v = nu[0];
        s2 = s * s; s4 = s2 * s2;
        v_s2 = v * s2; v_s4 = v * s4;
        trigamma_v_half = R::trigamma(0.5 * v);
        trigamma_v1_half = R::trigamma(0.5 * (v + 1.0));
    }

    for(int i = 0; i < n; i++) {
        if (!all_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            s = sigma_is_scalar ? sigma[0] : sigma[i];
            v = nu_is_scalar ? nu[0] : nu[i];
            s2 = s * s; s4 = s2 * s2;
            v_s2 = v * s2; v_s4 = v * s4;
            trigamma_v_half = R::trigamma(0.5 * v);
            trigamma_v1_half = R::trigamma(0.5 * (v + 1.0));
        }
        
        double res = y[i] - m;
        double res2 = res * res;
        double res4 = res2 * res2;
        double den = v_s2 + res2;
        double den2 = den * den;
        
        hess_mu_mu[i] = (v + 1.0) * (res2 - v_s2) / den2;
        hess_sigma_sigma[i] = (v * (v_s4 - (3.0 * v + 1.0) * s2 * res2 - res4)) / (s2 * den2);
        hess_nu_nu[i] = 0.25 * (-trigamma_v_half + trigamma_v1_half + (2.0 * (v_s4 + res4)) / (v * den2));
        
        hess_mu_sigma[i] = -(2.0 * v * (v + 1.0) * s * res) / den2;
        hess_mu_nu[i] = (res * (res2 - s2)) / den2;
        hess_sigma_nu[i] = (res4 - s2 * res2) / (s * den2);
    }
    
    return List::create(
        Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("nu_nu") = hess_nu_nu,
        Named("mu_sigma") = hess_mu_sigma, Named("mu_nu") = hess_mu_nu, Named("sigma_nu") = hess_sigma_nu
    );
}

// [[Rcpp::export]]
List student_t_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu) {
    int n = y.size();
    NumericVector hess_mu_mu(n), hess_sigma_sigma(n), hess_nu_nu(n);
    NumericVector hess_mu_sigma(n), hess_mu_nu(n), hess_sigma_nu(n);
    
    bool sigma_is_scalar = (sigma.size() == 1);
    bool nu_is_scalar = (nu.size() == 1);
    bool both_scalar = sigma_is_scalar && nu_is_scalar;

    double s = 0, v = 0, s2 = 0;
    double trigamma_v_half = 0, trigamma_v1_half = 0;

    if (both_scalar) {
        s = sigma[0]; v = nu[0];
        s2 = s * s;
        trigamma_v_half = R::trigamma(0.5 * v);
        trigamma_v1_half = R::trigamma(0.5 * (v + 1.0));
    }

    for(int i = 0; i < n; i++) {
        if (!both_scalar) {
            s = sigma_is_scalar ? sigma[0] : sigma[i];
            v = nu_is_scalar ? nu[0] : nu[i];
            s2 = s * s;
            trigamma_v_half = R::trigamma(0.5 * v);
            trigamma_v1_half = R::trigamma(0.5 * (v + 1.0));
        }
        
        hess_mu_mu[i] = -(v + 1.0) / (s2 * (v + 3.0));
        hess_sigma_sigma[i] = -(2.0 * v) / (s2 * (v + 3.0));
        hess_nu_nu[i] = 0.25 * (trigamma_v1_half - trigamma_v_half) + (v + 5.0) / (2.0 * v * (v + 1.0) * (v + 3.0));
        
        hess_sigma_nu[i] = 2.0 / ((v + 1.0) * (v + 3.0) * s);
        // mu_sigma and mu_nu default initialized to 0.0, no need to assign
    }
    
    return List::create(
        Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("nu_nu") = hess_nu_nu,
        Named("mu_sigma") = hess_mu_sigma, Named("mu_nu") = hess_mu_nu, Named("sigma_nu") = hess_sigma_nu
    );
}