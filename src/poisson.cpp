#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// [[Rcpp::export]]
List poisson_gradient_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector grad_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    double m0 = 0;
    
    if (mu_is_scalar) m0 = mu[0];

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m0;
        if (!mu_is_scalar) {
            m = mu[i];
        }
        grad_mu[i] = (y[i] - m) / m;
    });
    
    return List::create(Named("mu") = grad_mu);
}

// [[Rcpp::export]]
List poisson_hessian_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    double m20 = 0;
    
    if (mu_is_scalar) {
        double m = mu[0];
        m20 = m * m;
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m2 = m20;
        if (!mu_is_scalar) {
            m2 = mu[i] * mu[i];
        }
        hess_mu_mu[i] = -y[i] / m2;
    });
    
    return List::create(Named("mu_mu") = hess_mu_mu);
}

// [[Rcpp::export]]
List poisson_expected_hessian_cpp(NumericVector y, NumericVector mu,
                        int threads = 1) {
    int n = y.size();
    NumericVector hess_mu_mu(n);
    
    bool mu_is_scalar = (mu.size() == 1);
    double val0 = 0;
    
    if (mu_is_scalar) {
        val0 = -1.0 / mu[0];
    }

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double val = val0;
        if (!mu_is_scalar) {
            val = -1.0 / mu[i];
        }
        hess_mu_mu[i] = val;
    });

    return List::create(Named("mu_mu") = hess_mu_mu);
}

// The derivatives of E_mm = -1/mu: 1/mu^2 at order 1, -2/mu^3 at order 2.
// [[Rcpp::export]]
List poisson_dexpected_cpp(NumericVector y, NumericVector mu, int order,
                           int threads = 1) {
    int n = y.size();
    NumericVector out(n);
    bool mu_is_scalar = (mu.size() == 1);
    const double *mp = mu.begin();
    double *o = out.begin();
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = mu_is_scalar ? mp[0] : mp[i];
        double inv = 1.0 / m;
        o[i] = (order == 1) ? inv * inv : -2.0 * inv * inv * inv;
    });
    if (order == 1) return List::create(Named("mu_mu_mu") = out);
    return List::create(Named("mu_mu_mu_mu") = out);
}