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

// The derivatives of the expected information in the parameters. With
// alpha = mu phi, beta = (1 - mu) phi, n = 1 - mu, and a_k = psi^(k)(alpha),
// b_k = psi^(k)(beta), c_k = psi^(k)(phi):
//
//   E_mm = -phi^2 (a1 + b1)
//   E_pp = c1 - mu^2 a1 - n^2 b1
//   E_mp = -phi (mu a1 - n b1)
//
// and d alpha/d mu = phi, d beta/d mu = -phi, d alpha/d phi = mu,
// d beta/d phi = n. Every component below is the ordinary derivative of those
// three; psigamma at orders 2 and 3 is admissible in a worker since d7_par.h
// installs the calling thread's floating-point environment.
// [[Rcpp::export]]
List beta_dexpected_cpp(NumericVector y, NumericVector mu, NumericVector phi,
                        int order, int threads = 1) {
    int n_obs = y.size();
    bool m_s = (mu.size() == 1), p_s = (phi.size() == 1);
    const double *mp = mu.begin(), *pp = phi.begin();
    if (order == 1) {
        NumericVector mm_m(n_obs), mm_p(n_obs), pp_m(n_obs), pp_p(n_obs),
                      mp_m(n_obs), mp_p(n_obs);
        double *o1 = mm_m.begin(), *o2 = mm_p.begin(), *o3 = pp_m.begin(),
               *o4 = pp_p.begin(), *o5 = mp_m.begin(), *o6 = mp_p.begin();
        d7::par_for(n_obs, threads, d7::kMinCostly, [&](std::size_t i) {
            double m = m_s ? mp[0] : mp[i];
            double p = p_s ? pp[0] : pp[i];
            double n = 1.0 - m, al = m * p, be = n * p;
            double a1 = R::trigamma(al), b1 = R::trigamma(be);
            double a2 = R::psigamma(al, 2), b2 = R::psigamma(be, 2);
            double c2 = R::psigamma(p, 2);
            o1[i] = -p * p * p * (a2 - b2);
            o2[i] = -2.0 * p * (a1 + b1) - p * p * (m * a2 + n * b2);
            o3[i] = -2.0 * m * a1 - m * m * p * a2 + 2.0 * n * b1 + n * n * p * b2;
            o4[i] = c2 - m * m * m * a2 - n * n * n * b2;
            o5[i] = -p * (a1 + b1) - p * p * (m * a2 + n * b2);
            o6[i] = -(m * a1 - n * b1) - p * (m * m * a2 - n * n * b2);
        });
        return List::create(
            Named("mu_mu_mu") = mm_m, Named("mu_mu_phi") = mm_p,
            Named("phi_phi_mu") = pp_m, Named("phi_phi_phi") = pp_p,
            Named("mu_phi_mu") = mp_m, Named("mu_phi_phi") = mp_p);
    }
    NumericVector v[9];
    for (int k = 0; k < 9; ++k) v[k] = NumericVector(n_obs);
    double *o[9];
    for (int k = 0; k < 9; ++k) o[k] = v[k].begin();
    d7::par_for(n_obs, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i];
        double p = p_s ? pp[0] : pp[i];
        double n = 1.0 - m, al = m * p, be = n * p;
        double a1 = R::trigamma(al), b1 = R::trigamma(be);
        double a2 = R::psigamma(al, 2), b2 = R::psigamma(be, 2);
        double a3 = R::psigamma(al, 3), b3 = R::psigamma(be, 3);
        double c3 = R::psigamma(p, 3);
        double p2 = p * p, p3 = p2 * p;
        // E_mm
        o[0][i] = -p2 * p2 * (a3 + b3);                                   // mu mu
        o[1][i] = -2.0 * (a1 + b1) - 4.0 * p * (m * a2 + n * b2)
                  - p2 * (m * m * a3 + n * n * b3);                       // phi phi
        o[2][i] = -3.0 * p2 * (a2 - b2) - p3 * (m * a3 - n * b3);         // mu phi
        // E_pp
        o[3][i] = -2.0 * a1 - 4.0 * m * p * a2 - m * m * p2 * a3
                  - 2.0 * b1 - 4.0 * n * p * b2 - n * n * p2 * b3;        // mu mu
        o[4][i] = c3 - m * m * m * m * a3 - n * n * n * n * b3;           // phi phi
        o[5][i] = -3.0 * m * m * a2 - m * m * m * p * a3
                  + 3.0 * n * n * b2 + n * n * n * p * b3;                // mu phi
        // E_mp
        o[6][i] = -2.0 * p2 * (a2 - b2) - p3 * (m * a3 - n * b3);         // mu mu
        o[7][i] = -2.0 * (m * m * a2 - n * n * b2)
                  - p * (m * m * m * a3 - n * n * n * b3);                // phi phi
        o[8][i] = -(a1 + b1) - 3.0 * p * (m * a2 + n * b2)
                  - p2 * (m * m * a3 + n * n * b3);                       // mu phi
    });
    return List::create(
        Named("mu_mu_mu_mu") = v[0], Named("mu_mu_phi_phi") = v[1],
        Named("mu_mu_mu_phi") = v[2],
        Named("phi_phi_mu_mu") = v[3], Named("phi_phi_phi_phi") = v[4],
        Named("phi_phi_mu_phi") = v[5],
        Named("mu_phi_mu_mu") = v[6], Named("mu_phi_phi_phi") = v[7],
        Named("mu_phi_mu_phi") = v[8]);
}