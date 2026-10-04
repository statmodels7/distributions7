#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include "pt_pseudohuber.h"
#include "pt_sqrt.h"
using namespace Rcpp;

// The pseudo-Huber score and Hessian in (mu, sigma, nu), one kernel each,
// from the component functions of pt_pseudohuber.h; the off-diagonal
// components are written here. The Bessel ratios depend on nu alone, so they
// are computed once when every parameter is a scalar.

namespace d7 {

// Resolved on first use, on the calling thread. Not a guarded static, for
// the reason vm_bessel() gives in vonmises.cpp.
static N7BesselK n7bk_ptr = nullptr;

N7BesselK bessel_k_fn() {
    if (n7bk_ptr == nullptr)
        n7bk_ptr = (N7BesselK) R_GetCCallable("numericals7", "n7_bessel_k");
    return n7bk_ptr;
}

}  // namespace d7

// [[Rcpp::export]]
List pseudohuber_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu) {
    int n = y.size();
    NumericVector grad_mu(n), grad_sigma(n), grad_nu(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    bool nu_is_scalar = (nu.size() == 1);
    bool all_scalar = mu_is_scalar && sigma_is_scalar && nu_is_scalar;

    double m = 0, s = 0, v = 0, s2 = 0;
    d7::PhNu P = {0, 0, 0, 0, 0, 0, 0, 0};

    if (all_scalar) {
        m = mu[0]; s = sigma[0]; v = nu[0];
        s2 = s * s;
        P = d7::ph_nu(v);
    }

    for(int i = 0; i < n; i++) {
        if (!all_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            s = sigma_is_scalar ? sigma[0] : sigma[i];
            v = nu_is_scalar ? nu[0] : nu[i];
            s2 = s * s;
            P = d7::ph_nu(v);
        }

        double res = y[i] - m;
        double res2 = res * res;
        double D = d7::sqrt_cr(v + res2 / s2);

        grad_mu[i] = d7::pseudohuber_score_mu(res, s2, D);
        grad_sigma[i] = d7::pseudohuber_score_sigma(res2, s, s2, D);
        grad_nu[i] = d7::pseudohuber_score_nu(v, D, P);
    }

    return List::create(Named("mu") = grad_mu, Named("sigma") = grad_sigma, Named("nu") = grad_nu);
}

// [[Rcpp::export]]
List pseudohuber_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma, NumericVector nu) {
    int n = y.size();
    NumericVector hess_mu_mu(n), hess_sigma_sigma(n), hess_nu_nu(n);
    NumericVector hess_mu_sigma(n), hess_mu_nu(n), hess_sigma_nu(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    bool nu_is_scalar = (nu.size() == 1);
    bool all_scalar = mu_is_scalar && sigma_is_scalar && nu_is_scalar;

    double m = 0, s = 0, v = 0, s2 = 0, s4 = 0;
    d7::PhNu P = {0, 0, 0, 0, 0, 0, 0, 0};

    if (all_scalar) {
        m = mu[0]; s = sigma[0]; v = nu[0];
        s2 = s * s; s4 = s2 * s2;
        P = d7::ph_nu(v);
    }

    for(int i = 0; i < n; i++) {
        if (!all_scalar) {
            m = mu_is_scalar ? mu[0] : mu[i];
            s = sigma_is_scalar ? sigma[0] : sigma[i];
            v = nu_is_scalar ? nu[0] : nu[i];
            s2 = s * s; s4 = s2 * s2;
            P = d7::ph_nu(v);
        }

        double res = y[i] - m;
        double res2 = res * res;
        double D = d7::sqrt_cr(v + res2 / s2);
        double D3 = D * D * D;

        hess_mu_mu[i] = d7::pseudohuber_hess_mu_mu(v, s2, D3);
        hess_sigma_sigma[i] = d7::pseudohuber_hess_sigma_sigma(res2, s2, s4, D, D3);
        hess_nu_nu[i] = d7::pseudohuber_hess_nu_nu(D3, P);

        hess_mu_sigma[i] = (-2.0 * v * s2 * res - res2 * res) / (s2 * std::pow(v * s2 + res2, 1.5));
        hess_mu_nu[i] = -res / (2.0 * s2 * D3);
        hess_sigma_nu[i] = -res2 / (2.0 * s2 * s * D3);
    }

    return List::create(
        Named("mu_mu") = hess_mu_mu, Named("sigma_sigma") = hess_sigma_sigma, Named("nu_nu") = hess_nu_nu,
        Named("mu_sigma") = hess_mu_sigma, Named("mu_nu") = hess_mu_nu, Named("sigma_nu") = hess_sigma_nu
    );
}
