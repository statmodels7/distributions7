#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
#include "pt_gumbel.h"
using namespace Rcpp;

// The Gumbel's score, Hessian and expected information, one kernel each,
// from the component functions of pt_gumbel.h; the off-diagonal components
// are written here. With z = (y - mu)/sigma and w = exp(-z),
//   l_ms = -(1 - w + z w)/sigma^2.

// [[Rcpp::export]]
List gumbel_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                         int threads = 1) {
    int n = y.size();
    NumericVector g_mu(n), g_s(n);
    bool m_s = (mu.size() == 1), s_s = (sigma.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *sp = sigma.begin();
    double *a = g_mu.begin(), *b = g_s.begin();
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i], s = s_s ? sp[0] : sp[i];
        double z = (yp[i] - m) / s, w = std::exp(-z);
        a[i] = d7::gumbel_score_mu(w, s);
        b[i] = d7::gumbel_score_sigma(z, w, s);
    });
    return List::create(Named("mu") = g_mu, Named("sigma") = g_s);
}

// [[Rcpp::export]]
List gumbel_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_ss(n), h_ms(n);
    bool m_s = (mu.size() == 1), s_s = (sigma.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *sp = sigma.begin();
    double *a = h_mm.begin(), *b = h_ss.begin(), *c = h_ms.begin();
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i], s = s_s ? sp[0] : sp[i];
        double z = (yp[i] - m) / s, w = std::exp(-z);
        a[i] = d7::gumbel_hess_mu_mu(w, s);
        b[i] = d7::gumbel_hess_sigma_sigma(z, w, s);
        c[i] = -(1.0 - w + z * w) / (s * s);
    });
    return List::create(Named("mu_mu") = h_mm, Named("sigma_sigma") = h_ss,
                        Named("mu_sigma") = h_ms);
}

// [[Rcpp::export]]
List gumbel_expected_hessian_cpp(NumericVector y, NumericVector mu,
                                 NumericVector sigma, int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_ss(n), h_ms(n);
    bool s_s = (sigma.size() == 1);
    const double *sp = sigma.begin();
    double *a = h_mm.begin(), *b = h_ss.begin(), *c = h_ms.begin();
    // the information depends on sigma alone: once when it is a scalar
    double a0 = 0.0, b0 = 0.0, c0 = 0.0;
    if (s_s) {
        double s = sp[0];
        a0 = d7::gumbel_expected_mu_mu(s);
        b0 = d7::gumbel_expected_sigma_sigma(s);
        c0 = (1.0 - d7::kEulerGamma) / (s * s);
    }
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        if (s_s) {
            a[i] = a0; b[i] = b0; c[i] = c0;
            return;
        }
        double s = sp[i];
        a[i] = d7::gumbel_expected_mu_mu(s);
        b[i] = d7::gumbel_expected_sigma_sigma(s);
        c[i] = (1.0 - d7::kEulerGamma) / (s * s);
    });
    return List::create(Named("mu_mu") = h_mm, Named("sigma_sigma") = h_ss,
                        Named("mu_sigma") = h_ms);
}
