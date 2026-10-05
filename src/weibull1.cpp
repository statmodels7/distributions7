#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
#include "pt_weibull1.h"
using namespace Rcpp;

// The Weibull's score, Hessian and expected information in (mu, sigma), one
// kernel each, from the component functions of pt_weibull1.h; the
// off-diagonal components are written here. With z = y/mu, lz = log z and
// u = z^sigma,
//   l_ms = (u - 1 + sigma u lz)/mu,   E_ms = (1 - gamma)/mu.

// [[Rcpp::export]]
List weibull1_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                           int threads = 1) {
    int n = y.size();
    NumericVector g_mu(n), g_s(n);
    bool m_s = (mu.size() == 1), s_s = (sigma.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *sp = sigma.begin();
    double *a = g_mu.begin(), *b = g_s.begin();
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i], s = s_s ? sp[0] : sp[i];
        double lz = std::log(yp[i] / m), u = std::exp(s * lz);
        a[i] = d7::weibull1_score_mu(u, m, s);
        b[i] = d7::weibull1_score_sigma(lz, u, s);
    });
    return List::create(Named("mu") = g_mu, Named("sigma") = g_s);
}

// [[Rcpp::export]]
List weibull1_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                          int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_ss(n), h_ms(n);
    bool m_s = (mu.size() == 1), s_s = (sigma.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *sp = sigma.begin();
    double *a = h_mm.begin(), *b = h_ss.begin(), *c = h_ms.begin();
    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i], s = s_s ? sp[0] : sp[i];
        double lz = std::log(yp[i] / m), u = std::exp(s * lz);
        a[i] = d7::weibull1_hess_mu_mu(u, m, s);
        b[i] = d7::weibull1_hess_sigma_sigma(lz, u, s);
        c[i] = (u - 1.0 + s * u * lz) / m;
    });
    return List::create(Named("mu_mu") = h_mm, Named("sigma_sigma") = h_ss,
                        Named("mu_sigma") = h_ms);
}

// [[Rcpp::export]]
List weibull1_expected_hessian_cpp(NumericVector y, NumericVector mu,
                                   NumericVector sigma, int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_ss(n), h_ms(n);
    bool m_s = (mu.size() == 1), s_s = (sigma.size() == 1);
    const double *mp = mu.begin(), *sp = sigma.begin();
    double *a = h_mm.begin(), *b = h_ss.begin(), *c = h_ms.begin();
    // the information depends on the parameters alone: once when both are
    // scalars
    bool both = m_s && s_s;
    double a0 = 0.0, b0 = 0.0, c0 = 0.0;
    if (both) {
        double m = mp[0], s = sp[0];
        a0 = d7::weibull1_expected_mu_mu(m, s);
        b0 = d7::weibull1_expected_sigma_sigma(s);
        c0 = (1.0 - d7::kEulerGamma) / m;
    }
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        if (both) {
            a[i] = a0; b[i] = b0; c[i] = c0;
            return;
        }
        double m = m_s ? mp[0] : mp[i], s = s_s ? sp[0] : sp[i];
        a[i] = d7::weibull1_expected_mu_mu(m, s);
        b[i] = d7::weibull1_expected_sigma_sigma(s);
        c[i] = (1.0 - d7::kEulerGamma) / m;
    });
    return List::create(Named("mu_mu") = h_mm, Named("sigma_sigma") = h_ss,
                        Named("mu_sigma") = h_ms);
}
