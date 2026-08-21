#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// Gaussian in the mean and the PRECISION, l = log(t)/2 - log(2 pi)/2 - t r^2/2,
// with r = y - mu. This is the flattest of the three parametrizations: the
// third derivatives are free of the data except through nothing at all, and
// the fourth has a single non-zero component.

// [[Rcpp::export]]
List gaussian3_gradient_cpp(NumericVector y, NumericVector mu, NumericVector tau,
                        int threads = 1) {
    int n = y.size();
    NumericVector g_mu(n), g_t(n);
    bool m_s = (mu.size() == 1), t_s = (tau.size() == 1);

    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = m_s ? mu[0] : mu[i];
        double t = t_s ? tau[0] : tau[i];
        double r = y[i] - m;
        g_mu[i] = t * r;
        g_t[i] = 0.5 / t - 0.5 * r * r;
    });
    return List::create(Named("mu") = g_mu, Named("tau") = g_t);
}

// [[Rcpp::export]]
List gaussian3_hessian_cpp(NumericVector y, NumericVector mu, NumericVector tau,
                        int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_mt(n), h_tt(n);
    bool m_s = (mu.size() == 1), t_s = (tau.size() == 1);

    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = m_s ? mu[0] : mu[i];
        double t = t_s ? tau[0] : tau[i];
        h_mm[i] = -t;
        h_mt[i] = y[i] - m;
        h_tt[i] = -0.5 / (t * t);
    });
    return List::create(Named("mu_mu") = h_mm, Named("mu_tau") = h_mt,
                        Named("tau_tau") = h_tt);
}

// [[Rcpp::export]]
List gaussian3_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector tau,
                        int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_mt(n), h_tt(n);
    bool t_s = (tau.size() == 1);

    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double t = t_s ? tau[0] : tau[i];
        h_mm[i] = -t;
        h_mt[i] = 0.0;
        h_tt[i] = -0.5 / (t * t);
    });
    return List::create(Named("mu_mu") = h_mm, Named("mu_tau") = h_mt,
                        Named("tau_tau") = h_tt);
}

// Every third derivative is free of y, so the observed and the expected ones
// coincide and one routine serves both.
// [[Rcpp::export]]
List gaussian3_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector tau,
                        int threads = 1) {
    int n = y.size();
    NumericVector a(n), b(n), c(n), d(n);
    bool t_s = (tau.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double t = t_s ? tau[0] : tau[i];
        a[i] = 0.0;
        b[i] = -1.0;
        c[i] = 0.0;
        d[i] = 1.0 / (t * t * t);
    });
    return List::create(Named("mu_mu_mu") = a, Named("mu_mu_tau") = b,
                        Named("mu_tau_tau") = c, Named("tau_tau_tau") = d);
}

// [[Rcpp::export]]
List gaussian3_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector tau,
                        int threads = 1) {
    int n = y.size();
    NumericVector a(n), b(n), c(n), d(n), e(n);
    bool t_s = (tau.size() == 1);

    d7::par_for(n, threads, d7::kMinMid, [&](std::size_t i) {
        double t = t_s ? tau[0] : tau[i];
        double t2 = t * t;
        a[i] = 0.0;
        b[i] = 0.0;
        c[i] = 0.0;
        d[i] = 0.0;
        e[i] = -3.0 / (t2 * t2);
    });
    return List::create(Named("mu_mu_mu_mu") = a, Named("mu_mu_mu_tau") = b,
                        Named("mu_mu_tau_tau") = c, Named("mu_tau_tau_tau") = d,
                        Named("tau_tau_tau_tau") = e);
}
