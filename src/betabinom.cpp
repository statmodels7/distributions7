#include <Rcpp.h>
#include "d7_par.h"
#include "psi_diff.h"
#include "pt_betabinom.h"
#include <limits>
using namespace Rcpp;
using d7::BBderiv;
using d7::BBmap;
using d7::bb_map;
using d7::bb_shape_derivs;
using d7::bb_log_mass;

// Beta-binomial with size n, written in the mean proportion mu and the
// dispersion sigma. The shape parameters are A = mu/sigma and
// B = (1-mu)/sigma, so S = A + B = 1/sigma, and
//   l = lchoose(n, y) + lgamma(y+A) + lgamma(n-y+B) - lgamma(n+S)
//                     - lgamma(A) - lgamma(B) + lgamma(S).
// The whole dependence on the parameters is through A and B, where every
// derivative is a difference of polygammas:
//   l_A = psi(y+A) - psi(A) - [psi(n+S) - psi(S)],
// and likewise for B, with l_AB carrying only the S part. The chain rule to
// (mu, sigma) is then two variables deep and written out below.
//
// Var(Y) = n mu (1-mu) [1 + (n-1) sigma/(1+sigma)], so sigma -> 0 recovers
// the binomial and the family is overdispersed for every sigma > 0.

// [[Rcpp::export]]
List betabinom_gradient_cpp(NumericVector y, NumericVector mu,
                            NumericVector sigma, double size,
                        int threads = 1) {
    int n = y.size();
    NumericVector g_mu(n), g_sigma(n);
    bool mu_s = (mu.size() == 1), si_s = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double s = si_s ? sigma[0] : sigma[i];
        BBmap mp = bb_map(m, s);
        d7::BBd1 d = d7::bb_shape_d1(y[i], size, mp.A, mp.B);
        g_mu[i]    = d7::betabinom1_score_mu(d, mp);
        g_sigma[i] = d7::betabinom1_score_sigma(d, mp);
    });
    return List::create(Named("mu") = g_mu, Named("sigma") = g_sigma);
}

// [[Rcpp::export]]
List betabinom_hessian_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma, double size,
                        int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_ms(n), h_ss(n);
    bool mu_s = (mu.size() == 1), si_s = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double s = si_s ? sigma[0] : sigma[i];
        BBmap mp = bb_map(m, s);
        BBderiv d = bb_shape_derivs(y[i], size, mp.A, mp.B);

        h_mm[i] = d7::betabinom1_hess_mu_mu(d, mp);
        h_ms[i] = d7::betabinom1_hess_mu_sigma(d, mp);
        h_ss[i] = d7::betabinom1_hess_sigma_sigma(d, mp);
    });
    return List::create(Named("mu_mu") = h_mm,
                        Named("mu_sigma") = h_ms,
                        Named("sigma_sigma") = h_ss);
}

// The expected information, by summing the observed Hessian against the mass
// over the finite support {0, ..., n}. The support being finite makes this an
// EXACT expectation rather than a quadrature, and the result does not depend
// on the data, only on the parameters.
// [[Rcpp::export]]
List betabinom_expected_hessian_cpp(NumericVector y, NumericVector mu,
                                    NumericVector sigma, double size,
                        int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_ms(n), h_ss(n);
    bool mu_s = (mu.size() == 1), si_s = (sigma.size() == 1);
    int N = (int) size;

    double last_m = R_NegInf, last_s = R_NegInf;
    double e_mm = 0, e_ms = 0, e_ss = 0;

    // this loop stays SEQUENTIAL: it memoizes an expensive expectation
    // across consecutive observations (last_m), which is state the
    // parallel decomposition may not share.
    (void) threads;
    for (int i = 0; i < n; i++) {
        double m = mu_s ? mu[0] : mu[i];
        double s = si_s ? sigma[0] : sigma[i];

        if (m != last_m || s != last_s) {
            BBmap mp = bb_map(m, s);
            e_mm = 0; e_ms = 0; e_ss = 0;
            for (int k = 0; k <= N; k++) {
                double p = std::exp(bb_log_mass(k, mp.A, mp.B, size));
                BBderiv d = bb_shape_derivs(k, size, mp.A, mp.B);
                e_mm += p * d7::betabinom1_hess_mu_mu(d, mp);
                e_ms += p * d7::betabinom1_hess_mu_sigma(d, mp);
                e_ss += p * d7::betabinom1_hess_sigma_sigma(d, mp);
            }
            last_m = m; last_s = s;
        }
        h_mm[i] = e_mm; h_ms[i] = e_ms; h_ss[i] = e_ss;
    }
    return List::create(Named("mu_mu") = h_mm,
                        Named("mu_sigma") = h_ms,
                        Named("sigma_sigma") = h_ss);
}

// [[Rcpp::export]]
NumericVector betabinom_logpmf_cpp(NumericVector y, NumericVector mu,
                                   NumericVector sigma, double size,
                        int threads = 1) {
    int n = y.size();
    NumericVector out(n);
    bool mu_s = (mu.size() == 1), si_s = (sigma.size() == 1);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double s = si_s ? sigma[0] : sigma[i];
        double A = m / s, B = (1.0 - m) / s;
        if (y[i] < 0 || y[i] > size || y[i] != std::floor(y[i])) {
            out[i] = R_NegInf;
        } else {
            out[i] = bb_log_mass(y[i], A, B, size);
        }
    });
    return out;
}

// The shape parametrization (alpha, beta): the score and the Hessian are the
// shapes' derivatives of bb_shape_derivs() themselves.
// [[Rcpp::export]]
List betabinom2_gradient_cpp(NumericVector y, NumericVector alpha,
                             NumericVector beta, double size, int threads = 1) {
    int n = y.size();
    NumericVector g_a(n), g_b(n);
    bool a_s = (alpha.size() == 1), b_s = (beta.size() == 1);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double a = a_s ? alpha[0] : alpha[i];
        double b = b_s ? beta[0] : beta[i];
        d7::BBd1 d = d7::bb_shape_d1_mixed(y[i], size, a, b);
        g_a[i] = d.lA;
        g_b[i] = d.lB;
    });
    return List::create(Named("alpha") = g_a, Named("beta") = g_b);
}

// [[Rcpp::export]]
List betabinom2_hessian_cpp(NumericVector y, NumericVector alpha,
                            NumericVector beta, double size, int threads = 1) {
    int n = y.size();
    NumericVector h_aa(n), h_bb(n), h_ab(n);
    bool a_s = (alpha.size() == 1), b_s = (beta.size() == 1);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double a = a_s ? alpha[0] : alpha[i];
        double b = b_s ? beta[0] : beta[i];
        d7::BBd2 d = d7::bb_shape_d2_mixed(y[i], size, a, b);
        h_aa[i] = d.lAA;
        h_bb[i] = d.lBB;
        h_ab[i] = d.lAB;
    });
    return List::create(Named("alpha_alpha") = h_aa, Named("beta_beta") = h_bb,
                        Named("alpha_beta") = h_ab);
}
