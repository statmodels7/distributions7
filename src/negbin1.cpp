#include <Rcpp.h>
#include "d7_par.h"
#include "psi_diff.h"
#include "pt_negbin1.h"
#include "pt_sqrt.h"
using namespace Rcpp;

// Negative binomial with a variance LINEAR in the mean: Var(Y) = mu (1 + theta),
// against the quadratic mu + mu^2/theta of negbin_distrib(). The two are
// different families rather than two parametrizations of one, and the
// difference shows in where the mean sits: here the size is r = mu/theta, so
// mu appears INSIDE the gamma functions, while in NB2 it stays outside them.
//
// With r = mu/theta and p = 1/(1+theta),
//   l = lgamma(y+r) - lgamma(r) - lgamma(y+1) - r log(1+theta)
//       + y log(theta) - y log(1+theta).
// Writing P = dl/dr and Q = dl/dtheta at fixed r,
//   P    = psi(y+r) - psi(r) - log(1+theta),
//   Q    = -r/(1+theta) + y/theta - y/(1+theta),
//   P_r  = psi'(y+r) - psi'(r),
//   P_th = -1/(1+theta)  (which is also Q_r, the mixed second derivative),
//   Q_th = r/(1+theta)^2 - y/theta^2 + y/(1+theta)^2,
// and the chain rule through r = mu/theta gives the rest.

// NB1parts, nb1_parts() and the support sums are in pt_negbin1.h.

// The three components of the expected hessian at one (mu, theta), from the
// two sums of nb1_sums0() (pt_negbin1.h, which records the derivation)
static inline void nb1_E_parts(double m, double t,
                               double &emm, double &emt, double &ett) {
    double A, W;
    nb1_sums0(m, t, &A, &W);
    emm = d7::negbin1_expected_mu_mu(A);
    emt = d7::negbin1_expected_mu_theta(m, t, W);
    ett = d7::negbin1_expected_theta_theta(t, W);
}

// [[Rcpp::export]]
NumericVector negbin1_logpmf_cpp(NumericVector y, NumericVector mu,
                                 NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector out(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        out[i] = R::dnbinom(y[i], m / t, 1.0 / (1.0 + t), 1);
    });
    return out;
}

// [[Rcpp::export]]
List negbin1_gradient_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector g_mu(n), g_th(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        NB1parts z = nb1_parts(y[i], m, t, 1);
        g_mu[i] = d7::negbin1_score_mu(z);
        g_th[i] = d7::negbin1_score_theta(z);
    });
    return List::create(Named("mu") = g_mu, Named("theta") = g_th);
}

// [[Rcpp::export]]
List negbin1_hessian_cpp(NumericVector y, NumericVector mu, NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_mt(n), h_tt(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        NB1parts z = nb1_parts(y[i], m, t, 2);
        h_mm[i] = d7::negbin1_hess_mu_mu(z);
        h_mt[i] = (z.Pr * z.rt + z.Pth) * z.rm + z.P * z.rmt;
        h_tt[i] = d7::negbin1_hess_theta_theta(z);
    });
    return List::create(Named("mu_mu") = h_mm,
                        Named("mu_theta") = h_mt,
                        Named("theta_theta") = h_tt);
}

// The expected Hessian, from the support sums of pt_negbin1.h.
//
// One observation is computed and written in full by one thread, as every
// other kernel here does. What this replaced was a sequential loop that
// memoized the expectation across CONSECUTIVE equal parameters: that fires
// wherever a design repeats a mean in adjacent rows, and nowhere at all
// under an offset, which gives every row its own -- measured on a design
// of 28800 cells with an offset, 28800 distinct means of 28800, so the memo never hit once
// and the loop could not be parallelized because of it.
// [[Rcpp::export]]
List negbin1_expected_hessian_cpp(NumericVector y, NumericVector mu,
                                  NumericVector theta,
                        int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_mt(n), h_tt(n);
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    bool both_scalar = mu_s && th_s;

    double emm0 = 0, emt0 = 0, ett0 = 0;
    if (both_scalar) nb1_E_parts(mu[0], theta[0], emm0, emt0, ett0);

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double emm = emm0, emt = emt0, ett = ett0;
        if (!both_scalar) {
            double m = mu_s ? mu[0] : mu[i];
            double t = th_s ? theta[0] : theta[i];
            nb1_E_parts(m, t, emm, emt, ett);
        }
        h_mm[i] = emm; h_mt[i] = emt; h_tt[i] = ett;
    });
    return List::create(Named("mu_mu") = h_mm,
                        Named("mu_theta") = h_mt,
                        Named("theta_theta") = h_tt);
}

// ---- the derivatives of the expected information ------------------------
//
// The components are those of nb1_E_parts(), -A, W/theta^2 and
// -W/(mu theta), differentiated with A and W: the first order reads
// nb1_sums1(), the second nb1_sums2() (pt_negbin1.h).

// [[Rcpp::export]]
List negbin1_dexpected1_cpp(NumericVector y, NumericVector mu,
                            NumericVector theta, int threads = 1) {
    int n = y.size();
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    const int K = 6;
    std::vector<NumericVector> v(K);
    std::vector<double*> p(K);
    for (int k = 0; k < K; ++k) { v[k] = NumericVector(n); p[k] = v[k].begin(); }
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        double A[3], W[3];
        nb1_sums1(m, t, A, W);
        p[0][i] = d7::negbin1_dexpected_mu_mu_mu(A);
        p[1][i] = d7::negbin1_dexpected_mu_mu_theta(A);
        p[2][i] = d7::negbin1_dexpected_theta_theta_mu(t, W);
        p[3][i] = d7::negbin1_dexpected_theta_theta_theta(t, W);
        p[4][i] = d7::negbin1_dexpected_mu_theta_mu(m, t, W);
        p[5][i] = d7::negbin1_dexpected_mu_theta_theta(m, t, W);
    });
    CharacterVector nm = CharacterVector::create(
        "mu_mu_mu", "mu_mu_theta", "theta_theta_mu",
        "theta_theta_theta", "mu_theta_mu", "mu_theta_theta");
    List out(K);
    for (int k = 0; k < K; ++k) out[k] = v[k];
    out.attr("names") = nm;
    return out;
}

// [[Rcpp::export]]
List negbin1_dexpected2_cpp(NumericVector y, NumericVector mu,
                            NumericVector theta, int threads = 1) {
    int n = y.size();
    bool mu_s = (mu.size() == 1), th_s = (theta.size() == 1);
    const int K = 9;
    std::vector<NumericVector> v(K);
    std::vector<double*> p(K);
    for (int k = 0; k < K; ++k) { v[k] = NumericVector(n); p[k] = v[k].begin(); }
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = mu_s ? mu[0] : mu[i];
        double t = th_s ? theta[0] : theta[i];
        // A and W as (value, d_mu, d_theta, d_mu,mu, d_mu,theta, d_theta,theta)
        double A[6], W[6];
        nb1_sums2(m, t, A, W);
        double it = 1.0 / t, it2 = it * it, it3 = it2 * it, it4 = it2 * it2;
        double im = 1.0 / m, im2 = im * im, im3 = im2 * im;
        // E_mm = -A
        p[0][i] = -A[3];
        p[1][i] = -A[5];
        p[2][i] = -A[4];
        // E_tt = W/theta^2
        p[3][i] = W[3] * it2;
        p[4][i] = W[5] * it2 - 4.0 * W[2] * it3 + 6.0 * W[0] * it4;
        p[5][i] = W[4] * it2 - 2.0 * W[1] * it3;
        // E_mt = -W/(mu theta)
        p[6][i] = -W[3] * im * it + 2.0 * W[1] * im2 * it - 2.0 * W[0] * im3 * it;
        p[7][i] = -W[5] * im * it + 2.0 * W[2] * im * it2 - 2.0 * W[0] * im * it3;
        p[8][i] = -W[4] * im * it + W[2] * im2 * it + W[1] * im * it2
                  - W[0] * im2 * it2;
    });
    CharacterVector nm = CharacterVector::create(
        "mu_mu_mu_mu", "mu_mu_theta_theta", "mu_mu_mu_theta",
        "theta_theta_mu_mu", "theta_theta_theta_theta",
        "theta_theta_mu_theta", "mu_theta_mu_mu",
        "mu_theta_theta_theta", "mu_theta_mu_theta");
    List out(K);
    for (int k = 0; k < K; ++k) out[k] = v[k];
    out.attr("names") = nm;
    return out;
}
