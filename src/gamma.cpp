#include <Rcpp.h>
#include "d7_par.h"
#include "pt_gamma2.h"
using namespace Rcpp;

// The vector kernels of the gamma in the mean and the variance, orders one
// and two. Each component is d7::gamma2_component() (pt_gamma2.h) at its
// multi-index; the order-two components are returned in the order of
// hess_names(): mu_mu, sigma2_sigma2, mu_sigma2.

// The components (i, j) listed in idx, of one order, observed or expected,
// into the columns of out.
static void gamma2_fill(const NumericVector& y, const NumericVector& mu,
                        const NumericVector& sigma2, int order,
                        const int (*idx)[2], int ncomp, bool expected,
                        std::vector<NumericVector>& out, int threads) {
    const int n = y.size();
    const bool m_s = (mu.size() == 1), v_s = (sigma2.size() == 1);
    const bool both = m_s && v_s;
    double Gd0[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    if (both) d7::gamma2_shape_rests(mu[0] * mu[0] / sigma2[0], order, Gd0, expected);
    d7::gamma2_terms(1, 0);   // builds the static table before the threads start
    std::vector<double*> po(ncomp);
    for (int c = 0; c < ncomp; ++c) po[c] = out[c].begin();
    const double *yp = y.begin(), *mp = mu.begin(), *vp = sigma2.begin();

    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i];
        double v = v_s ? vp[0] : vp[i];
        d7::Gamma2Tables T;
        for (int r = 0; r < 5; ++r) T.Gd[r] = Gd0[r];
        if (!both) d7::gamma2_shape_rests(m * m / v, order, T.Gd, expected);
        d7::gamma2_fill_tables(yp[i], m, v, order, expected, T);
        for (int c = 0; c < ncomp; ++c) {
            po[c][i] = d7::gamma2_component(idx[c][0], idx[c][1], T);
        }
    });
}

static const int G2_IDX1[2][2] = {{1, 0}, {0, 1}};
static const int G2_IDX2[3][2] = {{2, 0}, {0, 2}, {1, 1}};

// [[Rcpp::export]]
List gamma_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    int n = y.size();
    std::vector<NumericVector> o = {NumericVector(n), NumericVector(n)};
    gamma2_fill(y, mu, sigma2, 1, G2_IDX1, 2, false, o, threads);
    return List::create(Named("mu") = o[0], Named("sigma2") = o[1]);
}

// [[Rcpp::export]]
List gamma_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    int n = y.size();
    std::vector<NumericVector> o = {NumericVector(n), NumericVector(n),
                                    NumericVector(n)};
    gamma2_fill(y, mu, sigma2, 2, G2_IDX2, 3, false, o, threads);
    return List::create(Named("mu_mu") = o[0], Named("sigma2_sigma2") = o[1],
                        Named("mu_sigma2") = o[2]);
}

// [[Rcpp::export]]
List gamma_expected_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    int n = y.size();
    std::vector<NumericVector> o = {NumericVector(n), NumericVector(n),
                                    NumericVector(n)};
    gamma2_fill(y, mu, sigma2, 2, G2_IDX2, 3, true, o, threads);
    return List::create(Named("mu_mu") = o[0], Named("sigma2_sigma2") = o[1],
                        Named("mu_sigma2") = o[2]);
}

static const int G2_IDX3[4][2] = {{3, 0}, {2, 1}, {1, 2}, {0, 3}};
static const int G2_IDX4[5][2] = {{4, 0}, {3, 1}, {2, 2}, {1, 3}, {0, 4}};

static List gamma2_order3(NumericVector y, NumericVector mu,
                          NumericVector sigma2, bool expected, int threads) {
    int n = y.size();
    std::vector<NumericVector> o(4);
    for (auto& c : o) c = NumericVector(n);
    gamma2_fill(y, mu, sigma2, 3, G2_IDX3, 4, expected, o, threads);
    return List::create(Named("mu_mu_mu") = o[0], Named("mu_mu_sigma2") = o[1],
                        Named("mu_sigma2_sigma2") = o[2],
                        Named("sigma2_sigma2_sigma2") = o[3]);
}

static List gamma2_order4(NumericVector y, NumericVector mu,
                          NumericVector sigma2, bool expected, int threads) {
    int n = y.size();
    std::vector<NumericVector> o(5);
    for (auto& c : o) c = NumericVector(n);
    gamma2_fill(y, mu, sigma2, 4, G2_IDX4, 5, expected, o, threads);
    return List::create(Named("mu_mu_mu_mu") = o[0],
                        Named("mu_mu_mu_sigma2") = o[1],
                        Named("mu_mu_sigma2_sigma2") = o[2],
                        Named("mu_sigma2_sigma2_sigma2") = o[3],
                        Named("sigma2_sigma2_sigma2_sigma2") = o[4]);
}

// [[Rcpp::export]]
List gamma_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    return gamma2_order3(y, mu, sigma2, false, threads);
}

// [[Rcpp::export]]
List gamma_deriv3_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    return gamma2_order3(y, mu, sigma2, true, threads);
}

// [[Rcpp::export]]
List gamma_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    return gamma2_order4(y, mu, sigma2, false, threads);
}

// [[Rcpp::export]]
List gamma_deriv4_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma2,
                        int threads = 1) {
    return gamma2_order4(y, mu, sigma2, true, threads);
}
