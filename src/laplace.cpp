#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
#include "pt_laplace.h"
#include "pt_laplace2.h"
using namespace Rcpp;

// The two Laplace parametrizations: score, Hessian, expected information
// and its derivatives, one kernel each, from the component functions of
// pt_laplace.h and pt_laplace2.h; the off-diagonal components are written
// here. With r = y - mu,
//   laplace (mu, b):        l_mb = -sgn(r)/b^2,  E_mb = 0
//   laplace2 (mu, lambda):  l_ml = sgn(r),       E_ml = 0
// and every derivative of an expected information that is not on the
// diagonal is zero, E_mm and E_bb (E_ll) depending on the scale alone.

// --- laplace ------------------------------------------------------------------

// [[Rcpp::export]]
List laplace_gradient_cpp(NumericVector y, NumericVector mu, NumericVector b,
                          int threads = 1) {
    int n = y.size();
    NumericVector g_mu(n), g_b(n);
    bool m_s = (mu.size() == 1), b_s = (b.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *bp = b.begin();
    double *o1 = g_mu.begin(), *o2 = g_b.begin();
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i], bb = b_s ? bp[0] : bp[i];
        double r = yp[i] - m;
        o1[i] = d7::laplace_score_mu(d7::laplace_sign(r), bb);
        o2[i] = d7::laplace_score_sigma(std::fabs(r), bb);
    });
    return List::create(Named("mu") = g_mu, Named("sigma") = g_b);
}

// [[Rcpp::export]]
List laplace_hessian_cpp(NumericVector y, NumericVector mu, NumericVector b,
                         int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_bb(n), h_mb(n);
    bool m_s = (mu.size() == 1), b_s = (b.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *bp = b.begin();
    double *o1 = h_mm.begin(), *o2 = h_bb.begin(), *o3 = h_mb.begin();
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i], bb = b_s ? bp[0] : bp[i];
        double r = yp[i] - m;
        o1[i] = d7::laplace_hess_mu_mu();
        o2[i] = d7::laplace_hess_sigma_sigma(std::fabs(r), bb);
        o3[i] = -d7::laplace_sign(r) / (bb * bb);
    });
    return List::create(Named("mu_mu") = h_mm, Named("sigma_sigma") = h_bb,
                        Named("mu_sigma") = h_mb);
}

// [[Rcpp::export]]
List laplace_expected_hessian_cpp(NumericVector y, NumericVector mu,
                                  NumericVector b, int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_bb(n), h_mb(n);
    bool b_s = (b.size() == 1);
    const double *bp = b.begin();
    double *o1 = h_mm.begin(), *o2 = h_bb.begin();
    // once when the scale is a scalar
    double e10 = 0.0, e20 = 0.0;
    if (b_s) {
        e10 = d7::laplace_expected_mu_mu(bp[0]);
        e20 = d7::laplace_expected_sigma_sigma(bp[0]);
    }
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        if (b_s) {
            o1[i] = e10; o2[i] = e20;
            return;
        }
        double bb = bp[i];
        o1[i] = d7::laplace_expected_mu_mu(bb);
        o2[i] = d7::laplace_expected_sigma_sigma(bb);
    });
    return List::create(Named("mu_mu") = h_mm, Named("sigma_sigma") = h_bb,
                        Named("mu_sigma") = h_mb);
}

// d_b E_mm = d_b E_bb = 2/b^3 at order one, -6/b^4 at order two
// [[Rcpp::export]]
List laplace_dexpected1_cpp(NumericVector y, NumericVector mu, NumericVector b,
                            int threads = 1) {
    int n = y.size();
    NumericVector mm_b(n), bb_b(n);
    bool b_s = (b.size() == 1);
    const double *bp = b.begin();
    double *o1 = mm_b.begin(), *o2 = bb_b.begin();
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double bb = b_s ? bp[0] : bp[i];
        o1[i] = 2.0 / (bb * bb * bb);
        o2[i] = d7::laplace_dexpected_sigma_sigma_sigma(bb);
    });
    NumericVector zero(n);
    return List::create(
        Named("mu_mu_mu") = zero, Named("mu_mu_sigma") = mm_b,
        Named("sigma_sigma_mu") = clone(zero), Named("sigma_sigma_sigma") = bb_b,
        Named("mu_sigma_mu") = clone(zero), Named("mu_sigma_sigma") = clone(zero));
}

// [[Rcpp::export]]
List laplace_dexpected2_cpp(NumericVector y, NumericVector mu, NumericVector b,
                            int threads = 1) {
    int n = y.size();
    NumericVector mm_bb(n), bb_bb(n);
    bool b_s = (b.size() == 1);
    const double *bp = b.begin();
    double *o1 = mm_bb.begin(), *o2 = bb_bb.begin();
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double bb = b_s ? bp[0] : bp[i];
        double b2 = bb * bb;
        o1[i] = -6.0 / (b2 * b2);
        o2[i] = -6.0 / (b2 * b2);
    });
    NumericVector zero(n);
    return List::create(
        Named("mu_mu_mu_mu") = zero, Named("mu_mu_sigma_sigma") = mm_bb,
        Named("mu_mu_mu_sigma") = clone(zero),
        Named("sigma_sigma_mu_mu") = clone(zero),
        Named("sigma_sigma_sigma_sigma") = bb_bb,
        Named("sigma_sigma_mu_sigma") = clone(zero),
        Named("mu_sigma_mu_mu") = clone(zero),
        Named("mu_sigma_sigma_sigma") = clone(zero),
        Named("mu_sigma_mu_sigma") = clone(zero));
}

// --- laplace2 -----------------------------------------------------------------

// [[Rcpp::export]]
List laplace2_gradient_cpp(NumericVector y, NumericVector mu, NumericVector lambda,
                           int threads = 1) {
    int n = y.size();
    NumericVector g_mu(n), g_l(n);
    bool m_s = (mu.size() == 1), l_s = (lambda.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *lp = lambda.begin();
    double *o1 = g_mu.begin(), *o2 = g_l.begin();
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i], lam = l_s ? lp[0] : lp[i];
        double r = yp[i] - m;
        o1[i] = d7::laplace2_score_mu(d7::laplace_sign(r), lam);
        o2[i] = d7::laplace2_score_lambda(std::fabs(r), lam);
    });
    return List::create(Named("mu") = g_mu, Named("lambda") = g_l);
}

// [[Rcpp::export]]
List laplace2_hessian_cpp(NumericVector y, NumericVector mu, NumericVector lambda,
                          int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_ll(n), h_ml(n);
    bool m_s = (mu.size() == 1), l_s = (lambda.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *lp = lambda.begin();
    double *o1 = h_mm.begin(), *o2 = h_ll.begin(), *o3 = h_ml.begin();
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i], lam = l_s ? lp[0] : lp[i];
        o1[i] = d7::laplace2_hess_mu_mu();
        o2[i] = d7::laplace2_hess_lambda_lambda(lam);
        o3[i] = d7::laplace_sign(yp[i] - m);
    });
    return List::create(Named("mu_mu") = h_mm, Named("lambda_lambda") = h_ll,
                        Named("mu_lambda") = h_ml);
}

// [[Rcpp::export]]
List laplace2_expected_hessian_cpp(NumericVector y, NumericVector mu,
                                   NumericVector lambda, int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_ll(n), h_ml(n);
    bool l_s = (lambda.size() == 1);
    const double *lp = lambda.begin();
    double *o1 = h_mm.begin(), *o2 = h_ll.begin();
    // once when the rate is a scalar
    double e10 = 0.0, e20 = 0.0;
    if (l_s) {
        e10 = d7::laplace2_expected_mu_mu(lp[0]);
        e20 = d7::laplace2_expected_lambda_lambda(lp[0]);
    }
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        if (l_s) {
            o1[i] = e10; o2[i] = e20;
            return;
        }
        double lam = lp[i];
        o1[i] = d7::laplace2_expected_mu_mu(lam);
        o2[i] = d7::laplace2_expected_lambda_lambda(lam);
    });
    return List::create(Named("mu_mu") = h_mm, Named("lambda_lambda") = h_ll,
                        Named("mu_lambda") = h_ml);
}

// E_mm = -lambda^2, E_ll = -1/lambda^2: d_l E_mm = -2 lambda,
// d_l E_ll = 2/lambda^3; d_ll E_mm = -2, d_ll E_ll = -6/lambda^4
// [[Rcpp::export]]
List laplace2_dexpected1_cpp(NumericVector y, NumericVector mu,
                             NumericVector lambda, int threads = 1) {
    int n = y.size();
    NumericVector mm_l(n), ll_l(n);
    bool l_s = (lambda.size() == 1);
    const double *lp = lambda.begin();
    double *o1 = mm_l.begin(), *o2 = ll_l.begin();
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double lam = l_s ? lp[0] : lp[i];
        o1[i] = -2.0 * lam;
        o2[i] = d7::laplace2_dexpected_lambda_lambda_lambda(lam);
    });
    NumericVector zero(n);
    return List::create(
        Named("mu_mu_mu") = zero, Named("mu_mu_lambda") = mm_l,
        Named("lambda_lambda_mu") = clone(zero),
        Named("lambda_lambda_lambda") = ll_l,
        Named("mu_lambda_mu") = clone(zero), Named("mu_lambda_lambda") = clone(zero));
}

// [[Rcpp::export]]
List laplace2_dexpected2_cpp(NumericVector y, NumericVector mu,
                             NumericVector lambda, int threads = 1) {
    int n = y.size();
    NumericVector mm_ll(n), ll_ll(n);
    bool l_s = (lambda.size() == 1);
    const double *lp = lambda.begin();
    double *o1 = mm_ll.begin(), *o2 = ll_ll.begin();
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double lam = l_s ? lp[0] : lp[i];
        double l2 = lam * lam;
        o1[i] = -2.0;
        o2[i] = -6.0 / (l2 * l2);
    });
    NumericVector zero(n);
    return List::create(
        Named("mu_mu_mu_mu") = zero, Named("mu_mu_lambda_lambda") = mm_ll,
        Named("mu_mu_mu_lambda") = clone(zero),
        Named("lambda_lambda_mu_mu") = clone(zero),
        Named("lambda_lambda_lambda_lambda") = ll_ll,
        Named("lambda_lambda_mu_lambda") = clone(zero),
        Named("mu_lambda_mu_mu") = clone(zero),
        Named("mu_lambda_lambda_lambda") = clone(zero),
        Named("mu_lambda_mu_lambda") = clone(zero));
}
