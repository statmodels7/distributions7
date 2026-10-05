#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
#include "pt_beta2.h"
using namespace Rcpp;

// The beta's score, Hessian and expected information in its shapes, one
// kernel each, from the component functions of pt_beta2.h. The second
// derivatives are free of y, l_ab = psi'(alpha + beta), so the Hessian and
// the expected information are the same numbers; the polygammas are formed
// once when both shapes are scalars.

// [[Rcpp::export]]
List beta2_gradient_cpp(NumericVector y, NumericVector alpha, NumericVector beta,
                        int threads = 1) {
    int n = y.size();
    NumericVector g_a(n), g_b(n);
    bool a_s = (alpha.size() == 1), b_s = (beta.size() == 1);
    bool both = a_s && b_s;
    const double *yp = y.begin(), *ap = alpha.begin(), *bp = beta.begin();
    double a00 = 0.0, b00 = 0.0, s00 = 0.0;
    if (both) {
        a00 = R::digamma(ap[0]);
        b00 = R::digamma(bp[0]);
        s00 = R::digamma(ap[0] + bp[0]);
    }
    double *o1 = g_a.begin(), *o2 = g_b.begin();
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double a0 = a00, b0 = b00, s0 = s00;
        if (!both) {
            double a = a_s ? ap[0] : ap[i], b = b_s ? bp[0] : bp[i];
            a0 = R::digamma(a);
            b0 = R::digamma(b);
            s0 = R::digamma(a + b);
        }
        o1[i] = d7::beta2_score_alpha(std::log(yp[i]), a0, s0);
        o2[i] = d7::beta2_score_beta(std::log1p(-yp[i]), b0, s0);
    });
    return List::create(Named("alpha") = g_a, Named("beta") = g_b);
}

// the Hessian and the expected information share this body
static List beta2_second(int n, NumericVector alpha, NumericVector beta,
                         int threads, bool expected) {
    NumericVector h_aa(n), h_ab(n), h_bb(n);
    bool a_s = (alpha.size() == 1), b_s = (beta.size() == 1);
    bool both = a_s && b_s;
    const double *ap = alpha.begin(), *bp = beta.begin();
    double a10 = 0.0, b10 = 0.0, s10 = 0.0;
    if (both) {
        a10 = R::trigamma(ap[0]);
        b10 = R::trigamma(bp[0]);
        s10 = R::trigamma(ap[0] + bp[0]);
    }
    double *o1 = h_aa.begin(), *o2 = h_ab.begin(), *o3 = h_bb.begin();
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double a1 = a10, b1 = b10, s1 = s10;
        if (!both) {
            double a = a_s ? ap[0] : ap[i], b = b_s ? bp[0] : bp[i];
            a1 = R::trigamma(a);
            b1 = R::trigamma(b);
            s1 = R::trigamma(a + b);
        }
        o1[i] = expected ? d7::beta2_expected_alpha_alpha(a1, s1)
                         : d7::beta2_hess_alpha_alpha(a1, s1);
        o2[i] = s1;
        o3[i] = expected ? d7::beta2_expected_beta_beta(b1, s1)
                         : d7::beta2_hess_beta_beta(b1, s1);
    });
    return List::create(Named("alpha_alpha") = h_aa, Named("beta_beta") = h_bb,
                        Named("alpha_beta") = h_ab);
}

// [[Rcpp::export]]
List beta2_hessian_cpp(NumericVector y, NumericVector alpha, NumericVector beta,
                       int threads = 1) {
    return beta2_second(y.size(), alpha, beta, threads, false);
}

// [[Rcpp::export]]
List beta2_expected_hessian_cpp(NumericVector y, NumericVector alpha,
                                NumericVector beta, int threads = 1) {
    return beta2_second(y.size(), alpha, beta, threads, true);
}
