#include <Rcpp.h>
#include <cmath>
#include "pt_skewt.h"
using namespace Rcpp;

// The skew t's score and Hessian in (mu, sigma, alpha, nu), one kernel each,
// from the component functions of pt_skewt.h; the off-diagonal components
// are written here. With d = a + q b, dd = da + dq b^2 + q db and
// g = dq b z c + q e,
//   l_ms = (d + z dd)/sigma^2,   l_ma = -g/sigma,   l_sa = -z g/sigma,
// and the three mixed components in nu step the closed-form score in nu,
// one stencil on an analytic quantity.
//
// The t distribution function may signal a warning, which a worker thread
// must not do, so both kernels run on the calling thread whatever `threads`
// says.

namespace {

inline double par_at(const NumericVector& v, std::size_t i) {
    return v.size() == 1 ? v[0] : v[i];
}

}  // namespace

// [[Rcpp::export]]
List skewt_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        NumericVector alpha, NumericVector nu, int threads = 1) {
    int n = y.size();
    NumericVector g_m(n), g_s(n), g_a(n), g_n(n);
    for (int i = 0; i < n; ++i) {
        double m = par_at(mu, i), s = par_at(sigma, i), a = par_at(alpha, i),
               v = par_at(nu, i);
        const d7::SkewtPieces P = d7::skewt_pieces(y[i], m, s, a, v);
        double h = d7::skewt_nu_step(v);
        double lp[5];
        for (int k = 0; k < 5; ++k) {
            if (k == 2) { lp[k] = 0.0; continue; }
            lp[k] = d7::skewt_logpdf(y[i], m, s, a, v + (k - 2) * h);
        }
        g_m[i] = d7::skewt_score_mu(P, s);
        g_s[i] = d7::skewt_score_sigma(P, s);
        g_a[i] = d7::skewt_score_alpha(P);
        g_n[i] = d7::skewt_score_nu(lp, h);
    }
    return List::create(Named("mu") = g_m, Named("sigma") = g_s,
                        Named("alpha") = g_a, Named("nu") = g_n);
}

// [[Rcpp::export]]
List skewt_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                       NumericVector alpha, NumericVector nu, int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_ss(n), h_aa(n), h_nn(n), h_ms(n), h_ma(n),
                  h_mn(n), h_sa(n), h_sn(n), h_an(n);
    for (int i = 0; i < n; ++i) {
        double m = par_at(mu, i), s = par_at(sigma, i), a = par_at(alpha, i),
               v = par_at(nu, i);
        const d7::SkewtPieces P = d7::skewt_pieces(y[i], m, s, a, v);
        double d = P.a + P.q * P.b;
        double dd = P.da + P.dq * P.b * P.b + P.q * P.db;
        double s2 = s * s;
        double zc = P.z * P.c;
        double g = P.dq * P.b * zc + P.q * P.e;
        double h = d7::skewt_nu_step(v);
        double lp[5];
        double gm[5], gs[5], ga[5];
        for (int k = 0; k < 5; ++k) {
            double vk = v + (k - 2) * h;
            lp[k] = d7::skewt_logpdf(y[i], m, s, a, vk);
            if (k == 2) { gm[k] = gs[k] = ga[k] = 0.0; continue; }
            const d7::SkewtPieces Pk = d7::skewt_pieces(y[i], m, s, a, vk);
            gm[k] = d7::skewt_score_mu(Pk, s);
            gs[k] = d7::skewt_score_sigma(Pk, s);
            ga[k] = d7::skewt_score_alpha(Pk);
        }
        h_mm[i] = d7::skewt_hess_mu_mu(P, s);
        h_ss[i] = d7::skewt_hess_sigma_sigma(P, s);
        h_aa[i] = d7::skewt_hess_alpha_alpha(P);
        h_nn[i] = d7::skewt_hess_nu_nu(lp, h);
        h_ms[i] = (d + P.z * dd) / s2;
        h_ma[i] = -g / s;
        h_mn[i] = d7::skewt_score_nu(gm, h);
        h_sa[i] = -P.z * g / s;
        h_sn[i] = d7::skewt_score_nu(gs, h);
        h_an[i] = d7::skewt_score_nu(ga, h);
    }
    return List::create(Named("mu_mu") = h_mm, Named("sigma_sigma") = h_ss,
                        Named("alpha_alpha") = h_aa, Named("nu_nu") = h_nn,
                        Named("mu_sigma") = h_ms, Named("mu_alpha") = h_ma,
                        Named("mu_nu") = h_mn, Named("sigma_alpha") = h_sa,
                        Named("sigma_nu") = h_sn, Named("alpha_nu") = h_an);
}
