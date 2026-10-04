#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <cmath>
#include "d7_par.h"
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
// The t distribution function is numericals7's n7_pt(), resolved before the
// loop starts, so the kernels run on `threads` threads.

namespace d7 {

// Resolved on first use, on the calling thread. Not a guarded static, for
// the reason vm_bessel() gives in vonmises.cpp.
static N7Pt n7pt_ptr = nullptr;

N7Pt skewt_pt() {
    if (n7pt_ptr == nullptr)
        n7pt_ptr = (N7Pt) R_GetCCallable("numericals7", "n7_pt");
    return n7pt_ptr;
}

}  // namespace d7

// [[Rcpp::export]]
List skewt_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        NumericVector alpha, NumericVector nu, int threads = 1) {
    int n = y.size();
    NumericVector g_m(n), g_s(n), g_a(n), g_n(n);
    bool m_s = (mu.size() == 1), s_s = (sigma.size() == 1),
         a_s = (alpha.size() == 1), v_s = (nu.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *sp = sigma.begin(),
                 *ap = alpha.begin(), *vp = nu.begin();
    double *o1 = g_m.begin(), *o2 = g_s.begin(), *o3 = g_a.begin(),
           *o4 = g_n.begin();
    d7::skewt_pt();
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i], s = s_s ? sp[0] : sp[i],
               a = a_s ? ap[0] : ap[i], v = v_s ? vp[0] : vp[i];
        const d7::SkewtPieces P = d7::skewt_pieces(yp[i], m, s, a, v);
        double h = d7::skewt_nu_step(v);
        double lp[5];
        for (int k = 0; k < 5; ++k) {
            if (k == 2) { lp[k] = 0.0; continue; }
            lp[k] = d7::skewt_logpdf(yp[i], m, s, a, v + (k - 2) * h);
        }
        o1[i] = d7::skewt_score_mu(P, s);
        o2[i] = d7::skewt_score_sigma(P, s);
        o3[i] = d7::skewt_score_alpha(P);
        o4[i] = d7::skewt_score_nu(lp, h);
    });
    return List::create(Named("mu") = g_m, Named("sigma") = g_s,
                        Named("alpha") = g_a, Named("nu") = g_n);
}

// [[Rcpp::export]]
List skewt_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                       NumericVector alpha, NumericVector nu, int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_ss(n), h_aa(n), h_nn(n), h_ms(n), h_ma(n),
                  h_mn(n), h_sa(n), h_sn(n), h_an(n);
    bool m_s = (mu.size() == 1), s_s = (sigma.size() == 1),
         a_s = (alpha.size() == 1), v_s = (nu.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *sp = sigma.begin(),
                 *ap = alpha.begin(), *vp = nu.begin();
    double *o_mm = h_mm.begin(), *o_ss = h_ss.begin(), *o_aa = h_aa.begin(),
           *o_nn = h_nn.begin(), *o_ms = h_ms.begin(), *o_ma = h_ma.begin(),
           *o_mn = h_mn.begin(), *o_sa = h_sa.begin(), *o_sn = h_sn.begin(),
           *o_an = h_an.begin();
    d7::skewt_pt();
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i], s = s_s ? sp[0] : sp[i],
               a = a_s ? ap[0] : ap[i], v = v_s ? vp[0] : vp[i];
        const d7::SkewtPieces P = d7::skewt_pieces(yp[i], m, s, a, v);
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
            lp[k] = d7::skewt_logpdf(yp[i], m, s, a, vk);
            if (k == 2) { gm[k] = gs[k] = ga[k] = 0.0; continue; }
            const d7::SkewtPieces Pk = d7::skewt_pieces(yp[i], m, s, a, vk);
            gm[k] = d7::skewt_score_mu(Pk, s);
            gs[k] = d7::skewt_score_sigma(Pk, s);
            ga[k] = d7::skewt_score_alpha(Pk);
        }
        o_mm[i] = d7::skewt_hess_mu_mu(P, s);
        o_ss[i] = d7::skewt_hess_sigma_sigma(P, s);
        o_aa[i] = d7::skewt_hess_alpha_alpha(P);
        o_nn[i] = d7::skewt_hess_nu_nu(lp, h);
        o_ms[i] = (d + P.z * dd) / s2;
        o_ma[i] = -g / s;
        o_mn[i] = d7::skewt_score_nu(gm, h);
        o_sa[i] = -P.z * g / s;
        o_sn[i] = d7::skewt_score_nu(gs, h);
        o_an[i] = d7::skewt_score_nu(ga, h);
    });
    return List::create(Named("mu_mu") = h_mm, Named("sigma_sigma") = h_ss,
                        Named("alpha_alpha") = h_aa, Named("nu_nu") = h_nn,
                        Named("mu_sigma") = h_ms, Named("mu_alpha") = h_ma,
                        Named("mu_nu") = h_mn, Named("sigma_alpha") = h_sa,
                        Named("sigma_nu") = h_sn, Named("alpha_nu") = h_an);
}
