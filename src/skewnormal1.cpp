#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
#include "pt_skewnormal1.h"
using namespace Rcpp;

// The skew normal's score and Hessian in (mu, sigma, alpha), one kernel
// each, from the component functions of pt_skewnormal1.h; the off-diagonal
// components are written here. With z = (y - mu)/sigma, r the Mills ratio
// at alpha z and dr its derivative,
//   l_ms = (alpha^2 z dr - 2z + alpha r)/sigma^2,
//   l_ma = -(r + alpha z dr)/sigma,   l_sa = -(z r + alpha z^2 dr)/sigma.

// [[Rcpp::export]]
List skewnormal1_gradient_cpp(NumericVector y, NumericVector mu,
                              NumericVector sigma, NumericVector alpha,
                              int threads = 1) {
    int n = y.size();
    NumericVector g_m(n), g_s(n), g_a(n);
    bool m_s = (mu.size() == 1), s_s = (sigma.size() == 1),
         a_s = (alpha.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *sp = sigma.begin(),
                 *ap = alpha.begin();
    double *o1 = g_m.begin(), *o2 = g_s.begin(), *o3 = g_a.begin();
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i], s = s_s ? sp[0] : sp[i],
               a = a_s ? ap[0] : ap[i];
        double z = (yp[i] - m) / s;
        double r = d7::mills_ratio(a * z);
        o1[i] = d7::skewnormal1_score_mu(z, s, a, r);
        o2[i] = d7::skewnormal1_score_sigma(z, s, a, r);
        o3[i] = d7::skewnormal1_score_alpha(z, r);
    });
    return List::create(Named("mu") = g_m, Named("sigma") = g_s,
                        Named("alpha") = g_a);
}

// [[Rcpp::export]]
List skewnormal1_hessian_cpp(NumericVector y, NumericVector mu,
                             NumericVector sigma, NumericVector alpha,
                             int threads = 1) {
    int n = y.size();
    NumericVector h_mm(n), h_ss(n), h_aa(n), h_ms(n), h_ma(n), h_sa(n);
    bool m_s = (mu.size() == 1), s_s = (sigma.size() == 1),
         a_s = (alpha.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *sp = sigma.begin(),
                 *ap = alpha.begin();
    double *o1 = h_mm.begin(), *o2 = h_ss.begin(), *o3 = h_aa.begin(),
           *o4 = h_ms.begin(), *o5 = h_ma.begin(), *o6 = h_sa.begin();
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = m_s ? mp[0] : mp[i], s = s_s ? sp[0] : sp[i],
               a = a_s ? ap[0] : ap[i];
        double z = (yp[i] - m) / s;
        double t = a * z;
        double r = d7::mills_ratio(t);
        double dr = -r * (t + r);
        double s2 = s * s;
        o1[i] = d7::skewnormal1_hess_mu_mu(s, a, dr);
        o2[i] = d7::skewnormal1_hess_sigma_sigma(z, s, a, r, dr);
        o3[i] = d7::skewnormal1_hess_alpha_alpha(z, dr);
        o4[i] = (a * a * z * dr - 2.0 * z + a * r) / s2;
        o5[i] = -(r + a * z * dr) / s;
        o6[i] = -(z * r + a * z * z * dr) / s;
    });
    return List::create(Named("mu_mu") = h_mm, Named("sigma_sigma") = h_ss,
                        Named("alpha_alpha") = h_aa, Named("mu_sigma") = h_ms,
                        Named("mu_alpha") = h_ma, Named("sigma_alpha") = h_sa);
}
