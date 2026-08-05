#include <Rcpp.h>
using namespace Rcpp;

// Generalized gamma in Stacy's form, where the nesting is visible:
//   f(y) = p / (a^d Gamma(d/p)) * y^(d-1) * exp(-(y/a)^p),
//   l    = log p - d log a - lgamma(k) + (d-1) log y - w,
// with w = (y/a)^p, L = log(y/a) and k = d/p. Then p = 1 is the gamma with
// shape d and scale a, d = p is the Weibull with shape p and scale a, and
// p -> 0 approaches the lognormal.
//
// Every expectation is a moment of u = w ~ Gamma(k, 1):
//   E[u] = k, E[log u] = psi(k), E[(log u)^2] = psi(k)^2 + psi'(k),
//   E[u log u] = k psi(k+1), E[u (log u)^2] = k (psi(k+1)^2 + psi'(k+1)),
// which is what makes the expected information closed form.

struct GGparts { double w, L, k; };

static inline GGparts gg_parts(double y, double a, double d, double p) {
    GGparts z;
    z.L = std::log(y) - std::log(a);
    z.w = std::exp(p * z.L);
    z.k = d / p;
    return z;
}

// [[Rcpp::export]]
NumericVector gengamma_logpdf_cpp(NumericVector y, NumericVector a,
                                  NumericVector d, NumericVector p) {
    int n = y.size();
    NumericVector out(n);
    bool a_s = (a.size() == 1), d_s = (d.size() == 1), p_s = (p.size() == 1);
    for (int i = 0; i < n; i++) {
        double av = a_s ? a[0] : a[i], dv = d_s ? d[0] : d[i],
               pv = p_s ? p[0] : p[i];
        if (y[i] <= 0) { out[i] = R_NegInf; continue; }
        GGparts z = gg_parts(y[i], av, dv, pv);
        out[i] = std::log(pv) - dv * std::log(av) - R::lgammafn(z.k)
               + (dv - 1.0) * std::log(y[i]) - z.w;
    }
    return out;
}

// [[Rcpp::export]]
List gengamma_gradient_cpp(NumericVector y, NumericVector a, NumericVector d,
                           NumericVector p) {
    int n = y.size();
    NumericVector g_a(n), g_d(n), g_p(n);
    bool a_s = (a.size() == 1), d_s = (d.size() == 1), p_s = (p.size() == 1);
    for (int i = 0; i < n; i++) {
        double av = a_s ? a[0] : a[i], dv = d_s ? d[0] : d[i],
               pv = p_s ? p[0] : p[i];
        GGparts z = gg_parts(y[i], av, dv, pv);
        double psi = R::digamma(z.k);
        g_a[i] = (pv * z.w - dv) / av;
        g_d[i] = z.L - psi / pv;
        g_p[i] = 1.0 / pv + dv * psi / (pv * pv) - z.w * z.L;
    }
    return List::create(Named("a") = g_a, Named("d") = g_d, Named("p") = g_p);
}

// [[Rcpp::export]]
List gengamma_hessian_cpp(NumericVector y, NumericVector a, NumericVector d,
                          NumericVector p) {
    int n = y.size();
    NumericVector h_aa(n), h_ad(n), h_ap(n), h_dd(n), h_dp(n), h_pp(n);
    bool a_s = (a.size() == 1), d_s = (d.size() == 1), p_s = (p.size() == 1);
    for (int i = 0; i < n; i++) {
        double av = a_s ? a[0] : a[i], dv = d_s ? d[0] : d[i],
               pv = p_s ? p[0] : p[i];
        GGparts z = gg_parts(y[i], av, dv, pv);
        double psi = R::digamma(z.k), psi1 = R::trigamma(z.k);
        double p2 = pv * pv, p3 = p2 * pv, p4 = p3 * pv;

        h_aa[i] = (-p2 * z.w - pv * z.w + dv) / (av * av);
        h_ad[i] = -1.0 / av;
        h_ap[i] = (z.w + pv * z.w * z.L) / av;
        h_dd[i] = -psi1 / p2;
        h_dp[i] = psi / p2 + dv * psi1 / p3;
        h_pp[i] = -1.0 / p2 - 2.0 * dv * psi / p3 - dv * dv * psi1 / p4
                - z.w * z.L * z.L;
    }
    return List::create(Named("a_a") = h_aa, Named("a_d") = h_ad,
                        Named("a_p") = h_ap, Named("d_d") = h_dd,
                        Named("d_p") = h_dp, Named("p_p") = h_pp);
}

// [[Rcpp::export]]
List gengamma_expected_hessian_cpp(NumericVector y, NumericVector a,
                                   NumericVector d, NumericVector p) {
    int n = y.size();
    NumericVector h_aa(n), h_ad(n), h_ap(n), h_dd(n), h_dp(n), h_pp(n);
    bool a_s = (a.size() == 1), d_s = (d.size() == 1), p_s = (p.size() == 1);
    for (int i = 0; i < n; i++) {
        double av = a_s ? a[0] : a[i], dv = d_s ? d[0] : d[i],
               pv = p_s ? p[0] : p[i];
        double k = dv / pv;
        double psi = R::digamma(k), psi1 = R::trigamma(k);
        double psiA = R::digamma(k + 1.0), psi1A = R::trigamma(k + 1.0);
        double p2 = pv * pv, p3 = p2 * pv, p4 = p3 * pv;

        // E[w] = k, E[w L] = k psi(k+1)/p, E[w L^2] = k(psi(k+1)^2+psi'(k+1))/p^2
        h_aa[i] = -pv * dv / (av * av);
        h_ad[i] = -1.0 / av;
        h_ap[i] = k * (1.0 + psiA) / av;
        h_dd[i] = -psi1 / p2;
        h_dp[i] = psi / p2 + dv * psi1 / p3;
        h_pp[i] = -1.0 / p2 - 2.0 * dv * psi / p3 - dv * dv * psi1 / p4
                - k * (psiA * psiA + psi1A) / p2;
    }
    return List::create(Named("a_a") = h_aa, Named("a_d") = h_ad,
                        Named("a_p") = h_ap, Named("d_d") = h_dd,
                        Named("d_p") = h_dp, Named("p_p") = h_pp);
}
