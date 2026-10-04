#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// Generalized gamma in Stacy's form, where the nesting is visible:
//   f(y) = p / (a^d Gamma(d/p)) * y^(d-1) * exp(-(y/a)^p),
//   l    = log p - d log a - lgamma(k) + (d-1) log y - w,
// with w = (y/a)^p, L = log(y/a) and k = d/p. Then p = 1 is the gamma with
// shape d and scale a, d = p is the Weibull with shape p and scale a, and
// p -> 0 approaches the lognormal.
//
// The log-density is shared by gengamma1 and gengamma2 (through the scale a
// of its mean); the derivatives of gengamma1 are in gengamma1.cpp, generated
// by stabilita/gen_gengamma1.py.

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
                                  NumericVector d, NumericVector p,
                        int threads = 1) {
    int n = y.size();
    NumericVector out(n);
    bool a_s = (a.size() == 1), d_s = (d.size() == 1), p_s = (p.size() == 1);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double av = a_s ? a[0] : a[i], dv = d_s ? d[0] : d[i],
               pv = p_s ? p[0] : p[i];
        if (y[i] <= 0) { out[i] = R_NegInf; return; }
        GGparts z = gg_parts(y[i], av, dv, pv);
        out[i] = std::log(pv) - dv * std::log(av) - R::lgammafn(z.k)
               + (dv - 1.0) * std::log(y[i]) - z.w;
    });
    return out;
}
