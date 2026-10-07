#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <cmath>
#include <string>
#include <vector>
#include "d7_par.h"
#include "pt_skewt.h"
#include "pt_skewt_exact.h"
using namespace Rcpp;

// The skew t's derivatives of orders one to five in (mu, sigma, alpha, nu),
// one kernel per order, each the generated closed forms of
// pt_skewt_table.h at one observation (pt_skewt_exact.h, which records the
// derivation). Every component is exact in nu: the derivatives of the t
// distribution function in its degrees of freedom are integrals, taken by
// quadrature, and nothing is differenced.
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

namespace {

// the components of order K, in the order of pt_skewt_table.h: the
// multisets of the four parameters in lexicographic order, as deriv_names()
std::vector<std::string> skewt_names(int K) {
    static const char* P[4] = {"mu", "sigma", "alpha", "nu"};
    std::vector<std::string> out;
    std::vector<int> idx(K, 0);
    while (true) {
        std::string nm = P[idx[0]];
        for (int i = 1; i < K; ++i) nm += std::string("_") + P[idx[i]];
        out.push_back(nm);
        int k = K - 1;
        while (k >= 0 && idx[k] == 3) --k;
        if (k < 0) break;
        ++idx[k];
        for (int j = k + 1; j < K; ++j) idx[j] = idx[k];
    }
    return out;
}

template <int K>
List skewt_kernel(NumericVector y, NumericVector mu, NumericVector sigma,
                  NumericVector alpha, NumericVector nu, int threads,
                  const std::vector<int>& order) {
    const int n = y.size();
    const std::vector<std::string> nm = skewt_names(K);
    const int nc = nm.size();
    std::vector<NumericVector> v(nc);
    std::vector<double*> p(nc);
    for (int c = 0; c < nc; ++c) { v[c] = NumericVector(n); p[c] = v[c].begin(); }
    const bool m_s = (mu.size() == 1), s_s = (sigma.size() == 1),
               a_s = (alpha.size() == 1), v_s = (nu.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *sp = sigma.begin(),
                 *ap = alpha.begin(), *vp = nu.begin();
    d7::skewt_pt();
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        thread_local d7::SkewtConst cc;
        const double nv = v_s ? vp[0] : vp[i];
        cc.at(nv, K);
        double out[56];
        d7::skewt_obs<K>(yp[i], m_s ? mp[0] : mp[i], s_s ? sp[0] : sp[i],
                         a_s ? ap[0] : ap[i], nv, cc, out);
        for (int c = 0; c < nc; ++c) p[c][i] = out[c];
    });
    List res(nc);
    CharacterVector names(nc);
    for (int c = 0; c < nc; ++c) {
        const int src = order.empty() ? c : order[c];
        res[c] = v[src];
        names[c] = nm[src];
    }
    res.attr("names") = names;
    return res;
}

}  // namespace

// [[Rcpp::export]]
List skewt_gradient_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                        NumericVector alpha, NumericVector nu, int threads = 1) {
    return skewt_kernel<1>(y, mu, sigma, alpha, nu, threads, {});
}

// [[Rcpp::export]]
List skewt_hessian_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                       NumericVector alpha, NumericVector nu, int threads = 1) {
    // hess_names() order: the diagonal first
    return skewt_kernel<2>(y, mu, sigma, alpha, nu, threads,
                           {0, 4, 7, 9, 1, 2, 3, 5, 6, 8});
}

// [[Rcpp::export]]
List skewt_deriv3_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                      NumericVector alpha, NumericVector nu, int threads = 1) {
    return skewt_kernel<3>(y, mu, sigma, alpha, nu, threads, {});
}

// [[Rcpp::export]]
List skewt_deriv4_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                      NumericVector alpha, NumericVector nu, int threads = 1) {
    return skewt_kernel<4>(y, mu, sigma, alpha, nu, threads, {});
}

// [[Rcpp::export]]
List skewt_deriv5_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                      NumericVector alpha, NumericVector nu, int threads = 1) {
    return skewt_kernel<5>(y, mu, sigma, alpha, nu, threads, {});
}

// (d_m^j T_m(w))/T_m(w), j = 1..5, the ratios the kernels read, for tests
// [[Rcpp::export]]
NumericMatrix skewt_tdf_ratio_cpp(NumericVector w, double m) {
    const int n = w.size();
    NumericMatrix out(n, 5);
    d7::skewt_pt();
    double Cm[7];
    d7::skewt_cder(m, 5, Cm);
    for (int i = 0; i < n; ++i) {
        const double wi = w[i], w2 = wi * wi;
        const double q = std::exp(R::dt(wi, m, 1) - d7::skewt_pt()(wi, m, 1, 1));
        double ell[6], Pw[6], G[6] = {0, 0, 0, 0, 0, 0};
        d7::skewt_ell<5>(w2, m, 1.0 / m, 0.5 * (m + 1.0), std::log1p(w2 / m), Cm, ell);
        d7::skewt_bell<5>(ell, Pw);
        d7::skewt_R0<5>(wi, m, q, Cm, Pw, G);
        for (int j = 1; j <= 5; ++j) out(i, j - 1) = G[j] + Pw[j];
    }
    return out;
}
