#include <Rcpp.h>
#include <cmath>
#include <cstring>
#include <vector>
#include "d7_par.h"
#include "pt_loc_scale.h"
#include "pt_skewnormal1.h"
#include "pt_skewt.h"
#include "pt_pseudohuber.h"
#include "pt_pseudohuber2.h"
#include "pt_skewnormal2.h"
using namespace Rcpp;

// The exp-sinh rule of loc_scale_rule() (expected_loc_scale.R), built once
// when the library is loaded, the same arithmetic in the same order: on
// (0, Inf) the nodes exp(c0 sinh(t)) at t = i h, h = 1/32, c0 = pi/2, for
// |x| from 1e-16 to 1e60, mirrored onto the negative half.

namespace {

struct RuleData {
    std::vector<double> x, w;
    RuleData() {
        const double c0 = M_PI / 2.0;
        const double h = 1.0 / 32.0;
        const double tlo = -std::asinh(16.0 * std::log(10.0) / c0);
        const double thi = std::asinh(60.0 * std::log(10.0) / c0);
        const int lo = (int) std::ceil(tlo / h), hi = (int) std::floor(thi / h);
        const int m = hi - lo + 1;
        std::vector<double> xp(m), wp(m);
        for (int i = 0; i < m; ++i) {
            double t = (double) (lo + i) * h;
            xp[i] = std::exp(c0 * std::sinh(t));
            wp[i] = h * c0 * std::cosh(t) * xp[i];
        }
        x.resize(2 * m);
        w.resize(2 * m);
        for (int i = 0; i < m; ++i) {
            x[i] = -xp[m - 1 - i];
            w[i] = wp[m - 1 - i];
            x[m + i] = xp[i];
            w[m + i] = wp[i];
        }
    }
};

const RuleData rule_data;

thread_local d7::LocScaleSlot slots[4] = {
    {-1, 0, {0.0, 0.0}, {0.0, 0.0}}, {-1, 0, {0.0, 0.0}, {0.0, 0.0}},
    {-1, 0, {0.0, 0.0}, {0.0, 0.0}}, {-1, 0, {0.0, 0.0}, {0.0, 0.0}}};

}  // namespace

namespace d7 {

LocScaleRule loc_scale_rule() {
    LocScaleRule r = {(int) rule_data.x.size(), rule_data.x.data(),
                      rule_data.w.data()};
    return r;
}

LocScaleSlot& loc_scale_slot(int k) { return slots[k]; }

}  // namespace d7

// [[Rcpp::export]]
List loc_scale_rule_cpp() {
    return List::create(Named("x") = NumericVector(rule_data.x.begin(),
                                                   rule_data.x.end()),
                        Named("w") = NumericVector(rule_data.w.begin(),
                                                   rule_data.w.end()));
}

// The standardized diagonal for each row of U, the distinct shapes: column k
// of `info` is I_kk and of `dinfo` D_kk (left zero at order 0). `fam` is the
// class name; the families are those loc_scale_compiled() names.
// [[Rcpp::export]]
List loc_scale_diag_cpp(std::string fam, NumericMatrix U, int order,
                        int threads = 1) {
    int tag = -1, p = 0;
    if (fam == "SkewNormal1Distrib") { tag = 0; p = 3; }
    else if (fam == "SkewTDistrib") { tag = 1; p = 4; d7::skewt_pt(); }
    else if (fam == "PseudoHuberDistrib") { tag = 2; p = 3; d7::bessel_k_fn(); }
    else if (fam == "PseudoHuber2Distrib") { tag = 3; p = 3; d7::bessel_k_fn(); }
    else if (fam == "SkewNormal2Distrib") { tag = 4; p = 3; }
    else Rcpp::stop("no compiled quadrature for this family");
    const int m = U.nrow(), ns = U.ncol();
    NumericMatrix info(m, p), dinfo(m, p);
    const double* up = U.begin();
    double *ip = info.begin(), *dp = dinfo.begin();
    const bool want_d = order >= 1;
    d7::par_for((std::size_t) m * p, threads, d7::kMinCostly,
                [&](std::size_t task) {
        const int r = (int) (task % m), k = (int) (task / m);
        double shape[2] = {0.0, 0.0};
        for (int j = 0; j < ns && j < 2; ++j) shape[j] = up[r + (std::size_t) j * m];
        double out[2] = {0.0, 0.0};
        switch (tag) {
        case 0: d7::skewnormal1_quad_diag(k, shape, out, want_d); break;
        case 1: d7::skewt_quad_diag(k, shape, out, want_d); break;
        case 2: d7::pseudohuber_quad_diag(k, shape, out, want_d); break;
        case 3: d7::pseudohuber2_quad_diag(k, shape, out, want_d); break;
        case 4: d7::skewnormal2_quad_diag(k, shape, out, want_d); break;
        }
        ip[r + (std::size_t) k * m] = out[0];
        dp[r + (std::size_t) k * m] = want_d ? out[1] : 0.0;
    });
    return List::create(Named("info") = info, Named("dinfo") = dinfo);
}
