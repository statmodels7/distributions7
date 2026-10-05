#ifndef D7_PT_LOC_SCALE_H
#define D7_PT_LOC_SCALE_H

#include <Rcpp.h>
#include <cmath>

// The diagonal of the expected information of a location-scale family, and
// the derivative of each diagonal entry in its own parameter, by the exp-sinh
// rule of loc_scale_rule() (expected_loc_scale.R). At location 0 and scale 1,
//   I_kk = sum_j w_j f(x_j) l_kk(x_j),
//   D_kk = sum_j w_j f(x_j) (l_kkk(x_j) + l_kk(x_j) l_k(x_j)),
// and at (mu, sigma) each is divided by sigma^m, m the number of its indices
// on mu or sigma: 2 and 3 for k a location or a scale, 0 for a shape. The
// sums are accumulated in long double in node order, as R's colSums() does,
// a node with f = 0 contributes nothing, and a non-finite sum is NA.
//
// Each family supplies <fam>_quad_diag(k, shape, out, want_d), which loops
// over the nodes calling its own component functions. The vector path
// (loc_scale_diag_cpp) and the scalar registry both call it, the registry
// through loc_scale_cached(), which keeps the last standardized pair per
// parameter index on each thread: one quadrature per distinct shape.

namespace d7 {

struct LocScaleRule {
    int n;
    const double* x;
    const double* w;
};

// the nodes and weights, built when the library is loaded
LocScaleRule loc_scale_rule();

struct LocScaleSum {
    long double s0 = 0.0L, s1 = 0.0L;
    void add(double fw, double g, double H, double T3) {
        if (fw == 0.0) return;
        s0 += H * fw;
        s1 += (T3 + H * g) * fw;
    }
    void add0(double fw, double H) {
        if (fw == 0.0) return;
        s0 += H * fw;
    }
    static double finish(long double s) {
        double r = (double) s;
        return std::isfinite(r) ? r : NA_REAL;
    }
    void result(double* out) const {
        out[0] = finish(s0);
        out[1] = finish(s1);
    }
};

// the standardized pair to the pair at scale sigma, k the parameter index;
// R_pow() is what R's sigma^3 evaluates (powl() on 64-bit Windows)
inline void loc_scale_unstandardize(int k, double sigma, const double* std,
                                    double* out) {
    if (k < 2) {
        out[0] = std[0] / (sigma * sigma);
        out[1] = std[1] / R_pow(sigma, 3.0);
    } else {
        out[0] = std[0];
        out[1] = std[1];
    }
}

struct LocScaleSlot {
    int tag;
    int ns;
    double shape[2];
    double val[2];
};

// the slot of parameter index k (k < 4) on the calling thread
LocScaleSlot& loc_scale_slot(int k);

// the standardized pair of family `tag` at `shape`, from the thread's cache
// when the last call for this k was the same family at the same shape
template <class F>
inline void loc_scale_cached(int tag, int k, const double* shape, int ns,
                             F compute, double* std) {
    LocScaleSlot& s = loc_scale_slot(k);
    bool hit = s.tag == tag && s.ns == ns;
    for (int j = 0; hit && j < ns; ++j) hit = s.shape[j] == shape[j];
    if (!hit) {
        compute(s.val);
        s.tag = tag;
        s.ns = ns;
        for (int j = 0; j < ns; ++j) s.shape[j] = shape[j];
    }
    std[0] = s.val[0];
    std[1] = s.val[1];
}

} // namespace d7

#endif
