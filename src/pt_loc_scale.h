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

// THE FULL MATRICES, for the joint registry entries. At location 0 and scale
// 1, with P parameters (P <= 4),
//   I_ab  = sum_j w_j f(x_j) l_ab(x_j),
//   D_abc = sum_j w_j f(x_j) (l_abc(x_j) + l_ab(x_j) l_c(x_j)),
// P x P and P x P x P column-major (I[a + P*b], D[a + P*b + P*P*c]), both
// triangles filled. These are the sums loc_scale_expected() forms for the
// entries off the diagonal (integrate_rows() over T3 + Hab * g), in long
// double and node order, skipping a node with f = 0; a non-finite sum is NA.
// Each family supplies <fam>_quad_full(shape, I, D), which loops over the
// nodes, fills g (P), H (P x P) and T3 (P x P x P) at each and calls add().
struct LocScaleFull {
    int P;
    long double s0[16], s1[64];
    explicit LocScaleFull(int p) : P(p) {
        for (int i = 0; i < 16; ++i) s0[i] = 0.0L;
        for (int i = 0; i < 64; ++i) s1[i] = 0.0L;
    }
    void add(double fw, const double* g, const double* H, const double* T3) {
        if (fw == 0.0) return;
        const int PP = P * P;
        for (int b = 0; b < P; ++b)
            for (int a = 0; a < P; ++a) {
                const double Hab = H[a + P * b];
                s0[a + P * b] += Hab * fw;
                for (int c = 0; c < P; ++c)
                    s1[a + P * b + PP * c] +=
                        (T3[a + P * b + PP * c] + Hab * g[c]) * fw;
            }
    }
    void result(double* I, double* D) const {
        const int PP = P * P;
        for (int i = 0; i < PP; ++i) I[i] = LocScaleSum::finish(s0[i]);
        for (int i = 0; i < PP * P; ++i) D[i] = LocScaleSum::finish(s1[i]);
    }
};

struct LocScaleFullSlot {
    int tag;
    int ns;
    double shape[2];
    double I[16];
    double D[64];
};

// the full slot of the calling thread
LocScaleFullSlot& loc_scale_full_slot();

// the standardized full pair of family `tag` at `shape`, from the thread's
// cache when the last call was the same family at the same shape;
// compute(I, D) fills both
template <class F>
inline const LocScaleFullSlot& loc_scale_full_cached(int tag,
                                                     const double* shape,
                                                     int ns, F compute) {
    LocScaleFullSlot& s = loc_scale_full_slot();
    bool hit = s.tag == tag && s.ns == ns;
    for (int j = 0; hit && j < ns; ++j) hit = s.shape[j] == shape[j];
    if (!hit) {
        compute(s.I, s.D);
        s.tag = tag;
        s.ns = ns;
        for (int j = 0; j < ns; ++j) s.shape[j] = shape[j];
    }
    return s;
}

// sigma^m, m the number of indices on the location or the scale (index 0
// or 1), as R's sig^m evaluates it
inline double loc_scale_sigma_pow(int m, double sigma) {
    if (m == 0) return 1.0;
    if (m == 1) return sigma;
    if (m == 2) return sigma * sigma;
    return R_pow(sigma, (double) m);
}

// the standardized full matrices to the expected Hessian at scale sigma
// (P x P) and its derivative (P x P x P), the joint entries' layout
inline void loc_scale_full_info(const LocScaleFullSlot& s, int P,
                                double sigma, double* out) {
    for (int b = 0; b < P; ++b)
        for (int a = 0; a < P; ++a) {
            const int m = (a < 2) + (b < 2);
            out[a + P * b] = s.I[a + P * b] / loc_scale_sigma_pow(m, sigma);
        }
}

inline void loc_scale_full_dinfo(const LocScaleFullSlot& s, int P,
                                 double sigma, double* out) {
    const int PP = P * P;
    for (int c = 0; c < P; ++c)
        for (int b = 0; b < P; ++b)
            for (int a = 0; a < P; ++a) {
                const int m = (a < 2) + (b < 2) + (c < 2);
                out[a + P * b + PP * c] =
                    s.D[a + P * b + PP * c] / loc_scale_sigma_pow(m, sigma);
            }
}

} // namespace d7

#endif
