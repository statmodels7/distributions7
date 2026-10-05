#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <cstring>
#include <cstdlib>
#include <string>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <vector>
#include "pt_gaussian1.h"
#include "pt_gamma1.h"
#include "pt_poisson.h"
#include "pt_negbin2.h"
#include "pt_beta1.h"
#include "pt_bernoulli.h"
#include "pt_binomial.h"
#include "pt_exponential.h"
#include "pt_geometric.h"
#include "pt_chisq.h"
#include "pt_cauchy.h"
#include "pt_logistic.h"
#include "pt_gaussian2.h"
#include "pt_gaussian3.h"
#include "pt_lognormal1.h"
#include "pt_invgauss1.h"
#include "pt_invgauss2.h"
#include "pt_gamma2.h"
#include "pt_gpd.h"
#include "pt_vonmises.h"
#include "pt_weibull3.h"
#include "pt_lognormal2.h"
#include "pt_student_t1.h"
#include "pt_student_t2.h"
#include "pt_gengamma1.h"
#include "pt_gengamma2.h"
#include "pt_gumbel.h"
#include "pt_laplace.h"
#include "pt_laplace2.h"
#include "pt_weibull1.h"
#include "pt_beta2.h"
#include "pt_enet.h"
#include "pt_negbin1.h"
#include "pt_skewnormal1.h"
#include "pt_skewt.h"
#include "pt_pseudohuber.h"
#include "pt_pseudohuber2.h"
#include "pt_skewnormal2.h"
#include "pt_pig.h"
#include "pt_betabinom.h"
#include "pt_logpdf.h"
#include "pt_wrappers.h"
#include "pt_transform.h"
#include "pt_trunc_rule.h"
#include "pt_cdf.h"
#include "pt_loc_scale.h"
#include "pt_sqrt.h"

// The scalar C entry points of the fast route piano_parallel.txt section 2a
// describes: the score and the second derivative of the log-density in ONE
// parameter at ONE observation, on the parameter scale, which is what a
// score-driven filter reads at every step of its recursion, and the
// expected second derivative with its derivative in the same parameter,
// which a filter driven by a scaled score reads. A consumer resolves them
// once with R_GetCCallable and its loop then calls plain function pointers,
// touching no R API -- the precondition for any thread near that loop.
//
// The families covered are keyed by their S7 CLASS name, which is
// unambiguous where a distrib_name is shared across parametrizations, and
// an unknown name answers -1: the consumer keeps its R callbacks, so
// coverage is a speed property and never a correctness one.
//
// This file holds no formula. Each family's header (pt_<family>.h) carries
// one function per component and order, which the family's vector kernels
// call as well, and two routers by parameter index that compute only the
// special functions that index needs; a case below calls a router. The
// fast route is therefore bit-identical to the callback route by
// construction, and the twin tests in test-ccallable.R assert it.

namespace {

// the id of a family is its position here
const char* const d7_scalar_classes[] = {
    "Gaussian1Distrib",     // 0
    "Gamma1Distrib",        // 1
    "PoissonDistrib",       // 2
    "NegBin2Distrib",       // 3
    "Beta1Distrib",         // 4
    "BernoulliDistrib",     // 5
    "BinomialDistrib",      // 6
    "ExponentialDistrib",   // 7
    "GeometricDistrib",     // 8
    "ChisqDistrib",         // 9
    "CauchyDistrib",        // 10
    "LogisticDistrib",      // 11
    "Gaussian2Distrib",     // 12
    "Gaussian3Distrib",     // 13
    "Lognormal1Distrib",    // 14
    "InvGauss1Distrib",     // 15
    "InvGauss2Distrib",     // 16
    "Gamma2Distrib",        // 17
    "GPDDistrib",           // 18
    "VonMises1Distrib",     // 19
    "VonMises2Distrib",     // 20
    "Weibull3Distrib",      // 21
    "Lognormal2Distrib",    // 22
    "StudentT1Distrib",     // 23
    "StudentT2Distrib",     // 24
    "GenGamma1Distrib",     // 25
    "GenGamma2Distrib",     // 26
    "GumbelDistrib",        // 27
    "LaplaceDistrib",       // 28
    "Laplace2Distrib",      // 29
    "Weibull1Distrib",      // 30
    "Beta2Distrib",         // 31
    "EnetDistrib",          // 32
    "NegBin1Distrib",       // 33
    "SkewNormal1Distrib",   // 34
    "SkewTDistrib",         // 35
    "PseudoHuberDistrib",   // 36
    "PseudoHuber2Distrib",  // 37
    "SkewNormal2Distrib",   // 38
    "Pig1Distrib",          // 39
    "Pig2Distrib",          // 40
    "BetaBinom1Distrib",    // 41
    "BetaBinom2Distrib"     // 42
};

const int d7_n_scalar_classes =
    sizeof(d7_scalar_classes) / sizeof(d7_scalar_classes[0]);

// the number of parameters and of constants of each family, by id: what a
// wrapper reads to rebuild the inner family's vector
const int d7_n_params[] = {2, 2, 1, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 2, 3, 2, 3, 4, 3, 3, 3, 2, 2, 2, 2};
const int d7_n_constants[] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1};

} // namespace

namespace d7 {

// Resolved on first use, on the calling thread. Not a guarded static, for
// the reason vm_bessel() gives in vonmises.cpp.
static N7LogBesselI n7lbi_ptr = nullptr;

N7LogBesselI vm_log_i0_fn() {
    if (n7lbi_ptr == nullptr)
        n7lbi_ptr = (N7LogBesselI) R_GetCCallable("numericals7",
                                                  "n7_log_bessel_i");
    return n7lbi_ptr;
}

}  // namespace d7

namespace {

// ---- the wrappers ------------------------------------------------------------
//
// A wrapped family is named "<WrapperClass>[:aux]|<inner name>" by
// distrib_scalar_route(), the inner name being a family's class or, for a
// wrapper of a wrapper, another such name. aux is what the wrapper needs
// besides its constants (for fixed() the mask of the fixed parameters, for
// transformation() the code of the transformer in pt_transform.h). The
// first d7_scalar_id() call on a name enters a node {wrapper, aux, inner id}
// in a fixed table and the id is kWrapBase plus the node's position; a
// later call on the same name returns the same id. th is the wrapper's
// parameters, then the inner family's constants, then the wrapper's own; k
// counts the wrapper's parameters from zero. Each entry rebuilds the inner
// family's vector and reads the inner family's entries, which for a wrapped
// inner family are a node's entries in turn; the composition formulas are
// in pt_wrappers.h.

const int kWrapBase = 100000;

struct WrapName { const char* name; int code; };
const WrapName d7_wrappers[] = {
    {"FixedContinuousDistrib", 1},
    {"FixedDiscreteDistrib", 1},
    {"ZeroInflatedDistrib", 2},
    {"ZeroAdjustedDiscreteDistrib", 3},
    {"ZeroAdjustedContinuousDistrib", 4},
    {"FoldedDistrib", 5},
    {"TransformedDistrib", 6},
    {"TruncatedDiscreteDistrib", 7},
    {"TruncatedContinuousDistrib", 8}
};
const int d7_n_wrappers = sizeof(d7_wrappers) / sizeof(d7_wrappers[0]);

// The node table. It is written only by d7_scalar_id(), which a consumer
// calls on its own thread before any loop, and a slot is complete before
// its id is returned; the entries only read it. Nodes are never removed, so
// the table bounds the number of distinct wrapped names in a session.
struct WrapNode { int code, aux, inner, np, nc; };
const int kMaxNodes = 1024;
WrapNode d7_nodes[kMaxNodes];
int d7_n_nodes = 0;

// the length of the vector one observation's entries read, bounded by the
// rebuilding buffers below
const int kMaxVec = 16;

bool node_valid(int id) {
    return id >= kWrapBase && id - kWrapBase < d7_n_nodes;
}
const WrapNode& node_of(int id) { return d7_nodes[id - kWrapBase]; }

// the number of parameters and of constants of a family or a node
int id_np(int id) { return id >= kWrapBase ? node_of(id).np : d7_n_params[id]; }
int id_nc(int id) { return id >= kWrapBase ? node_of(id).nc : d7_n_constants[id]; }

// the family at the bottom of a chain of wrappers
int base_of(int id) {
    while (id >= kWrapBase) id = node_of(id).inner;
    return id;
}

int popcount(int m) {
    int c = 0;
    for (; m; m >>= 1) c += m & 1;
    return c;
}

}  // namespace

extern "C" {
int d7_scalar_id(const char* cls);
void d7_score_curv(int id, int k, double y, const double* th, double* out);
void d7_info_dinfo(int id, int k, double y, const double* th, double* out);
double d7_logpdf(int id, double y, const double* th);
}

namespace {

// the families with a compiled distribution function (pt_cdf.h)
bool quad_cdf_family(int id);

bool cont_has_cdf(int id) {
    if (quad_cdf_family(id)) return true;
    switch (id) {
    case 0: case 10: case 11: case 27: case 28: case 12: case 13: case 29: case 14: case 7: case 30: case 18: case 21: case 22: case 15: case 16: case 32: return true;
    default: return false;
    }
}

// whether a continuous truncation can be routed over `inner`: a family with
// a compiled distribution function and at most four parameters, or a chain
// of fixed(), folded() and transformation() nodes over one
bool cont_trunc_ok(int inner) {
    while (inner >= kWrapBase) {
        const int code = node_of(inner).code;
        if (code != 1 && code != 5 && code != 6) return false;
        inner = node_of(inner).inner;
    }
    return d7_n_params[inner] <= 4 && cont_has_cdf(inner);
}

int wrap_id(const char* cls, const char* bar) {
    std::string head(cls, bar - cls);
    int aux = 0;
    std::size_t colon = head.find(':');
    if (colon != std::string::npos) {
        aux = std::atoi(head.c_str() + colon + 1);
        head = head.substr(0, colon);
    }
    int code = -1;
    for (int i = 0; i < d7_n_wrappers; ++i)
        if (head == d7_wrappers[i].name) code = d7_wrappers[i].code;
    if (code < 0) return -1;
    const int inner = d7_scalar_id(bar + 1);
    if (inner < 0 || aux < 0) return -1;
    const int P = id_np(inner), nc = id_nc(inner);
    int np = P, ncw = nc;
    if (code == 1) {
        if (aux == 0 || P > 30 || aux >= (1 << P)) return -1;
        np = P - popcount(aux);
        ncw = nc + popcount(aux);
    } else if (code == 6) {
        if (aux < 1 || aux > d7::kTransformN) return -1;
        ncw = nc + d7::transform_n_par(aux);
    } else if (aux != 0) {
        return -1;
    } else if (code >= 2 && code <= 4) {
        np = P + 1;
    } else if (code == 7 || code == 8) {
        ncw = nc + 2;
    }
    if (P + nc > kMaxVec || np + ncw > kMaxVec) return -1;
    if (code == 8 && !cont_trunc_ok(inner)) return -1;
    for (int i = 0; i < d7_n_nodes; ++i) {
        const WrapNode& w = d7_nodes[i];
        if (w.code == code && w.aux == aux && w.inner == inner) return kWrapBase + i;
    }
    if (d7_n_nodes >= kMaxNodes) return -1;
    d7_nodes[d7_n_nodes] = WrapNode{code, aux, inner, np, ncw};
    return kWrapBase + d7_n_nodes++;
}

// fixed(): the inner vector from the free values, the inner constants and
// the fixed values, and the inner index of free parameter k
int fixed_inner(int inner, int mask, int k, const double* th, double* full) {
    const int P = id_np(inner), nc = id_nc(inner);
    int nf = 0;
    for (int j = 0; j < P; ++j) nf += !((mask >> j) & 1);
    int fi = 0, xi = 0, jk = -1;
    for (int j = 0; j < P; ++j) {
        if ((mask >> j) & 1) {
            full[j] = th[nf + nc + xi++];
        } else {
            if (fi == k) jk = j;
            full[j] = th[fi++];
        }
    }
    for (int c = 0; c < nc; ++c) full[P + c] = th[nf + c];
    return jk;
}

// the zero wrappers: the parent's vector (its parameters, then its
// constants) and the wrapper's probability, which follows the parameters
// a center and a scale of the parent, from its parameters, where the
// parent is a location-scale family on the real line: what folded()'s
// quadrature centers and scales its rule by; false elsewhere. Through
// fixed() the parent's vector is rebuilt first.
bool center_scale(int inner, const double* th, double* c, double* s) {
    if (inner >= kWrapBase) {
        const WrapNode& w = node_of(inner);
        if (w.code != 1) return false;
        double full[kMaxVec];
        fixed_inner(w.inner, w.aux, 0, th, full);
        return center_scale(w.inner, full, c, s);
    }
    switch (inner) {
    case 0: case 10: case 11: case 23: case 24: case 27: case 28:
    case 34: case 35: case 36: case 37: case 38:
        *c = th[0]; *s = th[1]; return true;
    case 12: *c = th[0]; *s = d7::sqrt_cr(th[1]); return true;
    case 13: *c = th[0]; *s = 1.0 / d7::sqrt_cr(th[1]); return true;
    case 29: *c = th[0]; *s = 1.0 / th[1]; return true;
    case 32: {
        double a = th[1] * th[2], cc = th[1] * (1.0 - th[2]);
        *c = th[0]; *s = 1.0 / (a + d7::sqrt_cr(cc)); return true;
    }
    default: return false;
    }
}

// folded(): the parent's (k, k) pieces at y and -y, and the folded score
// and second derivative from them
void fold_score_curv(int inner, int k, double y, const double* full,
                     double* out) {
    double ghp[2], ghm[2];
    d7_score_curv(inner, k, y, full, ghp);
    d7_score_curv(inner, k, -y, full, ghm);
    const d7::FoldW W = d7::fold_w(d7_logpdf(inner, y, full),
                                   d7_logpdf(inner, -y, full));
    out[0] = d7::fold_score(W.w, ghp[0], ghm[0]);
    out[1] = d7::fold_curv(W.w, ghp[0], ghm[0], ghp[1], ghm[1]);
}

// The nodes and weights of folded()'s quadrature over y in (0, Inf), for a
// parent centered at c with scale s: a tanh-sinh rule on [0, |c|], whose
// nodes cluster at both ends -- at zero, where the weight f(y)/(f(y) +
// f(-y)) turns from 1/2 to 1 over a width of order s^2/|c|, and at |c|,
// where the mass is and where a parent with a kink at its center puts it --
// and the exp-sinh rule of pt_loc_scale.h scaled by s on [|c|, Inf). Both
// have step 1/32; the tanh-sinh nodes run until their distance to an end
// falls below 1e-17 of the interval.
void fold_rule(double c, double s, std::vector<double>& y,
               std::vector<double>& w) {
    y.clear();
    w.clear();
    const double a = std::fabs(c);
    const double h = 1.0 / 32.0, hp = M_PI / 2.0;
    if (a > 0) {
        const int K = (int) std::floor(3.2 / h);
        for (int i = -K; i <= K; ++i) {
            const double t = i * h;
            const double u = hp * std::sinh(t);
            const double e = std::exp(-2.0 * std::fabs(u));
            const double lo = a * e / (1.0 + e);       // distance to the near end
            const double yy = (u < 0) ? lo : a - lo;
            const double ww = h * hp * std::cosh(t) * a * 2.0 * e /
                ((1.0 + e) * (1.0 + e));
            if (!(yy > 0) || !(yy < a)) continue;
            y.push_back(yy);
            w.push_back(ww);
        }
    }
    const d7::LocScaleRule R = d7::loc_scale_rule();
    for (int j = R.n / 2; j < R.n; ++j) {
        y.push_back(a + s * R.x[j]);
        w.push_back(s * R.w[j]);
    }
}

// whether parameter k of a parent is the center of a family with a kink
// there (laplace, laplace2, enet), through fixed()
bool kink_center(int inner, int k) {
    while (inner >= kWrapBase) {
        const WrapNode& w = node_of(inner);
        if (w.code != 1) return false;
        const int P = id_np(w.inner);
        int fi = 0, j = -1;
        for (int i = 0; i < P && j < 0; ++i)
            if (!((w.aux >> i) & 1) && fi++ == k) j = i;
        if (j < 0) return false;
        k = j;
        inner = w.inner;
    }
    return k == 0 && (inner == 28 || inner == 29 || inner == 32);
}

// folded()'s (k, k) expected second derivative and its derivative in k, as
// integrals over y of the folded density by fold_rule(); by the second
// Bartlett identity E[l_kk] = -E[l_k^2] and d_k E[l_kk] =
// -E[2 l_kk l_k + l_k^3]. The sums are accumulated in long double in node
// order, as folded_expected() in R accumulates them.
void fold_info_dinfo(int inner, int k, const double* full, double* out) {
    double c, s;
    if (!center_scale(inner, full, &c, &s)) {
        out[0] = out[1] = R_NaN;
        return;
    }
    std::vector<double> ys, ws;
    fold_rule(c, s, ys, ws);
    long double a0 = 0.0L, a1 = 0.0L;
    for (std::size_t j = 0; j < ys.size(); ++j) {
        const double y = ys[j];
        const d7::FoldW Wt = d7::fold_w(d7_logpdf(inner, y, full),
                                        d7_logpdf(inner, -y, full));
        const double fw = Wt.L * ws[j];
        if (fw == 0.0) continue;
        double ghp[2], ghm[2];
        d7_score_curv(inner, k, y, full, ghp);
        d7_score_curv(inner, k, -y, full, ghm);
        const double g = d7::fold_score(Wt.w, ghp[0], ghm[0]);
        const double h = d7::fold_curv(Wt.w, ghp[0], ghm[0], ghp[1], ghm[1]);
        a0 += (g * g) * fw;
        a1 += (h * g + g * h + g * g * g) * fw;
    }
    const double r0 = (double) a0;
    double r1 = (double) a1;
    // a parent with a kink at its center (laplace, laplace2, enet): the
    // folded score in the center jumps at y* = |c|, a point that moves with
    // c, and d_c E[l_cc] carries the boundary term L(y*) [G(y*-) - G(y*+)]
    // sign(c), G = l_c^2, the one-sided scores read four units in the last
    // place either side of y*
    if (kink_center(inner, k) && c != 0) {
        const double a = std::fabs(c);
        const d7::FoldW W0 = d7::fold_w(d7_logpdf(inner, a, full),
                                        d7_logpdf(inner, -a, full));
        double lo[2], hi[2];
        fold_score_curv(inner, k, a * (1 + (-4) * DBL_EPSILON), full, lo);
        fold_score_curv(inner, k, a * (1 + 4 * DBL_EPSILON), full, hi);
        r1 = r1 + W0.L * (lo[0] * lo[0] - hi[0] * hi[0]) * (c > 0 ? 1.0 : -1.0);
    }
    out[0] = R_FINITE(r0) ? -r0 : NA_REAL;
    out[1] = R_FINITE(r1) ? -r1 : NA_REAL;
}

// transformation(): the parent's vector is the wrapper's own (the same
// parameters and the parent's constants), and the transformer's parameters
// follow it
const double* transform_par(int inner, const double* th) {
    return th + id_np(inner) + id_nc(inner);
}

double zero_inner(int inner, const double* th, double* full) {
    const int P = id_np(inner), nc = id_nc(inner);
    for (int j = 0; j < P; ++j) full[j] = th[j];
    for (int c = 0; c < nc; ++c) full[P + c] = th[P + 1 + c];
    return th[P];
}

// ---- truncated(), discrete parents ------------------------------------------
//
// th is the parent's vector followed by the truncation points lower and
// upper, both included in the support. The retained mass and its
// derivatives, and the truncated expectations the information reads, are
// finite sums of the parent's own entries over support points:
//   Z = sum f,  Z_k = sum f s_k,  Z_kk = sum f (l_kk + s_k^2),
//   S = sum f s_k^2,  D = sum f (s_k^3 + l_kk s_k + s_k l_kk),
// with f = exp(log f). When the retained set is finite the sums run over it.
// When it is not (no upper point on an unbounded support), Z, Z_k and Z_kk
// for the density, the score and the curvature run over the points removed
// below `lower`, as the complements Z = 1 - sum f, Z_k = -sum f s_k and
// Z_kk = -sum f (...); a complement loses digits when the removed mass is
// most of the total, so above a removed mass of one half, and always for
// the information, the sums run over the retained points as a series,
// stopped beyond the mode once a term weighted by (1 + y)^6 falls below
// 1e-20 of the sum. (The information's complement would read the parent's
// E[s_k^2] and its derivative, whose error the subtraction amplifies.) The
// sums accumulate in long double in increasing y, as trunc_route_parts() in
// R accumulates them over the same points, which trunc_rule_cpp() hands it.

const double kTruncSeriesTol = 1e-20;
const double kTruncSeriesMax = 1e7;

// the largest support point of a discrete family or chain: the size of the
// binomial and the beta-binomials, one for the Bernoulli, infinite
// otherwise
double support_max(int id, const double* th) {
    if (id >= kWrapBase) {
        const WrapNode& w = node_of(id);
        double full[kMaxVec];
        if (w.code == 1) {
            fixed_inner(w.inner, w.aux, 0, th, full);
            return support_max(w.inner, full);
        }
        if (w.code == 2 || w.code == 3) {
            zero_inner(w.inner, th, full);
            return support_max(w.inner, full);
        }
        if (w.code == 7) {
            const double up = th[id_np(w.inner) + id_nc(w.inner) + 1];
            return std::min(up, support_max(w.inner, th));
        }
        return R_PosInf;
    }
    switch (id) {
    case 5: return 1.0;
    case 6: return th[1];
    case 41: case 42: return th[2];
    default: return R_PosInf;
    }
}

struct TruncRule { int branch; double y0, y1; };

// which points a truncated family's sums run over: branch 0 the finite
// retained set [y0, y1], branch 1 the removed points [y0, y1] below lower,
// branch 2 the retained points from y0 to the series' stop y1; y1 < y0 for
// an empty set and NaN when the series does not stop; `info` excludes
// branch 1
TruncRule trunc_disc_rule(int inner, const double* thi, double lo, double up,
                          bool info) {
    TruncRule r;
    const double a = (lo > 0) ? std::ceil(lo) : 0.0;
    const double b = std::min(std::floor(up), support_max(inner, thi));
    if (R_FINITE(b)) {
        r.branch = 0; r.y0 = a; r.y1 = b;
        return r;
    }
    if (!info) {
        long double removed = 0.0L;
        for (double y = 0; y < a; ++y)
            removed += std::exp(d7_logpdf(inner, y, thi));
        if ((double) removed <= 0.5) {
            r.branch = 1; r.y0 = 0; r.y1 = a - 1;
            return r;
        }
    }
    r.branch = 2; r.y0 = a;
    long double z = 0.0L;
    double fprev = -1.0;
    for (double y = a;; ++y) {
        const double f = std::exp(d7_logpdf(inner, y, thi));
        z += f;
        if (y > a && f <= fprev &&
            f * std::pow(1 + y, 6.0) <= kTruncSeriesTol * (double) z) {
            r.y1 = y;
            return r;
        }
        fprev = f;
        if (y - a > kTruncSeriesMax) {
            r.y1 = R_NaN;
            return r;
        }
    }
}

struct TruncSums { double z, zk, zkk, s, d; };

// the retained mass and the sums of parameter k, over the rule's points;
// `info` asks for S and D as well
TruncSums trunc_disc_sums(int inner, int k, const double* thi,
                          double lo, double up, bool info) {
    TruncSums t;
    const TruncRule r = trunc_disc_rule(inner, thi, lo, up, info);
    if (ISNAN(r.y1)) {
        t.z = t.zk = t.zkk = t.s = t.d = R_NaN;
        return t;
    }
    long double F = 0.0L, G = 0.0L, H = 0.0L, S = 0.0L, D = 0.0L;
    for (double x = r.y0; x <= r.y1; ++x) {
        const double f = std::exp(d7_logpdf(inner, x, thi));
        double sc[2];
        d7_score_curv(inner, k, x, thi, sc);
        const double g = sc[0], l = sc[1];
        F += f;
        G += f * g;
        H += f * (l + g * g);
        if (info) {
            S += f * (g * g);
            D += f * (g * g * g + l * g + g * l);
        }
    }
    if (r.branch == 1) {
        t.z = 1 - (double) F;
        t.zk = -(double) G;
        t.zkk = -(double) H;
        t.s = t.d = R_NaN;
    } else {
        t.z = (double) F;
        t.zk = (double) G;
        t.zkk = (double) H;
        t.s = (double) S;
        t.d = (double) D;
    }
    return t;
}

const double* trunc_points(int inner, const double* th) {
    return th + id_np(inner) + id_nc(inner);
}

void trunc_score_curv(int inner, int k, double y, const double* th, double* out) {
    const double* tp = trunc_points(inner, th);
    const TruncSums t = trunc_disc_sums(inner, k, th, tp[0], tp[1], false);
    double sc[2];
    d7_score_curv(inner, k, y, th, sc);
    const double m = t.zk / t.z, M = t.zkk / t.z;
    out[0] = sc[0] - m;
    out[1] = sc[1] - M + m * m;
}

void trunc_info_dinfo(int inner, int k, double y, const double* th, double* out) {
    const double* tp = trunc_points(inner, th);
    const TruncSums t = trunc_disc_sums(inner, k, th, tp[0], tp[1], true);
    const double m = t.zk / t.z, M = t.zkk / t.z;
    const double dm = M - m * m;
    out[0] = -(t.s / t.z - m * m);
    out[1] = -(t.d / t.z - t.s * t.zk / (t.z * t.z) - (dm * m + m * dm));
}

double trunc_logpdf(int inner, double y, const double* th) {
    const double* tp = trunc_points(inner, th);
    const TruncSums t = trunc_disc_sums(inner, 0, th, tp[0], tp[1], false);
    const double ld = d7_logpdf(inner, y, th) - std::log(t.z);
    return (y < tp[0] || y > tp[1]) ? R_NegInf : ld;
}

// ---- truncated(), continuous parents ----------------------------------------
//
// The retained mass is Z = F(b) - F(a), or S(a) - S(b) in the survival
// function when F(a) > 1/2, so that a lower point in the upper tail keeps
// its digits; an end at or beyond the parent's support contributes 0 (or 1).
// Its derivatives are those of F at the two ends, from the compiled
// distribution functions of pt_cdf.h. The information's sums
//   S = int f s_k^2,  D = int f (s_k^3 + l_kk s_k + s_k l_kk)
// over [a, b] are taken by the rule of pt_trunc_rule.h, at the parent's
// center and scale and with its kinks as cuts; a node where f w is zero
// contributes nothing. trunc_cont_ends_cpp() and trunc_cont_rule_cpp() hand
// R the same ends and the same nodes.

// A continuous family seen through a chain of fixed() nodes (cont_trunc_ok):
// the family at the bottom, its full vector, and the position in that vector
// of each free parameter. The distribution function, its support, center,
// scale and kinks are the bottom family's at the full vector, and a
// derivative in free parameter k is the bottom family's in parameter at[k].
struct ContView { int base, np; int at[4]; double full[kMaxVec]; };

void cont_view(int id, const double* th, ContView& v) {
    if (id < kWrapBase) {
        v.base = id;
        v.np = d7_n_params[id];
        const int n = v.np + d7_n_constants[id];
        for (int j = 0; j < n; ++j) v.full[j] = th[j];
        for (int k = 0; k < v.np; ++k) v.at[k] = k;
        return;
    }
    const WrapNode& w = node_of(id);
    double full[kMaxVec];
    int jk[4] = {0, 0, 0, 0};
    fixed_inner(w.inner, w.aux, 0, th, full);
    for (int k = 0; k < w.np && k < 4; ++k)
        jk[k] = fixed_inner(w.inner, w.aux, k, th, full);
    ContView in;
    cont_view(w.inner, full, in);
    v = in;
    v.np = w.np;
    for (int k = 0; k < w.np && k < 4; ++k) v.at[k] = in.at[jk[k]];
}

int hess_pos(int P, int i, int j);

// The two beta families read the response only through log y and
// log(1 - y). On the rule's logit variable v both are exact, -log1p(e^-v)
// and -log1p(e^v) (written so that the exponential never overflows), where
// y itself rounds to 1 once v passes about 37 and the mass beyond it, of
// order eps^b, was lost with the node: measured on a beta2 truncated below
// 0.3, the information was off by up to 1.5e-2 at b = 0.24 and by 0.48 at
// b = 0.1. The rule hands these families the two logarithms
// (mapped_rule()), and they are evaluated from them here.
void logit_logs(double v, double* ly, double* l1y) {
    if (v > 0) {
        const double t = std::log1p(std::exp(-v));
        *ly = -t;
        *l1y = -v - t;
    } else {
        const double t = std::log1p(std::exp(v));
        *ly = v - t;
        *l1y = -t;
    }
}

// one fixed() node: the inner family's vector, and the inner index of each
// free parameter; returns the inner id
int fixed_level(int id, const double* th, double* full, int* at) {
    const WrapNode& w = node_of(id);
    fixed_inner(w.inner, w.aux, 0, th, full);
    for (int k = 0; k < w.np && k < 4; ++k)
        at[k] = fixed_inner(w.inner, w.aux, k, th, full);
    return w.inner;
}

// whether a family, or a fixed() chain over one, is a beta (a folded or
// transformed beta is not: its response is not the beta variable)
bool logs_family(int id) {
    while (id >= kWrapBase && node_of(id).code == 1) id = node_of(id).inner;
    return id == 4 || id == 31;
}

// the log-density of a beta or a fixed() chain over one, its score in each
// free parameter (g) and its Hessian in hess_names() order (h), from log y
// and log(1 - y)
double beta_logs_parts(int id, double ly, double l1y, const double* th,
                       double* g, double* h) {
    ContView v;
    cont_view(id, th, v);
    const double* f = v.full;
    double a, b, G[2], H[3];
    if (v.base == 4) {
        const double m = f[0], p = f[1];
        a = m * p;
        b = (1.0 - m) * p;
        const double a0 = R::digamma(a), b0 = R::digamma(b);
        const double a1 = R::trigamma(a), b1 = R::trigamma(b);
        G[0] = d7::beta1_score_mu(p, ly, l1y, a0, b0);
        G[1] = d7::beta1_score_phi(m, ly, l1y, R::digamma(p), a0, b0);
        H[0] = d7::beta1_hess_mu_mu(p, a1, b1);
        H[1] = d7::beta1_hess_phi_phi(m, R::trigamma(p), a1, b1);
        H[2] = -a0 + b0 + ly - l1y + p * (-a1 * m + b1 * (1 - m));
    } else {
        a = f[0];
        b = f[1];
        const double s0 = R::digamma(a + b), s1 = R::trigamma(a + b);
        G[0] = d7::beta2_score_alpha(ly, R::digamma(a), s0);
        G[1] = d7::beta2_score_beta(l1y, R::digamma(b), s0);
        H[0] = d7::beta2_hess_alpha_alpha(R::trigamma(a), s1);
        H[1] = d7::beta2_hess_beta_beta(R::trigamma(b), s1);
        H[2] = s1;
    }
    if (g != nullptr)
        for (int k = 0; k < v.np; ++k) g[k] = G[v.at[k]];
    if (h != nullptr)
        for (int i = 0; i < v.np; ++i)
            for (int j = i; j < v.np; ++j)
                h[hess_pos(v.np, i, j)] = H[hess_pos(2, v.at[i], v.at[j])];
    return (a - 1.0) * ly + (b - 1.0) * l1y - R::lbeta(a, b);
}

// the support of a continuous family at its parameters
void cont_support(int id, const double* th, double* lo, double* hi) {
    if (id >= kWrapBase) {
        const WrapNode& w = node_of(id);
        if (w.code == 1) {
            double full[kMaxVec];
            int at[4];
            cont_support(fixed_level(id, th, full, at), full, lo, hi);
            return;
        }
        double a, b;
        cont_support(w.inner, th, &a, &b);
        if (w.code == 5) {
            *lo = 0.0;
            *hi = std::max(std::fabs(a), std::fabs(b));
            return;
        }
        // transformation(): the image of the parent's support
        const double* tp = transform_par(w.inner, th);
        const double ga = d7::transform_fwd(w.aux, a, tp);
        const double gb = d7::transform_fwd(w.aux, b, tp);
        *lo = std::min(ga, gb);
        *hi = std::max(ga, gb);
        return;
    }
    *lo = R_NegInf;
    *hi = R_PosInf;
    switch (id) {
    case 1: case 7: case 9: case 14: case 15: case 16: case 17: case 21:
    case 22: case 25: case 26: case 30:
        *lo = 0.0; return;
    case 18:
        *lo = 0.0;
        if (th[1] < 0) *hi = th[0] / -th[1];
        return;
    case 4: case 31:
        *lo = 0.0; *hi = 1.0; return;
    case 19: case 20:
        *lo = -M_PI; *hi = M_PI; return;
    default: return;
    }
}

// a center and a scale for the rule, for every continuous family
bool cont_center_scale(int id, const double* th, double* c, double* s) {
    if (id >= kWrapBase) {
        const WrapNode& w = node_of(id);
        if (w.code == 1) {
            double full[kMaxVec];
            int at[4];
            return cont_center_scale(fixed_level(id, th, full, at), full, c, s);
        }
        double c0, s0;
        if (!cont_center_scale(w.inner, th, &c0, &s0)) return false;
        if (w.code == 5) {
            *c = std::fabs(c0);
            *s = s0;
            return true;
        }
        // transformation(): the parent's center carried forward, and its
        // scale times |dy/dx| there
        const double* tp = transform_par(w.inner, th);
        *c = d7::transform_fwd(w.aux, c0, tp);
        *s = s0 * std::exp(-d7::transform_log_jac(w.aux, *c, tp));
        return R_FINITE(*c) && R_FINITE(*s) && *s > 0;
    }
    if (center_scale(id, th, c, s)) return true;
    switch (id) {
    case 1: *c = th[0]; *s = th[0] * d7::sqrt_cr(th[1]); return true;
    case 7: *c = th[0]; *s = th[0]; return true;
    case 9: *c = th[0]; *s = d7::sqrt_cr(2 * th[0]); return true;
    case 14: *c = std::exp(th[0]); *s = *c * d7::sqrt_cr(th[1]); return true;
    case 15: *c = th[0]; *s = d7::sqrt_cr(th[1] * th[0]) * th[0]; return true;
    case 16: *c = th[0]; *s = d7::sqrt_cr(th[0] / th[1]) * th[0]; return true;
    case 17: *c = th[0]; *s = d7::sqrt_cr(th[1]); return true;
    case 18: *c = 0.0; *s = th[0]; return true;
    case 19: *c = th[0]; *s = std::min(M_PI, 1 / d7::sqrt_cr(th[1])); return true;
    case 20: *c = th[0];
        *s = std::min(M_PI, d7::sqrt_cr(-2 * std::log(th[1]))); return true;
    case 21: case 26: *c = th[0]; *s = th[0]; return true;
    case 22: {
        const double r = th[1] / (th[0] * th[0]);
        *c = th[0] / d7::sqrt_cr(1 + r);
        *s = *c * d7::sqrt_cr(std::log1p(r));
        return true;
    }
    case 25: case 30: *c = th[0]; *s = th[0]; return true;
    case 4: {
        *c = th[0];
        *s = d7::sqrt_cr(th[0] * (1 - th[0]) / (1 + th[1]));
        return true;
    }
    case 31: {
        const double ab = th[0] + th[1];
        *c = th[0] / ab;
        *s = d7::sqrt_cr(th[0] * th[1] / (ab * ab * (ab + 1)));
        return true;
    }
    default: return false;
    }
}

// the kinks of the density in the response
int cont_kinks(int id, const double* th, double* kk) {
    if (id >= kWrapBase) {
        const WrapNode& w = node_of(id);
        if (w.code == 1) {
            double full[kMaxVec];
            int at[4];
            return cont_kinks(fixed_level(id, th, full, at), full, kk);
        }
        const int nk = cont_kinks(w.inner, th, kk);
        const double* tp = transform_par(w.inner, th);
        for (int i = 0; i < nk; ++i)
            kk[i] = (w.code == 5) ? std::fabs(kk[i]) : d7::transform_fwd(w.aux, kk[i], tp);
        return nk;
    }
    if (id == 28 || id == 29 || id == 32) { kk[0] = th[0]; return 1; }
    return 0;
}

void mapped_rule(int id, const double* th, double a, double b,
                 const double* kk, int nk, std::vector<double>& y,
                 std::vector<double>& w, std::vector<double>* ly = nullptr,
                 std::vector<double>* l1y = nullptr);

// ---- distribution functions by quadrature ----------------------------------
//
// For the families whose distribution function has no closed derivative in
// a shape parameter, the derivatives are integrals of the density's own:
//   F_k(q) = int_{lo}^{q} f s_k,  F_kl(q) = int_{lo}^{q} f (l_kl + s_k s_l),
// taken by the rule of pt_trunc_rule.h, or as minus the same integrals over
// [q, hi) when q lies above the family's center, so that neither side is a
// difference of two numbers near one. A location-scale family integrates
// only its shape parameters: with z = (q - mu)/sigma and l_y = -s_mu,
// F_mu = -f, F_sigma = -z f, F_mumu = -f s_mu, F_sigmasigma =
// f (-z^2 s_mu + 2z/sigma), F_musigma = f (-z s_mu + 1/sigma), and for a
// shape k, F_muk = -f s_k and F_sigmak = -z f s_k, all at q. The
// distribution function itself is R's where it has one (pt, pgamma, pchisq,
// pbeta) and the integral of f otherwise.

// the off-diagonal second derivative of the log-density in (k, l), k < l
double hess_off(int id, int k, int l, double y, const double* th) {
    switch (id) {
    case 1: {  // gamma1 (mu, phi)
        const double mu = th[0], ph = th[1];
        return -(y - mu) / (ph * ph * mu * mu);
    }
    case 17: {  // gamma2 (mu, v)
        const double mu = th[0], v = th[1], kk = mu * mu / v;
        return (-2 * mu * (std::log(mu / v) + std::log(y) - R::digamma(kk)) -
                3 * mu + y) / (v * v) + 2 * mu * mu * mu * R::trigamma(kk) / (v * v * v);
    }
    case 4: {  // beta1 (mu, phi)
        const double mu = th[0], ph = th[1], a = mu * ph, b = (1 - mu) * ph;
        return -R::digamma(a) + R::digamma(b) + std::log(y) - std::log1p(-y) +
            ph * (-R::trigamma(a) * mu + R::trigamma(b) * (1 - mu));
    }
    case 31:  // beta2 (alpha, beta)
        return R::trigamma(th[0] + th[1]);
    case 19:  // vonmises1 (mu, kappa)
        return std::sin(y - th[0]);
    case 20: {  // vonmises2 (mu, rho)
        double kk[2];
        d7::vm2_kappa(d7::vm_bessel(), th[1], 1, kk);
        return std::sin(y - th[0]) * kk[1];
    }
    case 25: case 26: {  // gengamma1 (a, d, p), gengamma2 (m, d, p)
        const double d = th[1], p = th[2], u = d / p, w = (d + 1) / p;
        double a = th[0];
        if (id == 26) a = std::exp(std::log(th[0]) + R::lgammafn(u) - R::lgammafn(w));
        const double lya = std::log(y / a), t = std::exp(p * lya);
        const double La = -d / a + p * t / a;
        const double Laa = (d - p * (p + 1) * t) / (a * a);
        const double Lad = -1 / a, Lap = t * (1 + p * lya) / a;
        const double Ldp = R::trigamma(u) * d / (p * p * p) + R::digamma(u) / (p * p);
        if (id == 25) {
            if (k == 0 && l == 1) return Lad;
            if (k == 0 && l == 2) return Lap;
            return Ldp;
        }
        // a = m e^g, g = lgamma(u) - lgamma(w): A = log a
        const double pu = R::digamma(u), pw = R::digamma(w);
        const double qu = R::trigamma(u), qw = R::trigamma(w);
        const double A1[3] = {1 / th[0], (pu - pw) / p,
                              (-d * pu + (d + 1) * pw) / (p * p)};
        double A2;  // A_kl
        if (k == 0) A2 = 0.0;
        else A2 = (-qu * d / (p * p) + qw * (d + 1) / (p * p)) / p - (pu - pw) / (p * p);
        const double ak = a * A1[k], al = a * A1[l];
        const double akl = a * (A2 + A1[k] * A1[l]);
        // direct mixed partials L_a(theta) for theta in (d, p)
        const double Lak = (k == 1) ? Lad : (k == 2 ? Lap : 0.0);
        const double Lal = (l == 1) ? Lad : (l == 2 ? Lap : 0.0);
        const double Lkl = (k == 1 && l == 2) ? Ldp : 0.0;
        return Laa * ak * al + La * akl + Lak * al + Lal * ak + Lkl;
    }
    case 35: {  // skewt (mu, sigma, alpha, nu): alpha-nu by the family's own
        // stencil in nu on the alpha score
        if (!(k == 2 && l == 3)) return R_NaN;
        const double h = d7::skewt_nu_step(th[3]);
        double acc = 0.0;
        for (int j = 0; j < 5; ++j) {
            if (d7::kSkewtW1[j] == 0.0) continue;
            const d7::SkewtPieces P =
                d7::skewt_pieces(y, th[0], th[1], th[2], th[3] + (j - 2) * h);
            acc = acc + d7::kSkewtW1[j] * d7::skewt_score_alpha(P);
        }
        return acc / h;
    }
    default: return R_NaN;
    }
}

// the integrated parameters of a family (bit j for parameter j) and whether
// it is location-scale in its first two
int quad_cdf_mask(int id, bool* ls) {
    *ls = false;
    switch (id) {
    case 23: case 24: case 34: case 36: case 37: case 38:
        *ls = true; return 4;
    case 35: *ls = true; return 12;
    case 9: return 1;
    case 1: case 17: case 4: case 31: case 19: case 20: return 3;
    case 25: case 26: return 7;
    default: return 0;
    }
}

bool quad_cdf_family(int id) { bool ls; return quad_cdf_mask(id, &ls) != 0; }

// R's closed distribution function where the family has one
bool closed_cdf(int id, double q, const double* th, bool lower, double* F) {
    switch (id) {
    case 23: *F = R::pt((q - th[0]) / th[1], th[2], lower, 0); return true;
    case 24: *F = R::pt((q - th[0]) / (th[1] * std::sqrt(1 - 2 / th[2])), th[2],
                        lower, 0); return true;
    case 1: *F = R::pgamma(q, 1 / th[1], th[1] * th[0], lower, 0); return true;
    case 17: *F = R::pgamma(q, th[0] * th[0] / th[1], th[1] / th[0], lower, 0);
        return true;
    case 9: *F = R::pchisq(q, th[0], lower, 0); return true;
    case 4: *F = R::pbeta(q, th[0] * th[1], (1 - th[0]) * th[1], lower, 0);
        return true;
    case 31: *F = R::pbeta(q, th[0], th[1], lower, 0); return true;
    case 25: case 26: {
        const double d = th[1], p = th[2];
        double a = th[0];
        if (id == 26) a = std::exp(std::log(th[0]) + R::lgammafn(d / p) -
                                   R::lgammafn((d + 1) / p));
        *F = R::pgamma(std::pow(std::max(q, 0.0) / a, p), d / p, 1.0, lower, 0);
        return true;
    }
    default: return false;
    }
}

// the integrals over the side of q: I0 = int f, Ik = int f s_k and
// Ikl = int f (l_kl + s_k s_l) over the parameters of `mask`, with
// `upper` true for [q, hi)
void quad_cdf_sums(int id, double q, const double* th, int mask, bool upper,
                   bool want_h, double* I0, double* Ik, double* Ikl) {
    const int P = d7_n_params[id];
    double slo, shi, c, s;
    cont_support(id, th, &slo, &shi);
    cont_center_scale(id, th, &c, &s);
    std::vector<double> ys, ws;
    const double a = upper ? std::max(q, slo) : slo;
    const double b = upper ? shi : std::min(q, shi);
    const bool logs = logs_family(id);
    std::vector<double> lys, l1ys;
    mapped_rule(id, th, a, b, nullptr, 0, ys, ws, logs ? &lys : nullptr,
                logs ? &l1ys : nullptr);
    long double S0 = 0.0L, S1[4] = {0, 0, 0, 0}, S2[16] = {0};
    for (std::size_t j = 0; j < ys.size(); ++j) {
        const double y = ys[j];
        if (logs) {
            double g[2], h[3];
            const double fw = std::exp(beta_logs_parts(id, lys[j], l1ys[j], th,
                                                       g, h)) * ws[j];
            if (fw == 0.0) continue;
            S0 += fw;
            for (int k = 0; k < 2; ++k) {
                if (!((mask >> k) & 1)) continue;
                S1[k] += fw * g[k];
            }
            if (!want_h) continue;
            for (int k = 0; k < 2; ++k) {
                if (!((mask >> k) & 1)) continue;
                S2[k * 4 + k] += fw * (h[k] + g[k] * g[k]);
                for (int m = k + 1; m < 2; ++m) {
                    if (!((mask >> m) & 1)) continue;
                    S2[k * 4 + m] += fw * (h[hess_pos(2, k, m)] + g[k] * g[m]);
                }
            }
            continue;
        }
        const double fw = std::exp(d7_logpdf(id, y, th)) * ws[j];
        if (fw == 0.0) continue;
        S0 += fw;
        double g[4], l[4];
        for (int k = 0; k < P; ++k) {
            if (!((mask >> k) & 1)) continue;
            double sc[2];
            d7_score_curv(id, k, y, th, sc);
            g[k] = sc[0]; l[k] = sc[1];
            S1[k] += fw * g[k];
        }
        if (!want_h) continue;
        for (int k = 0; k < P; ++k) {
            if (!((mask >> k) & 1)) continue;
            S2[k * 4 + k] += fw * (l[k] + g[k] * g[k]);
            for (int m = k + 1; m < P; ++m) {
                if (!((mask >> m) & 1)) continue;
                S2[k * 4 + m] += fw * (hess_off(id, k, m, y, th) + g[k] * g[m]);
            }
        }
    }
    *I0 = (double) S0;
    for (int k = 0; k < 4; ++k) Ik[k] = (double) S1[k];
    for (int k = 0; k < 16; ++k) Ikl[k] = (double) S2[k];
}

// the side of q to integrate over: above the center, the upper side,
// except on (0, 1) when the second shape is below one, where the mass
// within one unit in the last place of 1 is not representable in y (it is
// eps^b of the total) and the lower side is used instead
bool quad_side_upper(int id, double q, const double* th) {
    double c, s;
    cont_center_scale(id, th, &c, &s);
    if (id == 4 && (1 - th[0]) * th[1] < 1) return false;
    if (id == 31 && th[1] < 1) return false;
    return q > c;
}

double quad_cdf(int id, double q, const double* th, bool lower) {
    double F;
    if (closed_cdf(id, q, th, lower, &F)) return F;
    const bool up = quad_side_upper(id, q, th);
    double I0, Ik[4], Ikl[16];
    quad_cdf_sums(id, q, th, 0, up, false, &I0, Ik, Ikl);
    // I0 is F on the lower side and S on the upper
    if (up) return lower ? 1 - I0 : I0;
    return lower ? I0 : 1 - I0;
}

// the gradient (want_h false) or the Hessian, in hess_names() order
void quad_cdf_derivs(int id, double q, const double* th, bool want_h,
                     double* out) {
    const int P = d7_n_params[id];
    bool ls;
    const int mask = quad_cdf_mask(id, &ls);
    const bool up = quad_side_upper(id, q, th);
    double I0, Ik[4], Ikl[16];
    quad_cdf_sums(id, q, th, mask, up, want_h, &I0, Ik, Ikl);
    const double sgn = up ? -1.0 : 1.0;
    // the point quantities a location-scale family reads at q
    double f = 0.0, z = 0.0, sq[4] = {0, 0, 0, 0};
    if (ls) {
        f = std::exp(d7_logpdf(id, q, th));
        z = (q - th[0]) / th[1];
        for (int k = 0; k < P; ++k) {
            if (k != 0 && !((mask >> k) & 1)) continue;
            double sc[2];
            d7_score_curv(id, k, q, th, sc);
            sq[k] = sc[0];
        }
    }
    if (!want_h) {
        for (int k = 0; k < P; ++k) {
            if ((mask >> k) & 1) out[k] = sgn * Ik[k];
            else if (k == 0) out[k] = -f;
            else out[k] = -z * f;
        }
        return;
    }
    const double sg = ls ? th[1] : 1.0;
    for (int i = 0; i < P; ++i) {
        for (int j = i; j < P; ++j) {
            const bool mi = (mask >> i) & 1, mj = (mask >> j) & 1;
            double v;
            if (mi && mj) {
                v = sgn * Ikl[i * 4 + j];
            } else if (mi || mj) {
                // a location or scale with a shape: -f s_k or -z f s_k
                const int sh = mi ? i : j, lsj = mi ? j : i;
                v = (lsj == 0) ? -f * sq[sh] : -z * f * sq[sh];
            } else if (i == 0 && j == 0) {
                v = -f * sq[0];
            } else if (i == 1 && j == 1) {
                v = f * (-z * z * sq[0] + 2 * z / sg);
            } else {
                v = f * (-z * sq[0] + 1 / sg);
            }
            out[hess_pos(P, i, j)] = v;
        }
    }
}

// The rule over [a, b] in the variable where the family's density is
// regular. On (0, Inf) a density may behave like y^(k-1) at zero, whose
// mass a tanh-sinh rule truncated at |t| = 3.2 misses for a small k; in
// v = log y it is a tail decaying like e^(k v), which the exp-sinh half
// reaches. On (0, 1) the same holds at both ends in v = logit(y). The nodes
// are mapped back to y and the weights multiplied by dy/dv; the center and
// the scale go over as v(c) and the scale's relative size, kept within
// [0.05, 20].
int cont_map(int id, const double* th) {
    if (id >= kWrapBase) {
        const WrapNode& w = node_of(id);
        if (w.code == 1) {
            double full[kMaxVec];
            int at[4];
            return cont_map(fixed_level(id, th, full, at), full);
        }
        // a folded or transformed variable: by its support, (0, Inf) in
        // log y and (0, 1) in logit y
        double lo, hi;
        cont_support(id, th, &lo, &hi);
        if (lo == 0.0 && hi == R_PosInf) return 1;
        if (lo == 0.0 && hi == 1.0) return 2;
        return 0;
    }
    switch (id) {
    case 1: case 7: case 9: case 14: case 15: case 16: case 17: case 18:
    case 21: case 22: case 25: case 26: case 30:
        return 1;
    case 4: case 31:
        return 2;
    default:
        return 0;
    }
}

void mapped_rule(int id, const double* th, double a, double b,
                 const double* kk, int nk, std::vector<double>& y,
                 std::vector<double>& w, std::vector<double>* ly,
                 std::vector<double>* l1y) {
    if (ly != nullptr) ly->clear();
    if (l1y != nullptr) l1y->clear();
    double c, s;
    cont_center_scale(id, th, &c, &s);
    const int mp = cont_map(id, th);
    if (mp == 0) {
        d7::trunc_rule(a, b, c, s, kk, nk, y, w);
        return;
    }
    auto fwd = [mp](double x) {
        if (mp == 1) return (x > 0) ? std::log(x) : R_NegInf;
        if (!(x > 0)) return R_NegInf;
        if (!(x < 1)) return R_PosInf;
        return std::log(x) - std::log1p(-x);
    };
    const double cv = fwd(c);
    double sv = (mp == 1) ? s / c : s / (c * (1 - c));
    sv = std::min(20.0, std::max(0.05, sv));
    double kv[4];
    for (int i = 0; i < nk && i < 4; ++i) kv[i] = fwd(kk[i]);
    d7::trunc_rule(fwd(a), fwd(b), cv, sv, kv, nk, y, w);
    // a node whose image leaves the open support, or whose weight is not
    // finite, is dropped: the density's tail has no mass there. On (0, 1),
    // when the caller takes the logarithms, a node whose y rounds to 0 or 1
    // is kept, since log y and log(1 - y) still place it, and dy/dv =
    // y (1 - y) is formed from them
    const bool logs = (mp == 2 && ly != nullptr && l1y != nullptr);
    std::size_t k = 0;
    for (std::size_t j = 0; j < y.size(); ++j) {
        const double v = y[j];
        double yy, ww;
        if (mp == 1) {
            yy = std::exp(v);
            ww = w[j] * yy;
            if (!(yy > 0) || !R_FINITE(yy)) continue;
        } else if (logs) {
            double a1, b1;
            logit_logs(v, &a1, &b1);
            yy = std::exp(a1);
            ww = w[j] * std::exp(a1 + b1);
            if (!R_FINITE(a1) || !R_FINITE(b1)) continue;
            if (!R_FINITE(ww) || !(ww > 0)) continue;
            ly->push_back(a1);
            l1y->push_back(b1);
        } else {
            yy = 1 / (1 + std::exp(-v));
            ww = w[j] * (yy * (1 - yy));
            if (!(yy > 0) || !(yy < 1)) continue;
        }
        if (!R_FINITE(ww) || !(ww > 0)) continue;
        y[k] = yy;
        w[k] = ww;
        ++k;
    }
    y.resize(k);
    w.resize(k);
}

// the compiled distribution function and its derivatives, by family; false
// where the family has none
bool cont_cdf(int id, double q, const double* th, bool lower, double* F) {
    if (id >= kWrapBase) {
        const WrapNode& w = node_of(id);
        if (w.code == 1) {
            double full[kMaxVec];
            int at[4];
            return cont_cdf(fixed_level(id, th, full, at), q, full, lower, F);
        }
        if (w.code == 5) {
            // F(q) - F(-q) below, S(q) + F(-q) above
            if (!(q > 0)) {
                *F = lower ? 0.0 : 1.0;
                return true;
            }
            double a, b;
            if (!cont_cdf(w.inner, q, th, lower, &a)) return false;
            if (!cont_cdf(w.inner, -q, th, true, &b)) return false;
            *F = lower ? a - b : a + b;
            return true;
        }
        const double* tp = transform_par(w.inner, th);
        const bool dec = d7::transform_decreasing(w.aux, tp);
        return cont_cdf(w.inner, d7::transform_inv(w.aux, q, tp), th,
                        dec ? !lower : lower, F);
    }
    switch (id) {
    case 0: *F = d7::gaussian1_cdf(q, th, lower); return true;
    case 32: *F = d7::enet_cdf(q, th, lower); return true;
    case 18: *F = d7::gpd_cdf(q, th, lower); return true;
    case 21: *F = d7::weibull3_cdf(q, th, lower); return true;
    case 22: *F = d7::lognormal2_cdf(q, th, lower); return true;
    case 15: *F = d7::invgauss1_cdf(q, th, lower); return true;
    case 16: *F = d7::invgauss2_cdf(q, th, lower); return true;
    case 10: *F = d7::cauchy_cdf(q, th, lower); return true;
    case 11: *F = d7::logistic_cdf(q, th, lower); return true;
    case 27: *F = d7::gumbel_cdf(q, th, lower); return true;
    case 28: *F = d7::laplace_cdf(q, th, lower); return true;
    case 12: *F = d7::gaussian2_cdf(q, th, lower); return true;
    case 13: *F = d7::gaussian3_cdf(q, th, lower); return true;
    case 29: *F = d7::laplace2_cdf(q, th, lower); return true;
    case 14: *F = d7::lognormal1_cdf(q, th, lower); return true;
    case 7: *F = d7::exponential_cdf(q, th, lower); return true;
    case 30: *F = d7::weibull1_cdf(q, th, lower); return true;
    default:
        if (!quad_cdf_family(id)) return false;
        *F = quad_cdf(id, q, th, lower);
        return true;
    }
}
bool cont_cdf_grad(int id, double q, const double* th, double* g) {
    if (id >= kWrapBase) {
        const WrapNode& w = node_of(id);
        const int np = w.np;
        if (w.code == 1) {
            double full[kMaxVec], gf[4];
            int at[4];
            const int inner = fixed_level(id, th, full, at);
            if (!cont_cdf_grad(inner, q, full, gf)) return false;
            for (int k = 0; k < np; ++k) g[k] = gf[at[k]];
            return true;
        }
        if (w.code == 5) {
            if (!(q > 0)) {
                for (int k = 0; k < np; ++k) g[k] = 0.0;
                return true;
            }
            double gp[4], gm[4];
            if (!cont_cdf_grad(w.inner, q, th, gp)) return false;
            if (!cont_cdf_grad(w.inner, -q, th, gm)) return false;
            for (int k = 0; k < np; ++k) g[k] = gp[k] - gm[k];
            return true;
        }
        const double* tp = transform_par(w.inner, th);
        const double sg = d7::transform_decreasing(w.aux, tp) ? -1.0 : 1.0;
        double gi[4];
        if (!cont_cdf_grad(w.inner, d7::transform_inv(w.aux, q, tp), th, gi))
            return false;
        for (int k = 0; k < np; ++k) g[k] = sg * gi[k];
        return true;
    }
    switch (id) {
    case 0: d7::gaussian1_cdf_grad(q, th, g); return true;
    case 32: d7::enet_cdf_grad(q, th, g); return true;
    case 18: d7::gpd_cdf_grad(q, th, g); return true;
    case 21: d7::weibull3_cdf_grad(q, th, g); return true;
    case 22: d7::lognormal2_cdf_grad(q, th, g); return true;
    case 15: d7::invgauss1_cdf_grad(q, th, g); return true;
    case 16: d7::invgauss2_cdf_grad(q, th, g); return true;
    case 10: d7::cauchy_cdf_grad(q, th, g); return true;
    case 11: d7::logistic_cdf_grad(q, th, g); return true;
    case 27: d7::gumbel_cdf_grad(q, th, g); return true;
    case 28: d7::laplace_cdf_grad(q, th, g); return true;
    case 12: d7::gaussian2_cdf_grad(q, th, g); return true;
    case 13: d7::gaussian3_cdf_grad(q, th, g); return true;
    case 29: d7::laplace2_cdf_grad(q, th, g); return true;
    case 14: d7::lognormal1_cdf_grad(q, th, g); return true;
    case 7: d7::exponential_cdf_grad(q, th, g); return true;
    case 30: d7::weibull1_cdf_grad(q, th, g); return true;
    default:
        if (!quad_cdf_family(id)) return false;
        quad_cdf_derivs(id, q, th, false, g);
        return true;
    }
}
bool cont_cdf_hess(int id, double q, const double* th, double* h) {
    if (id >= kWrapBase) {
        const WrapNode& w = node_of(id);
        const int np = w.np, np2 = np * (np + 1) / 2;
        if (w.code == 1) {
            double full[kMaxVec], hf[10];
            int at[4];
            const int inner = fixed_level(id, th, full, at);
            if (!cont_cdf_hess(inner, q, full, hf)) return false;
            const int P = id_np(inner);
            for (int i = 0; i < np; ++i)
                for (int j = i; j < np; ++j)
                    h[hess_pos(np, i, j)] = hf[hess_pos(P, at[i], at[j])];
            return true;
        }
        if (w.code == 5) {
            if (!(q > 0)) {
                for (int k = 0; k < np2; ++k) h[k] = 0.0;
                return true;
            }
            double hp[10], hm[10];
            if (!cont_cdf_hess(w.inner, q, th, hp)) return false;
            if (!cont_cdf_hess(w.inner, -q, th, hm)) return false;
            for (int k = 0; k < np2; ++k) h[k] = hp[k] - hm[k];
            return true;
        }
        const double* tp = transform_par(w.inner, th);
        const double sg = d7::transform_decreasing(w.aux, tp) ? -1.0 : 1.0;
        double hi[10];
        if (!cont_cdf_hess(w.inner, d7::transform_inv(w.aux, q, tp), th, hi))
            return false;
        for (int k = 0; k < np2; ++k) h[k] = sg * hi[k];
        return true;
    }
    switch (id) {
    case 0: d7::gaussian1_cdf_hess(q, th, h); return true;
    case 32: d7::enet_cdf_hess(q, th, h); return true;
    case 18: d7::gpd_cdf_hess(q, th, h); return true;
    case 21: d7::weibull3_cdf_hess(q, th, h); return true;
    case 22: d7::lognormal2_cdf_hess(q, th, h); return true;
    case 15: d7::invgauss1_cdf_hess(q, th, h); return true;
    case 16: d7::invgauss2_cdf_hess(q, th, h); return true;
    case 10: d7::cauchy_cdf_hess(q, th, h); return true;
    case 11: d7::logistic_cdf_hess(q, th, h); return true;
    case 27: d7::gumbel_cdf_hess(q, th, h); return true;
    case 28: d7::laplace_cdf_hess(q, th, h); return true;
    case 12: d7::gaussian2_cdf_hess(q, th, h); return true;
    case 13: d7::gaussian3_cdf_hess(q, th, h); return true;
    case 29: d7::laplace2_cdf_hess(q, th, h); return true;
    case 14: d7::lognormal1_cdf_hess(q, th, h); return true;
    case 7: d7::exponential_cdf_hess(q, th, h); return true;
    case 30: d7::weibull1_cdf_hess(q, th, h); return true;
    default:
        if (!quad_cdf_family(id)) return false;
        quad_cdf_derivs(id, q, th, true, h);
        return true;
    }
}

// the position of pair (i, j) in hess_names() order
int hess_pos(int P, int i, int j) {
    if (i == j) return i;
    if (i > j) std::swap(i, j);
    int pos = P;
    for (int r = 0; r < i; ++r) pos += P - 1 - r;
    return pos + (j - i - 1);
}

struct TruncEnds { double z; double g[4]; double h[10]; bool ok; };

// the retained mass and, as `level` asks (0, 1 or 2), its first and
// second derivatives in every parameter, at one observation; inner a base
// family or a chain of fixed() nodes over one. The density needs the mass alone, and for a family whose
// derivatives are quadratures the derivatives cost far more than the mass.
TruncEnds trunc_cont_ends(int inner, const double* th, double lo, double up,
                          int level) {
    TruncEnds e;
    const int P = id_np(inner), np2 = P * (P + 1) / 2;
    double slo, shi;
    cont_support(inner, th, &slo, &shi);
    const bool lo_in = lo > slo, up_in = up < shi;
    double Fa = 0.0, gl[4] = {0, 0, 0, 0}, gu[4] = {0, 0, 0, 0};
    double hl[10] = {0}, hu[10] = {0};
    e.ok = true;
    if (lo_in) {
        e.ok = e.ok && cont_cdf(inner, lo, th, true, &Fa);
        if (level >= 1) e.ok = e.ok && cont_cdf_grad(inner, lo, th, gl);
        if (level >= 2) e.ok = e.ok && cont_cdf_hess(inner, lo, th, hl);
    }
    if (up_in) {
        if (level >= 1) e.ok = e.ok && cont_cdf_grad(inner, up, th, gu);
        if (level >= 2) e.ok = e.ok && cont_cdf_hess(inner, up, th, hu);
    }
    if (!e.ok) return e;
    if (Fa > 0.5) {
        double Sa = 1.0, Sb = 0.0;
        if (lo_in) cont_cdf(inner, lo, th, false, &Sa);
        if (up_in) cont_cdf(inner, up, th, false, &Sb);
        e.z = Sa - Sb;
    } else {
        double Fb = 1.0;
        if (up_in) cont_cdf(inner, up, th, true, &Fb);
        e.z = Fb - Fa;
    }
    for (int i = 0; i < P; ++i) e.g[i] = gu[i] - gl[i];
    for (int i = 0; i < np2; ++i) e.h[i] = hu[i] - hl[i];
    return e;
}

// the nodes and weights of the information's rule at one observation
// at one observation, and for a beta the logarithms its entries read
void trunc_cont_nodes(int inner, const double* th, double lo, double up,
                      std::vector<double>& y, std::vector<double>& w,
                      std::vector<double>* ly = nullptr,
                      std::vector<double>* l1y = nullptr) {
    double slo, shi, c, s, kk[4];
    cont_support(inner, th, &slo, &shi);
    if (!cont_center_scale(inner, th, &c, &s)) {
        y.clear();
        w.clear();
        if (ly != nullptr) ly->clear();
        if (l1y != nullptr) l1y->clear();
        return;
    }
    const int nk = cont_kinks(inner, th, kk);
    mapped_rule(inner, th, std::max(lo, slo), std::min(up, shi), kk, nk, y, w,
                ly, l1y);
}

void trunc_cont_score_curv(int inner, int k, double y, const double* th,
                           double* out) {
    const double* tp = trunc_points(inner, th);
    const TruncEnds e = trunc_cont_ends(inner, th, tp[0], tp[1], 2);
    if (!e.ok) { out[0] = out[1] = R_NaN; return; }
    const int P = id_np(inner);
    double sc[2];
    d7_score_curv(inner, k, y, th, sc);
    const double m = e.g[k] / e.z, M = e.h[hess_pos(P, k, k)] / e.z;
    out[0] = sc[0] - m;
    out[1] = sc[1] - M + m * m;
}

void trunc_cont_info_dinfo(int inner, int k, double y, const double* th,
                           double* out) {
    const double* tp = trunc_points(inner, th);
    const TruncEnds e = trunc_cont_ends(inner, th, tp[0], tp[1], 2);
    if (!e.ok) { out[0] = out[1] = R_NaN; return; }
    const int P = id_np(inner);
    const bool logs = logs_family(inner);
    std::vector<double> ys, ws, lys, l1ys;
    trunc_cont_nodes(inner, th, tp[0], tp[1], ys, ws, logs ? &lys : nullptr,
                     logs ? &l1ys : nullptr);
    if (ys.empty()) { out[0] = out[1] = R_NaN; return; }
    long double S = 0.0L, D = 0.0L;
    for (std::size_t j = 0; j < ys.size(); ++j) {
        double sc[2], lf;
        if (logs) {
            double gg[2], hh[3];
            lf = beta_logs_parts(inner, lys[j], l1ys[j], th, gg, hh);
            sc[0] = gg[k];
            sc[1] = hh[k];
        } else {
            lf = d7_logpdf(inner, ys[j], th);
        }
        const double fw = std::exp(lf) * ws[j];
        if (fw == 0.0) continue;
        if (!logs) d7_score_curv(inner, k, ys[j], th, sc);
        const double g = sc[0], l = sc[1];
        S += fw * (g * g);
        D += fw * (g * g * g + l * g + g * l);
    }
    const double z = e.z, zk = e.g[k], zkk = e.h[hess_pos(P, k, k)];
    const double m = zk / z, M = zkk / z, dm = M - m * m;
    const double Sd = (double) S, Dd = (double) D;
    out[0] = -(Sd / z - m * m);
    out[1] = -(Dd / z - Sd * zk / (z * z) - (dm * m + m * dm));
}

double trunc_cont_logpdf(int inner, double y, const double* th) {
    const double* tp = trunc_points(inner, th);
    const TruncEnds e = trunc_cont_ends(inner, th, tp[0], tp[1], 0);
    if (!e.ok) return R_NaN;
    const double ld = d7_logpdf(inner, y, th) - std::log(e.z);
    return (y < tp[0] || y > tp[1]) ? R_NegInf : ld;
}

void wrap_score_curv(int id, int k, double y, const double* th, double* out) {
    const WrapNode& w = node_of(id);
    const int code = w.code, aux = w.aux, inner = w.inner;
    double full[kMaxVec];
    if (code == 1) {
        const int j = fixed_inner(inner, aux, k, th, full);
        d7_score_curv(inner, j, y, full, out);
        return;
    }
    if (code == 5) {
        fold_score_curv(inner, k, y, th, out);
        return;
    }
    if (code == 6) {
        d7_score_curv(inner, k, d7::transform_inv(aux, y, transform_par(inner, th)),
                      th, out);
        return;
    }
    if (code == 7) {
        trunc_score_curv(inner, k, y, th, out);
        return;
    }
    if (code == 8) {
        trunc_cont_score_curv(inner, k, y, th, out);
        return;
    }
    if (code >= 2 && code <= 4) {
        const int P = id_np(inner);
        const double z = zero_inner(inner, th, full);
        const double f0 = (code == 4) ? 0.0 : std::exp(d7_logpdf(inner, 0.0, full));
        if (k == P) {
            if (code == 2) {
                out[0] = d7::zi_score_zi(y, z, f0);
                out[1] = d7::zi_curv_zi(y, z, f0);
            } else {
                out[0] = d7::za_score_za(y, z);
                out[1] = d7::za_curv_za(y, z);
            }
            return;
        }
        double gh[2], sh0[2] = {0.0, 0.0};
        d7_score_curv(inner, k, y, full, gh);
        if (code != 4) d7_score_curv(inner, k, 0.0, full, sh0);
        if (code == 2) {
            out[0] = d7::zi_score_parent(y, z, f0, gh[0]);
            out[1] = d7::zi_curv_parent(y, z, f0, sh0[0], sh0[1], gh[1]);
        } else if (code == 3) {
            out[0] = d7::zad_score_parent(y, f0, sh0[0], gh[0]);
            out[1] = d7::zad_curv_parent(y, f0, sh0[0], sh0[1], gh[1]);
        } else {
            out[0] = d7::zac_score_parent(y, gh[0]);
            out[1] = d7::zac_curv_parent(y, gh[1]);
        }
        return;
    }
    out[0] = R_NaN; out[1] = R_NaN;
}

void wrap_info_dinfo(int id, int k, double y, const double* th, double* out) {
    const WrapNode& w = node_of(id);
    const int code = w.code, aux = w.aux, inner = w.inner;
    double full[kMaxVec];
    if (code == 1) {
        const int j = fixed_inner(inner, aux, k, th, full);
        d7_info_dinfo(inner, j, y, full, out);
        return;
    }
    if (code == 5) {
        fold_info_dinfo(inner, k, th, out);
        return;
    }
    if (code == 6) {
        d7_info_dinfo(inner, k, d7::transform_inv(aux, y, transform_par(inner, th)),
                      th, out);
        return;
    }
    if (code == 7) {
        trunc_info_dinfo(inner, k, y, th, out);
        return;
    }
    if (code == 8) {
        trunc_cont_info_dinfo(inner, k, y, th, out);
        return;
    }
    if (code >= 2 && code <= 4) {
        const int P = id_np(inner);
        const double z = zero_inner(inner, th, full);
        const double f0 = (code == 4) ? 0.0 : std::exp(d7_logpdf(inner, 0.0, full));
        if (k == P) {
            if (code == 2) {
                out[0] = d7::zi_info_zi(z, f0);
                out[1] = d7::zi_dinfo_zi(z, f0);
            } else {
                out[0] = d7::za_info_za(z);
                out[1] = d7::za_dinfo_za(z);
            }
            return;
        }
        double ed[2], sh0[2] = {0.0, 0.0};
        d7_info_dinfo(inner, k, y, full, ed);
        if (code != 4) d7_score_curv(inner, k, 0.0, full, sh0);
        if (code == 2) {
            out[0] = d7::zi_info_parent(z, f0, sh0[0], sh0[1], ed[0]);
            out[1] = d7::zi_dinfo_parent(z, f0, sh0[0], sh0[1], ed[1]);
        } else if (code == 3) {
            out[0] = d7::zad_info_parent(z, f0, sh0[0], sh0[1], ed[0]);
            out[1] = d7::zad_dinfo_parent(z, f0, sh0[0], sh0[1], ed[0], ed[1]);
        } else {
            out[0] = d7::zac_info_parent(z, ed[0]);
            out[1] = d7::zac_dinfo_parent(z, ed[1]);
        }
        return;
    }
    out[0] = R_NaN; out[1] = R_NaN;
}

double wrap_logpdf(int id, double y, const double* th) {
    const WrapNode& w = node_of(id);
    const int code = w.code, aux = w.aux, inner = w.inner;
    double full[kMaxVec];
    if (code == 1) {
        fixed_inner(inner, aux, 0, th, full);
        return d7_logpdf(inner, y, full);
    }
    if (code == 5) {
        const d7::FoldW W = d7::fold_w(d7_logpdf(inner, y, th),
                                       d7_logpdf(inner, -y, th));
        return d7::fold_logpdf(y, W);
    }
    if (code == 6) {
        const double* tp = transform_par(inner, th);
        return d7::transform_logpdf(
            d7_logpdf(inner, d7::transform_inv(aux, y, tp), th),
            d7::transform_log_jac(aux, y, tp));
    }
    if (code == 7) return trunc_logpdf(inner, y, th);
    if (code == 8) return trunc_cont_logpdf(inner, y, th);
    if (code >= 2 && code <= 4) {
        const double z = zero_inner(inner, th, full);
        const double lf = d7_logpdf(inner, y, full);
        if (code == 2) return d7::zi_logpdf(y, z, lf);
        if (code == 4) return d7::zac_logpdf(y, z, lf);
        const double f0 = std::exp(d7_logpdf(inner, 0.0, full));
        return d7::zad_logpdf(y, z, f0, lf);
    }
    return R_NaN;
}

}  // namespace

extern "C" {

// Called once by the consumer on its own thread, before any loop: the von
// Mises families resolve numericals7's Bessel-ratio entry points here, and
// the skew t numericals7's t distribution function and the pseudo-Huber
// families its Bessel K, so that no lookup into R happens later from a
// worker.
int d7_scalar_id(const char* cls) {
    const char* bar = std::strchr(cls, '|');
    if (bar != nullptr) return wrap_id(cls, bar);
    for (int i = 0; i < d7_n_scalar_classes; ++i) {
        if (std::strcmp(cls, d7_scalar_classes[i]) == 0) {
            if (std::strncmp(cls, "VonMises", 8) == 0) {
                d7::vm_bessel();
                d7::vm_log_i0_fn();
            }
            if (std::strcmp(cls, "SkewTDistrib") == 0) d7::skewt_pt();
            if (std::strncmp(cls, "PseudoHuber", 11) == 0) d7::bessel_k_fn();
            return i;
        }
    }
    return -1;
}

// k is the 0-based index of the parameter among the family's own, th the
// full parameter vector at this observation; out[0] the score component,
// out[1] the (k, k) second derivative, both on the parameter scale
void d7_score_curv(int id, int k, double y, const double* th, double* out) {
    if (id >= kWrapBase) {
        if (node_valid(id)) wrap_score_curv(id, k, y, th, out);
        else { out[0] = R_NaN; out[1] = R_NaN; }
        return;
    }
    switch (id) {
    case 0: d7::gaussian1_score_curv(k, y, th, out); break;
    case 1: d7::gamma1_score_curv(k, y, th, out); break;
    case 2: d7::poisson_score_curv(k, y, th, out); break;
    case 3: d7::negbin2_score_curv(k, y, th, out); break;
    case 4: d7::beta1_score_curv(k, y, th, out); break;
    case 5: d7::bernoulli_score_curv(k, y, th, out); break;
    case 6: d7::binomial_score_curv(k, y, th, out); break;
    case 7: d7::exponential_score_curv(k, y, th, out); break;
    case 8: d7::geometric_score_curv(k, y, th, out); break;
    case 9: d7::chisq_score_curv(k, y, th, out); break;
    case 10: d7::cauchy_score_curv(k, y, th, out); break;
    case 11: d7::logistic_score_curv(k, y, th, out); break;
    case 12: d7::gaussian2_score_curv(k, y, th, out); break;
    case 13: d7::gaussian3_score_curv(k, y, th, out); break;
    case 14: d7::lognormal1_score_curv(k, y, th, out); break;
    case 15: d7::invgauss1_score_curv(k, y, th, out); break;
    case 16: d7::invgauss2_score_curv(k, y, th, out); break;
    case 17: d7::gamma2_score_curv(k, y, th, out); break;
    case 18: d7::gpd_score_curv(k, y, th, out); break;
    case 19: d7::vonmises1_score_curv(k, y, th, out); break;
    case 20: d7::vonmises2_score_curv(k, y, th, out); break;
    case 21: d7::weibull3_score_curv(k, y, th, out); break;
    case 22: d7::lognormal2_score_curv(k, y, th, out); break;
    case 23: d7::student_t1_score_curv(k, y, th, out); break;
    case 24: d7::student_t2_score_curv(k, y, th, out); break;
    case 25: d7::gengamma1_score_curv(k, y, th, out); break;
    case 26: d7::gengamma2_score_curv(k, y, th, out); break;
    case 27: d7::gumbel_score_curv(k, y, th, out); break;
    case 28: d7::laplace_score_curv(k, y, th, out); break;
    case 29: d7::laplace2_score_curv(k, y, th, out); break;
    case 30: d7::weibull1_score_curv(k, y, th, out); break;
    case 31: d7::beta2_score_curv(k, y, th, out); break;
    case 32: d7::enet_score_curv(k, y, th, out); break;
    case 33: d7::negbin1_score_curv(k, y, th, out); break;
    case 34: d7::skewnormal1_score_curv(k, y, th, out); break;
    case 35: d7::skewt_score_curv(k, y, th, out); break;
    case 36: d7::pseudohuber_score_curv(k, y, th, out); break;
    case 37: d7::pseudohuber2_score_curv(k, y, th, out); break;
    case 38: d7::skewnormal2_score_curv(k, y, th, out); break;
    case 39: d7::pig1_score_curv(k, y, th, out); break;
    case 40: d7::pig2_score_curv(k, y, th, out); break;
    case 41: d7::betabinom1_score_curv(k, y, th, out); break;
    case 42: d7::betabinom2_score_curv(k, y, th, out); break;
    default:
        out[0] = R_NaN; out[1] = R_NaN;
    }
}

// 1 when the family's entries never reach the R API, so that a consumer may
// call them from a worker thread, 0 otherwise; -1 for an unknown id. The
// audit behind the 1s: the entries call R's digamma, trigamma and psigamma
// (orders below 100), whose failures set errno without a warning;
// lgammafn and dt at positive arguments, where nmath's warnings (precision
// near a negative integer, lgammacor's underflow past 3.7e306, reached only
// above 4.9e6) cannot fire; dnorm and pnorm, whose only warning is the
// silent ME_DOMAIN; numericals7's Bessel ratio, which is plain C; and
// numericals7's n7_pt and n7_bessel_k, R's pt() and bessel_k() compiled
// without their warnings and, for the latter, without R's allocator. The
// quadrature families (pt_loc_scale.h) keep their cache in thread-local
// storage. A family whose entries evaluate a function that may warn (R's
// pt, pbeta, lchoose, R's Bessel K) answers 0. A transformed family adds
// R_pow and nmath's plogis, qlogis and dlogis, which never warn.
int d7_scalar_thread_safe(int id) {
    if (id >= kWrapBase) {
        if (!node_valid(id)) return -1;
        id = base_of(id);
    }
    if (id < 0 || id >= d7_n_scalar_classes) return -1;
    return 1;
}

// out[0] the (k, k) expected second derivative E[l_kk], out[1] its
// derivative in the same parameter, both on the parameter scale
void d7_info_dinfo(int id, int k, double y, const double* th, double* out) {
    if (id >= kWrapBase) {
        if (node_valid(id)) wrap_info_dinfo(id, k, y, th, out);
        else { out[0] = R_NaN; out[1] = R_NaN; }
        return;
    }
    switch (id) {
    case 0: d7::gaussian1_info_dinfo(k, y, th, out); break;
    case 1: d7::gamma1_info_dinfo(k, y, th, out); break;
    case 2: d7::poisson_info_dinfo(k, y, th, out); break;
    case 3: d7::negbin2_info_dinfo(k, y, th, out); break;
    case 4: d7::beta1_info_dinfo(k, y, th, out); break;
    case 5: d7::bernoulli_info_dinfo(k, y, th, out); break;
    case 6: d7::binomial_info_dinfo(k, y, th, out); break;
    case 7: d7::exponential_info_dinfo(k, y, th, out); break;
    case 8: d7::geometric_info_dinfo(k, y, th, out); break;
    case 9: d7::chisq_info_dinfo(k, y, th, out); break;
    case 10: d7::cauchy_info_dinfo(k, y, th, out); break;
    case 11: d7::logistic_info_dinfo(k, y, th, out); break;
    case 12: d7::gaussian2_info_dinfo(k, y, th, out); break;
    case 13: d7::gaussian3_info_dinfo(k, y, th, out); break;
    case 14: d7::lognormal1_info_dinfo(k, y, th, out); break;
    case 15: d7::invgauss1_info_dinfo(k, y, th, out); break;
    case 16: d7::invgauss2_info_dinfo(k, y, th, out); break;
    case 17: d7::gamma2_info_dinfo(k, y, th, out); break;
    case 18: d7::gpd_info_dinfo(k, y, th, out); break;
    case 19: d7::vonmises1_info_dinfo(k, y, th, out); break;
    case 20: d7::vonmises2_info_dinfo(k, y, th, out); break;
    case 21: d7::weibull3_info_dinfo(k, y, th, out); break;
    case 22: d7::lognormal2_info_dinfo(k, y, th, out); break;
    case 23: d7::student_t1_info_dinfo(k, y, th, out); break;
    case 24: d7::student_t2_info_dinfo(k, y, th, out); break;
    case 25: d7::gengamma1_info_dinfo(k, y, th, out); break;
    case 26: d7::gengamma2_info_dinfo(k, y, th, out); break;
    case 27: d7::gumbel_info_dinfo(k, y, th, out); break;
    case 28: d7::laplace_info_dinfo(k, y, th, out); break;
    case 29: d7::laplace2_info_dinfo(k, y, th, out); break;
    case 30: d7::weibull1_info_dinfo(k, y, th, out); break;
    case 31: d7::beta2_info_dinfo(k, y, th, out); break;
    case 32: d7::enet_info_dinfo(k, y, th, out); break;
    case 33: d7::negbin1_info_dinfo(k, y, th, out); break;
    case 34: d7::skewnormal1_info_dinfo(k, y, th, out); break;
    case 35: d7::skewt_info_dinfo(k, y, th, out); break;
    case 36: d7::pseudohuber_info_dinfo(k, y, th, out); break;
    case 37: d7::pseudohuber2_info_dinfo(k, y, th, out); break;
    case 38: d7::skewnormal2_info_dinfo(k, y, th, out); break;
    case 39: d7::pig1_info_dinfo(k, y, th, out); break;
    case 40: d7::pig2_info_dinfo(k, y, th, out); break;
    case 41: d7::betabinom1_info_dinfo(k, y, th, out); break;
    case 42: d7::betabinom2_info_dinfo(k, y, th, out); break;
    default:
        out[0] = R_NaN; out[1] = R_NaN;
    }
}


// the log-density of one observation, th as for d7_score_curv(); NaN for an
// unknown id
double d7_logpdf(int id, double y, const double* th) {
    if (id >= kWrapBase) return node_valid(id) ? wrap_logpdf(id, y, th) : R_NaN;
    switch (id) {
    case 0: return d7::gaussian1_logpdf(y, th);
    case 1: return d7::gamma1_logpdf(y, th);
    case 2: return d7::poisson_logpdf(y, th);
    case 3: return d7::negbin2_logpdf(y, th);
    case 4: return d7::beta1_logpdf(y, th);
    case 5: return d7::bernoulli_logpdf(y, th);
    case 6: return d7::binomial_logpdf(y, th);
    case 7: return d7::exponential_logpdf(y, th);
    case 8: return d7::geometric_logpdf(y, th);
    case 9: return d7::chisq_logpdf(y, th);
    case 10: return d7::cauchy_logpdf(y, th);
    case 11: return d7::logistic_logpdf(y, th);
    case 12: return d7::gaussian2_logpdf(y, th);
    case 13: return d7::gaussian3_logpdf(y, th);
    case 14: return d7::lognormal1_logpdf(y, th);
    case 15: return d7::invgauss1_logpdf(y, th);
    case 16: return d7::invgauss2_logpdf(y, th);
    case 17: return d7::gamma2_logpdf(y, th);
    case 18: return d7::gpd_logpdf(y, th);
    case 19: return d7::vonmises1_logpdf(y, th);
    case 20: return d7::vonmises2_logpdf(y, th);
    case 21: return d7::weibull3_logpdf(y, th);
    case 22: return d7::lognormal2_logpdf(y, th);
    case 23: return d7::student_t1_logpdf(y, th);
    case 24: return d7::student_t2_logpdf(y, th);
    case 25: return d7::gengamma1_logpdf(y, th);
    case 26: return d7::gengamma2_logpdf(y, th);
    case 27: return d7::gumbel_logpdf(y, th);
    case 28: return d7::laplace_logpdf(y, th);
    case 29: return d7::laplace2_logpdf(y, th);
    case 30: return d7::weibull1_logpdf(y, th);
    case 31: return d7::beta2_logpdf(y, th);
    case 32: return d7::enet_logpdf(y, th);
    case 33: return d7::negbin1_logpdf(y, th);
    case 34: return d7::skewnormal1_logpdf_th(y, th);
    case 35: return d7::skewt_logpdf_th(y, th);
    case 36: return d7::pseudohuber_logpdf_th(y, th);
    case 37: return d7::pseudohuber2_logpdf_th(y, th);
    case 38: return d7::skewnormal2_logpdf_th(y, th);
    case 39: return d7::pig1_logpdf(y, th);
    case 40: return d7::pig2_logpdf(y, th);
    case 41: return d7::betabinom1_logpdf(y, th);
    case 42: return d7::betabinom2_logpdf(y, th);
    default: return R_NaN;
    }
}

} // extern "C"

// exposed to this package's own tests: the twin comparison against
// distrib_gradient()/distrib_hessian() lives where the formulas do
// [[Rcpp::export]]
Rcpp::List d7_scalar_probe(std::string cls, int k, Rcpp::NumericVector y,
                           Rcpp::NumericMatrix theta) {
    int id = d7_scalar_id(cls.c_str());
    int n = y.size(), np = theta.ncol();
    Rcpp::NumericVector g(n), h(n);
    std::vector<double> th(np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        double out[2];
        d7_score_curv(id, k - 1, y[i], th.data(), out);
        g[i] = out[0]; h[i] = out[1];
    }
    return Rcpp::List::create(Rcpp::_["id"] = id, Rcpp::_["score"] = g,
                              Rcpp::_["curvature"] = h);
}

// the twin of d7_scalar_probe() for the expected information
// [[Rcpp::export]]
Rcpp::List d7_info_probe(std::string cls, int k, Rcpp::NumericVector y,
                         Rcpp::NumericMatrix theta) {
    int id = d7_scalar_id(cls.c_str());
    int n = y.size(), np = theta.ncol();
    Rcpp::NumericVector e(n), de(n);
    std::vector<double> th(np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        double out[2];
        d7_info_dinfo(id, k - 1, y[i], th.data(), out);
        e[i] = out[0]; de[i] = out[1];
    }
    return Rcpp::List::create(Rcpp::_["id"] = id, Rcpp::_["expected"] = e,
                              Rcpp::_["dexpected"] = de);
}

// d7_scalar_thread_safe() by class name, for the tests
// the twin of d7_scalar_probe() for the log-density
// [[Rcpp::export]]
Rcpp::List d7_logpdf_probe(std::string cls, Rcpp::NumericVector y,
                           Rcpp::NumericMatrix theta) {
    int id = d7_scalar_id(cls.c_str());
    int n = y.size(), np = theta.ncol();
    Rcpp::NumericVector lp(n);
    std::vector<double> th(np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        lp[i] = d7_logpdf(id, y[i], th.data());
    }
    return Rcpp::List::create(Rcpp::_["id"] = id, Rcpp::_["logpdf"] = lp);
}

// the center and scale folded()'s quadrature reads, one row per observation,
// NA where the family has none
// [[Rcpp::export]]
Rcpp::NumericMatrix d7_center_scale_probe(std::string cls,
                                          Rcpp::NumericMatrix theta) {
    int id = d7_scalar_id(cls.c_str());
    int n = theta.nrow(), np = theta.ncol();
    Rcpp::NumericMatrix out(n, 2);
    std::vector<double> th(np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        double c = NA_REAL, s = NA_REAL;
        if (id < 0 || (id >= kWrapBase && !node_valid(id)) ||
            !center_scale(id, th.data(), &c, &s)) {
            c = NA_REAL;
            s = NA_REAL;
        }
        out(i, 0) = c;
        out(i, 1) = s;
    }
    return out;
}

// trunc_disc_rule() for trunc_route_parts() in R, one row of theta per
// observation (the truncated family's parameters, the parent's constants,
// lower and upper), as the route of the truncated family names them; `info`
// as for trunc_disc_rule()
// [[Rcpp::export]]
Rcpp::List trunc_rule_cpp(std::string cls, Rcpp::NumericMatrix theta,
                          bool info) {
    const int id = d7_scalar_id(cls.c_str());
    const int n = theta.nrow(), np = theta.ncol();
    Rcpp::IntegerVector branch(n);
    Rcpp::NumericVector y0(n), y1(n);
    if (!node_valid(id) || node_of(id).code != 7)
        Rcpp::stop("'%s' is not the route of a truncated family.", cls);
    const int inner = node_of(id).inner;
    std::vector<double> th(np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        const double* tp = trunc_points(inner, th.data());
        const TruncRule r = trunc_disc_rule(inner, th.data(), tp[0], tp[1], info);
        branch[i] = r.branch; y0[i] = r.y0; y1[i] = r.y1;
    }
    return Rcpp::List::create(Rcpp::_["branch"] = branch, Rcpp::_["y0"] = y0,
                              Rcpp::_["y1"] = y1);
}

// sums of x by group in index order, accumulated in long double as the
// registry's loops accumulate them; g counts groups from one
// [[Rcpp::export]]
Rcpp::NumericVector ld_group_sum(Rcpp::NumericVector x, Rcpp::IntegerVector g,
                                 int n) {
    std::vector<long double> acc(n, 0.0L);
    for (R_xlen_t i = 0; i < x.size(); ++i) acc[g[i] - 1] += x[i];
    Rcpp::NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = (double) acc[i];
    return out;
}

// the rule of pt_trunc_rule.h, for the tests
// [[Rcpp::export]]
Rcpp::List trunc_cont_rule_raw_cpp(double a, double b, double c, double s,
                                   Rcpp::NumericVector kinks,
                                   double h = 0.0625, int nd = 6) {
    std::vector<double> y, w;
    d7::trunc_rule(a, b, c, s, kinks.begin(), kinks.size(), y, w, h, nd);
    return Rcpp::List::create(Rcpp::_["y"] = Rcpp::wrap(y),
                              Rcpp::_["w"] = Rcpp::wrap(w));
}

// trunc_cont_ends() for R, one row of theta per observation as for
// trunc_rule_cpp(): the retained mass, its gradient (one column per
// parameter) and its Hessian (one column per pair, in hess_names() order)
// [[Rcpp::export]]
Rcpp::List trunc_cont_ends_cpp(std::string cls, Rcpp::NumericMatrix theta,
                               int level = 2) {
    const int id = d7_scalar_id(cls.c_str());
    if (!node_valid(id) || node_of(id).code != 8)
        Rcpp::stop("'%s' is not the route of a truncated continuous family.", cls);
    const int inner = node_of(id).inner, P = id_np(inner);
    const int n = theta.nrow(), np = theta.ncol(), np2 = P * (P + 1) / 2;
    Rcpp::NumericVector z(n);
    Rcpp::NumericMatrix g(n, P), h(n, np2);
    std::vector<double> th(np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        const double* tp = trunc_points(inner, th.data());
        const TruncEnds e = trunc_cont_ends(inner, th.data(), tp[0], tp[1], level);
        z[i] = e.ok ? e.z : NA_REAL;
        for (int j = 0; j < P; ++j) g(i, j) = e.g[j];
        for (int j = 0; j < np2; ++j) h(i, j) = e.h[j];
    }
    return Rcpp::List::create(Rcpp::_["z"] = z, Rcpp::_["g"] = g,
                              Rcpp::_["h"] = h);
}

// the rule's nodes for R: observation index (from one), node and weight
// [[Rcpp::export]]
Rcpp::List trunc_cont_rule_cpp(std::string cls, Rcpp::NumericMatrix theta) {
    const int id = d7_scalar_id(cls.c_str());
    if (!node_valid(id) || node_of(id).code != 8)
        Rcpp::stop("'%s' is not the route of a truncated continuous family.", cls);
    const int inner = node_of(id).inner;
    const int n = theta.nrow(), np = theta.ncol();
    const bool logs = logs_family(inner);
    std::vector<int> idx;
    std::vector<double> yy, ww, ys, ws, th(np), la, lb, las, lbs;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        const double* tp = trunc_points(inner, th.data());
        trunc_cont_nodes(inner, th.data(), tp[0], tp[1], ys, ws,
                         logs ? &las : nullptr, logs ? &lbs : nullptr);
        for (std::size_t j = 0; j < ys.size(); ++j) {
            idx.push_back(i + 1);
            yy.push_back(ys[j]);
            ww.push_back(ws[j]);
            if (logs) {
                la.push_back(las[j]);
                lb.push_back(lbs[j]);
            }
        }
    }
    Rcpp::List out = Rcpp::List::create(Rcpp::_["idx"] = Rcpp::wrap(idx),
                                        Rcpp::_["y"] = Rcpp::wrap(yy),
                                        Rcpp::_["w"] = Rcpp::wrap(ww));
    if (logs) {
        out["ly"] = Rcpp::wrap(la);
        out["l1y"] = Rcpp::wrap(lb);
    }
    return out;
}

// beta_logs_parts() for R at the rule's nodes, one row of theta per node:
// the log-density, the score (one column per free parameter) and the
// Hessian (one column per pair, in hess_names() order)
// [[Rcpp::export]]
Rcpp::List beta_logs_parts_cpp(std::string cls, Rcpp::NumericVector ly,
                               Rcpp::NumericVector l1y,
                               Rcpp::NumericMatrix theta) {
    const int id = d7_scalar_id(cls.c_str());
    if (id < 0 || !logs_family(id))
        Rcpp::stop("'%s' is not the route of a beta family.", cls);
    const int P = id_np(id), np2 = P * (P + 1) / 2;
    const int n = ly.size(), np = theta.ncol();
    Rcpp::NumericVector lf(n);
    Rcpp::NumericMatrix g(n, P), h(n, np2);
    std::vector<double> th(np);
    double gg[2], hh[3];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        lf[i] = beta_logs_parts(id, ly[i], l1y[i], th.data(), gg, hh);
        for (int j = 0; j < P; ++j) g(i, j) = gg[j];
        for (int j = 0; j < np2; ++j) h(i, j) = hh[j];
    }
    return Rcpp::List::create(Rcpp::_["logpdf"] = lf, Rcpp::_["g"] = g,
                              Rcpp::_["h"] = h);
}

// the compiled distribution function and its derivatives, by class name,
// one row of theta per point: what the R methods of the families without
// closed forms in R read
// [[Rcpp::export]]
Rcpp::NumericVector d7_cdf_cpp(std::string cls, Rcpp::NumericVector q,
                               Rcpp::NumericMatrix theta, bool lower) {
    const int id = d7_scalar_id(cls.c_str());
    const int n = q.size(), np = theta.ncol();
    Rcpp::NumericVector out(n);
    std::vector<double> th(np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        double F;
        out[i] = (id >= 0 && id < kWrapBase &&
                  cont_cdf(id, q[i], th.data(), lower, &F)) ? F : NA_REAL;
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix d7_cdf_grad_cpp(std::string cls, Rcpp::NumericVector q,
                                    Rcpp::NumericMatrix theta) {
    const int id = d7_scalar_id(cls.c_str());
    if (id < 0 || id >= kWrapBase) Rcpp::stop("no compiled cdf for '%s'", cls);
    const int n = q.size(), np = theta.ncol(), P = d7_n_params[id];
    Rcpp::NumericMatrix out(n, P);
    std::vector<double> th(np);
    double g[4];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        if (!cont_cdf_grad(id, q[i], th.data(), g)) Rcpp::stop("no compiled cdf");
        for (int j = 0; j < P; ++j) out(i, j) = g[j];
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix d7_cdf_hess_cpp(std::string cls, Rcpp::NumericVector q,
                                    Rcpp::NumericMatrix theta) {
    const int id = d7_scalar_id(cls.c_str());
    if (id < 0 || id >= kWrapBase) Rcpp::stop("no compiled cdf for '%s'", cls);
    const int n = q.size(), np = theta.ncol(), P = d7_n_params[id];
    const int np2 = P * (P + 1) / 2;
    Rcpp::NumericMatrix out(n, np2);
    std::vector<double> th(np);
    double h[10];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        if (!cont_cdf_hess(id, q[i], th.data(), h)) Rcpp::stop("no compiled cdf");
        for (int j = 0; j < np2; ++j) out(i, j) = h[j];
    }
    return out;
}

// fold_rule() for folded_expected() in R, which reads the same nodes
// [[Rcpp::export]]
Rcpp::List fold_rule_cpp(double c, double s) {
    std::vector<double> y, w;
    fold_rule(c, s, y, w);
    return Rcpp::List::create(Rcpp::_["y"] = Rcpp::wrap(y),
                              Rcpp::_["w"] = Rcpp::wrap(w));
}

// [[Rcpp::export]]
int d7_scalar_thread_safe_probe(std::string cls) {
    return d7_scalar_thread_safe(d7_scalar_id(cls.c_str()));
}

// the classes the registry covers, for the twin tests
// [[Rcpp::export]]
Rcpp::CharacterVector d7_scalar_classes_covered() {
    Rcpp::CharacterVector out(d7_n_scalar_classes);
    for (int i = 0; i < d7_n_scalar_classes; ++i) out[i] = d7_scalar_classes[i];
    return out;
}

// [[Rcpp::init]]
void d7_register_ccallable(DllInfo* dll) {
    R_RegisterCCallable("distributions7", "d7_scalar_id",
                        (DL_FUNC) d7_scalar_id);
    R_RegisterCCallable("distributions7", "d7_score_curv",
                        (DL_FUNC) d7_score_curv);
    R_RegisterCCallable("distributions7", "d7_info_dinfo",
                        (DL_FUNC) d7_info_dinfo);
    R_RegisterCCallable("distributions7", "d7_scalar_thread_safe",
                        (DL_FUNC) d7_scalar_thread_safe);
    R_RegisterCCallable("distributions7", "d7_logpdf",
                        (DL_FUNC) d7_logpdf);
}
