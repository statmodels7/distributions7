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
bool cont_has_cdf(int id) {
    switch (id) {
    case 0: return true;
    default: return false;
    }
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
    if (code == 8 && (inner >= kWrapBase || P > 4 || !cont_has_cdf(inner)))
        return -1;
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

// the support of a continuous family at its parameters
void cont_support(int id, const double* th, double* lo, double* hi) {
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
    if (id == 28 || id == 29 || id == 32) { kk[0] = th[0]; return 1; }
    return 0;
}

// the compiled distribution function and its derivatives, by family; false
// where the family has none
bool cont_cdf(int id, double q, const double* th, bool lower, double* F) {
    switch (id) {
    case 0: *F = d7::gaussian1_cdf(q, th, lower); return true;
    default: return false;
    }
}
bool cont_cdf_grad(int id, double q, const double* th, double* g) {
    switch (id) {
    case 0: d7::gaussian1_cdf_grad(q, th, g); return true;
    default: return false;
    }
}
bool cont_cdf_hess(int id, double q, const double* th, double* h) {
    switch (id) {
    case 0: d7::gaussian1_cdf_hess(q, th, h); return true;
    default: return false;
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

// the retained mass and its first and second derivatives in every
// parameter, at one observation; inner a base family
TruncEnds trunc_cont_ends(int inner, const double* th, double lo, double up) {
    TruncEnds e;
    const int P = d7_n_params[inner], np2 = P * (P + 1) / 2;
    double slo, shi;
    cont_support(inner, th, &slo, &shi);
    const bool lo_in = lo > slo, up_in = up < shi;
    double Fa = 0.0, gl[4] = {0, 0, 0, 0}, gu[4] = {0, 0, 0, 0};
    double hl[10] = {0}, hu[10] = {0};
    e.ok = true;
    if (lo_in) e.ok = e.ok && cont_cdf(inner, lo, th, true, &Fa) &&
                   cont_cdf_grad(inner, lo, th, gl) && cont_cdf_hess(inner, lo, th, hl);
    if (up_in) e.ok = e.ok && cont_cdf_grad(inner, up, th, gu) &&
                   cont_cdf_hess(inner, up, th, hu);
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
void trunc_cont_nodes(int inner, const double* th, double lo, double up,
                      std::vector<double>& y, std::vector<double>& w) {
    double slo, shi, c, s, kk[4];
    cont_support(inner, th, &slo, &shi);
    if (!cont_center_scale(inner, th, &c, &s)) { y.clear(); w.clear(); return; }
    const int nk = cont_kinks(inner, th, kk);
    d7::trunc_rule(std::max(lo, slo), std::min(up, shi), c, s, kk, nk, y, w);
}

void trunc_cont_score_curv(int inner, int k, double y, const double* th,
                           double* out) {
    const double* tp = trunc_points(inner, th);
    const TruncEnds e = trunc_cont_ends(inner, th, tp[0], tp[1]);
    if (!e.ok) { out[0] = out[1] = R_NaN; return; }
    const int P = d7_n_params[inner];
    double sc[2];
    d7_score_curv(inner, k, y, th, sc);
    const double m = e.g[k] / e.z, M = e.h[hess_pos(P, k, k)] / e.z;
    out[0] = sc[0] - m;
    out[1] = sc[1] - M + m * m;
}

void trunc_cont_info_dinfo(int inner, int k, double y, const double* th,
                           double* out) {
    const double* tp = trunc_points(inner, th);
    const TruncEnds e = trunc_cont_ends(inner, th, tp[0], tp[1]);
    if (!e.ok) { out[0] = out[1] = R_NaN; return; }
    const int P = d7_n_params[inner];
    std::vector<double> ys, ws;
    trunc_cont_nodes(inner, th, tp[0], tp[1], ys, ws);
    long double S = 0.0L, D = 0.0L;
    for (std::size_t j = 0; j < ys.size(); ++j) {
        const double fw = std::exp(d7_logpdf(inner, ys[j], th)) * ws[j];
        if (fw == 0.0) continue;
        double sc[2];
        d7_score_curv(inner, k, ys[j], th, sc);
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
    const TruncEnds e = trunc_cont_ends(inner, th, tp[0], tp[1]);
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
Rcpp::List trunc_cont_ends_cpp(std::string cls, Rcpp::NumericMatrix theta) {
    const int id = d7_scalar_id(cls.c_str());
    if (!node_valid(id) || node_of(id).code != 8)
        Rcpp::stop("'%s' is not the route of a truncated continuous family.", cls);
    const int inner = node_of(id).inner, P = d7_n_params[inner];
    const int n = theta.nrow(), np = theta.ncol(), np2 = P * (P + 1) / 2;
    Rcpp::NumericVector z(n);
    Rcpp::NumericMatrix g(n, P), h(n, np2);
    std::vector<double> th(np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        const double* tp = trunc_points(inner, th.data());
        const TruncEnds e = trunc_cont_ends(inner, th.data(), tp[0], tp[1]);
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
    std::vector<int> idx;
    std::vector<double> yy, ww, ys, ws, th(np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        const double* tp = trunc_points(inner, th.data());
        trunc_cont_nodes(inner, th.data(), tp[0], tp[1], ys, ws);
        for (std::size_t j = 0; j < ys.size(); ++j) {
            idx.push_back(i + 1);
            yy.push_back(ys[j]);
            ww.push_back(ws[j]);
        }
    }
    return Rcpp::List::create(Rcpp::_["idx"] = Rcpp::wrap(idx),
                              Rcpp::_["y"] = Rcpp::wrap(yy),
                              Rcpp::_["w"] = Rcpp::wrap(ww));
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
