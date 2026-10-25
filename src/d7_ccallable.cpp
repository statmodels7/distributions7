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
int d7_joint_covered(int id);
void d7_joint_grad(int id, double y, const double* th, double* out);
void d7_joint_hess(int id, double y, const double* th, double* out);
void d7_joint_info(int id, double y, const double* th, double* out);
void d7_joint_dinfo(int id, double y, const double* th, double* out);
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

// the degrees of freedom of a parent whose density has the polynomial tail
// |z|^-(1+nu) log^k|z| (the skew t and student_t1, through fixed()), or NaN
// for any other parent: fold_rule() adds the anchor nodes of the skew t's
// tail correction for it
double fold_tail_nu(int inner, const double* th) {
    if (inner >= kWrapBase) {
        const WrapNode& w = node_of(inner);
        if (w.code != 1) return R_NaN;
        double full[kMaxVec];
        fixed_inner(w.inner, w.aux, 0, th, full);
        return fold_tail_nu(w.inner, full);
    }
    if (inner == 35) return th[3];
    if (inner == 23) return th[2];
    return R_NaN;
}

// The nodes and weights of folded()'s quadrature over y in (0, Inf), for a
// parent centered at c with scale s: a tanh-sinh rule on [0, |c|], whose
// nodes cluster at both ends -- at zero, where the weight f(y)/(f(y) +
// f(-y)) turns from 1/2 to 1 over a width of order s^2/|c|, and at |c|,
// where the mass is and where a parent with a kink at its center puts it --
// and the exp-sinh rule of pt_loc_scale.h scaled by s on [|c|, Inf). Both
// have step 1/32; the tanh-sinh nodes run until their distance to an end
// falls below 1e-17 of the interval. For a parent with the tail of the skew t
// at nu < 1 (`nu` finite) the nine anchor nodes of pt_skewt.h are appended on
// the positive side, at a + s z_i with the weights s Omega_i(nu), which carry
// the integral past the last exp-sinh node (|z| = 1e60), where the tail
// |z|^-(1+nu) is not negligible.
void fold_rule(double c, double s, double nu, std::vector<double>& y,
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
    if (R_FINITE(nu) && nu < 1.0) {
        double om[d7::kSkewtAnchors];
        d7::skewt_tail_omega(nu, om);
        for (int i = 0; i < d7::kSkewtAnchors; ++i) {
            y.push_back(a + s * d7::skewt_anchor(i));
            w.push_back(s * om[i]);
        }
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
    fold_rule(c, s, fold_tail_nu(inner, full), ys, ws);
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
    case 25: {
        // gengamma1: (y/a)^p is gamma(d/p), whose median is near d/p, so
        // the mass sits at a (d/p)^(1/p) and not at the scale a
        *c = th[0] * std::pow(th[1] / th[2], 1 / th[2]);
        *s = *c;
        if (R_FINITE(*c) && *c > 0) return true;
        *c = th[0]; *s = th[0];
        return true;
    }
    case 30: *c = th[0]; *s = th[0]; return true;
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
                 std::vector<double>* l1y = nullptr, double sv_max = 20.0);
int cont_map(int id, const double* th);

// ---- distribution functions by quadrature ----------------------------------
//
// For the families whose distribution function has no closed derivative in
// a shape parameter, the derivatives are integrals of the density's own:
//   F_k(q) = int_{lo}^{q} f s_k,  F_kl(q) = int_{lo}^{q} f (l_kl + s_k s_l),
// taken by the rule of pt_trunc_rule.h, or as minus the same integrals over
// [q, hi) when that side holds less mass (quad_side_upper()), so that
// neither side is a difference of two numbers near one. A location-scale family integrates
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
    case 35: {  // skewt (mu, sigma, alpha, nu): alpha-nu, exact
        if (!(k == 2 && l == 3)) return R_NaN;
        thread_local d7::SkewtConst cc;
        cc.at(th[3], 2);
        double d2[10];
        d7::skewt_obs<2>(y, th[0], th[1], th[2], th[3], cc, d2);
        return d2[8];   // alpha_nu in the order of pt_skewt_table.h
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

// The rule of the cdf quadrature in log y is cut at kCdfHeadCut, below
// which a family on (0, Inf) whose density behaves like y^(alpha - 1) is
// integrated in closed form, and its panels are at most kCdfSvMax wide in
// log y, so that the cutoff of a small-shape gamma's density (of width
// one in log y) is resolved.
const double kCdfHeadCut = 1e-290;
const double kCdfSvMax = 1.0;

// The power-law head of a family integrated in log y: on (0, y0], y0 the
// smaller of q and kCdfHeadCut, with L = log y, the log-density is
// c + (alpha - 1) L and every derivative of it in the parameters is linear
// in L, to rounding. Its mass is int_{-Inf}^{L0} e^{c + alpha u} g(u) du =
// E sum_j (-1)^j g^(j)(L0) / alpha^(j+1), E = f(y0) y0, for g a polynomial
// in u. The family's functions are read at ya = kCdfHeadCut and at 15 and
// 30 units of L below it, all normal doubles, and extended linearly to L0
// (q itself may be subnormal); the head is refused when the two slopes of
// the log-density or of a score disagree, and it is used where its mass
// is at least 1e-17 of F(q) (otherwise the rule from zero is exact and
// unchanged).
struct CdfHead {
    double ya, yb, dL, La, L0, y0, alpha, logE;
};

bool cdf_head(int id, double q, const double* th, int mask, CdfHead* h) {
    if (id >= kWrapBase || cont_map(id, th) != 1 || !(q > 0)) return false;
    const double ya = kCdfHeadCut, La = std::log(ya);
    const double yb = std::exp(La - 15), yc = std::exp(La - 30);
    const double d1 = La - std::log(yb), d2 = std::log(yb) - std::log(yc);
    const double y0 = std::min(q, ya), L0 = std::log(y0);
    const double la = d7_logpdf(id, ya, th), lb = d7_logpdf(id, yb, th);
    if (!R_FINITE(la) || !R_FINITE(lb)) return false;
    const double a1 = (la - lb) / d1;
    // the head's mass, E / alpha, against F(q): the closed distribution
    // function cannot give F(y0) where (y0/a)^p underflows
    const double logE = la + a1 * (L0 - La) + L0;
    double Fq;
    if (!closed_cdf(id, q, th, true, &Fq)) return false;
    if (!(1 + a1 > 0) || !(std::exp(logE) / (1 + a1) > 1e-17 * Fq))
        return false;
    const double lc = d7_logpdf(id, yc, th);
    if (!R_FINITE(lc)) return false;
    const double a2 = (lb - lc) / d2;
    if (std::fabs(a1 - a2) > 1e-9 * (1 + std::fabs(a1))) return false;
    const int P = d7_n_params[id];
    for (int k = 0; k < P; ++k) {
        if (!((mask >> k) & 1)) continue;
        double s0[2], s1[2], s2[2];
        d7_score_curv(id, k, ya, th, s0);
        d7_score_curv(id, k, yb, th, s1);
        d7_score_curv(id, k, yc, th, s2);
        const double b1 = (s0[0] - s1[0]) / d1, b2 = (s1[0] - s2[0]) / d2;
        if (std::fabs(b1 - b2) >
            1e-9 * (1 + std::fabs(b1) + std::fabs(s0[0]) / 15)) return false;
    }
    h->ya = ya;
    h->yb = yb;
    h->dL = d1;
    h->La = La;
    h->L0 = L0;
    h->y0 = y0;
    h->alpha = 1 + a1;
    h->logE = logE;
    return R_FINITE(h->alpha) && h->alpha > 0 && R_FINITE(h->logE);
}

// A cut of the cdf rule where the density of a family on (0, Inf) begins
// its fast decay, at the scale theta of exp(-(y/theta)^p): from there the
// decay has width 1/p in log y, while the panels that double away from the
// center reach it several units wide when the center lies far below it (a
// small shape). The cut is placed when it lies more than one unit of log y
// above the center; the panels then restart from it.
int cdf_edge_cut(int id, const double* th, double* kk) {
    double e;
    switch (id) {
    case 1: e = th[1] * th[0]; break;           // gamma1: phi mu
    case 17: e = th[1] / th[0]; break;          // gamma2: sigma2 / mu
    case 9: e = 2.0; break;                     // chisq
    case 25: e = th[0]; break;                  // gengamma1: a
    case 26:                                    // gengamma2: a from the mean
        e = std::exp(std::log(th[0]) + R::lgammafn(th[1] / th[2]) -
                     R::lgammafn((th[1] + 1) / th[2]));
        break;
    default: return 0;
    }
    double c, s;
    cont_center_scale(id, th, &c, &s);
    if (!R_FINITE(e) || !(e > 0) || !(c > 0) || !(std::log(e / c) > 1))
        return 0;
    kk[0] = e;
    return 1;
}

// the head's contributions to the sums of quad_cdf_sums(); a value at L0
// is the value at La extended along its slope
void cdf_head_sums(int id, const double* th, int mask, bool want_h,
                   const CdfHead& hd, long double* S0, long double* S1,
                   long double* S2) {
    const int P = d7_n_params[id];
    const double a = hd.alpha, E = std::exp(hd.logE), dl = hd.L0 - hd.La;
    // g(L0) / a - g'(L0) / a^2 + g''(L0) / a^3
    auto mom = [a](double g, double g1, double g2) {
        return g / a - g1 / (a * a) + g2 / (a * a * a);
    };
    *S0 += E / a;
    double g0[4] = {0, 0, 0, 0}, l0[4] = {0, 0, 0, 0}, bg[4] = {0, 0, 0, 0},
        bl[4] = {0, 0, 0, 0};
    for (int k = 0; k < P; ++k) {
        if (!((mask >> k) & 1)) continue;
        double s0[2], s1[2];
        d7_score_curv(id, k, hd.ya, th, s0);
        d7_score_curv(id, k, hd.yb, th, s1);
        bg[k] = (s0[0] - s1[0]) / hd.dL;   // per unit of L
        bl[k] = (s0[1] - s1[1]) / hd.dL;
        g0[k] = s0[0] + bg[k] * dl;
        l0[k] = s0[1] + bl[k] * dl;
        S1[k] += E * mom(g0[k], bg[k], 0.0);
    }
    if (!want_h) return;
    for (int k = 0; k < P; ++k) {
        if (!((mask >> k) & 1)) continue;
        S2[k * 4 + k] += E * mom(l0[k] + g0[k] * g0[k],
                                 bl[k] + 2 * bg[k] * g0[k],
                                 2 * bg[k] * bg[k]);
        for (int m = k + 1; m < P; ++m) {
            if (!((mask >> m) & 1)) continue;
            const double ha = hess_off(id, k, m, hd.ya, th);
            const double bh = (ha - hess_off(id, k, m, hd.yb, th)) / hd.dL;
            S2[k * 4 + m] += E * mom(ha + bh * dl + g0[k] * g0[m],
                                     bh + bg[k] * g0[m] + g0[k] * bg[m],
                                     2 * bg[k] * bg[m]);
        }
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
    const bool logs = logs_family(id);
    CdfHead hd;
    const bool head = !upper && !logs && cdf_head(id, q, th, mask, &hd);
    const double a = upper ? std::max(q, slo) : (head ? hd.y0 : slo);
    const double b = upper ? shi : std::min(q, shi);
    std::vector<double> lys, l1ys;
    double kk[1];
    const int nk = cdf_edge_cut(id, th, kk);
    if (b > a)
        mapped_rule(id, th, a, b, kk, nk, ys, ws, logs ? &lys : nullptr,
                    logs ? &l1ys : nullptr, kCdfSvMax);
    long double S0 = 0.0L, S1[4] = {0, 0, 0, 0}, S2[16] = {0};
    if (head) cdf_head_sums(id, th, mask, want_h, hd, &S0, S1, S2);
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

// the side of q to integrate over: the one whose mass is at most one half,
// so that a derivative small beside one is not the sum of terms of order
// one (the center of a skewed family can lie far in its lower tail: the
// generalized gamma's is its scale a). The mass is R's distribution
// function where the family has one, else the integral over the side of
// the center, flipped when it exceeds one half. On (0, 1), when the second
// shape is below one, the lower side is used always: the mass within one
// unit in the last place of 1 is not representable in y (it is eps^b of
// the total)
// the first choice of side, and whether it is final (`fixed`); when it is
// not, the caller flips it if the integral over it exceeds one half
bool quad_side_first(int id, double q, const double* th, bool* fixed) {
    *fixed = true;
    if (id == 4 && (1 - th[0]) * th[1] < 1) return false;
    if (id == 31 && th[1] < 1) return false;
    double F;
    if (closed_cdf(id, q, th, true, &F)) return F > 0.5;
    *fixed = false;
    double c, s;
    cont_center_scale(id, th, &c, &s);
    return q > c;
}

bool quad_side_upper(int id, double q, const double* th) {
    bool fixed;
    const bool up = quad_side_first(id, q, th, &fixed);
    if (fixed) return up;
    double I0, Ik[4], Ikl[16];
    quad_cdf_sums(id, q, th, 0, up, false, &I0, Ik, Ikl);
    return (I0 > 0.5) ? !up : up;
}

double quad_cdf(int id, double q, const double* th, bool lower) {
    double F;
    if (closed_cdf(id, q, th, lower, &F)) return F;
    bool fixed;
    bool up = quad_side_first(id, q, th, &fixed);
    double I0, Ik[4], Ikl[16];
    quad_cdf_sums(id, q, th, 0, up, false, &I0, Ik, Ikl);
    if (!fixed && I0 > 0.5) {
        up = !up;
        quad_cdf_sums(id, q, th, 0, up, false, &I0, Ik, Ikl);
    }
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
    bool fixed;
    bool up = quad_side_first(id, q, th, &fixed);
    double I0, Ik[4], Ikl[16];
    quad_cdf_sums(id, q, th, mask, up, want_h, &I0, Ik, Ikl);
    if (!fixed && I0 > 0.5) {
        up = !up;
        quad_cdf_sums(id, q, th, mask, up, want_h, &I0, Ik, Ikl);
    }
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
                 std::vector<double>* l1y, double sv_max) {
    if (ly != nullptr) ly->clear();
    if (l1y != nullptr) l1y->clear();
    double c, s;
    cont_center_scale(id, th, &c, &s);
    const int mp = cont_map(id, th);
    if (mp == 0) {
        d7::trunc_rule(a, b, c, s, kk, nk, y, w, d7::kTsStep, d7::kDoublings,
                       fold_tail_nu(id, th));
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
    sv = std::min(mp == 1 ? sv_max : 20.0, std::max(0.05, sv));
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
        const double lf0 = (code == 4) ? R_NegInf : d7_logpdf(inner, 0.0, full);
        const double f0 = std::exp(lf0), q = -std::expm1(lf0);
        if (k == P) {
            if (code == 2) {
                out[0] = d7::zi_score_zi(y, z, f0, q);
                out[1] = d7::zi_curv_zi(y, z, f0, q);
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
            out[0] = d7::zad_score_parent(y, f0, q, sh0[0], gh[0]);
            out[1] = d7::zad_curv_parent(y, f0, q, sh0[0], sh0[1], gh[1]);
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
        const double lf0 = (code == 4) ? R_NegInf : d7_logpdf(inner, 0.0, full);
        const double f0 = std::exp(lf0), q = -std::expm1(lf0);
        if (k == P) {
            if (code == 2) {
                out[0] = d7::zi_info_zi(z, f0, q);
                out[1] = d7::zi_dinfo_zi(z, f0, q);
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
            out[0] = d7::zad_info_parent(z, f0, q, sh0[0], sh0[1], ed[0]);
            out[1] = d7::zad_dinfo_parent(z, f0, q, sh0[0], sh0[1], ed[0], ed[1]);
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
        return d7::zad_logpdf(y, z, -std::expm1(d7_logpdf(inner, 0.0, full)), lf);
    }
    return R_NaN;
}

// ---- the joint entries of the wrappers --------------------------------------
//
// A wrapper's joint entries are composed from its inner family's joint
// entries, which for a wrapped inner family are a node's in turn, with the
// layout of d7_joint_grad() and the others; the expressions are those of the
// wrappers' R methods, extended from the (k, k) entries above to every pair.
// fixed() reindexes, transformation() reads the inner entries at the inverse
// of the transform, the zero wrappers are closed forms in the parent's
// entries at the observation and at zero, truncated() sums or integrates the
// parent's gradient and Hessian over the retained support, and folded()
// integrates the folded ones by fold_rule(). A node is covered to the order
// its inner family is (wrap_joint_covered()), except that the information of
// a truncated or folded family reads only the inner gradient and Hessian, and
// that folded()'s information needs a center and a scale. kJointMax bounds
// the number of parameters of the node and of its inner family, which sizes
// the buffers of the inner family's entries.

const int kJointMax = 8;
const int kJointMax2 = kJointMax * kJointMax;
const int kJointMax3 = kJointMax2 * kJointMax;

// whether center_scale() has a center and a scale for `inner`, which folded()'s
// quadrature needs: a location-scale family on the real line, through fixed()
bool fold_cs_family(int inner) {
    while (inner >= kWrapBase) {
        if (node_of(inner).code != 1) return false;
        inner = node_of(inner).inner;
    }
    switch (inner) {
    case 0: case 10: case 11: case 12: case 13: case 23: case 24: case 27:
    case 28: case 29: case 32: case 34: case 35: case 36: case 37: case 38:
        return true;
    default: return false;
    }
}

int wrap_joint_covered(int id) {
    const WrapNode& w = node_of(id);
    if (w.np > kJointMax || id_np(w.inner) > kJointMax) return 0;
    const int ci = std::max(d7_joint_covered(w.inner), 0);
    switch (w.code) {
    case 1: case 2: case 3: case 4: case 6: return ci;
    case 5: return ci < 2 ? 0 : (fold_cs_family(w.inner) ? 4 : 2);
    case 7: case 8: return ci < 2 ? 0 : 4;
    default: return 0;
    }
}

// fixed(): the positions of the free parameters in the inner vector, and
// their number
int fixed_free(int inner, int mask, int* at) {
    const int P = id_np(inner);
    int n = 0;
    for (int j = 0; j < P; ++j)
        if (!((mask >> j) & 1)) at[n++] = j;
    return n;
}

// the zero wrappers: the parent's score and Hessian at zero
void zero_parts0(int inner, const double* full, double* s0, double* H0) {
    d7_joint_grad(inner, 0.0, full, s0);
    d7_joint_hess(inner, 0.0, full, H0);
}

// folded(): the folded gradient and, where h is not null, Hessian at y, from
// the parent's joint entries at y and -y, the Hessian in the form of
// fold_curv(); a preimage whose weight is zero contributes zero, as in
// fold_ratio(). Returns the weights.
d7::FoldW fold_joint_gh(int inner, double y, const double* th, int P,
                        double* g, double* h) {
    const d7::FoldW W = d7::fold_w(d7_logpdf(inner, y, th),
                                   d7_logpdf(inner, -y, th));
    const double w = W.w, wm = 1 - w;
    double gp[kJointMax], gm[kJointMax];
    d7_joint_grad(inner, y, th, gp);
    d7_joint_grad(inner, -y, th, gm);
    for (int a = 0; a < P; ++a)
        g[a] = (w == 0 ? 0.0 : w * gp[a]) + (wm == 0 ? 0.0 : wm * gm[a]);
    if (h == nullptr) return W;
    double hp[kJointMax2], hm[kJointMax2];
    d7_joint_hess(inner, y, th, hp);
    d7_joint_hess(inner, -y, th, hm);
    for (int b = 0; b < P; ++b) {
        for (int a = 0; a < P; ++a) {
            const double rp = (w == 0) ? 0.0 : w * hp[a + P * b];
            const double rm = (wm == 0) ? 0.0 : wm * hm[a + P * b];
            const double c = (w == 0 || wm == 0) ? 0.0
                : w * wm * (gp[a] - gm[a]) * (gp[b] - gm[b]);
            h[a + P * b] = rp + rm + c;
        }
    }
    return W;
}

// folded(): the expected Hessian (into e, d null) or its derivative (into d),
// as fold_info_dinfo() takes the (k, k) ones: -E[l_a l_b] and
// -E[l_ac l_b + l_a l_bc + l_a l_b l_c] by fold_rule(), with the boundary
// term of a kink at the center added to the derivative in the center
void fold_joint_info(int inner, const double* th, int P, double* e, double* d) {
    const int PP = P * P, PPP = PP * P;
    double c, s;
    if (!center_scale(inner, th, &c, &s)) {
        if (d == nullptr) for (int j = 0; j < PP; ++j) e[j] = R_NaN;
        else for (int j = 0; j < PPP; ++j) d[j] = R_NaN;
        return;
    }
    std::vector<double> ys, ws;
    fold_rule(c, s, fold_tail_nu(inner, th), ys, ws);
    long double S[kJointMax2], D[kJointMax3];
    for (int j = 0; j < PP; ++j) S[j] = 0.0L;
    if (d != nullptr) for (int j = 0; j < PPP; ++j) D[j] = 0.0L;
    double g[kJointMax], h[kJointMax2];
    for (std::size_t j = 0; j < ys.size(); ++j) {
        const double y = ys[j];
        const d7::FoldW Wt = d7::fold_w(d7_logpdf(inner, y, th),
                                        d7_logpdf(inner, -y, th));
        const double fw = Wt.L * ws[j];
        if (fw == 0.0) continue;
        fold_joint_gh(inner, y, th, P, g, d == nullptr ? nullptr : h);
        if (d == nullptr) {
            for (int b = 0; b < P; ++b)
                for (int a = 0; a < P; ++a) S[a + P * b] += (g[a] * g[b]) * fw;
            continue;
        }
        for (int k = 0; k < P; ++k)
            for (int b = 0; b < P; ++b)
                for (int a = 0; a < P; ++a)
                    D[a + P * b + PP * k] += (h[a + P * k] * g[b] + g[a] * h[b + P * k] +
                                              g[a] * g[b] * g[k]) * fw;
    }
    if (d == nullptr) {
        for (int j = 0; j < PP; ++j) {
            const double r = (double) S[j];
            e[j] = R_FINITE(r) ? -r : NA_REAL;
        }
        return;
    }
    // the parameter at a kink of the parent's center (laplace, laplace2,
    // enet), whose folded score jumps at y* = |c|: d_c E[l_a l_b] carries
    // L(y*) [G(y*-) - G(y*+)] sign(c), G = l_a l_b, the one-sided scores
    // read four units in the last place either side of y*
    int kp = -1;
    for (int k = 0; k < P && kp < 0; ++k)
        if (kink_center(inner, k)) kp = k;
    double lo[kJointMax], hi[kJointMax], L0 = 0.0;
    const bool kink = kp >= 0 && c != 0;
    if (kink) {
        const double a0 = std::fabs(c);
        L0 = d7::fold_w(d7_logpdf(inner, a0, th), d7_logpdf(inner, -a0, th)).L;
        fold_joint_gh(inner, a0 * (1 + (-4) * DBL_EPSILON), th, P, lo, nullptr);
        fold_joint_gh(inner, a0 * (1 + 4 * DBL_EPSILON), th, P, hi, nullptr);
    }
    for (int k = 0; k < P; ++k) {
        for (int b = 0; b < P; ++b) {
            for (int a = 0; a < P; ++a) {
                double r = (double) D[a + P * b + PP * k];
                if (kink && k == kp)
                    r = r + L0 * (lo[a] * lo[b] - hi[a] * hi[b]) * (c > 0 ? 1.0 : -1.0);
                d[a + P * b + PP * k] = R_FINITE(r) ? -r : NA_REAL;
            }
        }
    }
}

// truncated(): the retained mass Z and its derivatives Z_a and Z_ab, and the
// information's sums S_ab = sum f s_a s_b and
// D_abc = sum f (s_a s_b s_c + l_ac s_b + s_a l_bc), over every parameter
struct TruncJ {
    double z, zi[kJointMax], zij[kJointMax2], s[kJointMax2], d[kJointMax3];
};

// a discrete parent: trunc_disc_sums() over every parameter, to the level
// asked (1 the mass and its first derivatives, 2 the second, 3 S, 4 D); the
// levels of the information run over the points that trunc_disc_rule() gives
// the information
bool trunc_disc_joint(int inner, const double* th, int level, TruncJ& t) {
    const int P = id_np(inner), PP = P * P;
    const double* tp = trunc_points(inner, th);
    const TruncRule r = trunc_disc_rule(inner, th, tp[0], tp[1], level >= 3);
    if (ISNAN(r.y1)) return false;
    long double F = 0.0L, G[kJointMax], H[kJointMax2], S[kJointMax2], D[kJointMax3];
    for (int a = 0; a < P; ++a) G[a] = 0.0L;
    for (int j = 0; j < PP; ++j) H[j] = S[j] = 0.0L;
    for (int j = 0; j < PP * P; ++j) D[j] = 0.0L;
    double g[kJointMax], h[kJointMax2];
    for (double x = r.y0; x <= r.y1; ++x) {
        const double f = std::exp(d7_logpdf(inner, x, th));
        d7_joint_grad(inner, x, th, g);
        if (level >= 2) d7_joint_hess(inner, x, th, h);
        F += f;
        for (int a = 0; a < P; ++a) G[a] += f * g[a];
        if (level >= 2)
            for (int b = 0; b < P; ++b)
                for (int a = 0; a < P; ++a)
                    H[a + P * b] += f * (h[a + P * b] + g[a] * g[b]);
        if (level >= 3)
            for (int b = 0; b < P; ++b)
                for (int a = 0; a < P; ++a) S[a + P * b] += f * (g[a] * g[b]);
        if (level >= 4)
            for (int c = 0; c < P; ++c)
                for (int b = 0; b < P; ++b)
                    for (int a = 0; a < P; ++a)
                        D[a + P * b + PP * c] += f * (g[a] * g[b] * g[c] +
                                                      h[a + P * c] * g[b] +
                                                      g[a] * h[b + P * c]);
    }
    const double sg = (r.branch == 1) ? -1.0 : 1.0;
    t.z = (r.branch == 1) ? 1 - (double) F : (double) F;
    for (int a = 0; a < P; ++a) t.zi[a] = sg * (double) G[a];
    for (int j = 0; j < PP; ++j) {
        t.zij[j] = sg * (double) H[j];
        t.s[j] = (double) S[j];
    }
    for (int j = 0; j < PP * P; ++j) t.d[j] = (double) D[j];
    return true;
}

// a continuous parent: the mass and its derivatives from trunc_cont_ends() at
// the level asked (as for trunc_disc_joint()), and S and D by the rule of
// trunc_cont_nodes(), a beta reading log y and log(1 - y)
bool trunc_cont_joint(int inner, const double* th, int level, TruncJ& t) {
    const int P = id_np(inner), PP = P * P;
    const double* tp = trunc_points(inner, th);
    const TruncEnds e = trunc_cont_ends(inner, th, tp[0], tp[1], std::min(level, 2));
    if (!e.ok) return false;
    t.z = e.z;
    for (int a = 0; a < P; ++a) t.zi[a] = e.g[a];
    if (level >= 2)
        for (int b = 0; b < P; ++b)
            for (int a = 0; a < P; ++a) t.zij[a + P * b] = e.h[hess_pos(P, a, b)];
    if (level < 3) return true;
    const bool logs = logs_family(inner);
    std::vector<double> ys, ws, lys, l1ys;
    trunc_cont_nodes(inner, th, tp[0], tp[1], ys, ws, logs ? &lys : nullptr,
                     logs ? &l1ys : nullptr);
    if (ys.empty()) return false;
    long double S[kJointMax2], D[kJointMax3];
    for (int j = 0; j < PP; ++j) S[j] = 0.0L;
    for (int j = 0; j < PP * P; ++j) D[j] = 0.0L;
    double g[kJointMax], h[kJointMax2];
    for (std::size_t j = 0; j < ys.size(); ++j) {
        double lf;
        if (logs) {
            double gg[2], hh[3];
            lf = beta_logs_parts(inner, lys[j], l1ys[j], th, gg, hh);
            for (int a = 0; a < P; ++a) g[a] = gg[a];
            for (int b = 0; b < P; ++b)
                for (int a = 0; a < P; ++a) h[a + P * b] = hh[hess_pos(P, a, b)];
        } else {
            lf = d7_logpdf(inner, ys[j], th);
        }
        const double fw = std::exp(lf) * ws[j];
        if (fw == 0.0) continue;
        if (!logs) {
            d7_joint_grad(inner, ys[j], th, g);
            if (level >= 4) d7_joint_hess(inner, ys[j], th, h);
        }
        for (int b = 0; b < P; ++b)
            for (int a = 0; a < P; ++a) S[a + P * b] += fw * (g[a] * g[b]);
        if (level >= 4)
            for (int c = 0; c < P; ++c)
                for (int b = 0; b < P; ++b)
                    for (int a = 0; a < P; ++a)
                        D[a + P * b + PP * c] += fw * (g[a] * g[b] * g[c] +
                                                       h[a + P * c] * g[b] +
                                                       g[a] * h[b + P * c]);
    }
    for (int j = 0; j < PP; ++j) t.s[j] = (double) S[j];
    for (int j = 0; j < PP * P; ++j) t.d[j] = (double) D[j];
    return true;
}

bool trunc_joint(int code, int inner, const double* th, int level, TruncJ& t) {
    return (code == 7) ? trunc_disc_joint(inner, th, level, t)
                       : trunc_cont_joint(inner, th, level, t);
}

void wrap_joint_grad(int id, double y, const double* th, double* out) {
    const WrapNode& w = node_of(id);
    const int inner = w.inner, Pi = id_np(inner);
    double full[kMaxVec], g[kJointMax];
    switch (w.code) {
    case 1: {
        int at[kJointMax];
        fixed_inner(inner, w.aux, 0, th, full);
        const int n = fixed_free(inner, w.aux, at);
        d7_joint_grad(inner, y, full, g);
        for (int a = 0; a < n; ++a) out[a] = g[at[a]];
        return;
    }
    case 2: {
        const double z = zero_inner(inner, th, full);
        const double lf0 = d7_logpdf(inner, 0.0, full);
        const double f0 = std::exp(lf0);
        d7_joint_grad(inner, y, full, g);
        for (int a = 0; a < Pi; ++a) out[a] = d7::zi_score_parent(y, z, f0, g[a]);
        out[Pi] = d7::zi_score_zi(y, z, f0, -std::expm1(lf0));
        return;
    }
    case 3: case 4: {
        const double z = zero_inner(inner, th, full);
        if (y == 0) {
            for (int a = 0; a < Pi; ++a) out[a] = 0.0;
        } else {
            d7_joint_grad(inner, y, full, g);
            if (w.code == 3) {
                double s0[kJointMax];
                const double lf0 = d7_logpdf(inner, 0.0, full);
                const double f0 = std::exp(lf0), q = -std::expm1(lf0);
                d7_joint_grad(inner, 0.0, full, s0);
                for (int a = 0; a < Pi; ++a)
                    out[a] = d7::zad_score_parent(y, f0, q, s0[a], g[a]);
            } else {
                for (int a = 0; a < Pi; ++a) out[a] = g[a];
            }
        }
        out[Pi] = d7::za_score_za(y, z);
        return;
    }
    case 5: {
        fold_joint_gh(inner, y, th, Pi, out, nullptr);
        return;
    }
    case 6: {
        d7_joint_grad(inner, d7::transform_inv(w.aux, y, transform_par(inner, th)),
                      th, out);
        return;
    }
    case 7: case 8: {
        TruncJ t;
        if (!trunc_joint(w.code, inner, th, 1, t)) {
            for (int a = 0; a < Pi; ++a) out[a] = R_NaN;
            return;
        }
        d7_joint_grad(inner, y, th, g);
        for (int a = 0; a < Pi; ++a) out[a] = g[a] - t.zi[a] / t.z;
        return;
    }
    }
}

void wrap_joint_hess(int id, double y, const double* th, double* out) {
    const WrapNode& w = node_of(id);
    const int inner = w.inner, Pi = id_np(inner);
    double full[kMaxVec], h[kJointMax2];
    switch (w.code) {
    case 1: {
        int at[kJointMax];
        fixed_inner(inner, w.aux, 0, th, full);
        const int n = fixed_free(inner, w.aux, at);
        d7_joint_hess(inner, y, full, h);
        for (int b = 0; b < n; ++b)
            for (int a = 0; a < n; ++a) out[a + n * b] = h[at[a] + Pi * at[b]];
        return;
    }
    case 2: {
        const int P = Pi + 1;
        const double z = zero_inner(inner, th, full);
        const double lf0 = d7_logpdf(inner, 0.0, full);
        const double f0 = std::exp(lf0);
        d7_joint_hess(inner, y, full, h);
        if (y != 0) {
            for (int b = 0; b < Pi; ++b) {
                for (int a = 0; a < Pi; ++a) out[a + P * b] = h[a + Pi * b];
                out[b + P * Pi] = out[Pi + P * b] = 0.0;
            }
        } else {
            double s0[kJointMax];
            d7_joint_grad(inner, 0.0, full, s0);
            for (int b = 0; b < Pi; ++b) {
                for (int a = 0; a < Pi; ++a)
                    out[a + P * b] = d7::zi_hess_pair0(z, f0, s0[a], s0[b],
                                                       h[a + Pi * b]);
                out[b + P * Pi] = out[Pi + P * b] = d7::zi_hess_cross0(z, f0, s0[b]);
            }
        }
        out[Pi + P * Pi] = d7::zi_curv_zi(y, z, f0, -std::expm1(lf0));
        return;
    }
    case 3: case 4: {
        const int P = Pi + 1;
        const double z = zero_inner(inner, th, full);
        for (int j = 0; j < P * P; ++j) out[j] = 0.0;
        if (y != 0) {
            d7_joint_hess(inner, y, full, h);
            if (w.code == 3) {
                double s0[kJointMax], H0[kJointMax2];
                const double lf0 = d7_logpdf(inner, 0.0, full);
                const double f0 = std::exp(lf0), q = -std::expm1(lf0);
                zero_parts0(inner, full, s0, H0);
                for (int b = 0; b < Pi; ++b)
                    for (int a = 0; a < Pi; ++a)
                        out[a + P * b] = h[a + Pi * b] +
                            d7::zad_hess_correction_pair(f0, q, s0[a], s0[b],
                                                         H0[a + Pi * b]);
            } else {
                for (int b = 0; b < Pi; ++b)
                    for (int a = 0; a < Pi; ++a) out[a + P * b] = h[a + Pi * b];
            }
        }
        out[Pi + P * Pi] = d7::za_curv_za(y, z);
        return;
    }
    case 5: {
        double g[kJointMax];
        fold_joint_gh(inner, y, th, Pi, g, out);
        return;
    }
    case 6: {
        d7_joint_hess(inner, d7::transform_inv(w.aux, y, transform_par(inner, th)),
                      th, out);
        return;
    }
    case 7: case 8: {
        TruncJ t;
        if (!trunc_joint(w.code, inner, th, 2, t)) {
            for (int j = 0; j < Pi * Pi; ++j) out[j] = R_NaN;
            return;
        }
        d7_joint_hess(inner, y, th, h);
        for (int b = 0; b < Pi; ++b) {
            for (int a = 0; a < Pi; ++a) {
                const double ma = t.zi[a] / t.z, mb = t.zi[b] / t.z;
                out[a + Pi * b] = h[a + Pi * b] - t.zij[a + Pi * b] / t.z + ma * mb;
            }
        }
        return;
    }
    }
}

void wrap_joint_info(int id, double y, const double* th, double* out) {
    const WrapNode& w = node_of(id);
    const int inner = w.inner, Pi = id_np(inner);
    double full[kMaxVec], e[kJointMax2];
    switch (w.code) {
    case 1: {
        int at[kJointMax];
        fixed_inner(inner, w.aux, 0, th, full);
        const int n = fixed_free(inner, w.aux, at);
        d7_joint_info(inner, y, full, e);
        for (int b = 0; b < n; ++b)
            for (int a = 0; a < n; ++a) out[a + n * b] = e[at[a] + Pi * at[b]];
        return;
    }
    case 2: case 3: case 4: {
        const int P = Pi + 1;
        const double z = zero_inner(inner, th, full);
        d7_joint_info(inner, y, full, e);
        if (w.code == 4) {
            for (int b = 0; b < Pi; ++b) {
                for (int a = 0; a < Pi; ++a)
                    out[a + P * b] = d7::zac_info_parent(z, e[a + Pi * b]);
                out[b + P * Pi] = out[Pi + P * b] = 0.0;
            }
            out[Pi + P * Pi] = d7::za_info_za(z);
            return;
        }
        const double lf0 = d7_logpdf(inner, 0.0, full);
        const double f0 = std::exp(lf0), q = -std::expm1(lf0);
        double s0[kJointMax], H0[kJointMax2];
        zero_parts0(inner, full, s0, H0);
        for (int b = 0; b < Pi; ++b) {
            for (int a = 0; a < Pi; ++a) {
                out[a + P * b] = (w.code == 2)
                    ? d7::zi_info_pair(z, f0, s0[a], s0[b], H0[a + Pi * b], e[a + Pi * b])
                    : d7::zad_info_pair(z, f0, q, s0[a], s0[b], H0[a + Pi * b],
                                        e[a + Pi * b]);
            }
            out[b + P * Pi] = out[Pi + P * b] =
                (w.code == 2) ? d7::zi_info_cross(z, f0, s0[b]) : 0.0;
        }
        out[Pi + P * Pi] = (w.code == 2) ? d7::zi_info_zi(z, f0, q) : d7::za_info_za(z);
        return;
    }
    case 5: {
        fold_joint_info(inner, th, Pi, out, nullptr);
        return;
    }
    case 6: {
        d7_joint_info(inner, d7::transform_inv(w.aux, y, transform_par(inner, th)),
                      th, out);
        return;
    }
    case 7: case 8: {
        TruncJ t;
        if (!trunc_joint(w.code, inner, th, 3, t)) {
            for (int j = 0; j < Pi * Pi; ++j) out[j] = R_NaN;
            return;
        }
        for (int b = 0; b < Pi; ++b) {
            for (int a = 0; a < Pi; ++a) {
                const double ma = t.zi[a] / t.z, mb = t.zi[b] / t.z;
                out[a + Pi * b] = -(t.s[a + Pi * b] / t.z - ma * mb);
            }
        }
        return;
    }
    }
}

void wrap_joint_dinfo(int id, double y, const double* th, double* out) {
    const WrapNode& w = node_of(id);
    const int inner = w.inner, Pi = id_np(inner);
    double full[kMaxVec], d[kJointMax3];
    switch (w.code) {
    case 1: {
        int at[kJointMax];
        fixed_inner(inner, w.aux, 0, th, full);
        const int n = fixed_free(inner, w.aux, at);
        d7_joint_dinfo(inner, y, full, d);
        for (int c = 0; c < n; ++c)
            for (int b = 0; b < n; ++b)
                for (int a = 0; a < n; ++a)
                    out[a + n * b + n * n * c] =
                        d[at[a] + Pi * at[b] + Pi * Pi * at[c]];
        return;
    }
    case 2: case 3: case 4: {
        const int P = Pi + 1, m = Pi;
        const double z = zero_inner(inner, th, full);
        double e[kJointMax2], s0[kJointMax], H0[kJointMax2];
        d7_joint_info(inner, y, full, e);
        d7_joint_dinfo(inner, y, full, d);
        const double lf0 = (w.code == 4) ? R_NegInf : d7_logpdf(inner, 0.0, full);
        const double f0 = std::exp(lf0), q = -std::expm1(lf0);
        if (w.code != 4) zero_parts0(inner, full, s0, H0);
        for (int c = 0; c < P; ++c) {
            for (int b = 0; b < P; ++b) {
                for (int a = 0; a < P; ++a) {
                    const int i = std::min(a, b), j = std::max(a, b);
                    double v;
                    if (w.code == 2) {
                        if (j < m) {
                            v = (c < m)
                                ? d7::zi_dinfo_ppp(z, f0, s0[i], s0[j], s0[c],
                                                   H0[i + Pi * c], H0[j + Pi * c],
                                                   d[i + Pi * j + Pi * Pi * c])
                                : d7::zi_dinfo_ppz(z, f0, q, s0[i], s0[j], e[i + Pi * j]);
                        } else if (i < m) {
                            v = (c < m) ? d7::zi_dinfo_pzp(z, f0, s0[i], s0[c], H0[i + Pi * c])
                                        : d7::zi_dinfo_pzz(z, f0, q, s0[i]);
                        } else {
                            v = (c < m) ? d7::zi_dinfo_zzp(z, f0, q, s0[c])
                                        : d7::zi_dinfo_zi(z, f0, q);
                        }
                    } else if (j < m) {
                        if (c < m) {
                            v = (w.code == 3)
                                ? d7::zad_dinfo_ppp(z, f0, q, s0[i], s0[j], s0[c],
                                                    H0[i + Pi * c], H0[j + Pi * c],
                                                    e[i + Pi * j],
                                                    d[i + Pi * j + Pi * Pi * c])
                                : d7::zac_dinfo_parent(z, d[i + Pi * j + Pi * Pi * c]);
                        } else {
                            v = (w.code == 3) ? d7::zad_dinfo_ppz(f0, q, s0[i], s0[j], e[i + Pi * j])
                                              : -e[i + Pi * j];
                        }
                    } else if (i < m || c < m) {
                        v = 0.0;
                    } else {
                        v = d7::za_dinfo_za(z);
                    }
                    out[a + P * b + P * P * c] = v;
                }
            }
        }
        return;
    }
    case 5: {
        fold_joint_info(inner, th, Pi, nullptr, out);
        return;
    }
    case 6: {
        d7_joint_dinfo(inner, d7::transform_inv(w.aux, y, transform_par(inner, th)),
                       th, out);
        return;
    }
    case 7: case 8: {
        const int PP = Pi * Pi;
        TruncJ t;
        if (!trunc_joint(w.code, inner, th, 4, t)) {
            for (int j = 0; j < PP * Pi; ++j) out[j] = R_NaN;
            return;
        }
        // d_c of -(S_ab / Z - m_a m_b), with d_c m_a = M_ac - m_a m_c
        double mm[kJointMax], dm[kJointMax2];
        for (int a = 0; a < Pi; ++a) mm[a] = t.zi[a] / t.z;
        for (int c = 0; c < Pi; ++c)
            for (int a = 0; a < Pi; ++a)
                dm[a + Pi * c] = t.zij[a + Pi * c] / t.z - mm[a] * mm[c];
        for (int c = 0; c < Pi; ++c) {
            for (int b = 0; b < Pi; ++b) {
                for (int a = 0; a < Pi; ++a) {
                    out[a + Pi * b + PP * c] =
                        -(t.d[a + Pi * b + PP * c] / t.z -
                          t.s[a + Pi * b] * t.zi[c] / (t.z * t.z) -
                          (dm[a + Pi * c] * mm[b] + mm[a] * dm[b + Pi * c]));
                }
            }
        }
        return;
    }
    }
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

// The joint registry: every component at one observation, one entry per
// order, for a consumer that reads the derivatives in several parameters at
// once (a joint score-driven filter). th as for d7_score_curv(); the
// Hessian and the expected information fill a full P x P matrix and the
// derivative of the expected information a P x P x P array, column-major,
// out[a + P b + P^2 c] being the derivative of E[l_ab] in parameter c, all on
// the parameter scale. A family is covered once its components are written
// out in its header, and a family with one parameter through its scalar
// entries. d7_joint_covered() answers how many of the four entries the family
// has, in order: 4 where all four are exact, 2 where only the gradient and
// the Hessian are (the expected information being approximated for that
// family, which serves an unscaled filter only), 0 for the other families,
// and -1 for an unknown id. A wrapper answers its inner family's coverage
// where its composition is written (wrap_joint_covered()) and 0 elsewhere.
// An entry beyond the family's coverage fills `out` with NaN.
int d7_joint_covered(int id) {
    if (id >= kWrapBase) return node_valid(id) ? wrap_joint_covered(id) : -1;
    if (id < 0 || id >= d7_n_scalar_classes) return -1;
    switch (id) {
    case 0: case 1: case 3: case 4: case 10: case 11: case 12: case 13:
    case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 21:
    case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29:
    case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37:
    case 38: case 39: case 40: case 41: case 42:
        return 4;
    default: return d7_n_params[id] == 1 ? 4 : 0;
    }
}

static void fill_nan(double* out, int m) {
    for (int j = 0; j < m; ++j) out[j] = R_NaN;
}

// the number of parameters of a family that the joint registry knows, which
// sizes the NaN written for an entry beyond its coverage; 0 otherwise
static int joint_np(int id) {
    return d7_joint_covered(id) > 0 ? id_np(id) : 0;
}

// a family with one parameter, whose joint entries are its scalar ones
static bool joint_scalar(int id) {
    return id >= 0 && id < d7_n_scalar_classes && d7_n_params[id] == 1;
}

// the wrappers and the families with one parameter, ahead of the switch on
// the family: order 0 to 3 for the gradient, the Hessian, the expected
// Hessian and its derivative; true when `out` is written
static bool joint_route(int id, int order, double y, const double* th,
                        double* out) {
    if (id >= kWrapBase) {
        const int P = joint_np(id);
        if (d7_joint_covered(id) < (order <= 1 ? 2 : 4)) {
            fill_nan(out, order == 0 ? P : order == 3 ? P * P * P : P * P);
        } else if (order == 0) {
            wrap_joint_grad(id, y, th, out);
        } else if (order == 1) {
            wrap_joint_hess(id, y, th, out);
        } else if (order == 2) {
            wrap_joint_info(id, y, th, out);
        } else {
            wrap_joint_dinfo(id, y, th, out);
        }
        return true;
    }
    if (!joint_scalar(id)) return false;
    double o[2];
    if (order <= 1) d7_score_curv(id, 0, y, th, o);
    else d7_info_dinfo(id, 0, y, th, o);
    out[0] = o[order % 2];
    return true;
}

void d7_joint_grad(int id, double y, const double* th, double* out) {
    if (joint_route(id, 0, y, th, out)) return;
    switch (d7_joint_covered(id) >= 2 ? id : -1) {
    case 0: d7::gaussian1_joint_grad(y, th, out); break;
    case 1: d7::gamma1_joint_grad(y, th, out); break;
    case 3: d7::negbin2_joint_grad(y, th, out); break;
    case 4: d7::beta1_joint_grad(y, th, out); break;
    case 10: d7::cauchy_joint_grad(y, th, out); break;
    case 11: d7::logistic_joint_grad(y, th, out); break;
    case 12: d7::gaussian2_joint_grad(y, th, out); break;
    case 13: d7::gaussian3_joint_grad(y, th, out); break;
    case 14: d7::lognormal1_joint_grad(y, th, out); break;
    case 15: d7::invgauss1_joint_grad(y, th, out); break;
    case 16: d7::invgauss2_joint_grad(y, th, out); break;
    case 17: d7::gamma2_joint_grad(y, th, out); break;
    case 18: d7::gpd_joint_grad(y, th, out); break;
    case 19: d7::vonmises1_joint_grad(y, th, out); break;
    case 20: d7::vonmises2_joint_grad(y, th, out); break;
    case 21: d7::weibull3_joint_grad(y, th, out); break;
    case 22: d7::lognormal2_joint_grad(y, th, out); break;
    case 23: d7::student_t1_joint_grad(y, th, out); break;
    case 24: d7::student_t2_joint_grad(y, th, out); break;
    case 25: d7::gengamma1_joint_grad(y, th, out); break;
    case 26: d7::gengamma2_joint_grad(y, th, out); break;
    case 27: d7::gumbel_joint_grad(y, th, out); break;
    case 28: d7::laplace_joint_grad(y, th, out); break;
    case 29: d7::laplace2_joint_grad(y, th, out); break;
    case 30: d7::weibull1_joint_grad(y, th, out); break;
    case 31: d7::beta2_joint_grad(y, th, out); break;
    case 32: d7::enet_joint_grad(y, th, out); break;
    case 33: d7::negbin1_joint_grad(y, th, out); break;
    case 34: d7::skewnormal1_joint_grad(y, th, out); break;
    case 35: d7::skewt_joint_grad(y, th, out); break;
    case 36: d7::pseudohuber_joint_grad(y, th, out); break;
    case 37: d7::pseudohuber2_joint_grad(y, th, out); break;
    case 38: d7::skewnormal2_joint_grad(y, th, out); break;
    case 39: d7::pig1_joint_grad(y, th, out); break;
    case 40: d7::pig2_joint_grad(y, th, out); break;
    case 41: d7::betabinom1_joint_grad(y, th, out); break;
    case 42: d7::betabinom2_joint_grad(y, th, out); break;
    default: fill_nan(out, joint_np(id));
    }
}

void d7_joint_hess(int id, double y, const double* th, double* out) {
    if (joint_route(id, 1, y, th, out)) return;
    switch (d7_joint_covered(id) >= 2 ? id : -1) {
    case 0: d7::gaussian1_joint_hess(y, th, out); break;
    case 1: d7::gamma1_joint_hess(y, th, out); break;
    case 3: d7::negbin2_joint_hess(y, th, out); break;
    case 4: d7::beta1_joint_hess(y, th, out); break;
    case 10: d7::cauchy_joint_hess(y, th, out); break;
    case 11: d7::logistic_joint_hess(y, th, out); break;
    case 12: d7::gaussian2_joint_hess(y, th, out); break;
    case 13: d7::gaussian3_joint_hess(y, th, out); break;
    case 14: d7::lognormal1_joint_hess(y, th, out); break;
    case 15: d7::invgauss1_joint_hess(y, th, out); break;
    case 16: d7::invgauss2_joint_hess(y, th, out); break;
    case 17: d7::gamma2_joint_hess(y, th, out); break;
    case 18: d7::gpd_joint_hess(y, th, out); break;
    case 19: d7::vonmises1_joint_hess(y, th, out); break;
    case 20: d7::vonmises2_joint_hess(y, th, out); break;
    case 21: d7::weibull3_joint_hess(y, th, out); break;
    case 22: d7::lognormal2_joint_hess(y, th, out); break;
    case 23: d7::student_t1_joint_hess(y, th, out); break;
    case 24: d7::student_t2_joint_hess(y, th, out); break;
    case 25: d7::gengamma1_joint_hess(y, th, out); break;
    case 26: d7::gengamma2_joint_hess(y, th, out); break;
    case 27: d7::gumbel_joint_hess(y, th, out); break;
    case 28: d7::laplace_joint_hess(y, th, out); break;
    case 29: d7::laplace2_joint_hess(y, th, out); break;
    case 30: d7::weibull1_joint_hess(y, th, out); break;
    case 31: d7::beta2_joint_hess(y, th, out); break;
    case 32: d7::enet_joint_hess(y, th, out); break;
    case 33: d7::negbin1_joint_hess(y, th, out); break;
    case 34: d7::skewnormal1_joint_hess(y, th, out); break;
    case 35: d7::skewt_joint_hess(y, th, out); break;
    case 36: d7::pseudohuber_joint_hess(y, th, out); break;
    case 37: d7::pseudohuber2_joint_hess(y, th, out); break;
    case 38: d7::skewnormal2_joint_hess(y, th, out); break;
    case 39: d7::pig1_joint_hess(y, th, out); break;
    case 40: d7::pig2_joint_hess(y, th, out); break;
    case 41: d7::betabinom1_joint_hess(y, th, out); break;
    case 42: d7::betabinom2_joint_hess(y, th, out); break;
    default: fill_nan(out, joint_np(id) * joint_np(id));
    }
}

void d7_joint_info(int id, double y, const double* th, double* out) {
    if (joint_route(id, 2, y, th, out)) return;
    switch (d7_joint_covered(id) >= 4 ? id : -1) {
    case 0: d7::gaussian1_joint_info(y, th, out); break;
    case 1: d7::gamma1_joint_info(y, th, out); break;
    case 3: d7::negbin2_joint_info(y, th, out); break;
    case 4: d7::beta1_joint_info(y, th, out); break;
    case 10: d7::cauchy_joint_info(y, th, out); break;
    case 11: d7::logistic_joint_info(y, th, out); break;
    case 12: d7::gaussian2_joint_info(y, th, out); break;
    case 13: d7::gaussian3_joint_info(y, th, out); break;
    case 14: d7::lognormal1_joint_info(y, th, out); break;
    case 15: d7::invgauss1_joint_info(y, th, out); break;
    case 16: d7::invgauss2_joint_info(y, th, out); break;
    case 17: d7::gamma2_joint_info(y, th, out); break;
    case 18: d7::gpd_joint_info(y, th, out); break;
    case 19: d7::vonmises1_joint_info(y, th, out); break;
    case 20: d7::vonmises2_joint_info(y, th, out); break;
    case 21: d7::weibull3_joint_info(y, th, out); break;
    case 22: d7::lognormal2_joint_info(y, th, out); break;
    case 23: d7::student_t1_joint_info(y, th, out); break;
    case 24: d7::student_t2_joint_info(y, th, out); break;
    case 25: d7::gengamma1_joint_info(y, th, out); break;
    case 26: d7::gengamma2_joint_info(y, th, out); break;
    case 27: d7::gumbel_joint_info(y, th, out); break;
    case 28: d7::laplace_joint_info(y, th, out); break;
    case 29: d7::laplace2_joint_info(y, th, out); break;
    case 30: d7::weibull1_joint_info(y, th, out); break;
    case 31: d7::beta2_joint_info(y, th, out); break;
    case 32: d7::enet_joint_info(y, th, out); break;
    case 33: d7::negbin1_joint_info(y, th, out); break;
    case 39: d7::pig1_joint_info(y, th, out); break;
    case 40: d7::pig2_joint_info(y, th, out); break;
    case 41: d7::betabinom1_joint_info(y, th, out); break;
    case 42: d7::betabinom2_joint_info(y, th, out); break;
    case 34: d7::skewnormal1_joint_info(y, th, out); break;
    case 35: d7::skewt_joint_info(y, th, out); break;
    case 36: d7::pseudohuber_joint_info(y, th, out); break;
    case 37: d7::pseudohuber2_joint_info(y, th, out); break;
    case 38: d7::skewnormal2_joint_info(y, th, out); break;
    default: fill_nan(out, joint_np(id) * joint_np(id));
    }
}

void d7_joint_dinfo(int id, double y, const double* th, double* out) {
    if (joint_route(id, 3, y, th, out)) return;
    const int P = joint_np(id);
    switch (d7_joint_covered(id) >= 4 ? id : -1) {
    case 0: d7::gaussian1_joint_dinfo(y, th, out); break;
    case 1: d7::gamma1_joint_dinfo(y, th, out); break;
    case 3: d7::negbin2_joint_dinfo(y, th, out); break;
    case 4: d7::beta1_joint_dinfo(y, th, out); break;
    case 10: d7::cauchy_joint_dinfo(y, th, out); break;
    case 11: d7::logistic_joint_dinfo(y, th, out); break;
    case 12: d7::gaussian2_joint_dinfo(y, th, out); break;
    case 13: d7::gaussian3_joint_dinfo(y, th, out); break;
    case 14: d7::lognormal1_joint_dinfo(y, th, out); break;
    case 15: d7::invgauss1_joint_dinfo(y, th, out); break;
    case 16: d7::invgauss2_joint_dinfo(y, th, out); break;
    case 17: d7::gamma2_joint_dinfo(y, th, out); break;
    case 18: d7::gpd_joint_dinfo(y, th, out); break;
    case 19: d7::vonmises1_joint_dinfo(y, th, out); break;
    case 20: d7::vonmises2_joint_dinfo(y, th, out); break;
    case 21: d7::weibull3_joint_dinfo(y, th, out); break;
    case 22: d7::lognormal2_joint_dinfo(y, th, out); break;
    case 23: d7::student_t1_joint_dinfo(y, th, out); break;
    case 24: d7::student_t2_joint_dinfo(y, th, out); break;
    case 25: d7::gengamma1_joint_dinfo(y, th, out); break;
    case 26: d7::gengamma2_joint_dinfo(y, th, out); break;
    case 27: d7::gumbel_joint_dinfo(y, th, out); break;
    case 28: d7::laplace_joint_dinfo(y, th, out); break;
    case 29: d7::laplace2_joint_dinfo(y, th, out); break;
    case 30: d7::weibull1_joint_dinfo(y, th, out); break;
    case 31: d7::beta2_joint_dinfo(y, th, out); break;
    case 32: d7::enet_joint_dinfo(y, th, out); break;
    case 33: d7::negbin1_joint_dinfo(y, th, out); break;
    case 39: d7::pig1_joint_dinfo(y, th, out); break;
    case 40: d7::pig2_joint_dinfo(y, th, out); break;
    case 41: d7::betabinom1_joint_dinfo(y, th, out); break;
    case 42: d7::betabinom2_joint_dinfo(y, th, out); break;
    case 34: d7::skewnormal1_joint_dinfo(y, th, out); break;
    case 35: d7::skewt_joint_dinfo(y, th, out); break;
    case 36: d7::pseudohuber_joint_dinfo(y, th, out); break;
    case 37: d7::pseudohuber2_joint_dinfo(y, th, out); break;
    case 38: d7::skewnormal2_joint_dinfo(y, th, out); break;
    default: fill_nan(out, P * P * P);
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

// the twin of d7_scalar_probe() for the joint registry: the four orders at
// every observation, as an n x P matrix, an n x P x P and an n x P x P x P
// array; `covered` is d7_joint_covered(), and the orders beyond it are NaN
// [[Rcpp::export]]
Rcpp::List d7_joint_probe(std::string cls, Rcpp::NumericVector y,
                          Rcpp::NumericMatrix theta) {
    int id = d7_scalar_id(cls.c_str());
    int cov = d7_joint_covered(id);
    int n = y.size(), np = theta.ncol();
    int P = cov >= 2 ? id_np(id) : 0;
    Rcpp::NumericMatrix g(n, P);
    Rcpp::NumericVector h(n * P * P), e(n * P * P), de(n * P * P * P);
    std::vector<double> th(np), og(P), oh(P * P), oe(P * P), ode(P * P * P);
    for (int i = 0; i < n && P > 0; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        d7_joint_grad(id, y[i], th.data(), og.data());
        d7_joint_hess(id, y[i], th.data(), oh.data());
        d7_joint_info(id, y[i], th.data(), oe.data());
        d7_joint_dinfo(id, y[i], th.data(), ode.data());
        for (int a = 0; a < P; ++a) g(i, a) = og[a];
        for (int j = 0; j < P * P; ++j) {
            h[i + n * j] = oh[j];
            e[i + n * j] = oe[j];
        }
        for (int j = 0; j < P * P * P; ++j) de[i + n * j] = ode[j];
    }
    h.attr("dim") = Rcpp::IntegerVector::create(n, P, P);
    e.attr("dim") = Rcpp::IntegerVector::create(n, P, P);
    de.attr("dim") = Rcpp::IntegerVector::create(n, P, P, P);
    return Rcpp::List::create(Rcpp::_["id"] = id, Rcpp::_["covered"] = cov,
                              Rcpp::_["grad"] = g, Rcpp::_["hess"] = h,
                              Rcpp::_["info"] = e, Rcpp::_["dinfo"] = de);
}

// the classes the joint registry covers, for the twin tests
// [[Rcpp::export]]
Rcpp::CharacterVector d7_joint_classes_covered() {
    std::vector<std::string> out;
    for (int i = 0; i < d7_n_scalar_classes; ++i) {
        if (d7_joint_covered(i) > 0) out.push_back(d7_scalar_classes[i]);
    }
    return Rcpp::wrap(out);
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

// the tail degrees of freedom fold_rule() reads, one row per observation, NA
// where the parent has no such tail
// [[Rcpp::export]]
Rcpp::NumericVector d7_fold_tail_probe(std::string cls,
                                       Rcpp::NumericMatrix theta) {
    int id = d7_scalar_id(cls.c_str());
    int n = theta.nrow(), np = theta.ncol();
    Rcpp::NumericVector out(n, NA_REAL);
    std::vector<double> th(np);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        if (id < 0 || (id >= kWrapBase && !node_valid(id))) continue;
        const double v = fold_tail_nu(id, th.data());
        if (R_FINITE(v)) out[i] = v;
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

// the rule of quad_cdf_sums() for R, by class name, one row of theta per
// point: the side each point's integrals run over (TRUE for [q, hi)), and
// the nodes and weights with the point's index (from one); what the third
// and fourth derivatives of the distribution function integrate over
// [[Rcpp::export]]
Rcpp::List cdf_rule_cpp(std::string cls, Rcpp::NumericVector q,
                        Rcpp::NumericMatrix theta) {
    const int id = d7_scalar_id(cls.c_str());
    if (id < 0 || id >= kWrapBase || !quad_cdf_family(id))
        Rcpp::stop("no compiled cdf quadrature for '%s'", cls);
    const int n = q.size(), np = theta.ncol();
    bool ls;
    const int mask = quad_cdf_mask(id, &ls);
    Rcpp::LogicalVector upper(n);
    std::vector<int> idx, hidx;
    std::vector<double> yy, ww, ys, ws, th(np), hy0, hy1, hdl, hoff, ha, hle;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < np; ++j) th[j] = theta(i, j);
        const bool up = quad_side_upper(id, q[i], th.data());
        upper[i] = up;
        double slo, shi;
        cont_support(id, th.data(), &slo, &shi);
        CdfHead hd;
        const bool head = !up && !logs_family(id) &&
            cdf_head(id, q[i], th.data(), mask, &hd);
        if (head) {
            hidx.push_back(i + 1);
            hy0.push_back(hd.ya);
            hy1.push_back(hd.yb);
            hdl.push_back(hd.dL);
            hoff.push_back(hd.L0 - hd.La);
            ha.push_back(hd.alpha);
            hle.push_back(hd.logE);
        }
        const double a = up ? std::max(q[i], slo) : (head ? hd.y0 : slo);
        const double b = up ? shi : std::min(q[i], shi);
        ys.clear();
        ws.clear();
        double kk[1];
        const int nk = cdf_edge_cut(id, th.data(), kk);
        if (b > a)
            mapped_rule(id, th.data(), a, b, kk, nk, ys, ws, nullptr,
                        nullptr, kCdfSvMax);
        for (std::size_t j = 0; j < ys.size(); ++j) {
            idx.push_back(i + 1);
            yy.push_back(ys[j]);
            ww.push_back(ws[j]);
        }
    }
    return Rcpp::List::create(Rcpp::_["idx"] = Rcpp::wrap(idx),
                              Rcpp::_["y"] = Rcpp::wrap(yy),
                              Rcpp::_["w"] = Rcpp::wrap(ww),
                              Rcpp::_["upper"] = upper,
                              Rcpp::_["head_idx"] = Rcpp::wrap(hidx),
                              Rcpp::_["head_y0"] = Rcpp::wrap(hy0),
                              Rcpp::_["head_y1"] = Rcpp::wrap(hy1),
                              Rcpp::_["head_dL"] = Rcpp::wrap(hdl),
                              Rcpp::_["head_off"] = Rcpp::wrap(hoff),
                              Rcpp::_["head_alpha"] = Rcpp::wrap(ha),
                              Rcpp::_["head_logE"] = Rcpp::wrap(hle));
}

// fold_rule() for folded_expected() in R, which reads the same nodes; `nu` is
// the degrees of freedom of a skew t or student_t1 parent
// (d7_fold_tail_probe()), NA for any
// other parent
// [[Rcpp::export]]
Rcpp::List fold_rule_cpp(double c, double s, double nu) {
    std::vector<double> y, w;
    fold_rule(c, s, nu, y, w);
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
    R_RegisterCCallable("distributions7", "d7_joint_covered",
                        (DL_FUNC) d7_joint_covered);
    R_RegisterCCallable("distributions7", "d7_joint_grad",
                        (DL_FUNC) d7_joint_grad);
    R_RegisterCCallable("distributions7", "d7_joint_hess",
                        (DL_FUNC) d7_joint_hess);
    R_RegisterCCallable("distributions7", "d7_joint_info",
                        (DL_FUNC) d7_joint_info);
    R_RegisterCCallable("distributions7", "d7_joint_dinfo",
                        (DL_FUNC) d7_joint_dinfo);
}
