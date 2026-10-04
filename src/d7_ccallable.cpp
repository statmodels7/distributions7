#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <cstring>
#include <cstdlib>
#include <string>
#include <cmath>
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
// A wrapped family is named "<WrapperClass>[:aux]|<InnerClass>" by
// distrib_scalar_route() and identified by kWrapBase * w + 1000 * aux + inner,
// w the wrapper's code below and aux what the wrapper needs besides its
// constants (for fixed() the mask of the fixed parameters). th is the
// wrapper's parameters, then the inner family's constants, then the
// wrapper's own; k counts the wrapper's parameters from zero. Each entry
// rebuilds the inner family's vector and reads the inner family's entries;
// the composition formulas are in pt_wrappers.h.

const int kWrapBase = 100000;

struct WrapName { const char* name; int code; };
const WrapName d7_wrappers[] = {
    {"FixedContinuousDistrib", 1},
    {"FixedDiscreteDistrib", 1},
    {"ZeroInflatedDistrib", 2},
    {"ZeroAdjustedDiscreteDistrib", 3},
    {"ZeroAdjustedContinuousDistrib", 4}
};
const int d7_n_wrappers = sizeof(d7_wrappers) / sizeof(d7_wrappers[0]);

}  // namespace

extern "C" {
int d7_scalar_id(const char* cls);
void d7_score_curv(int id, int k, double y, const double* th, double* out);
void d7_info_dinfo(int id, int k, double y, const double* th, double* out);
double d7_logpdf(int id, double y, const double* th);
}

namespace {

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
    if (code < 0 || std::strchr(bar + 1, '|') != nullptr) return -1;
    int inner = d7_scalar_id(bar + 1);
    if (inner < 0 || aux < 0 || aux >= 100) return -1;
    if (code == 1) {
        const int P = d7_n_params[inner];
        if (aux == 0 || aux >= (1 << P)) return -1;
    } else if (aux != 0) {
        return -1;
    }
    return kWrapBase * code + 1000 * aux + inner;
}

// fixed(): the inner vector from the free values, the inner constants and
// the fixed values, and the inner index of free parameter k
int fixed_inner(int inner, int mask, int k, const double* th, double* full) {
    const int P = d7_n_params[inner], nc = d7_n_constants[inner];
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
double zero_inner(int inner, const double* th, double* full) {
    const int P = d7_n_params[inner], nc = d7_n_constants[inner];
    for (int j = 0; j < P; ++j) full[j] = th[j];
    for (int c = 0; c < nc; ++c) full[P + c] = th[P + 1 + c];
    return th[P];
}

void wrap_score_curv(int id, int k, double y, const double* th, double* out) {
    const int code = id / kWrapBase, aux = (id % kWrapBase) / 1000,
        inner = id % 1000;
    double full[16];
    if (code == 1) {
        const int j = fixed_inner(inner, aux, k, th, full);
        d7_score_curv(inner, j, y, full, out);
        return;
    }
    if (code >= 2 && code <= 4) {
        const int P = d7_n_params[inner];
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
    const int code = id / kWrapBase, aux = (id % kWrapBase) / 1000,
        inner = id % 1000;
    double full[16];
    if (code == 1) {
        const int j = fixed_inner(inner, aux, k, th, full);
        d7_info_dinfo(inner, j, y, full, out);
        return;
    }
    if (code >= 2 && code <= 4) {
        const int P = d7_n_params[inner];
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
    const int code = id / kWrapBase, aux = (id % kWrapBase) / 1000,
        inner = id % 1000;
    double full[16];
    if (code == 1) {
        fixed_inner(inner, aux, 0, th, full);
        return d7_logpdf(inner, y, full);
    }
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
        wrap_score_curv(id, k, y, th, out);
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
// pt, pbeta, lchoose at a non-integer, R's Bessel K) answers 0.
int d7_scalar_thread_safe(int id) {
    if (id >= kWrapBase) id = id % 1000;
    if (id < 0 || id >= d7_n_scalar_classes) return -1;
    return 1;
}

// out[0] the (k, k) expected second derivative E[l_kk], out[1] its
// derivative in the same parameter, both on the parameter scale
void d7_info_dinfo(int id, int k, double y, const double* th, double* out) {
    if (id >= kWrapBase) {
        wrap_info_dinfo(id, k, y, th, out);
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
    if (id >= kWrapBase) return wrap_logpdf(id, y, th);
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
