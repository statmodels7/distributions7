#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include <cstring>
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
    "NegBin1Distrib"        // 33
};

const int d7_n_scalar_classes =
    sizeof(d7_scalar_classes) / sizeof(d7_scalar_classes[0]);

} // namespace

extern "C" {

// Called once by the consumer on its own thread, before any loop: the von
// Mises families resolve numericals7's Bessel-ratio entry points here, so
// that no lookup into R happens later from a worker.
int d7_scalar_id(const char* cls) {
    for (int i = 0; i < d7_n_scalar_classes; ++i) {
        if (std::strcmp(cls, d7_scalar_classes[i]) == 0) {
            if (std::strncmp(cls, "VonMises", 8) == 0) d7::vm_bessel();
            return i;
        }
    }
    return -1;
}

// k is the 0-based index of the parameter among the family's own, th the
// full parameter vector at this observation; out[0] the score component,
// out[1] the (k, k) second derivative, both on the parameter scale
void d7_score_curv(int id, int k, double y, const double* th, double* out) {
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
// silent ME_DOMAIN; and numericals7's Bessel ratio, which is plain C. A
// family whose entries evaluate a function that may warn (pt, pbeta,
// lchoose at a non-integer, the Bessel K) answers 0.
int d7_scalar_thread_safe(int id) {
    if (id < 0 || id >= d7_n_scalar_classes) return -1;
    return 1;
}

// out[0] the (k, k) expected second derivative E[l_kk], out[1] its
// derivative in the same parameter, both on the parameter scale
void d7_info_dinfo(int id, int k, double y, const double* th, double* out) {
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
    default:
        out[0] = R_NaN; out[1] = R_NaN;
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
}
