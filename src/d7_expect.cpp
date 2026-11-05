// The expected Hessian of a family by simulation and by integration, as the
// R route computes it: expected_by_mc() (R/expected_derivatives.R), and
// expected_by_integrate() / expected_by_bartlett() over expectation()
// (R/numerical_functions.R), at order 2 and on the parameter scale. The
// entry d7_expect() returns the same doubles as the R route and, for the
// simulation, consumes R's random stream as the R route does.
//
// The family's pieces are the headers pt_expect_*.h (pt_expect_common.h
// describes them); the derivatives at a draw or a node are the parts route
// of the registry, identical to the family's vector kernels; the quadrature,
// the series, Brent's zero finder and R's mean() are numericals7's
// n7_expect.h.

// R evaluates these expressions one operation at a time; clang (arm64
// macOS) would otherwise fuse a*b + c into one rounding
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
#include "pt_expect_common.h"
#include "pt_expect_simple.h"
#include "pt_expect_invgauss1.h"
#include "pt_expect_invgauss2.h"
#include "pt_expect_vonmises1.h"
#include "pt_expect_vonmises2.h"
#include "pt_expect_weibull3.h"
#include "pt_expect_enet.h"
#include "pt_expect_skewnormal1.h"
#include "pt_expect_skewnormal2.h"
#include "pt_expect_skewt.h"
#include "pt_expect_pseudohuber.h"
#include "pt_expect_pseudohuber2.h"
#include "pt_expect_pig1.h"
#include "pt_expect_pig2.h"
#include "pt_expect_betabinom1.h"
#include "pt_expect_betabinom2.h"
#include "pt_parts.h"
#include "pt_expect_qfallback.h"
#include <cstring>
#include <unordered_map>
#include <string>
#include <algorithm>

extern "C" {
int d7_parts_covered(int id);
void d7_parts(int id, double y, const double* th, int level, double* buf);
void d7_parts_grad(int id, const double* buf, double* out);
void d7_parts_hess(int id, const double* buf, double* out);
}

namespace {

// ------------------------------------------------------------ dispatch
struct Fam {
    int (*rng_ok)(const double*);
    void (*rng)(const double*, int, double*);
    void (*pdf)(const double*, int, const double*, int, double*);
    int (*quantile)(const double*, int, const double*, int, double*);
    bool discrete;
};

#define D7X_FAM(f, disc) Fam{d7x::f##_rng_ok, d7x::f##_rng, d7x::f##_pdf, d7x::f##_quantile, disc}

const Fam* fam_of(int id) {
    static const Fam tab[] = {
        D7X_FAM(gaussian1, false),  D7X_FAM(gamma1, false),     D7X_FAM(poisson, true),
        D7X_FAM(negbin2, true),     D7X_FAM(beta1, false),      D7X_FAM(bernoulli, true),
        D7X_FAM(binomial, true),    D7X_FAM(exponential, false), D7X_FAM(geometric, true),
        D7X_FAM(chisq, false),      D7X_FAM(cauchy, false),     D7X_FAM(logistic, false),
        D7X_FAM(gaussian2, false),  D7X_FAM(gaussian3, false),  D7X_FAM(lognormal1, false),
        D7X_FAM(invgauss1, false),  D7X_FAM(invgauss2, false),  D7X_FAM(gamma2, false),
        D7X_FAM(gpd, false),        D7X_FAM(vonmises1, false),  D7X_FAM(vonmises2, false),
        D7X_FAM(weibull3, false),   D7X_FAM(lognormal2, false), D7X_FAM(student_t1, false),
        D7X_FAM(student_t2, false), D7X_FAM(gengamma1, false),  D7X_FAM(gengamma2, false),
        D7X_FAM(gumbel, false),     D7X_FAM(laplace, false),    D7X_FAM(laplace2, false),
        D7X_FAM(weibull1, false),   D7X_FAM(beta2, false),      D7X_FAM(enet, false),
        D7X_FAM(negbin1, true),     D7X_FAM(skewnormal1, false), D7X_FAM(skewt, false),
        D7X_FAM(pseudohuber, false), D7X_FAM(pseudohuber2, false), D7X_FAM(skewnormal2, false),
        D7X_FAM(pig1, true),        D7X_FAM(pig2, true),        D7X_FAM(betabinom1, true),
        D7X_FAM(betabinom2, true)};
    if (id < 0 || id >= (int) (sizeof(tab) / sizeof(tab[0]))) return nullptr;
    return &tab[id];
}
#undef D7X_FAM

// the families with a parameter at which the log density has a kink
// (param_smoothness() FALSE): laplace, laplace2, enet, all in mu, where the
// kink's argument is y - mu (kink_decomposition()@v)
bool kinked(int id) { return id == 28 || id == 29 || id == 32; }

// the test that R evaluates over a whole vector, at one setting
int vb_of(int id, const double* th) {
    if (id == 21) return d7x::weibull3_vb(th);
    if (id == 42) return d7x::betabinom2_vb(th);
    return 0;
}

// ---------------------------------------------- derivatives at points
struct Deriv {
    int id, P;
    std::vector<double> buf, h;
    Deriv(int id_, int P_) : id(id_), P(P_), buf(d7::kPartsMax), h((std::size_t) P_ * P_) {}
    // the full P x P Hessian (column-major) at (y, th)
    const double* hess(double y, const double* th) {
        d7_parts(id, y, th, 2, buf.data());
        d7_parts_hess(id, buf.data(), h.data());
        return h.data();
    }
    // the gradient at (y, th), in h[0..P-1]
    const double* grad(double y, const double* th) {
        d7_parts(id, y, th, 2, buf.data());
        d7_parts_grad(id, buf.data(), h.data());
        return h.data();
    }
};

// the order-2 components in hess_names() order: the diagonal, then a < b row
// by row
std::vector<std::pair<int, int>> hess_pairs(int P) {
    std::vector<std::pair<int, int>> out;
    for (int a = 0; a < P; ++a) out.push_back({a, a});
    for (int a = 0; a < P; ++a)
        for (int b = a + 1; b < P; ++b) out.push_back({a, b});
    return out;
}

// the value of the identity's integrand at one point from the derivatives:
// the Hessian component, or -(0 + (1 * g_a) * g_b) of score_identity_terms()
// (the sign is applied by the caller for the integral, as R does)
inline double prod_score(double ga, double gb) {
    double pv = 1;
    pv = pv * ga;
    pv = pv * gb;
    return pv;
}

// ------------------------------------------------------------- the MC
// expected_by_mc() at one setting: the K means, in hess_names() order
void mc_setting(const Fam& F, Deriv& D, bool score, const double* th, int nsim,
                double* comp) {
    const int P = D.P;
    const auto pr = hess_pairs(P);
    const int K = (int) pr.size();
    std::vector<double> ys(nsim);
    F.rng(th, nsim, ys.data());
    std::vector<double> v((std::size_t) K * nsim);
    for (int i = 0; i < nsim; ++i) {
        const double* d = score ? D.grad(ys[i], th) : D.hess(ys[i], th);
        for (int k = 0; k < K; ++k) {
            const int a = pr[k].first, b = pr[k].second;
            double val;
            if (score) {
                double acc = 0;
                acc = acc + prod_score(d[a], d[b]);
                val = -acc;
            } else {
                val = d[a + P * b];
            }
            v[(std::size_t) k * nsim + i] = val;
        }
    }
    for (int k = 0; k < K; ++k) comp[k] = n7::r_mean(v.data() + (std::size_t) k * nsim, nsim);
}

// --------------------------------------------------------- integration
struct Tol { double qa, qr; int mp; double sa, sr; int mt; };

// R's seq(from, to, length.out = 129)
void seq129(double from, double to, std::vector<double>& out) {
    out.resize(129);
    if (from == to) { std::fill(out.begin(), out.end(), from); return; }
    const double del = to - from;
    out[0] = from;
    for (int k = 1; k <= 127; ++k) out[k] = from + k * (del / 128);
    out[128] = to;
}

// sort(unique(x)) of R for finite doubles
void sort_unique(std::vector<double>& x) {
    std::vector<double> u;
    for (double v : x) {
        bool dup = false;
        for (double w : u) if (w == v) { dup = true; break; }
        if (!dup) u.push_back(v);
    }
    std::stable_sort(u.begin(), u.end());
    x = u;
}

// .kink_roots() of R/kink_decomposition.R with v(y) = y - mu, and the
// polish of .kink_root_polish()
void kink_roots(double mu, double lower, double upper, std::vector<double>& out) {
    out.clear();
    const double w = upper - lower;
    std::vector<double> s1, s2;
    seq129(lower, upper, s1);
    seq129(lower + w / 4, upper - w / 4, s2);
    std::vector<double> y(s1);
    y.insert(y.end(), s2.begin(), s2.end());
    sort_unique(y);
    const int n = (int) y.size();
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i) v[i] = y[i] - mu;
    for (int i = 0; i < n; ++i) if (R_FINITE(v[i]) && v[i] == 0) out.push_back(y[i]);
    auto f = [mu](double x) { return x - mu; };
    for (int i = 0; i + 1 < n; ++i) {
        if (!(R_FINITE(v[i]) && R_FINITE(v[i + 1]) && v[i] * v[i + 1] < 0)) continue;
        double r;
        if (n7::uniroot(f, y[i], y[i + 1], v[i], v[i + 1], 1e-14, 1000, &r) != 0) continue;
        if (!R_FINITE(r)) continue;
        // .kink_root_polish(f, a, b, fa, fb, r)
        double a = y[i], b = y[i + 1], fa = v[i], fb = v[i + 1];
        double fr = f(r);
        double res = r;
        if (R_FINITE(fr) && fr != 0) {
            if (d7x::r_sign(fr) == d7x::r_sign(fa)) { a = r; fa = fr; } else { b = r; fb = fr; }
            bool done = false;
            for (;;) {
                const double m = a + (b - a) / 2;
                if (m <= a || m >= b) break;
                const double fm = f(m);
                if (!R_FINITE(fm)) { res = r; done = true; break; }
                if (fm == 0) { res = m; done = true; break; }
                if (d7x::r_sign(fm) == d7x::r_sign(fa)) { a = m; fa = fm; } else { b = m; fb = fm; }
            }
            if (!done) res = std::fabs(fa) <= std::fabs(fb) ? a : b;
        }
        if (R_FINITE(res)) out.push_back(res);
    }
    sort_unique(out);
}

// the knots of expectation() for a continuous family at one setting;
// false where the quantile cannot be reproduced
bool cont_knots(int id, const Fam& F, const double* th, int vb, bool score,
                const double* bnd, std::vector<double>& knots) {
    static const double pr[3] = {0.1, 0.5, 0.9};
    double q[3];
    if (F.quantile(pr, 3, th, vb, q) != 0) return false;
    std::vector<double> kj(q, q + 3);
    if (kinked(id) && score) {
        std::vector<double> qj;
        for (double v : kj) if (R_FINITE(v)) qj.push_back(v);
        double sp = 1;
        if (qj.size() > 1) {
            const double rg = *std::max_element(qj.begin(), qj.end()) -
                              *std::min_element(qj.begin(), qj.end());
            if (rg > 0) sp = rg;
        }
        std::vector<double> q0(qj);
        q0.push_back(0);
        const double lo = R_FINITE(bnd[0]) ? bnd[0] : *std::min_element(q0.begin(), q0.end()) - 30 * sp;
        const double hi = R_FINITE(bnd[1]) ? bnd[1] : *std::max_element(q0.begin(), q0.end()) + 30 * sp;
        // kink_knots(): numeric(0) unless lower < upper, both finite
        if (R_FINITE(lo) && R_FINITE(hi) && lo < hi) {
            std::vector<double> kr;
            kink_roots(th[0], lo, hi, kr);
            kj.insert(kj.end(), kr.begin(), kr.end());
        }
    }
    // unique(kj[is.finite(kj) & kj > b1 & kj < b2]), then sorted with the bounds
    std::vector<double> in;
    for (double v : kj)
        if (R_FINITE(v) && v > bnd[0] && v < bnd[1]) {
            bool dup = false;
            for (double w : in) if (w == v) { dup = true; break; }
            if (!dup) in.push_back(v);
        }
    knots.assign(1, bnd[0]);
    knots.insert(knots.end(), in.begin(), in.end());
    knots.push_back(bnd[1]);
    std::stable_sort(knots.begin(), knots.end());
    return true;
}

// stats::integrate(g, lower, upper, rel.tol, abs.tol, subdivisions)$value:
// Rdqags with both limits finite, Rdqagi otherwise (inf = 1 above a finite
// lower limit, -1 below a finite upper one, 2 over the line, bound 0); 0, or
// 1 where integrate() signals an error
template <class G>
int r_integrate(G& g, double lower, double upper, double abs_tol, double rel_tol,
                int limit, double* value) {
    if (R_FINITE(lower) && R_FINITE(upper))
        return d7x::r_integrate_fin(g, lower, upper, abs_tol, rel_tol, limit, value);
    if (ISNAN(lower) || ISNAN(upper)) return 1;
    int inf;
    double bound;
    if (R_FINITE(lower)) { inf = 1; bound = lower; }
    else if (R_FINITE(upper)) { inf = -1; bound = upper; }
    else { inf = 2; bound = 0; }
    int lenw = 4 * limit;
    std::vector<int> iwork(limit);
    std::vector<double> work(lenw);
    double epsabs = abs_tol, epsrel = rel_tol, result = 0, abserr = 0;
    int neval = 0, ier = 0, last = 0;
    d7x::RdqCtx<G> ctx{&g, false};
    Rdqagi(d7x::rdq_trampoline<G>, &ctx, &bound, &inf, &epsabs, &epsrel, &result,
           &abserr, &neval, &ier, &limit, &lenw, &last, iwork.data(), work.data());
    if (ctx.bad || ier != 0) return 1;
    *value = result;
    return 0;
}

// expectation() of one component at one continuous setting; 0, or 2 where
// a sub-interval fails (R then runs its rescue)
int cont_component(const Fam& F, Deriv& D, const double* th, int vb, bool score,
                   int a, int b, const std::vector<double>& knots, const Tol& tol,
                   double* value) {
    const int P = D.P;
    std::vector<double> pdf;
    auto integrand = [&](const double* x, int nx, double* out) {
        pdf.resize(nx);
        F.pdf(x, nx, th, vb, pdf.data());
        for (int i = 0; i < nx; ++i) {
            double val_f;
            if (score) {
                const double* g = D.grad(x[i], th);
                val_f = prod_score(g[a], g[b]);
            } else {
                val_f = D.hess(x[i], th)[a + P * b];
            }
            out[i] = val_f * pdf[i];
        }
    };
    double tot = 0;
    bool failed = false;
    for (std::size_t k = 0; k + 1 < knots.size(); ++k) {
        double v;
        if (n7::quad_row(integrand, knots[k], knots[k + 1], tol.qa, tol.qr, tol.mp, &v) != 0) {
            failed = true;
            break;
        }
        tot = tot + v;
    }
    if (!failed) {
        *value = tot;
        return 0;
    }
    // expectation()'s rescue: the whole setting again by stats::integrate()
    // between consecutive knots, summed in long double; an error there is
    // the R route's error
    auto one = [&](double x) {
        double o;
        integrand(&x, 1, &o);
        return o;
    };
    const double rt = R_pow(DBL_EPSILON, 0.25);
    long double s = 0;
    for (std::size_t k = 0; k + 1 < knots.size(); ++k) {
        double v;
        if (r_integrate(one, knots[k], knots[k + 1], rt, rt, 100, &v) != 0) return 2;
        s += v;
    }
    *value = (double) s;
    return 0;
}

// expectation() of one component at one discrete setting
int disc_component(const Fam& F, Deriv& D, const double* th, int vb, bool score,
                   int a, int b, const double* bnd, const Tol& tol, double* value) {
    const int P = D.P;
    auto term_at = [&](double k) {
        double pdf;
        F.pdf(&k, 1, th, vb, &pdf);
        double val_f;
        if (score) {
            const double* g = D.grad(k, th);
            val_f = prod_score(g[a], g[b]);
        } else {
            val_f = D.hess(k, th)[a + P * b];
        }
        return val_f * pdf;
    };
    if (R_FINITE(bnd[0]) && R_FINITE(bnd[1])) {
        long double s = 0;
        for (double k = bnd[0]; k <= bnd[1]; k += 1) s += term_at(k);
        *value = (double) s;
        return 0;
    }
    if (!R_FINITE(bnd[0])) return 1;
    auto term = [&](const int* k, int nk, double* out) {
        for (int j = 0; j < nk; ++j) out[j] = term_at((double) k[j]);
    };
    double v;
    const int st = n7::series_row(term, (int) bnd[0], tol.sa, tol.sr, tol.mt, &v);
    if (st != 0) return 2;
    *value = v;
    return 0;
}

}  // namespace

extern "C" {

// The order-2 expected Hessian on the parameter scale at nobs rows of
// parameters: th is nobs x ncol, column-major, the P parameters then the
// family's constants. method: 0 expect_mc(identity = "hessian"), 1
// expect_mc(identity = "score"), 2 expect_integrate(identity = "hessian"),
// 3 expect_integrate(identity = "score"). bounds is the family's support,
// tol (quad_atol, quad_rtol, max_panels, series_atol, series_rtol,
// max_terms). out receives nobs full P x P matrices (column-major, one after
// the other). Returns 0; 1 where the call is not covered (decided before any
// draw); 2 where the quadrature, the series or a quantile fails, at which
// the R route runs its rescue or signals its error.
int d7_expect(int id, int method, int nobs, int ncol, const double* th,
              const double* bounds, double nsim, const double* tol_in, double* out) {
    const Fam* F = fam_of(id);
    if (!F || d7_parts_covered(id) < 2) return 1;
    const bool mc = method <= 1, score = method == 1 || method == 3;
    if (!score && kinked(id)) return 1;
    if (mc && !(nsim >= 1 && nsim == std::floor(nsim) && nsim <= INT_MAX)) return 1;

    // the distinct settings, by bit pattern, in the order of first appearance
    std::vector<int> first, map(nobs);
    {
        std::unordered_map<std::string, int> seen;
        std::string key((std::size_t) ncol * sizeof(double), '\0');
        for (int i = 0; i < nobs; ++i) {
            for (int c = 0; c < ncol; ++c) {
                const double v = th[i + (std::size_t) nobs * c];
                std::memcpy(&key[(std::size_t) c * sizeof(double)], &v, sizeof(double));
            }
            auto it = seen.find(key);
            if (it == seen.end()) {
                seen.emplace(key, (int) first.size());
                map[i] = (int) first.size();
                first.push_back(i);
            } else {
                map[i] = it->second;
            }
        }
    }
    const int nc = (int) first.size();
    std::vector<double> ths((std::size_t) nc * ncol);
    for (int j = 0; j < nc; ++j)
        for (int c = 0; c < ncol; ++c)
            ths[(std::size_t) j * ncol + c] = th[first[j] + (std::size_t) nobs * c];

    // whole-vector tests: covered only where every setting agrees
    int vb = 0;
    if (id == 21 || id == 42) {
        vb = vb_of(id, ths.data());
        for (int j = 1; j < nc; ++j)
            if (vb_of(id, ths.data() + (std::size_t) j * ncol) != vb) return 1;
    }

    // P from the registry's Hessian layout: ncol minus the constants
    int P = ncol;
    if (id == 6 || id == 41 || id == 42) P = ncol - 1;
    const auto pr = hess_pairs(P);
    const int K = (int) pr.size();
    std::vector<double> E((std::size_t) nc * K);
    Deriv D(id, P);

    if (mc) {
        for (int j = 0; j < nc; ++j)
            if (F->rng_ok(ths.data() + (std::size_t) j * ncol) != 0) return 1;
        GetRNGstate();
        for (int j = 0; j < nc; ++j)
            mc_setting(*F, D, score, ths.data() + (std::size_t) j * ncol, (int) nsim,
                       E.data() + (std::size_t) j * K);
        PutRNGstate();
    } else {
        const Tol tol{tol_in[0], tol_in[1], (int) tol_in[2], tol_in[3], tol_in[4], (int) tol_in[5]};
        for (int j = 0; j < nc; ++j) {
            const double* t = ths.data() + (std::size_t) j * ncol;
            std::vector<double> knots;
            if (!F->discrete && !cont_knots(id, *F, t, vb, score, bounds, knots)) return 2;
            for (int k = 0; k < K; ++k) {
                double v;
                const int st = F->discrete
                    ? disc_component(*F, D, t, vb, score, pr[k].first, pr[k].second, bounds, tol, &v)
                    : cont_component(*F, D, t, vb, score, pr[k].first, pr[k].second, knots, tol, &v);
                if (st != 0) return st;
                if (score) {
                    double acc = 0;
                    acc = acc + v;
                    v = -acc;
                }
                E[(std::size_t) j * K + k] = v;
            }
        }
    }
    for (int i = 0; i < nobs; ++i) {
        const double* e = E.data() + (std::size_t) map[i] * K;
        double* o = out + (std::size_t) i * P * P;
        for (int k = 0; k < K; ++k) {
            const int a = pr[k].first, b = pr[k].second;
            o[a + P * b] = e[k];
            o[b + P * a] = e[k];
        }
    }
    return 0;
}

}  // extern "C"

// The entry for the twin tests: the same call from R, returning the nobs x K
// components in hess_names() order, or NULL where the entry returns nonzero
// (the code as the attribute "code").
// [[Rcpp::export]]
SEXP d7_expect_probe(int id, int method, Rcpp::NumericMatrix th,
                     Rcpp::NumericVector bounds, double nsim,
                     Rcpp::NumericVector tol) {
    const int nobs = th.nrow(), ncol = th.ncol();
    int P = ncol;
    if (id == 6 || id == 41 || id == 42) P = ncol - 1;
    std::vector<double> out((std::size_t) nobs * P * P);
    const int st = d7_expect(id, method, nobs, ncol, th.begin(), bounds.begin(), nsim,
                             tol.begin(), out.data());
    if (st != 0) {
        Rcpp::IntegerVector r = Rcpp::IntegerVector::create(st);
        return r;
    }
    const auto pr = hess_pairs(P);
    Rcpp::NumericMatrix res(nobs, (int) pr.size());
    for (int i = 0; i < nobs; ++i)
        for (std::size_t k = 0; k < pr.size(); ++k)
            res(i, k) = out[(std::size_t) i * P * P + pr[k].first + P * pr[k].second];
    return res;
}

// The pieces of a family for their twin tests: what = 0 the density
// (log = FALSE) at x, 1 the quantile at the probabilities x (NULL where the
// family's quantile returns nonzero), 2 length(x) draws (x unused), 3 the
// whole-vector test vb; th is one setting, parameters then constants.
// [[Rcpp::export]]
SEXP d7_expect_piece_probe(int id, int what, Rcpp::NumericVector x,
                           Rcpp::NumericVector th, int vb = 0) {
    const Fam* F = fam_of(id);
    if (!F) Rcpp::stop("unknown family");
    const int n = x.size();
    Rcpp::NumericVector out(n);
    if (what == 0) {
        F->pdf(x.begin(), n, th.begin(), vb, out.begin());
    } else if (what == 1) {
        if (F->quantile(x.begin(), n, th.begin(), vb, out.begin()) != 0) return R_NilValue;
    } else if (what == 2) {
        if (F->rng_ok(th.begin()) != 0) return R_NilValue;
        GetRNGstate();
        F->rng(th.begin(), n, out.begin());
        PutRNGstate();
    } else {
        return Rcpp::wrap(vb_of(id, th.begin()));
    }
    return out;
}

