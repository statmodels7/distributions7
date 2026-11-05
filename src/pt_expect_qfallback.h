// The numerical quantile of a continuous family, as R computes it, at ONE
// parameter setting: a transcription of the S7 method
// distrib_quantile(<continuous_distrib>) (R/numerical_fallbacks.R:268-384)
// and of quantile_anchors() (:402-441), for lower.tail = TRUE, log.p = FALSE.
//
// A family supplies a type F with
//   double lb, ub;                    the support bounds (distrib@bounds)
//   bool   fail;                      set by a callback where R would signal
//                                     an error; the quantile then returns 1
//   double lpdf(double x);            distrib_pdf(distrib, x, th, log = TRUE)
//   double pdf(double x);             distrib_pdf(distrib, x, th)
//   double lcdf(double x, bool lower);distrib_cdf(distrib, x, th,
//                                       lower.tail = lower, log.p = TRUE)
//   double cdf(double x);             distrib_cdf(distrib, x, th)
//
// Every operation of the R method is elementwise except the call selections
// (`if (any(lo_i))`), so each probability is solved on its own; the anchor,
// its height and its probability are computed once, as R does for one
// distinct parameter setting. The R loops run over the unfinished elements
// with one common counter, which for an element is its own iteration count.
//
// Also here: r_pmin() and r_integrate_*(), the C form of stats::integrate()
// through R's Rdqags()/Rdqagi() (R_ext/Applic.h), which are what
// integrate() calls, with its error conditions reported as a return of 1.

#ifndef D7_PT_EXPECT_QFALLBACK_H
#define D7_PT_EXPECT_QFALLBACK_H

#include "pt_expect_common.h"
#include <R_ext/Applic.h>
#include <algorithm>
#include <type_traits>

namespace d7x {

// R's vectorized exp, expm1, log1p, sin, cos and tan are not the C functions:
// arithmetic.c (R 4.6.0, Rexp ... Rtan, :1300-1360) returns 1 + x, x or
// 1 - x*x*0.5 for a small argument. Every R-level call of these functions is
// transcribed with the functions below (Rmath routines are not affected).

// R's pmin(a, b) of two doubles: NaN if either is NaN
inline double r_pmin(double a, double b) {
    if (ISNAN(a)) return a;
    if (ISNAN(b)) return b;
    return a <= b ? a : b;
}

// ---- stats::integrate() ------------------------------------------------------

template <class G>
struct RdqCtx {
    G* g;
    bool bad;
};

// Rintfn() of R's integrate.c: the integrand is applied to the vector and a
// non-finite value is an error ("non-finite function value")
template <class G>
inline void rdq_trampoline(double* x, int n, void* ex) {
    RdqCtx<G>* c = (RdqCtx<G>*) ex;
    for (int i = 0; i < n; ++i) {
        double v = (*c->g)(x[i]);
        if (!R_FINITE(v)) c->bad = true;
        x[i] = v;
    }
}

// integrate(g, lower, upper, rel.tol = rel_tol, abs.tol = abs_tol,
//           subdivisions = limit)$value with both limits finite; returns 0
// with the value, or 1 where integrate() signals an error (invalid
// parameters, ier > 0, a non-finite function value). g(double) -> double.
template <class G>
inline int r_integrate_fin(G&& g, double lower, double upper, double abs_tol,
                           double rel_tol, int limit, double* value) {
    if (limit < 1 || (abs_tol <= 0 && rel_tol < std::max(50 * DBL_EPSILON, 5e-29)))
        return 1;
    if (ISNAN(lower) || ISNAN(upper)) return 1;
    if (!(R_FINITE(lower) && R_FINITE(upper))) return 1;
    int lenw = 4 * limit;
    std::vector<int> iwork(limit);
    std::vector<double> work(lenw);
    double a = lower, b = upper, epsabs = abs_tol, epsrel = rel_tol;
    double result = 0, abserr = 0;
    int neval = 0, ier = 0, last = 0;
    RdqCtx<typename std::remove_reference<G>::type> ctx{&g, false};
    Rdqags(rdq_trampoline<typename std::remove_reference<G>::type>, &ctx, &a,
           &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &limit,
           &lenw, &last, iwork.data(), work.data());
    if (ctx.bad || ier != 0) return 1;
    *value = result;
    return 0;
}

// the same with upper = Inf and a finite lower bound (integrate()'s dqagi
// branch with inf = 1)
template <class G>
inline int r_integrate_inf1(G&& g, double bound, double abs_tol,
                            double rel_tol, int limit, double* value) {
    if (limit < 1 || (abs_tol <= 0 && rel_tol < std::max(50 * DBL_EPSILON, 5e-29)))
        return 1;
    if (!R_FINITE(bound)) return 1;
    int lenw = 4 * limit;
    std::vector<int> iwork(limit);
    std::vector<double> work(lenw);
    double bd = bound, epsabs = abs_tol, epsrel = rel_tol;
    int inf = 1;
    double result = 0, abserr = 0;
    int neval = 0, ier = 0, last = 0;
    RdqCtx<typename std::remove_reference<G>::type> ctx{&g, false};
    Rdqagi(rdq_trampoline<typename std::remove_reference<G>::type>, &ctx, &bd,
           &inf, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &limit,
           &lenw, &last, iwork.data(), work.data());
    if (ctx.bad || ier != 0) return 1;
    *value = result;
    return 0;
}

// ---- the quantile fallback ---------------------------------------------------

// the variable of quantile_anchors()'s grid (:405-413)
inline double qf_to_y(double lb, double ub, double t) {
    const bool fl = R_FINITE(lb), fu = R_FINITE(ub);
    if (fl && fu) return lb + t * (ub - lb);
    if (fl) return lb + t / (1 - t);
    if (fu) return ub - (1 - t) / t;
    return r_tan(M_PI * (t - 0.5));
}

// quantile_anchors() for one setting (k = 1): 129-point grids refined up to
// 60 times; the grid is seq(lo, hi, length.out = 129L) of seq.default
template <class F>
inline double qf_anchor(F& F_) {
    const int G = 129;
    const double b1 = F_.lb, b2 = F_.ub;
    double lo = 1e-10, hi = 1 - 1e-10;
    bool act = true;
    std::vector<double> tt(G), v(G);
    for (int it = 1; it <= 60; ++it) {
        if (!act) break;
        if (lo == hi) {
            for (int k = 0; k < G; ++k) tt[k] = lo;
        } else {
            const double del = hi - lo;
            tt[0] = lo;
            tt[G - 1] = hi;
            for (int k = 1; k <= G - 2; ++k) tt[k] = lo + k * (del / (G - 1));
        }
        for (int k = 0; k < G; ++k) {
            v[k] = F_.lpdf(qf_to_y(b1, b2, tt[k]));
            if (ISNAN(v[k])) v[k] = R_NegInf;
        }
        int kk = 0;
        double large = v[0];
        for (int c = 1; c < G; ++c) if (large < v[c]) { large = v[c]; kk = c; }
        bool none = true;
        for (int k = 0; k < G; ++k) if (R_FINITE(v[k])) { none = false; break; }
        const double new_lo = tt[std::max(0, kk - 1)];
        const double new_hi = tt[std::min(G - 1, kk + 1)];
        if (!none) { lo = new_lo; hi = new_hi; }
        const double y_lo = qf_to_y(b1, b2, lo), y_hi = qf_to_y(b1, b2, hi);
        const double d = y_hi - y_lo;
        const bool fin = R_FINITE(d) &&
            d <= 1e-10 * r_pmax(r_pmax(1.0, std::fabs(y_lo)), std::fabs(y_hi));
        if (fin || none) act = false;
    }
    return qf_to_y(b1, b2, (lo + hi) / 2);
}

// distrib_quantile(distrib, p, th) for the numerical fallback; returns 0, or 1
// where F.fail was set or the anchor is not a number (R's behaviour is then an
// error or an NA index, which is not reproduced)
template <class F>
inline int qf_quantile(const double* p, int np, F& F_, double* out) {
    const double b1 = F_.lb, b2 = F_.ub;
    std::vector<double> llo(np), lup(np);
    std::vector<char> todo(np, 0);
    bool any = false;
    for (int i = 0; i < np; ++i) {
        // quantile_log_tails(p, TRUE, FALSE) (utility_functions.R:620-626)
        double lp = std::log(p[i]);
        if (!ISNAN(lp) && lp > 0) lp = R_NaN;
        double lq = r_log1p(-r_exp(lp));
        if (!ISNAN(lp) && lp > -std::log(2.0)) lq = std::log(-r_expm1(lp));
        llo[i] = lp;
        lup[i] = lq;
        const bool ok = !ISNAN(llo[i]) && !ISNAN(lup[i]) && llo[i] <= 0 && lup[i] <= 0;
        out[i] = R_NaN;
        if (ok && llo[i] == R_NegInf) out[i] = b1;
        if (ok && lup[i] == R_NegInf) out[i] = b2;
        if (ok && llo[i] > R_NegInf && lup[i] > R_NegInf) { todo[i] = 1; any = true; }
    }
    if (!any) return 0;

    const double m = qf_anchor(F_);
    const double f_u = F_.pdf(m);
    const double Fm = F_.cdf(m);
    if (F_.fail) return 1;
    if (ISNAN(m) || ISNAN(f_u) || ISNAN(Fm)) return 1;
    const double s = 1 / r_pmax(f_u, 1e-12);
    const double big = DBL_MAX / 4;
    const double eps34 = R_POW(DBL_EPSILON, 0.75);

    for (int i = 0; i < np; ++i) {
        if (!todo[i]) continue;
        const bool low = llo[i] <= lup[i];
        const double target = low ? llo[i] : lup[i];
        auto lcdf_at = [&](double x) { return F_.lcdf(x, low); };

        const bool left = low ? (target <= std::log(Fm)) : (target > r_log1p(-Fm));
        double lo = left ? b1 : m;
        double hi = left ? m : b2;
        bool beyond = false;

        // the bracket
        for (int side = 0; side < 2; ++side) {
            bool open = side == 0 ? (left && !R_FINITE(b1)) : (!left && !R_FINITE(b2));
            double step = s;
            for (int it = 1; it <= 150; ++it) {
                if (!open) break;
                const double x = side == 0 ? r_pmax(m - step, -big)
                                           : r_pmin(m + step, big);
                const double r = lcdf_at(x) - target;
                if (side == 0) lo = x; else hi = x;
                bool more;
                if (side == 0) more = (r > 0 && low) || (r < 0 && !low);
                else more = (r < 0 && low) || (r > 0 && !low);
                const bool edge = more && std::fabs(x) == big;
                if (edge) beyond = true;
                more = more && !edge;
                open = more;
                step = (it > 30 ? 1e3 : 2.0) * step;
            }
        }

        // the safeguarded Newton iteration
        double x = (R_FINITE(lo) && R_FINITE(hi)) ? (lo + hi) / 2
                                                  : (R_FINITE(lo) ? lo : hi);
        if (beyond) {
            out[i] = left ? R_NegInf : R_PosInf;
            continue;
        }
        for (int it = 1; it <= 100; ++it) {
            const double xa = x;
            const double lP = lcdf_at(xa);
            const double r = lP - target;
            const bool nar = ISNAN(r);
            const bool up = !nar && ((r > 0) == low);
            if (up) hi = xa;
            else if (!nar) lo = xa;
            const double slope = r_exp(F_.lpdf(xa) - lP) * (low ? 1.0 : -1.0);
            double xn = xa - r / slope;
            const bool bad = !R_FINITE(xn) || xn < lo || xn > hi;
            if (bad) xn = (lo + hi) / 2;
            const bool hit = !nar &&
                std::fabs(r) <= 4 * DBL_EPSILON * r_pmax(1.0, std::fabs(target));
            if (hit) xn = xa;
            const double tol = eps34 * r_pmax(1.0, std::fabs(xa));
            const bool done = hit || std::fabs(xn - xa) <= tol || (hi - lo) <= tol;
            x = xn;
            if (done) break;
        }
        out[i] = x;
    }
    return F_.fail ? 1 : 0;
}

}  // namespace d7x

#endif
