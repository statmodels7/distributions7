// vonmises1: generator, density and quantile for the compiled expected
// information by simulation and by integration (see pt_expect_common.h).
// th = (mu, kappa).
//
// The quantile is the numerical fallback of pt_expect_qfallback.h over this
// family's density and distribution function. The distribution function is
// transcribed from R/vonmises1_distrib.R: vm_cdf() (:989), vm_log_tail()
// (:1050), vm_log_piece() (:1074, which calls stats::integrate() through
// Rdqags) and vm_cdf_tail() (:1098), with numericals7's bessel_i_ratios()
// (R/special.R:288) and log_bessel_i() (compiled entry n7_log_bessel_i).

#ifndef D7_PT_EXPECT_VONMISES1_H
#define D7_PT_EXPECT_VONMISES1_H

#include "pt_expect_common.h"
#include "pt_expect_qfallback.h"
#include <R_ext/Rdynload.h>

namespace d7x {

// n7_log_bessel_i of numericals7, resolved on first use (not a guarded
// static: see vm_bessel() in vonmises.cpp)
typedef double (*vx_lbi_t)(double, double);
inline vx_lbi_t vx_lbi() {
    static vx_lbi_t p = nullptr;
    if (p == nullptr) p = (vx_lbi_t) R_GetCCallable("numericals7", "n7_log_bessel_i");
    return p;
}

// numericals7::bessel_i_ratios(kappa, m) for one finite positive kappa
// (R/special.R:288-312): the downward recurrence from n0, then cumprod, which
// R accumulates in long double
inline void vm1_ratios(double k, int m, std::vector<double>& R) {
    const double kmax = k;
    const double n0 = m + std::max(std::max(30.0, std::ceil(std::sqrt(40.0 * m))),
                                   std::ceil(kmax));
    std::vector<double> keep(m, 0.0);
    double rj = 0;
    for (long j = (long) n0; j >= 1; --j) {
        rj = 1 / (2 * (double) j / k + rj);
        if (j <= m) keep[j - 1] = rj;
    }
    R.assign(m, 0.0);
    if (m > 1) {
        long double prod = 1.0;
        for (int c = 0; c < m; ++c) {
            prod *= keep[c];
            R[c] = (double) prod;
        }
    } else {
        R[0] = keep[0];
    }
}

// The family at one setting: the density, the distribution function and the
// pieces they share. kappa is the concentration, which vonmises2 supplies as
// the inverse Bessel ratio of its rho.
struct Vm1Fam {
    double lb, ub;
    bool fail;
    double mu, k, log_i0;
    int m;
    std::vector<double> R;
    bool have_R;

    Vm1Fam(double mu_, double k_)
        : lb(-M_PI), ub(M_PI), fail(false), mu(mu_), k(k_), m(0), have_R(false) {
        log_i0 = vx_lbi()(k, 0);
    }

    // distrib_pdf(log = TRUE) (R/vonmises1_distrib.R:145-160)
    double lpdf_one(double y) const {
        double out = k * r_cos(y - mu) - std::log(2 * M_PI) - log_i0;
        if (y < -M_PI || y >= M_PI) out = R_NegInf;
        return out;
    }
    double lpdf(double y) { return lpdf_one(y); }
    double pdf(double y) { return r_exp(lpdf_one(y)); }

    // vm_cdf() for one point (:989-1017)
    double vm_cdf(double y) {
        double out = NA_REAL;
        bool ok = R_FINITE(y) && R_FINITE(mu) && R_FINITE(k) && k > 0;
        if (!ISNAN(y) && y < -M_PI) out = 0;
        if (!ISNAN(y) && y >= M_PI) out = 1;
        ok = ok && y >= -M_PI && y < M_PI;
        if (!ok) return r_pmin(r_pmax(out, 0.0), 1.0);
        if (!have_R) {
            m = std::max(20, (int) std::ceil(8.5 * std::sqrt(k)) + 10);
            vm1_ratios(k, m, R);
            have_R = true;
        }
        // outer() with "*" is a matrix product of one factor each
        long double s = 0.0;
        for (int c = 1; c <= m; ++c) {
            const double cj = (double) c;
            const double a1 = r_sin((y - mu) * cj);
            const double a2 = r_sin((M_PI + mu) * cj);
            s += (R[c - 1] * (a1 + a2)) / cj;
        }
        const double sd = (double) s;
        out = (y + M_PI) / (2 * M_PI) + sd / M_PI;
        return r_pmin(r_pmax(out, 0.0), 1.0);
    }

    // vm_log_piece() (:1074-1095)
    double vm_log_piece(double e, double dir, double L) {
        const double c0 = r_cos(e - mu);
        double a1 = k * std::fabs(r_sin(e - mu)), a2 = std::sqrt(k), a3 = 1e-300;
        double mx = (ISNAN(a1) || ISNAN(a2)) ? R_NaN : std::max(a1, a2);
        mx = ISNAN(mx) ? mx : std::max(mx, a3);
        const double h = 1 / mx;
        auto f = [&](double t) {
            const double dl = dir * h * t;
            return r_exp(-2 * k * r_sin(e - mu + dl / 2) * r_sin(dl / 2));
        };
        const double len = L / h;
        double tot = 0;
        double lo = 0;
        double hi = ISNAN(len) ? len : (len < 1 ? len : 1.0);
        for (;;) {
            double val = 0;
            if (r_integrate_fin(f, lo, hi, 1e-17 * tot, 1e-13, 2000, &val)) {
                fail = true;
                return R_NaN;
            }
            tot = tot + val;
            // `hi >= len || f(hi) * (len - hi) < 1e-17 * tot`: an NA in either
            // test is an error in R
            if (ISNAN(hi) || ISNAN(len)) {
                fail = true;
                return R_NaN;
            }
            if (hi >= len) break;
            const double t1 = f(hi) * (len - hi);
            if (ISNAN(t1) || ISNAN(1e-17 * tot)) {
                fail = true;
                return R_NaN;
            }
            if (t1 < 1e-17 * tot) break;
            lo = hi;
            hi = std::min(4 * hi, len);
        }
        return k * c0 + std::log(h * tot);
    }

    // vm_log_tail() (:1050-1067)
    double vm_log_tail(double q, bool lower) {
        const double a = lower ? -M_PI : q;
        const double b = lower ? q : M_PI;
        if (b <= a) return R_NegInf;
        const double anti = mu > 0 ? mu - M_PI : mu + M_PI;
        double lp[2];
        int nl;
        if (anti > a && anti < b) {
            lp[0] = vm_log_piece(a, 1, anti - a);
            lp[1] = vm_log_piece(b, -1, b - anti);
            nl = 2;
        } else if (r_cos(a - mu) >= r_cos(b - mu)) {
            lp[0] = vm_log_piece(a, 1, b - a);
            nl = 1;
        } else {
            lp[0] = vm_log_piece(b, -1, b - a);
            nl = 1;
        }
        double mx = lp[0];
        bool nan = ISNAN(lp[0]);
        if (nl == 2) { nan = nan || ISNAN(lp[1]); if (lp[1] > mx) mx = lp[1]; }
        if (nan) mx = R_NaN;
        long double sm = 0.0;
        for (int i = 0; i < nl; ++i) sm += r_exp(lp[i] - mx);
        return mx + std::log((double) sm) - vx_lbi()(k, 0) - std::log(2 * M_PI);
    }

    // vm_cdf_tail(q, mu, kappa, lower.tail, log.p) for one point (:1098-1111)
    double vm_cdf_tail(double q, bool lower, bool logp) {
        const double Fl = vm_cdf(q);
        double P = lower ? Fl : 1 - Fl;
        double lP = lower ? std::log(Fl) : r_log1p(-Fl);
        const bool small = !ISNAN(P) && P < 1e-6 && q > -M_PI && q < M_PI && k > 0;
        if (small) lP = vm_log_tail(q, lower);
        if (logp) return lP;
        if (small) P = r_exp(lP);
        return P;
    }

    double lcdf(double x, bool lower) { return vm_cdf_tail(x, lower, true); }
    double cdf(double x) { return vm_cdf_tail(x, true, false); }
};

// rng_ok: the generator has no case in which R signals an error.
inline int vonmises1_rng_ok(const double* th) {
    (void) th;
    return 0;
}

// The Best-Fisher generator of distrib_rng(VonMises1Distrib)
// (R/vonmises1_distrib.R:217-267), scalar kappa: a flat envelope draws
// runif(n, -pi, pi); otherwise rounds of three runif(m) vectors with
// m = 2 (n - accepted) + 16, accepted draws appended in index order (an NA
// ok appends NA, as ang[ok] does), the first n kept, then the wrap by %%.
inline void vm1_rng_kappa(double mu, double k, int n, double* out) {
    const double s = std::sqrt(1 + 4 * k * k);
    const double a = 1 + s;
    const double b = 2 * a * k / ((s + 1) * (a + std::sqrt(2 * a)));
    const double r = (1 + b * b) / (2 * b);
    const bool flat = !R_FINITE(r);
    std::vector<double> o;
    if (flat) {
        o.resize(n);
        for (int i = 0; i < n; ++i) o[i] = R::runif(-M_PI, M_PI);
    } else {
        while ((int) o.size() < n) {
            const int m = 2 * (n - (int) o.size()) + 16;
            std::vector<double> u1(m), u2(m), u3(m);
            for (int j = 0; j < m; ++j) u1[j] = R::runif(0, 1);
            for (int j = 0; j < m; ++j) u2[j] = R::runif(0, 1);
            for (int j = 0; j < m; ++j) u3[j] = R::runif(0, 1);
            for (int j = 0; j < m; ++j) {
                const double z = r_cos(M_PI * u1[j]);
                const double f = (1 + r * z) / (r + z);
                const double c0 = k * (r - f);
                // A | B in R's three-valued logic
                const double ta = c0 * (2 - c0) - u2[j];
                const double tb = std::log(c0 / u2[j]) + 1 - c0;
                const bool a_na = ISNAN(ta), b_na = ISNAN(tb);
                const bool a_t = !a_na && ta > 0, b_t = !b_na && tb >= 0;
                if (a_t || b_t) {
                    o.push_back(r_sign(u3[j] - 0.5) *
                                std::acos(r_pmin(r_pmax(f, -1.0), 1.0)));
                } else if (a_na || b_na) {
                    o.push_back(NA_REAL);
                }
            }
        }
    }
    for (int i = 0; i < n; ++i)
        out[i] = r_fmod(o[i] + mu + M_PI, 2 * M_PI) - M_PI;
}

inline void vonmises1_rng(const double* th, int n, double* out) {
    vm1_rng_kappa(th[0], th[1], n, out);
}

// distrib_pdf(log = FALSE) (R/vonmises1_distrib.R:145-160)
inline void vm1_pdf_kappa(const double* y, int n, double mu, double k, double* out) {
    const double log_i0 = vx_lbi()(k, 0);
    for (int i = 0; i < n; ++i) {
        double o = k * r_cos(y[i] - mu) - std::log(2 * M_PI) - log_i0;
        if (y[i] < -M_PI || y[i] >= M_PI) o = R_NegInf;
        out[i] = r_exp(o);
    }
}

inline void vonmises1_pdf(const double* y, int n, const double* th, int vb,
                          double* out) {
    (void) vb;
    vm1_pdf_kappa(y, n, th[0], th[1], out);
}

// distrib_quantile() by the numerical fallback (R/numerical_fallbacks.R:268)
inline int vonmises1_quantile(const double* p, int np, const double* th, int vb,
                              double* out) {
    (void) vb;
    Vm1Fam F(th[0], th[1]);
    return qf_quantile(p, np, F, out);
}

}  // namespace d7x

#endif
