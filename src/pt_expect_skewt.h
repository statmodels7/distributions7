// skewt: generator, density and quantile for the compiled expected
// information by simulation and by integration (see pt_expect_common.h).
// th = (mu, sigma, alpha, nu).
//
// The quantile is the numerical fallback of pt_expect_qfallback.h over this
// family's density and distribution function. The distribution function is
// compiled_cdf() (R/cdf_compiled.R:68-80, compiled_log_cdf() :107-122,
// compiled_far_tail() :89-91): the kernel is d7_cdf_cpp() of
// d7_ccallable.cpp, called here through its (external) C++ declaration with
// the class name R passes, and the far tails go through quad_log_tail_prob()
// (R/log_tail.R:403-454) with distrib_grad_y() (R/skewt_distrib.R:1092-1095
// over skewt_pieces() :6-48) and log_tail_nodes() (:458-471).
//
// d7_cdf_cpp is defined in d7_ccallable.cpp (global namespace) and exported
// to R by Rcpp, so any translation unit of the package links to it.

#ifndef D7_PT_EXPECT_SKEWT_H
#define D7_PT_EXPECT_SKEWT_H

#include "pt_expect_common.h"
#include "pt_expect_qfallback.h"

Rcpp::NumericVector d7_cdf_cpp(std::string cls, Rcpp::NumericVector q,
                               Rcpp::NumericMatrix theta, bool lower);

namespace d7x {

// log_tail_nodes() (R/log_tail.R:458-471): u = outer(nodes, h) + rep(a + h,
// each = 15) and w = outer(kronrod_weights, h) over the 16 intervals of br
inline const std::vector<double>& skt_tail_u() {
    static std::vector<double> u;
    if (u.empty()) {
        static const double br[17] = {0, 0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16,
                                      24, 32, 48, 64, 96, 128};
        const n7::GK15& r = n7::gk15();
        for (int k = 0; k < 16; ++k) {
            const double a = br[k];
            const double h = (br[k + 1] - br[k]) / 2;
            for (int j = 0; j < 15; ++j) u.push_back(r.node[j] * h + (a + h));
        }
    }
    return u;
}
inline const std::vector<double>& skt_tail_w() {
    static std::vector<double> w;
    if (w.empty()) {
        static const double br[17] = {0, 0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 16,
                                      24, 32, 48, 64, 96, 128};
        const n7::GK15& r = n7::gk15();
        for (int k = 0; k < 16; ++k) {
            const double h = (br[k + 1] - br[k]) / 2;
            for (int j = 0; j < 15; ++j) w.push_back(r.kw[j] * h);
        }
    }
    return w;
}

// the family at one setting
struct SktFam {
    double lb, ub;
    bool fail;
    double mu, sg, al, nu;

    SktFam(double mu_, double sg_, double al_, double nu_)
        : lb(R_NegInf), ub(R_PosInf), fail(false), mu(mu_), sg(sg_), al(al_), nu(nu_) {}

    // distrib_pdf(log = TRUE) (R/skewt_distrib.R:289-306)
    double lpdf(double y) {
        const double z = (y - mu) / sg;
        double w = al * z * std::sqrt((nu + 1) / (nu + z * z));
        if (std::fabs(z) > 1e20) {
            const double rs = std::fabs(z) * std::sqrt(1 + nu / (z * z));
            w = al * std::sqrt(nu + 1) * (z / rs);
        }
        return std::log(2.0) - std::log(sg) + R::dt(z, nu, 1) +
               R::pt(w, nu + 1, 1, 1);
    }
    double pdf(double y) { return r_exp(lpdf(y)); }

    // distrib_grad_y() = (p$a + p$q * p$b) / sigma, skewt_pieces() (:6-48)
    double grad_y(double y) {
        const double z = (y - mu) / sg;
        const double m = nu + 1;
        const double s = nu + z * z;
        double w = al * z * std::sqrt(m / s);
        double e = nu * std::sqrt(m) / R_POW(s, 1.5);
        double a = -m * z / s;
        if (std::fabs(z) > 1e20) {
            const double rs = std::fabs(z) * std::sqrt(1 + nu / (z * z));
            const double az = z / rs;
            const double irs = 1 / rs;
            const double b2 = nu * irs * irs;
            const double rm = std::sqrt(m);
            w = al * rm * az;
            a = -m * az * irs;
            e = rm * b2 * irs;
        }
        const double q = r_exp(R::dt(w, m, 1) - R::pt(w, m, 1, 1));
        const double b = al * e;
        return (a + q * b) / sg;
    }

    // d7_cdf_cpp("SkewTDistrib", q, theta, lower) for one point
    double kernel(double x, bool lower) {
        Rcpp::NumericVector qv(1);
        qv[0] = x;
        Rcpp::NumericMatrix tm(1, 4);
        tm(0, 0) = mu; tm(0, 1) = sg; tm(0, 2) = al; tm(0, 3) = nu;
        Rcpp::NumericVector o = d7_cdf_cpp("SkewTDistrib", qv, tm, lower);
        return o[0];
    }

    // quad_log_tail_prob(distrib, q, th, lower.tail) for one point
    // (R/log_tail.R:403-454, the infinite-bound branch)
    double quad_log_tail_prob(double q, bool lower) {
        const std::vector<double>& U = skt_tail_u();
        const std::vector<double>& W = skt_tail_w();
        const int K = (int) U.size();
        const double sgn = lower ? -1.0 : 1.0;
        const double gy = sgn * grad_y(q);
        double s = (R_FINITE(gy) && gy < 0) ? -1 / gy : R_PosInf;
        if (!R_FINITE(s)) s = 1;
        const double lfq = lpdf(q);
        long double rs = 0.0;
        for (int k = 0; k < K; ++k) {
            double Y = q + sgn * s * r_expm1(U[k]);
            const double jac = s * r_exp(U[k]);
            const bool ok = R_FINITE(Y) && (lower ? Y < q : Y > q) &&
                            Y > lb && Y < ub;
            if (!ok) Y = q;
            const double lfk = lpdf(Y);
            double Wt = r_exp(lfk - lfq) * jac * W[k];
            if (!ok || !R_FINITE(Wt)) Wt = 0;
            rs += Wt;
        }
        return lfq + std::log((double) rs);
    }

    // compiled_cdf(..., log.p = TRUE) (R/cdf_compiled.R:68-69, :107-122)
    double lcdf(double x, bool lower) {
        const double out = kernel(x, lower);
        double lo = std::log(out);
        if (!ISNAN(out) && out > 0.5) lo = r_log1p(-kernel(x, !lower));
        const bool far = !ISNAN(out) && out < r_exp(-5.0) && R_FINITE(x) &&
                         x > lb && x < ub;
        if (far) lo = quad_log_tail_prob(x, lower);
        return lo;
    }

    // compiled_cdf(..., log.p = FALSE) (:70-79)
    double cdf(double x) {
        double out = kernel(x, true);
        const bool far = !ISNAN(out) && out < r_exp(-5.0) && R_FINITE(x) &&
                         x > lb && x < ub;
        if (far) out = r_exp(quad_log_tail_prob(x, true));
        return out;
    }
};

// rng_ok: the generator has no case in which R signals an error.
inline int skewt_rng_ok(const double* th) {
    (void) th;
    return 0;
}

// distrib_rng(SkewTDistrib) (R/skewt_distrib.R:360-366): n normals, n
// normals, then n chi-squares
inline void skewt_rng(const double* th, int n, double* out) {
    const double alpha = th[2], nu = th[3];
    const double delta = alpha / std::sqrt(1 + alpha * alpha);
    std::vector<double> u(n), v(n), c(n);
    for (int i = 0; i < n; ++i) u[i] = R::rnorm(0.0, 1.0);
    for (int i = 0; i < n; ++i) v[i] = R::rnorm(0.0, 1.0);
    for (int i = 0; i < n; ++i) c[i] = R::rchisq(nu);
    for (int i = 0; i < n; ++i) {
        const double z = delta * std::fabs(u[i]) + std::sqrt(1 - delta * delta) * v[i];
        out[i] = th[0] + th[1] * z / std::sqrt(c[i] / nu);
    }
}

// distrib_pdf(log = FALSE) (R/skewt_distrib.R:289-306)
inline void skewt_pdf(const double* y, int n, const double* th, int vb, double* out) {
    (void) vb;
    SktFam F(th[0], th[1], th[2], th[3]);
    for (int i = 0; i < n; ++i) out[i] = F.pdf(y[i]);
}

inline int skewt_quantile(const double* p, int np, const double* th, int vb,
                          double* out) {
    (void) vb;
    SktFam F(th[0], th[1], th[2], th[3]);
    return qf_quantile(p, np, F, out);
}

}  // namespace d7x

#endif
