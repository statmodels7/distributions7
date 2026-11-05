// The generator, the density and the quantile of the inverse Gaussian family
// in (mean, dispersion) for the compiled expected information, value for
// value as the R methods compute them. The cores below (invgauss_*) are the
// statmod 1.5.2 functions rinvgauss, dinvgauss and qinvgauss for ONE setting
// (scalar mean and dispersion); pt_expect_invgauss2.h, pt_expect_pig1.h and
// pt_expect_pig2.h include this file and reuse them.
//
// Sources transcribed (printed from statmod 1.5.2):
//   statmod::rinvgauss, statmod::dinvgauss, statmod::qinvgauss,
//   statmod:::.pinvgauss, statmod:::.dinvgauss;
//   R/invgauss1_distrib.R:139 (pdf), :271 (quantile), :297
//   (invgauss_quantile()), :353 (rng); R/log_tail.R:510 (tail_quantile_refine).

#ifndef D7_PT_EXPECT_INVGAUSS1_H
#define D7_PT_EXPECT_INVGAUSS1_H

#include "pt_expect_common.h"
#include "pt_sqrt.h"

namespace d7x {

// 1 when mean and dispersion are both finite and positive, the cases of
// statmod::rinvgauss() that draw (R's `mu.ok & phi.ok`)
inline int invgauss_params_ok(double mu, double phi) {
    return (mu > 0 && R_FINITE(mu) && phi > 0 && R_FINITE(phi)) ? 1 : 0;
}

// statmod::rinvgauss(n, mean = mu, dispersion = phi), all elements ok:
//   Y <- rnorm(n)^2   (n draws), then runif(n) (n draws)
// The `any(bigphi)` switch of R is elementwise in value: both branches
// assign 1/Yphi to the elements with Yphi > 5e5 and the same expression to
// the others.
inline void invgauss_rng_core(double mu, double phi, int n, double* out) {
    std::vector<double> Y(n), Yphi(n), X1(n);
    for (int i = 0; i < n; ++i) {
        const double z = R::rnorm(0.0, 1.0);
        Y[i] = z * z;                                   // rnorm(n)^2
    }
    for (int i = 0; i < n; ++i) {
        Yphi[i] = (Y[i] * phi) * mu;                    // Y * phi * mu
        const bool bigphi = Yphi[i] > 5e+05;
        if (bigphi) X1[i] = 1 / Yphi[i];
        else X1[i] = 1 + (Yphi[i] / 2) * (1 - d7::sqrt_cr(1 + 4 / Yphi[i]));
    }
    for (int i = 0; i < n; ++i) {
        const double u = R::runif(0.0, 1.0);
        const bool firstroot = (u < 1 / (1 + X1[i]));
        const double r = firstroot ? X1[i] : 1 / X1[i];
        out[i] = mu * r;
    }
}

// statmod:::.pinvgauss(q, dispersion = phi, lower.tail = TRUE, log.p = TRUE),
// with mean = NULL (statmod .pinvgauss, lower.tail branch)
inline double invgauss_logp(double q, double phi) {
    const double pq = d7::sqrt_cr(phi * q);
    const double a = R::pnorm((q - 1) / pq, 0.0, 1.0, 1, 1);
    double b = 2 / phi + R::pnorm(-(q + 1) / pq, 0.0, 1.0, 1, 1);
    b = r_exp(b - a);
    return a + r_log1p(b);
}

// statmod:::.dinvgauss(x, dispersion = phi, log = FALSE), mean = NULL
inline double invgauss_dens(double x, double phi) {
    double t = (-std::log(phi)) - std::log(2 * M_PI);
    t = t - 3 * std::log(x);
    t = t - (((x - 1) * (x - 1)) / phi) / x;
    return r_exp(t / 2);
}

// statmod::dinvgauss(x, mean = mu, dispersion = phi, log = FALSE) at one
// point, parameters finite and positive. R evaluates .dinvgauss(log = TRUE)
// and then exp(); the special values of x (statmod's `any.special` branch,
// elementwise in value) give 0 for x <= 0 and x = Inf and NA for NaN.
inline double invgauss_pdf_one(double x, double mu, double phi) {
    if (ISNAN(x)) return NA_REAL;
    if (x < 0 || x == R_PosInf || (x == 0 && mu > 0 && phi < R_PosInf))
        return r_exp(R_NegInf);
    const double xs = x / mu;
    const double ph = phi * mu;
    double t = (-std::log(ph)) - std::log(2 * M_PI);
    t = t - 3 * std::log(xs);
    t = t - (((xs - 1) * (xs - 1)) / ph) / xs;
    double d = t / 2;
    d = d - std::log(mu);
    return r_exp(d);
}

// statmod::qinvgauss(p, mean = mu, dispersion = phi) for ONE p in (0, 1)
// with log(p) >= -5 (the points at which tail_quantile_refine() is a no-op),
// lower.tail = TRUE, log.p = FALSE, maxit = 200, tol = 1e-14; mu, phi finite
// and positive. Returns nonzero where R would take a path that is not
// reproduced (a NaN reaching the convergence test).
inline int invgauss_q_one(double p0, double mu, double phi0, double* out) {
    const double tol = 1e-14;
    const int maxit = 200;
    const double logp = std::log(p0);
    const double p = r_exp(logp);
    const double cv2 = phi0 * mu;
    if (!R_FINITE(cv2)) return 1;
    if (cv2 < 1e-08) {
        // qgamma(logp, shape = 1/cv2, scale = cv2 * mu, lower.tail, log.p = TRUE)
        *out = R::qgamma(logp, 1 / cv2, cv2 * mu, 1, 1);
        return 0;
    }
    const double phi = cv2;
    const double kappa = 1.5 * phi;
    double x = d7::sqrt_cr(1 + kappa * kappa) - kappa;
    if (kappa > 1000) {
        const double k1 = (1.0 / 2.0) / kappa;          // 1/2/kappa
        x = k1 * (1 - k1 * k1);
    }
    // small.left (logp < -11.51) cannot occur at log(p) >= -5
    if (logp > -1e-05) {
        const double alpha = 1 / phi;
        // qgamma(logp, shape = alpha, rate = alpha, ...): scale = 1/rate
        const double qgam = R::qgamma(logp, alpha, 1 / alpha, 1, 1);
        x = n7::r_pmax(x, qgam);
    }
    auto step = [&](double xx) -> double {
        const double logF = invgauss_logp(xx, phi);
        double dp, dlogp;
        dp = dlogp = logp - logF;
        if (ISNAN(dlogp)) return R_NaN;                 // see the report
        if (std::fabs(dlogp) < 1e-05)
            dp = r_exp(logp + r_log1p(-dlogp / 2)) * dlogp;
        else
            dp = p - r_exp(logF);
        return dp / invgauss_dens(xx, phi);
    };
    double dx = step(x);
    if (ISNAN(dx)) dx = 0;
    const double sdx = r_sign(dx);
    x = x + dx;
    bool i = std::fabs(dx) > tol;
    int iter = 0;
    while (i) {
        iter = iter + 1;
        if (iter > maxit) break;
        dx = step(x);
        if (ISNAN(dx) || dx * sdx < 0) dx = 0;
        x = x + dx;
        const double ratio = std::fabs(dx) / n7::r_pmax(x, 1.0);
        if (ISNAN(ratio)) return 2;
        i = ratio > tol;
    }
    *out = x * mu;
    return 0;
}

// invgauss_quantile(): statmod::qinvgauss on the non-far elements, then
// tail_quantile_refine(), which returns x unchanged when no element has
// log(p) < -5 (R/log_tail.R:515). Nonzero for any p not in (0, 1), for NaN
// and for log(p) < -5.
inline int invgauss_quantile_core(const double* p, int np, double mu,
                                  double phi, double* out) {
    if (!invgauss_params_ok(mu, phi)) return 1;
    for (int k = 0; k < np; ++k) {
        const double lt = std::log(p[k]);
        if (ISNAN(p[k]) || !(p[k] < 1) || !(lt >= -5.0)) return 1;
    }
    for (int k = 0; k < np; ++k)
        if (invgauss_q_one(p[k], mu, phi, &out[k])) return 1;
    return 0;
}

// th = (mu, phi)

// R/invgauss1_distrib.R:353: statmod::rinvgauss(n, mean = mu, dispersion = phi)
inline int invgauss1_rng_ok(const double* th) {
    return invgauss_params_ok(th[0], th[1]) ? 0 : 1;
}
inline void invgauss1_rng(const double* th, int n, double* out) {
    invgauss_rng_core(th[0], th[1], n, out);
}

// R/invgauss1_distrib.R:139: statmod::dinvgauss(y, mean = mu, dispersion =
// phi, log = FALSE). Parameters finite and positive (the caller checks them
// through invgauss1_rng_ok()).
inline void invgauss1_pdf(const double* y, int n, const double* th, int vb,
                          double* out) {
    (void) vb;
    for (int i = 0; i < n; ++i) out[i] = invgauss_pdf_one(y[i], th[0], th[1]);
}

// R/invgauss1_distrib.R:271 and :297
inline int invgauss1_quantile(const double* p, int np, const double* th,
                              int vb, double* out) {
    (void) vb;
    return invgauss_quantile_core(p, np, th[0], th[1], out);
}

}  // namespace d7x

#endif
