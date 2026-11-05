// The generator, the density (log = FALSE) and the quantile of the 28
// families whose R methods are base R or Rmath calls and closed forms, as
// pt_expect_common.h describes them. Each function cites the R method it
// transcribes; the argument transformations are R's, operation for
// operation (stats::rgamma(rate = r) passes the scale 1/r, stats::rexp(rate)
// passes 1/rate, x^2 is x * x).

#ifndef D7_PT_EXPECT_SIMPLE_H
#define D7_PT_EXPECT_SIMPLE_H

#include "pt_expect_common.h"
#include "pt_logpdf.h"

namespace d7x {

// ---------------------------------------------------------------- helpers
// R/utility_functions.R:536 log1mexp_neg()
inline double log1mexp_neg(double a) {
    double out = r_log1p(-r_exp(a));
    if (!ISNAN(a) && a > -std::log(2.0)) out = std::log(-r_expm1(a));
    return out;
}
// R/laplace_distrib.R:249 laplace_std_quantile(p, TRUE, FALSE), with
// quantile_log_tails() of R/utility_functions.R:620
inline double laplace_std_q(double p) {
    double lp = std::log(p);
    if (!ISNAN(lp) && lp > 0) lp = R_NaN;
    const double lq = log1mexp_neg(lp);
    double out = -(lq + std::log(2.0));
    if (!ISNAN(lp) && lp <= lq) out = lp + std::log(2.0);
    return out;
}
// R/log_tail.R:510 tail_quantile_refine() leaves x unchanged unless
// log(p) < log_tail_cut = -5; there the caller hands back
inline bool refine_noop(double p) {
    const double lt = std::log(p);
    return !(!ISNAN(lt) && lt < -5);
}
inline int no_quantile(const double*, int, const double*, int, double*) { return 1; }

#define D7X_RNG_ALWAYS(fam) inline int fam##_rng_ok(const double*) { return 0; }

// ---------------------------------------------------------- gaussian1 (0)
D7X_RNG_ALWAYS(gaussian1)
// R/gaussian1_distrib.R distrib_rng: stats::rnorm(n, theta[[1]], theta[[2]])
inline void gaussian1_rng(const double* th, int n, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::rnorm(th[0], th[1]);
}
// distrib_pdf: stats::dnorm(y, theta[[1]], theta[[2]])
inline void gaussian1_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::dnorm(y[i], th[0], th[1], 0);
}
// distrib_quantile: stats::qnorm(p, theta[[1]], theta[[2]])
inline int gaussian1_quantile(const double* p, int np, const double* th, int, double* out) {
    for (int i = 0; i < np; ++i) out[i] = R::qnorm(p[i], th[0], th[1], 1, 0);
    return 0;
}

// ------------------------------------------------------------- gamma1 (1)
// gamma1_shape_rate(): shape = 1/phi, rate = 1/(phi*mu); scale = 1/rate
D7X_RNG_ALWAYS(gamma1)
inline void gamma1_rng(const double* th, int n, double* out) {
    const double shape = 1 / th[1], rate = 1 / (th[1] * th[0]), scale = 1 / rate;
    for (int i = 0; i < n; ++i) out[i] = R::rgamma(shape, scale);
}
inline void gamma1_pdf(const double* y, int n, const double* th, int, double* out) {
    const double shape = 1 / th[1], rate = 1 / (th[1] * th[0]), scale = 1 / rate;
    for (int i = 0; i < n; ++i) out[i] = R::dgamma(y[i], shape, scale, 0);
}
inline int gamma1_quantile(const double* p, int np, const double* th, int, double* out) {
    const double shape = 1 / th[1], rate = 1 / (th[1] * th[0]), scale = 1 / rate;
    for (int i = 0; i < np; ++i) out[i] = R::qgamma(p[i], shape, scale, 1, 0);
    return 0;
}

// ------------------------------------------------------------ poisson (2)
D7X_RNG_ALWAYS(poisson)
inline void poisson_rng(const double* th, int n, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::rpois(th[0]);
}
inline void poisson_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::dpois(y[i], th[0], 0);
}
inline int poisson_quantile(const double* p, int np, const double* th, int vb, double* out) {
    return no_quantile(p, np, th, vb, out);
}

// ------------------------------------------------------------ negbin2 (3)
// stats::rnbinom(n, mu = theta[[1]], size = theta[[2]])
D7X_RNG_ALWAYS(negbin2)
inline void negbin2_rng(const double* th, int n, double* out) {
    for (int i = 0; i < n; ++i) out[i] = ::rnbinom_mu(th[1], th[0]);
}
inline void negbin2_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::dnbinom_mu(y[i], th[1], th[0], 0);
}
inline int negbin2_quantile(const double* p, int np, const double* th, int vb, double* out) {
    return no_quantile(p, np, th, vb, out);
}

// -------------------------------------------------------------- beta1 (4)
D7X_RNG_ALWAYS(beta1)
inline void beta1_rng(const double* th, int n, double* out) {
    const double a = th[0] * th[1], b = (1 - th[0]) * th[1];
    for (int i = 0; i < n; ++i) out[i] = R::rbeta(a, b);
}
inline void beta1_pdf(const double* y, int n, const double* th, int, double* out) {
    const double a = th[0] * th[1], b = (1 - th[0]) * th[1];
    for (int i = 0; i < n; ++i) out[i] = R::dbeta(y[i], a, b, 0);
}
inline int beta1_quantile(const double* p, int np, const double* th, int, double* out) {
    const double a = th[0] * th[1], b = (1 - th[0]) * th[1];
    for (int i = 0; i < np; ++i) out[i] = R::qbeta(p[i], a, b, 1, 0);
    return 0;
}

// ---------------------------------------------------------- bernoulli (5)
D7X_RNG_ALWAYS(bernoulli)
inline void bernoulli_rng(const double* th, int n, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::rbinom(1, th[0]);
}
inline void bernoulli_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::dbinom(y[i], 1, th[0], 0);
}
inline int bernoulli_quantile(const double* p, int np, const double* th, int vb, double* out) {
    return no_quantile(p, np, th, vb, out);
}

// ----------------------------------------------------------- binomial (6)
// th[1] is the size of this setting
D7X_RNG_ALWAYS(binomial)
inline void binomial_rng(const double* th, int n, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::rbinom(th[1], th[0]);
}
inline void binomial_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::dbinom(y[i], th[1], th[0], 0);
}
inline int binomial_quantile(const double* p, int np, const double* th, int vb, double* out) {
    return no_quantile(p, np, th, vb, out);
}

// -------------------------------------------------------- exponential (7)
// stats::rexp(n, rate = 1/theta[[1]]) passes the scale 1/rate
D7X_RNG_ALWAYS(exponential)
inline void exponential_rng(const double* th, int n, double* out) {
    const double rate = 1 / th[0], scale = 1 / rate;
    for (int i = 0; i < n; ++i) out[i] = R::rexp(scale);
}
inline void exponential_pdf(const double* y, int n, const double* th, int, double* out) {
    const double rate = 1 / th[0], scale = 1 / rate;
    for (int i = 0; i < n; ++i) out[i] = R::dexp(y[i], scale, 0);
}
inline int exponential_quantile(const double* p, int np, const double* th, int, double* out) {
    const double rate = 1 / th[0], scale = 1 / rate;
    for (int i = 0; i < np; ++i) out[i] = R::qexp(p[i], scale, 1, 0);
    return 0;
}

// ---------------------------------------------------------- geometric (8)
// geom_prob(mu) = 1/(1 + mu)
D7X_RNG_ALWAYS(geometric)
inline void geometric_rng(const double* th, int n, double* out) {
    const double pr = 1 / (1 + th[0]);
    for (int i = 0; i < n; ++i) out[i] = R::rgeom(pr);
}
inline void geometric_pdf(const double* y, int n, const double* th, int, double* out) {
    const double pr = 1 / (1 + th[0]);
    for (int i = 0; i < n; ++i) out[i] = R::dgeom(y[i], pr, 0);
}
inline int geometric_quantile(const double* p, int np, const double* th, int vb, double* out) {
    return no_quantile(p, np, th, vb, out);
}

// -------------------------------------------------------------- chisq (9)
D7X_RNG_ALWAYS(chisq)
inline void chisq_rng(const double* th, int n, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::rchisq(th[0]);
}
inline void chisq_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::dchisq(y[i], th[0], 0);
}
inline int chisq_quantile(const double* p, int np, const double* th, int, double* out) {
    for (int i = 0; i < np; ++i) out[i] = R::qchisq(p[i], th[0], 1, 0);
    return 0;
}

// ------------------------------------------------------------ cauchy (10)
D7X_RNG_ALWAYS(cauchy)
inline void cauchy_rng(const double* th, int n, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::rcauchy(th[0], th[1]);
}
// R/cauchy_distrib.R distrib_pdf: the log density by formula, then exp()
inline void cauchy_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) {
        const double z = (y[i] - th[0]) / th[1];
        const double az = std::fabs(z);
        const double l1a = 2 * std::log(az) + r_log1p(1 / (az * az));
        const double l1b = r_log1p(az * az);
        const double l1 = az > 1e150 ? l1a : l1b;
        const double log_d = -std::log(M_PI) - std::log(th[1]) - l1;
        out[i] = r_exp(log_d);
    }
}
inline int cauchy_quantile(const double* p, int np, const double* th, int, double* out) {
    for (int i = 0; i < np; ++i) out[i] = R::qcauchy(p[i], th[0], th[1], 1, 0);
    return 0;
}

// ---------------------------------------------------------- logistic (11)
D7X_RNG_ALWAYS(logistic)
inline void logistic_rng(const double* th, int n, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::rlogis(th[0], th[1]);
}
inline void logistic_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::dlogis(y[i], th[0], th[1], 0);
}
inline int logistic_quantile(const double* p, int np, const double* th, int, double* out) {
    for (int i = 0; i < np; ++i) out[i] = R::qlogis(p[i], th[0], th[1], 1, 0);
    return 0;
}

// --------------------------------------------------------- gaussian2 (12)
D7X_RNG_ALWAYS(gaussian2)
inline void gaussian2_rng(const double* th, int n, double* out) {
    const double sd = std::sqrt(th[1]);
    for (int i = 0; i < n; ++i) out[i] = R::rnorm(th[0], sd);
}
inline void gaussian2_pdf(const double* y, int n, const double* th, int, double* out) {
    const double sd = std::sqrt(th[1]);
    for (int i = 0; i < n; ++i) out[i] = R::dnorm(y[i], th[0], sd, 0);
}
inline int gaussian2_quantile(const double* p, int np, const double* th, int, double* out) {
    const double sd = std::sqrt(th[1]);
    for (int i = 0; i < np; ++i) out[i] = R::qnorm(p[i], th[0], sd, 1, 0);
    return 0;
}

// --------------------------------------------------------- gaussian3 (13)
D7X_RNG_ALWAYS(gaussian3)
inline void gaussian3_rng(const double* th, int n, double* out) {
    const double sd = 1 / std::sqrt(th[1]);
    for (int i = 0; i < n; ++i) out[i] = R::rnorm(th[0], sd);
}
inline void gaussian3_pdf(const double* y, int n, const double* th, int, double* out) {
    const double sd = 1 / std::sqrt(th[1]);
    for (int i = 0; i < n; ++i) out[i] = R::dnorm(y[i], th[0], sd, 0);
}
inline int gaussian3_quantile(const double* p, int np, const double* th, int, double* out) {
    const double sd = 1 / std::sqrt(th[1]);
    for (int i = 0; i < np; ++i) out[i] = R::qnorm(p[i], th[0], sd, 1, 0);
    return 0;
}

// -------------------------------------------------------- lognormal1 (14)
D7X_RNG_ALWAYS(lognormal1)
inline void lognormal1_rng(const double* th, int n, double* out) {
    const double sl = std::sqrt(th[1]);
    for (int i = 0; i < n; ++i) out[i] = R::rlnorm(th[0], sl);
}
inline void lognormal1_pdf(const double* y, int n, const double* th, int, double* out) {
    const double sl = std::sqrt(th[1]);
    for (int i = 0; i < n; ++i) out[i] = R::dlnorm(y[i], th[0], sl, 0);
}
inline int lognormal1_quantile(const double* p, int np, const double* th, int, double* out) {
    const double sl = std::sqrt(th[1]);
    for (int i = 0; i < np; ++i) out[i] = R::qlnorm(p[i], th[0], sl, 1, 0);
    return 0;
}

// ------------------------------------------------------------ gamma2 (17)
// shape = theta[[1]]^2 / theta[[2]], rate = theta[[1]] / theta[[2]]
D7X_RNG_ALWAYS(gamma2)
inline void gamma2_rng(const double* th, int n, double* out) {
    const double shape = (th[0] * th[0]) / th[1], rate = th[0] / th[1], scale = 1 / rate;
    for (int i = 0; i < n; ++i) out[i] = R::rgamma(shape, scale);
}
inline void gamma2_pdf(const double* y, int n, const double* th, int, double* out) {
    const double shape = (th[0] * th[0]) / th[1], rate = th[0] / th[1], scale = 1 / rate;
    for (int i = 0; i < n; ++i) out[i] = R::dgamma(y[i], shape, scale, 0);
}
inline int gamma2_quantile(const double* p, int np, const double* th, int, double* out) {
    const double shape = (th[0] * th[0]) / th[1], rate = th[0] / th[1], scale = 1 / rate;
    for (int i = 0; i < np; ++i) out[i] = R::qgamma(p[i], shape, scale, 1, 0);
    return 0;
}

// --------------------------------------------------------------- gpd (18)
// R/gpd_distrib.R distrib_quantile, lower tail, not log
inline int gpd_quantile(const double* p, int np, const double* th, int, double* out) {
    const double sv = th[0], xv = th[1];
    for (int i = 0; i < np; ++i) {
        double ls = r_log1p(-p[i]);
        if (!ISNAN(ls) && ls > 0) ls = R_NaN;
        double o = (-sv) * ls;
        if (!ISNAN(xv) && xv != 0) o = sv * r_expm1((-xv) * ls) / xv;
        out[i] = o;
    }
    return 0;
}
// distrib_rng: distrib_quantile(distrib, stats::runif(n), theta)
D7X_RNG_ALWAYS(gpd)
inline void gpd_rng(const double* th, int n, double* out) {
    std::vector<double> u(n);
    for (int i = 0; i < n; ++i) u[i] = R::runif(0.0, 1.0);
    gpd_quantile(u.data(), n, th, 0, out);
}
// distrib_pdf: exp(gpd_logpdf_cpp(y, sigma, xi)), whose body is d7::gpd_logpdf
inline void gpd_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) out[i] = r_exp(d7::gpd_logpdf(y[i], th));
}

// -------------------------------------------------------- lognormal2 (22)
// lognormal2_log_params(): s = log1p(v/m^2); meanlog = log(m) - s/2
struct Ln2 { double ml, sl; };
inline Ln2 lognormal2_lp(const double* th) {
    const double m = th[0];
    const double s = r_log1p(th[1] / (m * m));
    return {std::log(m) - s / 2, std::sqrt(s)};
}
D7X_RNG_ALWAYS(lognormal2)
inline void lognormal2_rng(const double* th, int n, double* out) {
    const Ln2 l = lognormal2_lp(th);
    for (int i = 0; i < n; ++i) out[i] = R::rlnorm(l.ml, l.sl);
}
inline void lognormal2_pdf(const double* y, int n, const double* th, int, double* out) {
    const Ln2 l = lognormal2_lp(th);
    for (int i = 0; i < n; ++i) out[i] = R::dlnorm(y[i], l.ml, l.sl, 0);
}
inline int lognormal2_quantile(const double* p, int np, const double* th, int, double* out) {
    const Ln2 l = lognormal2_lp(th);
    for (int i = 0; i < np; ++i) out[i] = R::qlnorm(p[i], l.ml, l.sl, 1, 0);
    return 0;
}

// -------------------------------------------------------- student_t1 (23)
D7X_RNG_ALWAYS(student_t1)
inline void student_t1_rng(const double* th, int n, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::rt(th[2]);
    for (int i = 0; i < n; ++i) out[i] = th[0] + th[1] * out[i];
}
inline void student_t1_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::dt((y[i] - th[0]) / th[1], th[2], 0) / th[1];
}
inline int student_t1_quantile(const double* p, int np, const double* th, int, double* out) {
    for (int i = 0; i < np; ++i) {
        if (!refine_noop(p[i])) return 1;
        out[i] = th[0] + th[1] * R::qt(p[i], th[2], 1, 0);
    }
    return 0;
}

// -------------------------------------------------------- student_t2 (24)
// student_t2_scale() = theta[[2]] * sqrt(1 - 2 / theta[[3]])
inline double student_t2_s0(const double* th) { return th[1] * std::sqrt(1 - 2 / th[2]); }
D7X_RNG_ALWAYS(student_t2)
inline void student_t2_rng(const double* th, int n, double* out) {
    const double s0 = student_t2_s0(th);
    for (int i = 0; i < n; ++i) out[i] = R::rt(th[2]);
    for (int i = 0; i < n; ++i) out[i] = th[0] + s0 * out[i];
}
inline void student_t2_pdf(const double* y, int n, const double* th, int, double* out) {
    const double s0 = student_t2_s0(th);
    for (int i = 0; i < n; ++i) out[i] = R::dt((y[i] - th[0]) / s0, th[2], 0) / s0;
}
inline int student_t2_quantile(const double* p, int np, const double* th, int, double* out) {
    const double s0 = student_t2_s0(th);
    for (int i = 0; i < np; ++i) {
        if (!refine_noop(p[i])) return 1;
        out[i] = th[0] + s0 * R::qt(p[i], th[2], 1, 0);
    }
    return 0;
}

// --------------------------------------------------------- gengamma1 (25)
// a * stats::rgamma(n, shape = d / p, rate = 1)^(1 / p)
D7X_RNG_ALWAYS(gengamma1)
inline void gengamma1_rng(const double* th, int n, double* out) {
    const double a = th[0], d = th[1], p = th[2];
    const double shape = d / p, scale = 1.0 / 1.0, e = 1 / p;
    for (int i = 0; i < n; ++i) out[i] = R::rgamma(shape, scale);
    for (int i = 0; i < n; ++i) out[i] = a * r_pow(out[i], e);
}
inline void gengamma1_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) out[i] = r_exp(d7::gengamma_logd(y[i], th[0], th[1], th[2]));
}
inline int gengamma1_quantile(const double* pr, int np, const double* th, int, double* out) {
    const double a = th[0], d = th[1], pw = th[2];
    for (int i = 0; i < np; ++i) {
        if (!refine_noop(pr[i])) return 1;
        const double g = R::qgamma(pr[i], d / pw, 1.0 / 1.0, 1, 0);
        out[i] = a * r_pow(g, 1 / pw);
    }
    return 0;
}

// --------------------------------------------------------- gengamma2 (26)
// gengamma2_scale() = exp(log(theta[[1]]) + lgamma(d / p) - lgamma((d + 1) / p))
inline double gengamma2_a(const double* th) {
    const double d = th[1], p = th[2];
    return r_exp(std::log(th[0]) + R::lgammafn(d / p) - R::lgammafn((d + 1) / p));
}
D7X_RNG_ALWAYS(gengamma2)
inline void gengamma2_rng(const double* th, int n, double* out) {
    const double pw = th[2];
    const double sc = gengamma2_a(th);
    const double shape = th[1] / pw, e = 1 / pw;
    for (int i = 0; i < n; ++i) out[i] = R::rgamma(shape, 1.0 / 1.0);
    for (int i = 0; i < n; ++i) out[i] = sc * r_pow(out[i], e);
}
inline void gengamma2_pdf(const double* y, int n, const double* th, int, double* out) {
    const double a = gengamma2_a(th);
    for (int i = 0; i < n; ++i) out[i] = r_exp(d7::gengamma_logd(y[i], a, th[1], th[2]));
}
inline int gengamma2_quantile(const double* pr, int np, const double* th, int, double* out) {
    const double pw = th[2];
    const double sc = gengamma2_a(th);
    for (int i = 0; i < np; ++i) {
        if (!refine_noop(pr[i])) return 1;
        out[i] = sc * r_pow(R::qgamma(pr[i], th[1] / pw, 1.0 / 1.0, 1, 0), 1 / pw);
    }
    return 0;
}

// ------------------------------------------------------------ gumbel (27)
// theta[[1]] - theta[[2]] * log(stats::rexp(n))
D7X_RNG_ALWAYS(gumbel)
inline void gumbel_rng(const double* th, int n, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::rexp(1.0 / 1.0);
    for (int i = 0; i < n; ++i) out[i] = th[0] - th[1] * std::log(out[i]);
}
inline void gumbel_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) {
        const double z = (y[i] - th[0]) / th[1];
        out[i] = r_exp(-std::log(th[1]) - z - r_exp(-z));
    }
}
inline int gumbel_quantile(const double* p, int np, const double* th, int, double* out) {
    for (int i = 0; i < np; ++i) {
        const double ll = std::log(-std::log(p[i]));
        out[i] = th[0] - th[1] * ll;
    }
    return 0;
}

// ----------------------------------------------------------- laplace (28)
inline int laplace_quantile(const double* p, int np, const double* th, int, double* out) {
    for (int i = 0; i < np; ++i) out[i] = th[0] + th[1] * laplace_std_q(p[i]);
    return 0;
}
D7X_RNG_ALWAYS(laplace)
inline void laplace_rng(const double* th, int n, double* out) {
    std::vector<double> u(n);
    for (int i = 0; i < n; ++i) u[i] = R::runif(0.0, 1.0);
    laplace_quantile(u.data(), n, th, 0, out);
}
inline void laplace_pdf(const double* y, int n, const double* th, int, double* out) {
    const double mu = th[0], b = th[1];
    for (int i = 0; i < n; ++i)
        out[i] = r_exp(-std::log(2.0) - std::log(b) - std::fabs(y[i] - mu) / b);
}

// ---------------------------------------------------------- laplace2 (29)
inline int laplace2_quantile(const double* p, int np, const double* th, int, double* out) {
    for (int i = 0; i < np; ++i) out[i] = th[0] + laplace_std_q(p[i]) / th[1];
    return 0;
}
D7X_RNG_ALWAYS(laplace2)
inline void laplace2_rng(const double* th, int n, double* out) {
    std::vector<double> u(n);
    for (int i = 0; i < n; ++i) u[i] = R::runif(0.0, 1.0);
    laplace2_quantile(u.data(), n, th, 0, out);
}
inline void laplace2_pdf(const double* y, int n, const double* th, int, double* out) {
    const double mu = th[0], lam = th[1];
    for (int i = 0; i < n; ++i)
        out[i] = r_exp(std::log(lam / 2) - lam * std::fabs(y[i] - mu));
}

// ---------------------------------------------------------- weibull1 (30)
// shape = theta[[2]], scale = theta[[1]]
D7X_RNG_ALWAYS(weibull1)
inline void weibull1_rng(const double* th, int n, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::rweibull(th[1], th[0]);
}
inline void weibull1_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::dweibull(y[i], th[1], th[0], 0);
}
inline int weibull1_quantile(const double* p, int np, const double* th, int, double* out) {
    for (int i = 0; i < np; ++i) out[i] = R::qweibull(p[i], th[1], th[0], 1, 0);
    return 0;
}

// ------------------------------------------------------------- beta2 (31)
D7X_RNG_ALWAYS(beta2)
inline void beta2_rng(const double* th, int n, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::rbeta(th[0], th[1]);
}
inline void beta2_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) out[i] = R::dbeta(y[i], th[0], th[1], 0);
}
inline int beta2_quantile(const double* p, int np, const double* th, int, double* out) {
    for (int i = 0; i < np; ++i) out[i] = R::qbeta(p[i], th[0], th[1], 1, 0);
    return 0;
}

// ----------------------------------------------------------- negbin1 (33)
// stats::rnbinom(n, size = mu / theta, mu = mu)
D7X_RNG_ALWAYS(negbin1)
inline void negbin1_rng(const double* th, int n, double* out) {
    const double size = th[0] / th[1];
    for (int i = 0; i < n; ++i) out[i] = ::rnbinom_mu(size, th[0]);
}
// exp(negbin1_logpmf_cpp(y, mu, theta)): dnbinom_mu(y, m / t, m, log = TRUE)
inline void negbin1_pdf(const double* y, int n, const double* th, int, double* out) {
    for (int i = 0; i < n; ++i) out[i] = r_exp(R::dnbinom_mu(y[i], th[0] / th[1], th[0], 1));
}
inline int negbin1_quantile(const double* p, int np, const double* th, int vb, double* out) {
    return no_quantile(p, np, th, vb, out);
}

#undef D7X_RNG_ALWAYS

}  // namespace d7x

#endif
