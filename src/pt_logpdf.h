#ifndef D7_PT_LOGPDF_H
#define D7_PT_LOGPDF_H

#include <Rcpp.h>
#include <cmath>
#include "pt_gpd.h"
#include "pt_vonmises.h"
#include "pt_skewnormal1.h"
#include "pt_skewnormal2.h"
#include "pt_skewt.h"
#include "pt_pseudohuber.h"
#include "pt_pseudohuber2.h"
#include "pt_pig.h"
#include "pt_betabinom.h"

// The log-density of one observation, one function per family, th being the
// family's parameters followed by its constants as the scalar registry
// reads them. Each writes the expression its distrib_pdf() method evaluates,
// with the same nmath functions where that method calls them, so that the
// two agree to the last bit; where the method is R code (the inverse
// Gaussian through statmod, the elastic net, the von Mises, the centered
// skew normal, the shape-parametrized beta-binomial) the expression is
// transcribed. A discrete family tests its support before calling nmath:
// nmath warns at a non-integer count, and a warning raised from a worker
// thread terminates the process.

namespace d7 {

inline bool lp_count(double y) {
    return R_FINITE(y) && y >= 0.0 && y == std::floor(y);
}

// log I_0, numericals7's n7_log_bessel_i: resolved by vm_log_i0_fn()
// (d7_ccallable.cpp) on the calling thread before any worker starts
typedef double (*N7LogBesselI)(double, double);
N7LogBesselI vm_log_i0_fn();

inline double gaussian1_logpdf(double y, const double* th) {
    return R::dnorm4(y, th[0], th[1], 1);
}
inline double gamma1_logpdf(double y, const double* th) {
    double shape = 1.0 / th[1], rate = 1.0 / (th[1] * th[0]);
    return R::dgamma(y, shape, 1.0 / rate, 1);
}
inline double poisson_logpdf(double y, const double* th) {
    if (!lp_count(y)) return R_NegInf;
    return R::dpois(y, th[0], 1);
}
inline double negbin2_logpdf(double y, const double* th) {
    if (!lp_count(y)) return R_NegInf;
    return Rf_dnbinom_mu(y, th[1], th[0], 1);
}
inline double beta1_logpdf(double y, const double* th) {
    return R::dbeta(y, th[0] * th[1], (1.0 - th[0]) * th[1], 1);
}
inline double bernoulli_logpdf(double y, const double* th) {
    if (!lp_count(y)) return R_NegInf;
    return R::dbinom(y, 1.0, th[0], 1);
}
inline double binomial_logpdf(double y, const double* th) {
    if (!lp_count(y)) return R_NegInf;
    return R::dbinom(y, th[1], th[0], 1);
}
inline double exponential_logpdf(double y, const double* th) {
    double rate = 1.0 / th[0];
    return R::dexp(y, 1.0 / rate, 1);
}
inline double geometric_logpdf(double y, const double* th) {
    if (!lp_count(y)) return R_NegInf;
    return R::dgeom(y, 1.0 / (1.0 + th[0]), 1);
}
inline double chisq_logpdf(double y, const double* th) {
    return R::dchisq(y, th[0], 1);
}
inline double cauchy_logpdf(double y, const double* th) {
    return R::dcauchy(y, th[0], th[1], 1);
}
inline double logistic_logpdf(double y, const double* th) {
    return R::dlogis(y, th[0], th[1], 1);
}
inline double gaussian2_logpdf(double y, const double* th) {
    return R::dnorm4(y, th[0], std::sqrt(th[1]), 1);
}
inline double gaussian3_logpdf(double y, const double* th) {
    return R::dnorm4(y, th[0], 1.0 / std::sqrt(th[1]), 1);
}
inline double lognormal1_logpdf(double y, const double* th) {
    return R::dlnorm(y, th[0], std::sqrt(th[1]), 1);
}

// statmod's .dinvgauss() at mean m and dispersion phi
inline double invgauss_logd(double y, double m, double phi) {
    if (ISNAN(y) || ISNAN(m) || ISNAN(phi)) return NA_REAL;
    if (!(y > 0.0) || y == R_PosInf) return R_NegInf;
    double x = y / m, d = phi * m;
    double v = (-std::log(d) - std::log(2.0 * M_PI) - 3.0 * std::log(x) -
                (x - 1.0) * (x - 1.0) / d / x) / 2.0;
    return v - std::log(m);
}
inline double invgauss1_logpdf(double y, const double* th) {
    return invgauss_logd(y, th[0], th[1]);
}
inline double invgauss2_logpdf(double y, const double* th) {
    return invgauss_logd(y, th[0], 1.0 / th[1]);
}
inline double gamma2_logpdf(double y, const double* th) {
    double shape = th[0] * th[0] / th[1], rate = th[0] / th[1];
    return R::dgamma(y, shape, 1.0 / rate, 1);
}
inline double gpd_logpdf(double y, const double* th) {
    const double sv = th[0], x = th[1];
    const double t = gpd_t(y, sv, x);
    if (y < 0.0 || !(t > 0.0)) return R_NegInf;
    const double z = y / sv, u = x * z;
    const double lt = gpd_logt(u, t);
    return -std::log(sv) - lt - z * gpd_phi(0, u, u / t, t, lt);
}
inline double vonmises_logd(double y, double mu, double k) {
    if (y < -M_PI || y >= M_PI) return R_NegInf;
    double log_i0 = vm_log_i0_fn()(k, 0.0);
    return k * std::cos(y - mu) - std::log(2.0 * M_PI) - log_i0;
}
inline double vonmises1_logpdf(double y, const double* th) {
    return vonmises_logd(y, th[0], th[1]);
}
inline double vonmises2_logpdf(double y, const double* th) {
    return vonmises_logd(y, th[0], vm_bessel().inv(th[1]));
}
inline double weibull3_logpdf(double y, const double* th) {
    double scale = std::exp(std::log(th[0]) - R::lgammafn(1.0 + 1.0 / th[1]));
    return R::dweibull(y, th[1], scale, 1);
}
inline double lognormal2_logpdf(double y, const double* th) {
    double m = th[0];
    double s = std::log1p(th[1] / (m * m));
    return R::dlnorm(y, std::log(m) - s / 2.0, std::sqrt(s), 1);
}
inline double student_t1_logpdf(double y, const double* th) {
    return R::dt((y - th[0]) / th[1], th[2], 1) - std::log(th[1]);
}
inline double student_t2_logpdf(double y, const double* th) {
    double s0 = th[1] * std::sqrt(1.0 - 2.0 / th[2]);
    return R::dt((y - th[0]) / s0, th[2], 1) - std::log(s0);
}
// gengamma_logpdf_cpp()'s body at scale a, shape d and power p
inline double gengamma_logd(double y, double a, double d, double p) {
    if (y <= 0) return R_NegInf;
    double L = std::log(y) - std::log(a);
    double w = std::exp(p * L);
    double k = d / p;
    return std::log(p) - d * std::log(a) - R::lgammafn(k)
        + (d - 1.0) * std::log(y) - w;
}
inline double gengamma1_logpdf(double y, const double* th) {
    return gengamma_logd(y, th[0], th[1], th[2]);
}
inline double gengamma2_logpdf(double y, const double* th) {
    double d = th[1], p = th[2];
    double a = std::exp(std::log(th[0]) + R::lgammafn(d / p) -
                        R::lgammafn((d + 1.0) / p));
    return gengamma_logd(y, a, d, p);
}
inline double gumbel_logpdf(double y, const double* th) {
    double z = (y - th[0]) / th[1];
    return -std::log(th[1]) - z - std::exp(-z);
}
inline double laplace_logpdf(double y, const double* th) {
    return -std::log(2.0 * th[1]) - std::fabs(y - th[0]) / th[1];
}
inline double laplace2_logpdf(double y, const double* th) {
    return std::log(th[1] / 2.0) - th[1] * std::fabs(y - th[0]);
}
inline double weibull1_logpdf(double y, const double* th) {
    return R::dweibull(y, th[1], th[0], 1);
}
inline double beta2_logpdf(double y, const double* th) {
    return R::dbeta(y, th[0], th[1], 1);
}
// .enet_logM(): log of the Mills-type normalizer at x
inline double enet_logM(double x) {
    if (x > 30.0) {
        double u = 1.0 / x, u2 = u * u;
        return -std::log(x) +
            std::log1p(u2 * (-1.0 + u2 * (3.0 + u2 * (-15.0 + u2 * 105.0))));
    }
    return R::pnorm5(-x, 0.0, 1.0, 1, 1) + x * x / 2.0 +
        0.5 * std::log(2.0 * M_PI);
}
inline double enet_logpdf(double y, const double* th) {
    double lam = th[1], al = th[2];
    double a = lam * al, c = lam * (1.0 - al);
    double x = a / std::sqrt(c);
    double z = y - th[0];
    double log_z = std::log(2.0) - 0.5 * std::log(c) + enet_logM(x);
    return -a * std::fabs(z) - c * (z * z) / 2.0 - log_z;
}
inline double negbin1_logpdf(double y, const double* th) {
    if (!lp_count(y)) return R_NegInf;
    return R::dnbinom(y, th[0] / th[1], 1.0 / (1.0 + th[1]), 1);
}
inline double skewnormal1_logpdf_th(double y, const double* th) {
    return skewnormal1_logpdf(y, th[0], th[1], th[2]);
}
inline double skewnormal2_logpdf_th(double y, const double* th) {
    return skewnormal2_logpdf(y, th[0], th[1], th[2]);
}
inline double skewt_logpdf_th(double y, const double* th) {
    return skewt_logpdf(y, th[0], th[1], th[2], th[3]);
}
inline double pseudohuber_logpdf_th(double y, const double* th) {
    return pseudohuber_logpdf(y, th[0], th[1], th[2], ph_nu(th[2]));
}
inline double pseudohuber2_logpdf_th(double y, const double* th) {
    return pseudohuber2_logpdf(y, th[0], th[1], th[2],
                               pseudohuber2_terms(th[2]));
}
inline double pig1_logpdf(double y, const double* th) {
    if (!pig_in_support(y)) return R_NegInf;
    double L[5] = {0, 0, 0, 0, 0};
    pig_logS_w(y, pig1_w(th[0], th[1]), 0, L);
    return pig1_value(y, th[0], th[1], L);
}
inline double pig2_logpdf(double y, const double* th) {
    if (!pig_in_support(y)) return R_NegInf;
    double p[5] = {0, 0, 0, 0, 0}, logS = 0.0, A1 = 0.0;
    psi_derivs(y, th[1], 0, p, &logS, &A1);
    return pig2_value(y, th[0], th[1], logS);
}
inline double betabinom1_logpdf(double y, const double* th) {
    double size = th[2];
    if (y < 0 || y > size || y != std::floor(y)) return R_NegInf;
    return bb_log_mass(y, th[0] / th[1], (1.0 - th[0]) / th[1], size);
}
inline double betabinom2_logpdf(double y, const double* th) {
    double size = th[2];
    if (y < 0 || y > size || y != std::floor(y)) return R_NegInf;
    return bb_log_mass(y, th[0], th[1], size);
}

} // namespace d7

#endif
