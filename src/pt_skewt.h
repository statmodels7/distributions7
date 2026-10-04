#ifndef D7_PT_SKEWT_H
#define D7_PT_SKEWT_H

#include <Rcpp.h>
#include <cmath>
#include <algorithm>

// Azzalini's skew t in (mu, sigma, alpha, nu): one function per component
// and order for the score and the diagonal of the Hessian, called by the
// vector kernels (skewt.cpp). With z = (y - mu)/sigma, m = nu + 1,
// s = nu + z^2, c = sqrt(m/s) and w = alpha z c,
//   l = log 2 - log sigma + log t_nu(z) + log T_{nu+1}(w).
// The components in (mu, sigma, alpha) are closed forms in the quantities
// skewt_pieces() returns (skewt_distrib.R's skewt_pieces(), written in C);
// those in nu are one five-point stencil on the log-density, at the step
// max(1e-3 |nu|, 1e-6) and with numericals7's weights at accuracy four, as
// the R methods take them. The caller evaluates the log-density at the
// stencil's nodes (skewt_logpdf) and passes the values.
//
// The t distribution function is numericals7's n7_pt(), R's pt() compiled
// without its warnings, resolved once through R_GetCCallable on the calling
// thread (skewt_pt()), so that the kernels and the registry may run these
// functions on worker threads. R::dt() is silent at positive degrees of
// freedom.

namespace d7 {

typedef double (*N7Pt)(double, double, int, int);

// n7_pt, resolved on first use; the first call must be on the calling thread
N7Pt skewt_pt();

struct SkewtPieces { double z, w, c, a, da, e, b, db, q, dq; };

inline SkewtPieces skewt_pieces(double y, double mu, double sigma,
                                double alpha, double nu) {
    SkewtPieces P;
    double z = (y - mu) / sigma;
    double m = nu + 1.0;
    double s = nu + z * z;
    double cc = std::sqrt(m / s);
    double w = alpha * z * cc;
    double q = std::exp(R::dt(w, m, 1) - skewt_pt()(w, m, 1, 1));
    double e = nu * std::sqrt(m) / std::pow(s, 1.5);
    P.z = z;
    P.w = w;
    P.c = cc;
    P.a = -m * z / s;
    P.da = -m * (nu - z * z) / (s * s);
    P.e = e;
    P.b = alpha * e;
    P.db = -3.0 * alpha * nu * std::sqrt(m) * z / std::pow(s, 2.5);
    P.q = q;
    P.dq = q * (-(m + 1.0) * w / (m + w * w) - q);
    return P;
}

inline double skewt_logpdf(double y, double mu, double sigma, double alpha,
                           double nu) {
    double z = (y - mu) / sigma;
    double w = alpha * z * std::sqrt((nu + 1.0) / (nu + z * z));
    return std::log(2.0) - std::log(sigma) + R::dt(z, nu, 1) +
        skewt_pt()(w, nu + 1.0, 1, 1);
}

// the step in nu, and the stencil weights on the offsets -2, -1, 0, 1, 2,
// exactly as numericals7::fd_weights() returns them
inline double skewt_nu_step(double nu) {
    return std::max(1e-3 * std::fabs(nu), 1e-6);
}

const double kSkewtW1[5] = {0x1.5555555555555p-4, -0x1.5555555555555p-1, 0.0,
                            0x1.5555555555555p-1, -0x1.5555555555555p-4};
const double kSkewtW2[5] = {-0x1.5555555555553p-4, 0x1.5555555555553p+0,
                            -0x1.3ffffffffffffp+1, 0x1.5555555555555p+0,
                            -0x1.5555555555555p-4};

inline double skewt_score_mu(const SkewtPieces& P, double sigma) {
    double d = P.a + P.q * P.b;
    return -d / sigma;
}

inline double skewt_score_sigma(const SkewtPieces& P, double sigma) {
    double d = P.a + P.q * P.b;
    return -(1.0 + P.z * d) / sigma;
}

inline double skewt_score_alpha(const SkewtPieces& P) {
    return P.q * P.z * P.c;
}

// lp[k] the log-density at nu + (k - 2) h, k = 0..4; lp[2] is not read
inline double skewt_score_nu(const double* lp, double h) {
    double acc = 0.0;
    for (int k = 0; k < 5; ++k) {
        if (kSkewtW1[k] == 0.0) continue;
        acc = acc + kSkewtW1[k] * lp[k];
    }
    return acc / h;
}

inline double skewt_hess_mu_mu(const SkewtPieces& P, double sigma) {
    double dd = P.da + P.dq * P.b * P.b + P.q * P.db;
    return dd / (sigma * sigma);
}

inline double skewt_hess_sigma_sigma(const SkewtPieces& P, double sigma) {
    double d = P.a + P.q * P.b;
    double dd = P.da + P.dq * P.b * P.b + P.q * P.db;
    return (1.0 + 2.0 * P.z * d + P.z * P.z * dd) / (sigma * sigma);
}

inline double skewt_hess_alpha_alpha(const SkewtPieces& P) {
    double zc = P.z * P.c;
    return P.dq * zc * zc;
}

inline double skewt_hess_nu_nu(const double* lp, double h) {
    double acc = 0.0;
    for (int k = 0; k < 5; ++k) acc = acc + kSkewtW2[k] * lp[k];
    return acc / (h * h);
}

} // namespace d7

#endif
