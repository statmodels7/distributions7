#include <Rcpp.h>
#include "d7_par.h"
using namespace Rcpp;

// Gaussian in the mean and the STANDARD DEVIATION. Every component is written
// in z = (y - mu)/sigma and inv = 1/sigma, and never in a positive power of
// the scale. The algebraically equivalent forms -- (res^2 - sigma^2)/sigma^3
// for the score, (sigma^2 - 3 res^2)/sigma^4 for the second derivative --
// overflow in the DENOMINATOR before the ratio does, so the score returns
// exactly 0 above sigma = 5.6e102 and NaN above 1.3e154 where its true value,
// -1/sigma on this scale and -1 on the link scale, stays representable. A
// score of zero is what a stopping rule reads as stationarity, so an
// optimizer that wanders out there is told it has arrived.
//
// Per-observation maps, so the loops run through d7::par_for: observation
// i's derivatives are written to slot i by one thread, and the result is
// bit-identical at any thread count.

// [[Rcpp::export]]
List gaussian_gradient_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma, int threads = 1) {
    int n = y.size();

    NumericVector grad_mu(n);
    NumericVector grad_sigma(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *sp = sigma.begin();
    double *gm = grad_mu.begin(), *gs = grad_sigma.begin();

    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = mu_is_scalar ? mp[0] : mp[i];
        double s = sigma_is_scalar ? sp[0] : sp[i];

        double inv = 1.0 / s;
        double z = (yp[i] - m) / s;

        gm[i] = z * inv;
        gs[i] = (z * z - 1.0) * inv;
    });

    return List::create(
        Named("mu") = grad_mu,
        Named("sigma") = grad_sigma
    );
}

// [[Rcpp::export]]
List gaussian_hessian_cpp(NumericVector y, NumericVector mu,
                          NumericVector sigma, int threads = 1) {
    int n = y.size();

    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma_sigma(n);
    NumericVector hess_mu_sigma(n);

    bool mu_is_scalar = (mu.size() == 1);
    bool sigma_is_scalar = (sigma.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *sp = sigma.begin();
    double *hmm = hess_mu_mu.begin(), *hss = hess_sigma_sigma.begin(),
           *hms = hess_mu_sigma.begin();

    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double m = mu_is_scalar ? mp[0] : mp[i];
        double s = sigma_is_scalar ? sp[0] : sp[i];

        double inv = 1.0 / s;
        double inv2 = inv * inv;
        double z = (yp[i] - m) / s;

        hmm[i] = -inv2;
        hss[i] = (1.0 - 3.0 * z * z) * inv2;
        hms[i] = -2.0 * z * inv2;
    });
    return List::create(
        Named("mu_mu") = hess_mu_mu,
        Named("sigma_sigma") = hess_sigma_sigma,
        Named("mu_sigma") = hess_mu_sigma
    );
}

// [[Rcpp::export]]
List gaussian_expected_hessian_cpp(NumericVector y, NumericVector mu,
                                   NumericVector sigma, int threads = 1) {
    int n = y.size();

    NumericVector hess_mu_mu(n);
    NumericVector hess_sigma_sigma(n);
    NumericVector hess_mu_sigma(n);

    bool sigma_is_scalar = (sigma.size() == 1);
    const double *sp = sigma.begin();
    double *hmm = hess_mu_mu.begin(), *hss = hess_sigma_sigma.begin(),
           *hms = hess_mu_sigma.begin();

    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double s = sigma_is_scalar ? sp[0] : sp[i];
        double inv = 1.0 / s;
        double inv2 = inv * inv;

        hmm[i] = -inv2;
        hss[i] = -2.0 * inv2;
        hms[i] = 0.0;
    });
    return List::create(
        Named("mu_mu") = hess_mu_mu,
        Named("sigma_sigma") = hess_sigma_sigma,
        Named("mu_sigma") = hess_mu_sigma
    );
}

// The derivatives of the expected information in the parameters. With
// E_mm = -1/sigma^2, E_ss = -2/sigma^2 and E_ms = 0, only the sigma
// derivatives survive:
//   d_s E_mm = 2/sigma^3,   d_s E_ss = 4/sigma^3,
//   d_ss E_mm = -6/sigma^4, d_ss E_ss = -12/sigma^4.
// Keys are "<ab>_<c>" at order 1 and "<ab>_<cd>" at order 2, with <ab> and
// <cd> spelled as hess_names() spells them.
// [[Rcpp::export]]
List gaussian_dexpected_cpp(NumericVector y, NumericVector mu,
                            NumericVector sigma, int order, int threads = 1) {
    int n = y.size();
    bool sigma_is_scalar = (sigma.size() == 1);
    const double *sp = sigma.begin();
    NumericVector zero(n);
    if (order == 1) {
        NumericVector mms(n), sss(n);
        double *a = mms.begin(), *b = sss.begin();
        d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
            double s = sigma_is_scalar ? sp[0] : sp[i];
            double inv = 1.0 / s, inv3 = inv * inv * inv;
            a[i] = 2.0 * inv3;
            b[i] = 4.0 * inv3;
        });
        return List::create(
            Named("mu_mu_mu") = zero, Named("mu_mu_sigma") = mms,
            Named("sigma_sigma_mu") = clone(zero), Named("sigma_sigma_sigma") = sss,
            Named("mu_sigma_mu") = clone(zero), Named("mu_sigma_sigma") = clone(zero));
    }
    NumericVector mmss(n), ssss(n);
    double *a = mmss.begin(), *b = ssss.begin();
    d7::par_for(n, threads, d7::kMinCheap, [&](std::size_t i) {
        double s = sigma_is_scalar ? sp[0] : sp[i];
        double inv = 1.0 / s, inv2 = inv * inv, inv4 = inv2 * inv2;
        a[i] = -6.0 * inv4;
        b[i] = -12.0 * inv4;
    });
    return List::create(
        Named("mu_mu_mu_mu") = zero, Named("mu_mu_sigma_sigma") = mmss,
        Named("mu_mu_mu_sigma") = clone(zero),
        Named("sigma_sigma_mu_mu") = clone(zero),
        Named("sigma_sigma_sigma_sigma") = ssss,
        Named("sigma_sigma_mu_sigma") = clone(zero),
        Named("mu_sigma_mu_mu") = clone(zero),
        Named("mu_sigma_sigma_sigma") = clone(zero),
        Named("mu_sigma_mu_sigma") = clone(zero));
}
