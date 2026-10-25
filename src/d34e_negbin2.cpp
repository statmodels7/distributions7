// The third and fourth derivatives of the expected Hessian of the negative
// binomial in its mean mu and size theta (negbin2), on the parameter scale.
//
// E_mm = -theta/(mu (theta + mu)) and E_mt = 0 are closed forms (the derivatives of
// E_mt are exact zeros). E_tt = E[l_theta,theta] is the sum over the support of
// negbin2_expected_theta_theta() (pt_negbin2.h); its derivatives of order three and
// four are the sums of pt_negbin2_d34e.h, which carry the mass derivatives and the
// summand derivatives along one walk of the support, as the order-two kernel
// nb_dE_ltt2() does, with the pure theta derivatives as running sums of merged
// increments. The form is E[l_ab] differentiated under the sum, so no derivative of
// the log-density above the second order is read. For mu < 1 the sums carry the
// mass divided by mu^k, and for theta >= 60 and mu/theta <= 0.65 a series in 1/theta
// with exact rational coefficients replaces the sum (which loses digits as theta
// grows); pt_negbin2_d34e.h records both.
//
// Each order is its own kernel: the order-3 kernel returns the twelve components of
// d3expected_names(), the order-4 kernel the fifteen of d4expected_names(), in that
// order (E_mm, E_tt, E_mt, each with the derivative indices mu^n, mu^(n-1) theta, ...,
// theta^n).

#include <Rcpp.h>
#include "d34e_run.h"
#include "pt_negbin2_d34e.h"
using namespace Rcpp;

// [[Rcpp::export]]
List negbin2_dexpected3_cpp(NumericVector y, NumericVector mu,
                            NumericVector theta, int threads = 1) {
    d34e::Par M(mu), T(theta);
    return d34e::run(y.size(), threads, d7::kMinCostly, 12,
                     [&](std::size_t i, d34e::Row o) {
        double e[4], t[4];
        nb2d34::emm3(M[i], T[i], e);
        nb2d34::ett3(M[i], T[i], t);
        for (int j = 0; j < 4; ++j) {
            o[j] = e[j];
            o[4 + j] = t[j];
            o[8 + j] = 0.0;
        }
    });
}

// [[Rcpp::export]]
List negbin2_dexpected4_cpp(NumericVector y, NumericVector mu,
                            NumericVector theta, int threads = 1) {
    d34e::Par M(mu), T(theta);
    return d34e::run(y.size(), threads, d7::kMinCostly, 15,
                     [&](std::size_t i, d34e::Row o) {
        double e[5], t[5];
        nb2d34::emm4(M[i], T[i], e);
        nb2d34::ett4(M[i], T[i], t);
        for (int j = 0; j < 5; ++j) {
            o[j] = e[j];
            o[5 + j] = t[j];
            o[10 + j] = 0.0;
        }
    });
}
