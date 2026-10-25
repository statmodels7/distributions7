// The third and fourth derivatives of the expected Hessian of the negative
// binomial with a variance linear in the mean (negbin1), on the parameter scale.
//
// The expected Hessian reduces to two scalar series (negbin1_dexpected2_cpp(),
// pt_negbin1.h): E_mm = -A, E_tt = G and E_mt = -theta G/mu. The kernels of this
// file carry A and G to the derivatives of order three and four, and E_mt without
// the division by mu where it matters, and assemble the components outside the loop;
// pt_negbin1_d34e.h records the construction. Three routes serve, as at orders 0 to 2:
// the sums over the support (D sums for theta below sums_cut(mu), W sums above, the
// mass written p = mu pb for mu < 0.5), the series in powers of 1/(1 + theta) for
// theta >= 1 and mu/theta <= 0.5, and the trapezoid rule on the integral of F for
// mu/theta > 0.5 (theta >= 20, or theta >= 0.3 with mu/theta > 20). Each order is its
// own kernel: the order-3 kernel returns the twelve components of d3expected_names(),
// the order-4 kernel the fifteen of d4expected_names(), in that order (E_mm, E_tt,
// E_mt, each with the derivative indices mu^n, mu^(n-1) theta, ..., theta^n).

#include <Rcpp.h>
#include "d34e_run.h"
#include "pt_negbin1_d34e.h"
using namespace Rcpp;

// [[Rcpp::export]]
List negbin1_dexpected3_cpp(NumericVector y, NumericVector mu,
                            NumericVector theta, int threads = 1) {
    d34e::Par M(mu), T(theta);
    return d34e::run(y.size(), threads, d7::kMinCostly, 12,
                     [&](std::size_t i, d34e::Row o) {
        double out[12];
        nb1d34::d3(M[i], T[i], out);
        for (int j = 0; j < 12; ++j) {
            o[j] = out[j];
        }
    });
}

// [[Rcpp::export]]
List negbin1_dexpected4_cpp(NumericVector y, NumericVector mu,
                            NumericVector theta, int threads = 1) {
    d34e::Par M(mu), T(theta);
    return d34e::run(y.size(), threads, d7::kMinCostly, 15,
                     [&](std::size_t i, d34e::Row o) {
        double out[15];
        nb1d34::d4(M[i], T[i], out);
        for (int j = 0; j < 15; ++j) {
            o[j] = out[j];
        }
    });
}
