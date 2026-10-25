#ifndef D7_PT_NEGBIN1_ASYM_H
#define D7_PT_NEGBIN1_ASYM_H

#include "pt_negbin1_series.h"

// THE NEGATIVE BINOMIAL WITH A LINEAR VARIANCE AT LARGE THETA: the quantities that nb1_sums0,
// nb1_sums1 and nb1_sums2 (pt_negbin1.h) return, A = -E[l_mu,mu] and G = E[l_theta,theta] with
// their derivatives in (mu, theta), computed without a sum over the support. The support sums need
// about theta terms at a ratio theta/(1 + theta), so a single row at theta = 1e8 costs seconds, and
// their cap of 1e6 terms makes every quantity wrong for theta >= 1e6 and mu/theta > 0.5.
//
// Two routes replace them, both in pt_negbin1_series.h, which records the derivations and the
// measured errors (against 60-digit values, relative to the largest component of the order):
//   nb1_asym_ok(): theta >= 1 and mu/theta <= 0.5, the series of the closed form of A in powers of
//     1/(1 + theta), carried termwise in (r, L) = (mu/theta, log1p(theta)) to the fourth derivatives
//     (1e-15 for theta >= 2 and r <= 0.3, 5e-13 at r = 0.5, as the series of this file's earlier
//     version, which stopped at r = 0.3 and theta = 2 and at the second derivatives);
//   the trapezoid rule on the integral of F = theta^2 A, for mu/theta > 0.5 and theta >= 20 (4e-14
//     or better up to mu/theta = 10).
// The outputs have the layouts of nb1_sums0/1/2: A and G as the value, then (value, d_mu, d_theta),
// then (value, d_mu, d_theta, d_mu,mu, d_mu,theta, d_theta,theta).

inline bool nb1_asym_ok(double mu, double th) {
    return nb1d34::series_ok(mu, th);
}

inline void nb1_asym_sums0(double mu, double th, double* A, double* G) {
    nb1d34::series_sums0(mu, th, A, G);
}

inline void nb1_asym_sums1(double mu, double th, double* A, double* G) {
    nb1d34::series_sums1(mu, th, A, G);
}

inline void nb1_asym_sums2(double mu, double th, double* A, double* G) {
    nb1d34::series_sums2(mu, th, A, G);
}

#endif
