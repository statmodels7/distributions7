// The third and fourth derivatives of the expected information of the gamma
// in its mean mu and its variance sigma2 (gamma2), on the parameter scale,
// in closed form. The expected information is
//   E_mm = -(1 + 4 G0)/v,   E_mv = 2 m G0/v^2,   E_vv = -a G0/v^2,
// with a = m^2/v and G0 = G(a) = a psi'(a) - 1, and G_j = d^j G/da^j. Writing
// psi'(a) = 1/a - r1(a), the shape functions are, with r_k the remainders of
// psi_diff.h (r_k = (-1)^(k+1) (k-1)!/a^k - psi^(k)(a), r_k' = r_(k+1)),
//   G_j = -(a r_(j+1) + j r_j),
// so G_0 .. G_3 need r_1 .. r_4 and G_4 needs r_5. r_4 and r_5 are written here
// by the construction of psi_diff.h: an asymptotic series from a = 50 and the
// direct difference below it. Each order is its own kernel. The derivatives in
// (m, v) are the total derivatives of the chain rule, derived symbolically
// (misure_d34e/agents/gamma2/derive.py) and printed by emit.py, so every
// component is a polynomial in m, 1/v and the G_j, with no cancellation
// between the shape functions and the prefactors at a large shape.
//
// The components are in the order of d3expected_names() and d4expected_names()
// for gamma2: for each pair (mu_mu, sigma2_sigma2, mu_sigma2), the multisets of
// the derivative indices in the order mu.mu.mu, mu.mu.sigma2, mu.sigma2.sigma2,
// sigma2.sigma2.sigma2 (order 3) and then the order-4 multisets, the same way.

#include <Rcpp.h>
#include "d7_par.h"
#include "psi_diff.h"
#include "d34e_run.h"
using namespace Rcpp;

namespace {

// r_4(a) = -6/a^4 - psi''''(a)
inline double d34_r4(double a) {
    if (a >= 50.0) {
        const double u = 1.0 / a, u2 = u * u;
        return u * u * u * u * u *
               (12.0 + u * (10.0 + u2 * (-7.0 + u2 * (12.0 + u2 * (-33.0)))));
    }
    return -6.0 / (a * a * a * a) - R::psigamma(a, 4);
}

// r_5(a) = 24/a^5 - psi^(5)(a), asymptotically
// -u^6 (60 + 60 u - 56 u^3 + 120 u^5 - 396 u^7 + 1820 u^9), u = 1/a
inline double d34_r5(double a) {
    if (a >= 50.0) {
        const double u = 1.0 / a, u2 = u * u;
        return -u * u * u * u * u * u *
               (60.0 + u * (60.0 + u2 * (-56.0 + u2 * (120.0 + u2 * (-396.0 +
                   u2 * 1820.0)))));
    }
    return 24.0 / (a * a * a * a * a) - R::psigamma(a, 5);
}

}  // namespace

// [[Rcpp::export]]
List gamma2_dexpected3_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma2, int threads = 1) {
    d34e::Par M(mu), V(sigma2);
    return d34e::run(y.size(), threads, d7::kMinCostly, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double m = M[i], v = V[i], a = m * m / v;
        const double r1 = d7::psi1_rest(a), r2 = d7::psi2_rest(a),
                     r3 = d7::psi3_rest(a), r4 = d34_r4(a);
        const double g0 = -a * r1, g1 = -(r1 + a * r2),
                     g2 = -(2.0 * r2 + a * r3), g3 = -(3.0 * r3 + a * r4);
        // below a = 1 the sigma2 components are read in q(a) = a psi'(1 + a)
        // and its derivatives, where the poles of the shape functions have
        // cancelled exactly; hq_j = q^(j)(a), hq_0 = q(a) - 1
        const bool sm = a < 1.0;
        double hq0 = 0.0, hq1 = 0.0, hq2 = 0.0, hq3 = 0.0;
        if (sm) {
            const double t1 = R::psigamma(1.0 + a, 1), t2 = R::psigamma(1.0 + a, 2),
                         t3 = R::psigamma(1.0 + a, 3), t4 = R::psigamma(1.0 + a, 4);
            hq0 = a * t1 - 1.0;
            hq1 = t1 + a * t2;
            hq2 = 2.0 * t2 + a * t3;
            hq3 = 3.0 * t3 + a * t4;
        }
        o[0] = (-48.0*g2*m/(v*v*v) - 32.0*g3*m*m*m/(v*v*v*v));
        o[1] = (16.0*g1/(v*v*v) + 16.0*g3*m*m*m*m/(v*v*v*v*v) + 56.0*g2*m*m/(v*v*v*v));
        o[2] = (-48.0*g1*m/(v*v*v*v) - 48.0*g2*m*m*m/(v*v*v*v*v) - 8.0*g3*m*m*m*m*m/(v*v*v*v*v*v));
        o[3] = (6.0/(v*v*v*v) + 24.0*g0/(v*v*v*v) + 4.0*g3*m*m*m*m*m*m/(v*v*v*v*v*v*v) + 36.0*g2*m*m*m*m/(v*v*v*v*v*v) + 72.0*g1*m*m/(v*v*v*v*v));
        o[4] = sm ? (-36.0*hq2*m*m*m/(v*v*v*v*v) - 24.0*hq1*m/(v*v*v*v) - 8.0*hq3*m*m*m*m*m/(v*v*v*v*v*v)) : (-36.0*g2*m*m*m/(v*v*v*v*v) - 24.0*g1*m/(v*v*v*v) - 8.0*g3*m*m*m*m*m/(v*v*v*v*v*v));
        o[5] = sm ? (6.0*hq0/(v*v*v*v) + 4.0*hq3*m*m*m*m*m*m/(v*v*v*v*v*v*v) + 30.0*hq2*m*m*m*m/(v*v*v*v*v*v) + 42.0*hq1*m*m/(v*v*v*v*v)) : (6.0*g0/(v*v*v*v) + 4.0*g3*m*m*m*m*m*m/(v*v*v*v*v*v*v) + 30.0*g2*m*m*m*m/(v*v*v*v*v*v) + 42.0*g1*m*m/(v*v*v*v*v));
        o[6] = sm ? (-56.0*hq1*m*m*m/(v*v*v*v*v*v) - 24.0*hq0*m/(v*v*v*v*v) - 22.0*hq2*m*m*m*m*m/(v*v*v*v*v*v*v) - 2.0*hq3*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v)) : (-56.0*g1*m*m*m/(v*v*v*v*v*v) - 24.0*g0*m/(v*v*v*v*v) - 22.0*g2*m*m*m*m*m/(v*v*v*v*v*v*v) - 2.0*g3*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v));
        o[7] = sm ? (24.0/(v*v*v*v*v) + hq3*m*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v) + 15.0*hq2*m*m*m*m*m*m/(v*v*v*v*v*v*v*v) + 60.0*hq0*m*m/(v*v*v*v*v*v) + 60.0*hq1*m*m*m*m/(v*v*v*v*v*v*v)) : (g3*m*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v) + 15.0*g2*m*m*m*m*m*m/(v*v*v*v*v*v*v*v) + 60.0*g0*m*m/(v*v*v*v*v*v) + 60.0*g1*m*m*m*m/(v*v*v*v*v*v*v));
        o[8] = (12.0*g1/(v*v*v) + 16.0*g3*m*m*m*m/(v*v*v*v*v) + 48.0*g2*m*m/(v*v*v*v));
        o[9] = (-44.0*g2*m*m*m/(v*v*v*v*v) - 36.0*g1*m/(v*v*v*v) - 8.0*g3*m*m*m*m*m/(v*v*v*v*v*v));
        o[10] = (12.0*g0/(v*v*v*v) + 4.0*g3*m*m*m*m*m*m/(v*v*v*v*v*v*v) + 34.0*g2*m*m*m*m/(v*v*v*v*v*v) + 60.0*g1*m*m/(v*v*v*v*v));
        o[11] = (-72.0*g1*m*m*m/(v*v*v*v*v*v) - 48.0*g0*m/(v*v*v*v*v) - 24.0*g2*m*m*m*m*m/(v*v*v*v*v*v*v) - 2.0*g3*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v));
    });
}

// [[Rcpp::export]]
List gamma2_dexpected4_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma2, int threads = 1) {
    d34e::Par M(mu), V(sigma2);
    return d34e::run(y.size(), threads, d7::kMinCostly, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double m = M[i], v = V[i], a = m * m / v;
        const double r1 = d7::psi1_rest(a), r2 = d7::psi2_rest(a),
                     r3 = d7::psi3_rest(a), r4 = d34_r4(a), r5 = d34_r5(a);
        const double g0 = -a * r1, g1 = -(r1 + a * r2),
                     g2 = -(2.0 * r2 + a * r3), g3 = -(3.0 * r3 + a * r4),
                     g4 = -(4.0 * r4 + a * r5);
        const bool sm = a < 1.0;
        double hq0 = 0.0, hq1 = 0.0, hq2 = 0.0, hq3 = 0.0, hq4 = 0.0;
        if (sm) {
            const double t1 = R::psigamma(1.0 + a, 1), t2 = R::psigamma(1.0 + a, 2),
                         t3 = R::psigamma(1.0 + a, 3), t4 = R::psigamma(1.0 + a, 4),
                         t5 = R::psigamma(1.0 + a, 5);
            hq0 = a * t1 - 1.0;
            hq1 = t1 + a * t2;
            hq2 = 2.0 * t2 + a * t3;
            hq3 = 3.0 * t3 + a * t4;
            hq4 = 4.0 * t4 + a * t5;
        }
        o[0] = (-48.0*g2/(v*v*v) - 192.0*g3*m*m/(v*v*v*v) - 64.0*g4*m*m*m*m/(v*v*v*v*v));
        o[1] = (32.0*g4*m*m*m*m*m/(v*v*v*v*v*v) + 144.0*g2*m/(v*v*v*v) + 176.0*g3*m*m*m/(v*v*v*v*v));
        o[2] = (-48.0*g1/(v*v*v*v) - 240.0*g2*m*m/(v*v*v*v*v) - 136.0*g3*m*m*m*m/(v*v*v*v*v*v) - 16.0*g4*m*m*m*m*m*m/(v*v*v*v*v*v*v));
        o[3] = (8.0*g4*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v) + 96.0*g3*m*m*m*m*m/(v*v*v*v*v*v*v) + 192.0*g1*m/(v*v*v*v*v) + 288.0*g2*m*m*m/(v*v*v*v*v*v));
        o[4] = (-24.0/(v*v*v*v*v) - 96.0*g0/(v*v*v*v*v) - 384.0*g1*m*m/(v*v*v*v*v*v) - 288.0*g2*m*m*m*m/(v*v*v*v*v*v*v) - 64.0*g3*m*m*m*m*m*m/(v*v*v*v*v*v*v*v) - 4.0*g4*m*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v));
        o[5] = sm ? (-24.0*hq1/(v*v*v*v) - 156.0*hq2*m*m/(v*v*v*v*v) - 112.0*hq3*m*m*m*m/(v*v*v*v*v*v) - 16.0*hq4*m*m*m*m*m*m/(v*v*v*v*v*v*v)) : (-24.0*g1/(v*v*v*v) - 156.0*g2*m*m/(v*v*v*v*v) - 112.0*g3*m*m*m*m/(v*v*v*v*v*v) - 16.0*g4*m*m*m*m*m*m/(v*v*v*v*v*v*v));
        o[6] = sm ? (8.0*hq4*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v) + 84.0*hq3*m*m*m*m*m/(v*v*v*v*v*v*v) + 96.0*hq1*m/(v*v*v*v*v) + 204.0*hq2*m*m*m/(v*v*v*v*v*v)) : (8.0*g4*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v) + 84.0*g3*m*m*m*m*m/(v*v*v*v*v*v*v) + 96.0*g1*m/(v*v*v*v*v) + 204.0*g2*m*m*m/(v*v*v*v*v*v));
        o[7] = sm ? (-24.0*hq0/(v*v*v*v*v) - 222.0*hq2*m*m*m*m/(v*v*v*v*v*v*v) - 216.0*hq1*m*m/(v*v*v*v*v*v) - 58.0*hq3*m*m*m*m*m*m/(v*v*v*v*v*v*v*v) - 4.0*hq4*m*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v)) : (-24.0*g0/(v*v*v*v*v) - 222.0*g2*m*m*m*m/(v*v*v*v*v*v*v) - 216.0*g1*m*m/(v*v*v*v*v*v) - 58.0*g3*m*m*m*m*m*m/(v*v*v*v*v*v*v*v) - 4.0*g4*m*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v));
        o[8] = sm ? (2.0*hq4*m*m*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v*v) + 38.0*hq3*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v) + 120.0*hq0*m/(v*v*v*v*v*v) + 210.0*hq2*m*m*m*m*m/(v*v*v*v*v*v*v*v) + 360.0*hq1*m*m*m/(v*v*v*v*v*v*v)) : (2.0*g4*m*m*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v*v) + 38.0*g3*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v) + 120.0*g0*m/(v*v*v*v*v*v) + 210.0*g2*m*m*m*m*m/(v*v*v*v*v*v*v*v) + 360.0*g1*m*m*m/(v*v*v*v*v*v*v));
        o[9] = sm ? (-120.0/(v*v*v*v*v*v) - hq4*m*m*m*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v*v*v) - 480.0*hq1*m*m*m*m/(v*v*v*v*v*v*v*v) - 360.0*hq0*m*m/(v*v*v*v*v*v*v) - 180.0*hq2*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v) - 24.0*hq3*m*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v*v)) : (-g4*m*m*m*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v*v*v) - 480.0*g1*m*m*m*m/(v*v*v*v*v*v*v*v) - 360.0*g0*m*m/(v*v*v*v*v*v*v) - 180.0*g2*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v) - 24.0*g3*m*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v*v));
        o[10] = (32.0*g4*m*m*m*m*m/(v*v*v*v*v*v) + 120.0*g2*m/(v*v*v*v) + 160.0*g3*m*m*m/(v*v*v*v*v));
        o[11] = (-36.0*g1/(v*v*v*v) - 204.0*g2*m*m/(v*v*v*v*v) - 128.0*g3*m*m*m*m/(v*v*v*v*v*v) - 16.0*g4*m*m*m*m*m*m/(v*v*v*v*v*v*v));
        o[12] = (8.0*g4*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v) + 92.0*g3*m*m*m*m*m/(v*v*v*v*v*v*v) + 144.0*g1*m/(v*v*v*v*v) + 256.0*g2*m*m*m/(v*v*v*v*v*v));
        o[13] = (-48.0*g0/(v*v*v*v*v) - 312.0*g1*m*m/(v*v*v*v*v*v) - 264.0*g2*m*m*m*m/(v*v*v*v*v*v*v) - 62.0*g3*m*m*m*m*m*m/(v*v*v*v*v*v*v*v) - 4.0*g4*m*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v));
        o[14] = (2.0*g4*m*m*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v*v) + 40.0*g3*m*m*m*m*m*m*m/(v*v*v*v*v*v*v*v*v) + 240.0*g0*m/(v*v*v*v*v*v) + 240.0*g2*m*m*m*m*m/(v*v*v*v*v*v*v*v) + 480.0*g1*m*m*m/(v*v*v*v*v*v*v));
    });
}
