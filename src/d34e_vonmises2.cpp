// The third and fourth derivatives of the expected information of the von Mises (mu, rho)
// family, on the parameter scale, in closed form (one kernel per order, vonmises2_dexpected3_cpp
// and vonmises2_dexpected4_cpp). The expected information does not depend on mu.
//
// With rho = A(kappa) and K(rho) = kappa(rho) the inverse, the components are
// E_mm = -rho K, so d^j E_mm = -(rho K^(j) + j K^(j-1)), and E_rr = -K', so
// d^j E_rr = -K^(j+1). The derivatives of the inverse follow the rule
// K^(n) = (1/A') d/dkappa K^(n-1) written in A', A'', ... at kappa (derive.py),
// with K^(1) = 1/A', K^(2) = -A''/A'^3, and so on. Order three reads K^(2) ..
// K^(4) and order four K^(2) .. K^(5). kappa = K(rho) comes from numericals7
// (n7_bessel_ratio_inverse), A^(1) .. A^(4) from n7_bessel_ratio_upto, and A^(5)
// from the Riccati equation, with A = rho, as in vonmises1.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "d7_par.h"
#include "d34e_run.h"
#include "pt_vonmises.h"
using namespace Rcpp;
using d7::N7Bessel;

namespace {

// A^(5)(kappa) for kappa < 1/2 from the Taylor series of A = I1/I0 at zero,
// A(k) = sum_n c_n k^(2n+1). The coefficients follow from I1 = A I0 with
// I0 = sum_m b_m k^(2m), I1 = sum_m e_m k^(2m+1), b_m = 1/(4^m m!^2),
// e_m = 1/(2^(2m+1) m! (m+1)!): c_N = e_N - sum_{n<N} c_n b_(N-n). The
// Riccati formula loses about 1/k^4 of its digits at small kappa (its terms
// are of order 1/k^5 against a result of order one), so below one half this
// series is used. Its radius is the first zero of I0, 2.40 i, so the twelve
// terms here give a relative error of order (0.5/2.40)^26, below 1e-17.
inline double vm_a5_series(double k) {
    const int N = 12;
    double b[N + 1], e[N + 1], c[N + 1], fac[N + 2];
    fac[0] = 1.0;
    for (int m = 1; m <= N + 1; ++m) fac[m] = fac[m - 1] * m;
    for (int m = 0; m <= N; ++m) {
        b[m] = 1.0 / (std::pow(4.0, m) * fac[m] * fac[m]);
        e[m] = 1.0 / (std::pow(2.0, 2 * m + 1) * fac[m] * fac[m + 1]);
    }
    for (int m = 0; m <= N; ++m) {
        double s = e[m];
        for (int n = 0; n < m; ++n) s -= c[n] * b[m - n];
        c[m] = s;  // b[0] = 1
    }
    // sum_{n >= 2} c_n (2n+1)(2n)(2n-1)(2n-2)(2n-3) k^(2n-4)
    double out = 0.0;
    for (int n = 2; n <= N; ++n) {
        const double poch = (2.0 * n + 1) * (2.0 * n) * (2.0 * n - 1) *
            (2.0 * n - 2) * (2.0 * n - 3);
        out += c[n] * poch * std::pow(k, 2 * n - 4);
    }
    return out;
}

}  // namespace

// [[Rcpp::export]]
List vonmises2_dexpected3_cpp(NumericVector y, NumericVector mu,
                         NumericVector rho, int threads = 1) {
    const std::size_t n = std::max({y.size(), mu.size(), rho.size()});
    const N7Bessel& B = d7::vm_bessel();
    const d34e::Par P(rho);
    return d34e::run(static_cast<int>(n), threads, d7::kMinCostly, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double rho = P[i];
        const double k = B.inv(rho);
        double a[5];
        B.upto(k, 4, a);
        const double A1 = a[1];
        const double A2 = a[2];
        const double A3 = a[3];
        const double A4 = a[4];
        // mu_mu_mu_mu_mu
        o[0] = 0.0;
        // mu_mu_mu_mu_rho
        o[1] = 0.0;
        // mu_mu_mu_rho_rho
        o[2] = 0.0;
        // mu_mu_rho_rho_rho
        o[3] = 3*A2/((A1)*(A1)*(A1)) - rho*(-A1*A3 + 3*((A2)*(A2)))/((A1)*(A1)*(A1)*(A1)*(A1));
        // rho_rho_mu_mu_mu
        o[4] = 0.0;
        // rho_rho_mu_mu_rho
        o[5] = 0.0;
        // rho_rho_mu_rho_rho
        o[6] = 0.0;
        // rho_rho_rho_rho_rho
        o[7] = -(-((A1)*(A1))*A4 + 6*A1*A2*A3 + A2*(4*A1*A3 - 15*((A2)*(A2))))/((A1)*(A1)*(A1)*(A1)*(A1)*(A1)*(A1));
        // mu_rho_mu_mu_mu
        o[8] = 0.0;
        // mu_rho_mu_mu_rho
        o[9] = 0.0;
        // mu_rho_mu_rho_rho
        o[10] = 0.0;
        // mu_rho_rho_rho_rho
        o[11] = 0.0;
    });
}

// [[Rcpp::export]]
List vonmises2_dexpected4_cpp(NumericVector y, NumericVector mu,
                         NumericVector rho, int threads = 1) {
    const std::size_t n = std::max({y.size(), mu.size(), rho.size()});
    const N7Bessel& B = d7::vm_bessel();
    const d34e::Par P(rho);
    return d34e::run(static_cast<int>(n), threads, d7::kMinCostly, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double rho = P[i];
        const double k = B.inv(rho);
        double a[5];
        B.upto(k, 4, a);
        const double A1 = a[1];
        const double A2 = a[2];
        const double A3 = a[3];
        const double A4 = a[4];
        const double A5 = (k < 0.5) ? vm_a5_series(k) : (-8*A1*A3 + 24*A1/((k)*(k)*(k)*(k)) - 6*((A2)*(A2)) - 12*A2/((k)*(k)*(k)) + 4*A3/((k)*(k)) - 2*A4*rho - A4/k - 24*rho/((k)*(k)*(k)*(k)*(k)));
        // mu_mu_mu_mu_mu_mu
        o[0] = 0.0;
        // mu_mu_mu_mu_mu_rho
        o[1] = 0.0;
        // mu_mu_mu_mu_rho_rho
        o[2] = 0.0;
        // mu_mu_mu_rho_rho_rho
        o[3] = 0.0;
        // mu_mu_rho_rho_rho_rho
        o[4] = -4*(-A1*A3 + 3*((A2)*(A2)))/((A1)*(A1)*(A1)*(A1)*(A1)) - rho*(-((A1)*(A1))*A4 + 6*A1*A2*A3 + A2*(4*A1*A3 - 15*((A2)*(A2))))/((A1)*(A1)*(A1)*(A1)*(A1)*(A1)*(A1));
        // rho_rho_mu_mu_mu_mu
        o[5] = 0.0;
        // rho_rho_mu_mu_mu_rho
        o[6] = 0.0;
        // rho_rho_mu_mu_rho_rho
        o[7] = 0.0;
        // rho_rho_mu_rho_rho_rho
        o[8] = 0.0;
        // rho_rho_rho_rho_rho_rho
        o[9] = -(-((A1)*(A1)*(A1))*A5 + 10*((A1)*(A1))*A2*A4 + 5*A1*A3*(2*A1*A3 - 9*((A2)*(A2))) + A2*(7*((A1)*(A1))*A4 - 42*A1*A2*A3 + 2*A1*(-A1*A4 + 5*A2*A3) - 7*A2*(4*A1*A3 - 15*((A2)*(A2)))))/((A1)*(A1)*(A1)*(A1)*(A1)*(A1)*(A1)*(A1)*(A1));
        // mu_rho_mu_mu_mu_mu
        o[10] = 0.0;
        // mu_rho_mu_mu_mu_rho
        o[11] = 0.0;
        // mu_rho_mu_mu_rho_rho
        o[12] = 0.0;
        // mu_rho_mu_rho_rho_rho
        o[13] = 0.0;
        // mu_rho_rho_rho_rho_rho
        o[14] = 0.0;
    });
}

