// The third and fourth derivatives of the expected information of the von Mises (mu, kappa)
// family, on the parameter scale, in closed form (one kernel per order, vonmises1_dexpected3_cpp
// and vonmises1_dexpected4_cpp). The expected information does not depend on mu.
//
// With A(k) = I1(k)/I0(k) and A^(j) its derivatives in kappa (A^(0) = A), the
// components are E_mm = -kappa A, so d^j E_mm = -(kappa A^(j) + j A^(j-1)),
// and E_kk = -A', so d^j E_kk = -A^(j+1). Order three reads A^(2) .. A^(4) and
// order four reads A^(3) .. A^(5). The derivatives A^(0) .. A^(4) come from
// numericals7 (n7_bessel_ratio_upto, orders up to four). A^(5) is the Riccati
// equation A' = 1 - A/kappa - A^2 differentiated four times,
//   A^(5) = -sum_j C(4,j) A^(j) (1/kappa)^(4-j) - sum_j C(4,j) A^(j) A^(4-j),
// with (1/kappa)^(m) = (-1)^m m! / kappa^(m+1), and it is formed here, so the
// fifth derivative is not taken from numericals7 (its upto() stops at four).

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
List vonmises1_dexpected3_cpp(NumericVector y, NumericVector mu,
                         NumericVector kappa, int threads = 1) {
    const std::size_t n = std::max({y.size(), mu.size(), kappa.size()});
    const N7Bessel& B = d7::vm_bessel();
    const d34e::Par P(kappa);
    return d34e::run(static_cast<int>(n), threads, d7::kMinCostly, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double k = P[i];
        double a[5];
        B.upto(k, 4, a);
        const double A2 = a[2];
        const double A3 = a[3];
        const double A4 = a[4];
        // mu_mu_mu_mu_mu
        o[0] = 0.0;
        // mu_mu_mu_mu_kappa
        o[1] = 0.0;
        // mu_mu_mu_kappa_kappa
        o[2] = 0.0;
        // mu_mu_kappa_kappa_kappa
        o[3] = -3*A2 - A3*k;
        // kappa_kappa_mu_mu_mu
        o[4] = 0.0;
        // kappa_kappa_mu_mu_kappa
        o[5] = 0.0;
        // kappa_kappa_mu_kappa_kappa
        o[6] = 0.0;
        // kappa_kappa_kappa_kappa_kappa
        o[7] = -A4;
        // mu_kappa_mu_mu_mu
        o[8] = 0.0;
        // mu_kappa_mu_mu_kappa
        o[9] = 0.0;
        // mu_kappa_mu_kappa_kappa
        o[10] = 0.0;
        // mu_kappa_kappa_kappa_kappa
        o[11] = 0.0;
    });
}

// [[Rcpp::export]]
List vonmises1_dexpected4_cpp(NumericVector y, NumericVector mu,
                         NumericVector kappa, int threads = 1) {
    const std::size_t n = std::max({y.size(), mu.size(), kappa.size()});
    const N7Bessel& B = d7::vm_bessel();
    const d34e::Par P(kappa);
    return d34e::run(static_cast<int>(n), threads, d7::kMinCostly, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double k = P[i];
        double a[5];
        B.upto(k, 4, a);
        const double A1 = a[1];
        const double A2 = a[2];
        const double A3 = a[3];
        const double A4 = a[4];
        const double A0 = a[0];
        const double A5 = (k < 0.5) ? vm_a5_series(k) : (-2*A0*A4 - 24*A0/((k)*(k)*(k)*(k)*(k)) - 8*A1*A3 + 24*A1/((k)*(k)*(k)*(k)) - 6*((A2)*(A2)) - 12*A2/((k)*(k)*(k)) + 4*A3/((k)*(k)) - A4/k);
        // mu_mu_mu_mu_mu_mu
        o[0] = 0.0;
        // mu_mu_mu_mu_mu_kappa
        o[1] = 0.0;
        // mu_mu_mu_mu_kappa_kappa
        o[2] = 0.0;
        // mu_mu_mu_kappa_kappa_kappa
        o[3] = 0.0;
        // mu_mu_kappa_kappa_kappa_kappa
        o[4] = -4*A3 - A4*k;
        // kappa_kappa_mu_mu_mu_mu
        o[5] = 0.0;
        // kappa_kappa_mu_mu_mu_kappa
        o[6] = 0.0;
        // kappa_kappa_mu_mu_kappa_kappa
        o[7] = 0.0;
        // kappa_kappa_mu_kappa_kappa_kappa
        o[8] = 0.0;
        // kappa_kappa_kappa_kappa_kappa_kappa
        o[9] = -A5;
        // mu_kappa_mu_mu_mu_mu
        o[10] = 0.0;
        // mu_kappa_mu_mu_mu_kappa
        o[11] = 0.0;
        // mu_kappa_mu_mu_kappa_kappa
        o[12] = 0.0;
        // mu_kappa_mu_kappa_kappa_kappa
        o[13] = 0.0;
        // mu_kappa_kappa_kappa_kappa_kappa
        o[14] = 0.0;
    });
}

