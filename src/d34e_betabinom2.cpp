// The third and fourth derivatives of the expected information of the
// beta-binomial in its shapes (alpha, beta), on the parameter scale, in closed
// form, for a fixed size N. The expected information is
//   E_pq = sum_{k=0}^{N} f_k l_pq(k),
// with f_k the beta-binomial mass and l_pq the second derivative of the
// log-mass. For a fixed observation the derivatives are taken under the sum:
// the k-th term's derivative of order T is written as exp(F) times a
// polynomial in the derivatives of the log-mass F, which is expanded by sympy
// (misure_d34e/agents/betabinom2/gen.py), so no term is a numerical
// derivative and no automatic differentiation is used.
//
// The derivatives of the log-mass are sums of negative powers at integer
// shifts, accumulated beside the mass: with S_m(x, n) = sum_{i<n} (x+i)^-(m+1)
// and the three power sums A = S(a, k), B = S(b, N-k), C = S(a+b, N), every
// derivative of the log-mass is a signed multiple of one of them, and no
// polygamma is called. The mass uses the product form
//   f_k = choose(N, k) prod_{i<k} (a+i)/(a+b+i) prod_{i<N-k} (b+i)/(a+b+k+i),
// whose factors are all in (0, 1], so it neither overflows nor takes a
// logarithm of a difference.
//
// Each order is its own kernel: the order-3 kernel returns the twelve
// components of d3expected_names(), the order-4 kernel the fifteen of
// d4expected_names(), in that order. The pair is the first two letters of a
// name and the derivative indices its last letters.

#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
#include "d34e_run.h"
using namespace Rcpp;

namespace {

// S[m] = sum_{i<n} (x+i)^-(m+1) for m = 0, ..., M-1.
inline void bb34_power_sums(double x, int n, double* S, int M) {
    for (int m = 0; m < M; ++m) S[m] = 0.0;
    for (int i = 0; i < n; ++i) {
        const double t = 1.0 / (x + i);
        double p = t;
        for (int m = 0; m < M; ++m) {
            S[m] += p;
            p *= t;
        }
    }
}

// The beta-binomial mass at k, from the product form in the header.
inline double bb34_mass(double a, double b, double s, int N, int k) {
    double f = 1.0;
    for (int j = 1; j <= k; ++j) f = f * static_cast<double>(N - k + j) / static_cast<double>(j);
    for (int j = 0; j < k; ++j) f *= (a + j) / (s + j);
    for (int j = 0; j < N - k; ++j) f *= (b + j) / (s + k + j);
    return f;
}

}  // namespace

// [[Rcpp::export]]
List betabinom2_dexpected3_cpp(NumericVector y, NumericVector alpha, NumericVector beta,
                    double size, int threads = 1) {
    d34e::Par Av(alpha), Bv(beta);
    const int N = static_cast<int>(size);
    return d34e::run(y.size(), threads, d7::kMinCostly, 12,
                     [&](std::size_t i, d34e::Row o) {
        const double a = Av[i], b = Bv[i];
        const double s = a + b;
        double C[5];
        bb34_power_sums(s, N, C, 5);
        double acc[12];
        for (int c = 0; c < 12; ++c) acc[c] = 0.0;
        for (int k = 0; k <= N; ++k) {
            const double f = bb34_mass(a, b, s, N, k);
            double A[5], B[5];
            bb34_power_sums(a, k, A, 5);
            bb34_power_sums(b, N - k, B, 5);
            const double L01 = 1.0 * (B[0] - C[0]);
            const double L02 = -1.0 * (B[1] - C[1]);
            const double L03 = 2.0 * (B[2] - C[2]);
            const double L04 = -6.0 * (B[3] - C[3]);
            const double L05 = 24.0 * (B[4] - C[4]);
            const double L10 = 1.0 * (A[0] - C[0]);
            const double L11 = 1.0 * C[1];
            const double L12 = -2.0 * C[2];
            const double L13 = 6.0 * C[3];
            const double L14 = -24.0 * C[4];
            const double L20 = -1.0 * (A[1] - C[1]);
            const double L21 = -2.0 * C[2];
            const double L22 = 6.0 * C[3];
            const double L23 = -24.0 * C[4];
            const double L30 = 2.0 * (A[2] - C[2]);
            const double L31 = 6.0 * C[3];
            const double L32 = -24.0 * C[4];
            const double L40 = -6.0 * (A[3] - C[3]);
            const double L41 = -24.0 * C[4];
            const double L50 = 24.0 * (A[4] - C[4]);
            acc[0] += f * (L50 + (L20*(L10*L10*L10)) + (3.0*L10*L40) + (3.0*L10*(L20*L20)) + (3.0*L30*(L10*L10)) + (4.0*L20*L30));
            acc[1] += f * (L41 + (L01*L40) + (L01*(L20*L20)) + (L21*(L10*L10)) + (2.0*L10*L31) + (2.0*L11*L30) + (2.0*L20*L21) + (L01*L20*(L10*L10)) + (2.0*L01*L10*L30) + (2.0*L10*L11*L20));
            acc[2] += f * (L32 + (L02*L30) + (L10*L22) + (L12*L20) + (L30*(L01*L01)) + (2.0*L01*L31) + (2.0*L11*L21) + (L02*L10*L20) + (L10*L20*(L01*L01)) + (2.0*L01*L10*L21) + (2.0*L01*L11*L20));
            acc[3] += f * (L23 + (L03*L20) + (L20*(L01*L01*L01)) + (3.0*L01*L22) + (3.0*L02*L21) + (3.0*L21*(L01*L01)) + (3.0*L01*L02*L20));
            acc[4] += f * (L32 + (L02*L30) + (L02*(L10*L10*L10)) + (3.0*L10*L22) + (3.0*L12*L20) + (3.0*L12*(L10*L10)) + (3.0*L02*L10*L20));
            acc[5] += f * (L23 + (L01*L22) + (L02*L21) + (L03*L20) + (L03*(L10*L10)) + (2.0*L10*L13) + (2.0*L11*L12) + (L01*L02*L20) + (L01*L02*(L10*L10)) + (2.0*L01*L10*L12) + (2.0*L02*L10*L11));
            acc[6] += f * (L14 + (L04*L10) + (L10*(L02*L02)) + (L12*(L01*L01)) + (2.0*L01*L13) + (2.0*L02*L12) + (2.0*L03*L11) + (L02*L10*(L01*L01)) + (2.0*L01*L02*L11) + (2.0*L01*L03*L10));
            acc[7] += f * (L05 + (L02*(L01*L01*L01)) + (3.0*L01*L04) + (3.0*L01*(L02*L02)) + (3.0*L03*(L01*L01)) + (4.0*L02*L03));
            acc[8] += f * (L41 + (L11*L30) + (L11*(L10*L10*L10)) + (3.0*L10*L31) + (3.0*L20*L21) + (3.0*L21*(L10*L10)) + (3.0*L10*L11*L20));
            acc[9] += f * (L32 + (L01*L31) + (L12*L20) + (L12*(L10*L10)) + (2.0*L10*L22) + (2.0*L10*(L11*L11)) + (3.0*L11*L21) + (L01*L11*L20) + (L01*L11*(L10*L10)) + (2.0*L01*L10*L21));
            acc[10] += f * (L23 + (L02*L21) + (L10*L13) + (L21*(L01*L01)) + (2.0*L01*L22) + (2.0*L01*(L11*L11)) + (3.0*L11*L12) + (L02*L10*L11) + (L10*L11*(L01*L01)) + (2.0*L01*L10*L12));
            acc[11] += f * (L14 + (L03*L11) + (L11*(L01*L01*L01)) + (3.0*L01*L13) + (3.0*L02*L12) + (3.0*L12*(L01*L01)) + (3.0*L01*L02*L11));
        }
        for (int c = 0; c < 12; ++c) o[c] = acc[c];
    });
}

// [[Rcpp::export]]
List betabinom2_dexpected4_cpp(NumericVector y, NumericVector alpha, NumericVector beta,
                    double size, int threads = 1) {
    d34e::Par Av(alpha), Bv(beta);
    const int N = static_cast<int>(size);
    return d34e::run(y.size(), threads, d7::kMinCostly, 15,
                     [&](std::size_t i, d34e::Row o) {
        const double a = Av[i], b = Bv[i];
        const double s = a + b;
        double C[6];
        bb34_power_sums(s, N, C, 6);
        double acc[15];
        for (int c = 0; c < 15; ++c) acc[c] = 0.0;
        for (int k = 0; k <= N; ++k) {
            const double f = bb34_mass(a, b, s, N, k);
            double A[6], B[6];
            bb34_power_sums(a, k, A, 6);
            bb34_power_sums(b, N - k, B, 6);
            const double L01 = 1.0 * (B[0] - C[0]);
            const double L02 = -1.0 * (B[1] - C[1]);
            const double L03 = 2.0 * (B[2] - C[2]);
            const double L04 = -6.0 * (B[3] - C[3]);
            const double L05 = 24.0 * (B[4] - C[4]);
            const double L06 = -120.0 * (B[5] - C[5]);
            const double L10 = 1.0 * (A[0] - C[0]);
            const double L11 = 1.0 * C[1];
            const double L12 = -2.0 * C[2];
            const double L13 = 6.0 * C[3];
            const double L14 = -24.0 * C[4];
            const double L15 = 120.0 * C[5];
            const double L20 = -1.0 * (A[1] - C[1]);
            const double L21 = -2.0 * C[2];
            const double L22 = 6.0 * C[3];
            const double L23 = -24.0 * C[4];
            const double L24 = 120.0 * C[5];
            const double L30 = 2.0 * (A[2] - C[2]);
            const double L31 = 6.0 * C[3];
            const double L32 = -24.0 * C[4];
            const double L33 = 120.0 * C[5];
            const double L40 = -6.0 * (A[3] - C[3]);
            const double L41 = -24.0 * C[4];
            const double L42 = 120.0 * C[5];
            const double L50 = 24.0 * (A[4] - C[4]);
            const double L51 = 120.0 * C[5];
            const double L60 = -120.0 * (A[5] - C[5]);
            acc[0] += f * (L60 + (3.0*(L20*L20*L20)) + (4.0*(L30*L30)) + (L20*(L10*L10*L10*L10)) + (4.0*L10*L50) + (4.0*L30*(L10*L10*L10)) + (6.0*L40*(L10*L10)) + (6.0*(L10*L10)*(L20*L20)) + (7.0*L20*L40) + (16.0*L10*L20*L30));
            acc[1] += f * (L51 + (L01*L50) + (L21*(L10*L10*L10)) + (3.0*L10*L41) + (3.0*L11*L40) + (3.0*L11*(L20*L20)) + (3.0*L31*(L10*L10)) + (4.0*L20*L31) + (4.0*L21*L30) + (L01*L20*(L10*L10*L10)) + (3.0*L01*L10*L40) + (3.0*L01*L10*(L20*L20)) + (3.0*L01*L30*(L10*L10)) + (3.0*L11*L20*(L10*L10)) + (4.0*L01*L20*L30) + (6.0*L10*L11*L30) + (6.0*L10*L20*L21));
            acc[2] += f * (L42 + (2.0*(L21*L21)) + (L02*L40) + (L02*(L20*L20)) + (L22*(L10*L10)) + (L40*(L01*L01)) + ((L01*L01)*(L20*L20)) + (2.0*L01*L41) + (2.0*L10*L32) + (2.0*L12*L30) + (2.0*L20*L22) + (2.0*L20*(L11*L11)) + (4.0*L11*L31) + (L02*L20*(L10*L10)) + (L20*(L01*L01)*(L10*L10)) + (2.0*L01*L21*(L10*L10)) + (2.0*L02*L10*L30) + (2.0*L10*L12*L20) + (2.0*L10*L30*(L01*L01)) + (4.0*L01*L10*L31) + (4.0*L01*L11*L30) + (4.0*L01*L20*L21) + (4.0*L10*L11*L21) + (4.0*L01*L10*L11*L20));
            acc[3] += f * (L33 + (L03*L30) + (L10*L23) + (L13*L20) + (L30*(L01*L01*L01)) + (3.0*L01*L32) + (3.0*L02*L31) + (3.0*L11*L22) + (3.0*L12*L21) + (3.0*L31*(L01*L01)) + (L03*L10*L20) + (L10*L20*(L01*L01*L01)) + (3.0*L01*L02*L30) + (3.0*L01*L10*L22) + (3.0*L01*L12*L20) + (3.0*L02*L10*L21) + (3.0*L02*L11*L20) + (3.0*L10*L21*(L01*L01)) + (3.0*L11*L20*(L01*L01)) + (6.0*L01*L11*L21) + (3.0*L01*L02*L10*L20));
            acc[4] += f * (L24 + (L04*L20) + (L20*(L01*L01*L01*L01)) + (3.0*L20*(L02*L02)) + (4.0*L01*L23) + (4.0*L03*L21) + (4.0*L21*(L01*L01*L01)) + (6.0*L02*L22) + (6.0*L22*(L01*L01)) + (4.0*L01*L03*L20) + (6.0*L02*L20*(L01*L01)) + (12.0*L01*L02*L21));
            acc[5] += f * (L42 + (L02*L40) + (L02*(L10*L10*L10*L10)) + (3.0*L02*(L20*L20)) + (4.0*L10*L32) + (4.0*L12*L30) + (4.0*L12*(L10*L10*L10)) + (6.0*L20*L22) + (6.0*L22*(L10*L10)) + (4.0*L02*L10*L30) + (6.0*L02*L20*(L10*L10)) + (12.0*L10*L12*L20));
            acc[6] += f * (L33 + (L01*L32) + (L02*L31) + (L03*L30) + (L03*(L10*L10*L10)) + (3.0*L10*L23) + (3.0*L11*L22) + (3.0*L12*L21) + (3.0*L13*L20) + (3.0*L13*(L10*L10)) + (L01*L02*L30) + (L01*L02*(L10*L10*L10)) + (3.0*L01*L10*L22) + (3.0*L01*L12*L20) + (3.0*L01*L12*(L10*L10)) + (3.0*L02*L10*L21) + (3.0*L02*L11*L20) + (3.0*L02*L11*(L10*L10)) + (3.0*L03*L10*L20) + (6.0*L10*L11*L12) + (3.0*L01*L02*L10*L20));
            acc[7] += f * (L24 + (2.0*(L12*L12)) + (L04*L20) + (L04*(L10*L10)) + (L20*(L02*L02)) + (L22*(L01*L01)) + ((L02*L02)*(L10*L10)) + (2.0*L01*L23) + (2.0*L02*L22) + (2.0*L02*(L11*L11)) + (2.0*L03*L21) + (2.0*L10*L14) + (4.0*L11*L13) + (L02*L20*(L01*L01)) + (L02*(L01*L01)*(L10*L10)) + (2.0*L01*L02*L21) + (2.0*L01*L03*L20) + (2.0*L01*L03*(L10*L10)) + (2.0*L10*L12*(L01*L01)) + (4.0*L01*L10*L13) + (4.0*L01*L11*L12) + (4.0*L02*L10*L12) + (4.0*L03*L10*L11) + (4.0*L01*L02*L10*L11));
            acc[8] += f * (L15 + (L05*L10) + (L12*(L01*L01*L01)) + (3.0*L01*L14) + (3.0*L04*L11) + (3.0*L11*(L02*L02)) + (3.0*L13*(L01*L01)) + (4.0*L02*L13) + (4.0*L03*L12) + (L02*L10*(L01*L01*L01)) + (3.0*L01*L04*L10) + (3.0*L01*L10*(L02*L02)) + (3.0*L02*L11*(L01*L01)) + (3.0*L03*L10*(L01*L01)) + (4.0*L02*L03*L10) + (6.0*L01*L02*L12) + (6.0*L01*L03*L11));
            acc[9] += f * (L06 + (3.0*(L02*L02*L02)) + (4.0*(L03*L03)) + (L02*(L01*L01*L01*L01)) + (4.0*L01*L05) + (4.0*L03*(L01*L01*L01)) + (6.0*L04*(L01*L01)) + (6.0*(L01*L01)*(L02*L02)) + (7.0*L02*L04) + (16.0*L01*L02*L03));
            acc[10] += f * (L51 + (L11*L40) + (L11*(L10*L10*L10*L10)) + (3.0*L11*(L20*L20)) + (4.0*L10*L41) + (4.0*L21*L30) + (4.0*L21*(L10*L10*L10)) + (6.0*L20*L31) + (6.0*L31*(L10*L10)) + (4.0*L10*L11*L30) + (6.0*L11*L20*(L10*L10)) + (12.0*L10*L20*L21));
            acc[11] += f * (L42 + (3.0*(L21*L21)) + (L01*L41) + (L12*L30) + (L12*(L10*L10*L10)) + (3.0*L10*L32) + (3.0*L20*L22) + (3.0*L20*(L11*L11)) + (3.0*L22*(L10*L10)) + (3.0*(L10*L10)*(L11*L11)) + (4.0*L11*L31) + (L01*L11*L30) + (L01*L11*(L10*L10*L10)) + (3.0*L01*L10*L31) + (3.0*L01*L20*L21) + (3.0*L01*L21*(L10*L10)) + (3.0*L10*L12*L20) + (9.0*L10*L11*L21) + (3.0*L01*L10*L11*L20));
            acc[12] += f * (L33 + (2.0*(L11*L11*L11)) + (L02*L31) + (L13*L20) + (L13*(L10*L10)) + (L31*(L01*L01)) + (2.0*L01*L32) + (2.0*L10*L23) + (4.0*L12*L21) + (5.0*L11*L22) + (L02*L11*L20) + (L02*L11*(L10*L10)) + (L11*L20*(L01*L01)) + (L11*(L01*L01)*(L10*L10)) + (2.0*L01*L12*L20) + (2.0*L01*L12*(L10*L10)) + (2.0*L02*L10*L21) + (2.0*L10*L21*(L01*L01)) + (4.0*L01*L10*L22) + (4.0*L01*L10*(L11*L11)) + (6.0*L01*L11*L21) + (6.0*L10*L11*L12));
            acc[13] += f * (L24 + (3.0*(L12*L12)) + (L03*L21) + (L10*L14) + (L21*(L01*L01*L01)) + (3.0*L01*L23) + (3.0*L02*L22) + (3.0*L02*(L11*L11)) + (3.0*L22*(L01*L01)) + (3.0*(L01*L01)*(L11*L11)) + (4.0*L11*L13) + (L03*L10*L11) + (L10*L11*(L01*L01*L01)) + (3.0*L01*L02*L21) + (3.0*L01*L10*L13) + (3.0*L02*L10*L12) + (3.0*L10*L12*(L01*L01)) + (9.0*L01*L11*L12) + (3.0*L01*L02*L10*L11));
            acc[14] += f * (L15 + (L04*L11) + (L11*(L01*L01*L01*L01)) + (3.0*L11*(L02*L02)) + (4.0*L01*L14) + (4.0*L03*L12) + (4.0*L12*(L01*L01*L01)) + (6.0*L02*L13) + (6.0*L13*(L01*L01)) + (4.0*L01*L03*L11) + (6.0*L02*L11*(L01*L01)) + (12.0*L01*L02*L12));
        }
        for (int c = 0; c < 15; ++c) o[c] = acc[c];
    });
}
