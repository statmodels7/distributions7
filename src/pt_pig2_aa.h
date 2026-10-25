// The alpha_alpha entry of the expected information of the Poisson-inverse Gaussian, family pig2
// (parameters mu, alpha), and its derivatives of orders 0 to 4 in (mu, alpha), from the series and the
// closed forms that do not sum over the support. Moved from src/d34e_pig2.cpp (orders 3 and 4) and
// extended to the orders 0, 1 and 2, which src/pig_hd.cpp uses for the expected information and its
// first and second derivatives.
//
// E_mu_mu = -1/mu + 1/sqrt(mu^2 + alpha^2) and E_mu_alpha = 0 for every (mu, alpha) (the parametrization
// is orthogonal), so only E_alpha_alpha needs more than a closed form. Its routes are
//
// (1) aa_series<K>: E_aa(mu, alpha) = sum_k sum_m c[k][m] mu^m alpha^-k, 4 <= k <= 70, 2 <= m <= k - 2
//     (pig2_aa_series.h), differentiated term by term. The series converges for t = (1 + mu)/alpha up to
//     about 0.5 once mu is of order 100 or more and is asymptotic for small mu. The sum is taken by
//     diagonals k and accepted when the last two diagonals are below 2e-11 of the largest component and
//     decreasing, or when three consecutive diagonals are below 1e-17 of it.
// (1b) aa_rows<K>: the same series summed by powers of 1/alpha at fixed u = mu/alpha (pig2_aa_rows.h for
//     K = 3, 4, pig2_aa_rows012.h for K = 0, 1, 2): each coefficient function is algebraic, a polynomial
//     over a power of 1 + u^2, accepted once alpha is above about 150 (mu >= 300), 100 (mu >= 100) or 70
//     (mu >= 30).
//
// Both return false when the sum is not accepted; the caller then sums over the support.
#ifndef DISTRIBUTIONS7_PT_PIG2_AA_H
#define DISTRIBUTIONS7_PT_PIG2_AA_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "pig2_aa_series.h"
#include "pig2_aa_rows.h"
#include "pig2_aa_rows012.h"

namespace pig2_aa {

// o[j] = d^K E_alpha_alpha / dmu^(K-j) dalpha^j, j = 0..K, from the series of
// pig2_aa_series.h differentiated term by term:
//   d^(i,j) [mu^m alpha^-k] = ff(m, i) (-1)^j rf(k, j) mu^(m-i) alpha^-(k+j),
// ff the falling and rf the rising factorial, i = K - j. With x = 1/alpha, u = mu x
// and z = 1/mu, the diagonal k of component j is
//   x^(K+2) u^(k-2-i) rf(k, j) (-1)^j sum_m c[k][m] ff(m, i) z^(k-2-m)     (mu >= 1),
//   x^(k+j) rf(k, j) (-1)^j mu^(m0-i) sum_m c[k][m] ff(m, i) mu^(m-m0)     (mu < 1),
// m running from m0 = max(2, i) to k - 2 and summed by Horner's rule in z <= 1 or
// mu < 1, so that no power mu^(k-2) is formed (mu^68 overflows at mu = 1e6).
// Returns false, o unspecified, when the series is not accepted (see the header).
static inline double aa_ff(int m, int i) {
    double r = 1.0;
    for (int q = 0; q < i; ++q) r *= (double) (m - q);
    return r;
}

template <int K>
static bool aa_series(double mu, double al, double* o) {
    using namespace d34e_pig2s;
    const double x = 1.0 / al;
    double xp[AA_KMAX + K + 4], up[AA_KMAX + 2], d[AA_KMAX + 1];
    xp[0] = 1.0;
    for (int n = 1; n < AA_KMAX + K + 4; ++n) xp[n] = xp[n - 1] * x;
    const double u = mu * x;
    up[0] = 1.0;
    for (int n = 1; n < AA_KMAX + 2; ++n) up[n] = up[n - 1] * u;
    const bool big = mu >= 1.0;
    const double z = big ? 1.0 / mu : 0.0;
    double S[K + 1];
    for (int j = 0; j <= K; ++j) S[j] = 0.0;
    for (int k = 0; k < 4; ++k) d[k] = 0.0;
    for (int k = 4; k <= AA_KMAX; ++k) {
        const double* c = AA_C + AA_OFF[k];              // c[m - 2], m = 2 .. k - 2
        double dmax = 0.0;
        for (int j = 0; j <= K; ++j) {
            const int i = K - j;
            const int m0 = i > 2 ? i : 2;
            double Dj = 0.0;
            if (m0 <= k - 2) {
                double rfk = 1.0;
                for (int q = 0; q < j; ++q) rfk *= (double) (k + q);
                double h = 0.0;
                if (big) {
                    for (int m = m0; m <= k - 2; ++m) h = h * z + c[m - 2] * aa_ff(m, i);
                    Dj = xp[K + 2] * up[k - 2 - i] * rfk * h;
                } else {
                    for (int m = k - 2; m >= m0; --m) h = h * mu + c[m - 2] * aa_ff(m, i);
                    for (int q = 0; q < m0 - i; ++q) h *= mu;
                    Dj = xp[k + j] * rfk * h;
                }
                if (j & 1) Dj = -Dj;
            }
            S[j] += Dj;
            dmax = std::max(dmax, std::fabs(Dj));
        }
        d[k] = dmax;
        if (k >= K + 8) {
            double sm = 0.0;
            for (int j = 0; j <= K; ++j) {
                if (!R_finite(S[j])) return false;
                sm = std::max(sm, std::fabs(S[j]));
            }
            if (sm > 0.0 && d[k] <= 1e-17 * sm && d[k - 1] <= 1e-17 * sm && d[k - 2] <= 1e-17 * sm) {
                for (int j = 0; j <= K; ++j) o[j] = S[j];
                return true;
            }
        }
    }
    double sm = 0.0;
    for (int j = 0; j <= K; ++j) {
        if (!R_finite(S[j])) return false;
        sm = std::max(sm, std::fabs(S[j]));
    }
    const double tail = d[AA_KMAX] + d[AA_KMAX - 1];
    const double prev = d[AA_KMAX - 2] + d[AA_KMAX - 3];
    if (!(tail <= 2e-11 * sm)) return false;
    if (!(tail <= 1e-15 * sm || tail <= 0.8 * prev)) return false;
    for (int j = 0; j <= K; ++j) o[j] = S[j];
    return true;
}

// The tables of the rows of order K: offsets, exponents s2 and coefficients.
template <int K> struct RowTab;
template <> struct RowTab<0> {
    static const int* off() { return d34e_pig2s::ROW0_OFF; }
    static const int* s2() { return d34e_pig2s::ROW0_S2; }
    static const double* c() { return d34e_pig2s::ROW0_C; }
};
template <> struct RowTab<1> {
    static const int* off() { return d34e_pig2s::ROW1_OFF; }
    static const int* s2() { return d34e_pig2s::ROW1_S2; }
    static const double* c() { return d34e_pig2s::ROW1_C; }
};
template <> struct RowTab<2> {
    static const int* off() { return d34e_pig2s::ROW2_OFF; }
    static const int* s2() { return d34e_pig2s::ROW2_S2; }
    static const double* c() { return d34e_pig2s::ROW2_C; }
};
template <> struct RowTab<3> {
    static const int* off() { return d34e_pig2s::ROW3_OFF; }
    static const int* s2() { return d34e_pig2s::ROW3_S2; }
    static const double* c() { return d34e_pig2s::ROW3_C; }
};
template <> struct RowTab<4> {
    static const int* off() { return d34e_pig2s::ROW4_OFF; }
    static const int* s2() { return d34e_pig2s::ROW4_S2; }
    static const double* c() { return d34e_pig2s::ROW4_C; }
};

// o[j] = d^K E_alpha_alpha / dmu^(K-j) dalpha^j, j = 0..K, from the same series summed
// the other way, by powers of 1/alpha at fixed u = mu/alpha (pig2_aa_rows.h):
//   E_ij = sum_{n=2}^{ROW_NMAX} alpha^-(n+K) sum_{parts} N(u) (1 + u^2)^(-s2/2),
// each N a polynomial of exact coefficients. For u <= 1 the polynomial is summed by
// Horner's rule in u; for u > 1 in v = 1/u, with N(u) = u^deg sum c_m v^(deg-m) and
// (1 + u^2)^(-s2/2) = u^-s2 (1 + v^2)^(-s2/2), so that no power of u beyond the
// bounded u^(deg-s2) is formed. The rows are the terms of an expansion in 1/alpha,
// asymptotic for small alpha and small mu. The acceptance test is that of aa_series: the
// last two rows below 2e-11 of the block's largest component and decreasing, or two
// consecutive rows below 1e-17 of it. Returns false, o unspecified, when the sum is not
// accepted.
template <int K>
static bool aa_rows(double mu, double al, double* o) {
    using namespace d34e_pig2s;
    const int* OFF = RowTab<K>::off();
    const int* S2 = RowTab<K>::s2();
    const double* C = RowTab<K>::c();
    const double x = 1.0 / al;
    const double u = mu * x;
    const bool big = u > 1.0;
    const double rho = 1.0 + u * u;
    const double v = big ? 1.0 / u : 0.0;
    const double w = 1.0 + v * v;
    double xn = 1.0;
    for (int q = 0; q < K + 1; ++q) xn *= x;           // x^(n+K) at n = 1
    double S[K + 1], d[ROW_NMAX + 1];
    for (int j = 0; j <= K; ++j) S[j] = 0.0;
    for (int n = 0; n <= ROW_NMAX; ++n) d[n] = 0.0;
    for (int n = 2; n <= ROW_NMAX; ++n) {
        xn *= x;
        double dm = 0.0;
        for (int j = 0; j <= K; ++j) {
            double tot = 0.0;
            for (int part = 0; part < 2; ++part) {
                const int idx = ((n - 2) * (K + 1) + j) * 2 + part;
                const double* c = C + OFF[idx];
                const int deg = OFF[idx + 1] - OFF[idx] - 1;
                const int s2 = S2[idx];
                double h = 0.0;
                if (!big) {
                    for (int m = deg; m >= 0; --m) h = h * u + c[m];
                    tot += h * std::pow(rho, -0.5 * s2);
                } else {
                    for (int m = 0; m <= deg; ++m) h = h * v + c[m];
                    tot += h * std::pow(u, (double) (deg - s2)) * std::pow(w, -0.5 * s2);
                }
            }
            tot *= xn;
            S[j] += tot;
            dm = std::max(dm, std::fabs(tot));
        }
        d[n] = dm;
        if (n >= 6) {
            double sm = 0.0;
            for (int j = 0; j <= K; ++j) {
                if (!R_finite(S[j])) return false;
                sm = std::max(sm, std::fabs(S[j]));
            }
            if (sm > 0.0 && d[n] <= 1e-17 * sm && d[n - 1] <= 1e-17 * sm) {
                for (int j = 0; j <= K; ++j) o[j] = S[j];
                return true;
            }
        }
    }
    double sm = 0.0;
    for (int j = 0; j <= K; ++j) {
        if (!R_finite(S[j])) return false;
        sm = std::max(sm, std::fabs(S[j]));
    }
    const double tail = d[ROW_NMAX] + d[ROW_NMAX - 1];
    const double prev = d[ROW_NMAX - 2] + d[ROW_NMAX - 3];
    if (!(tail <= 2e-11 * sm)) return false;
    if (!(tail <= 1e-15 * sm || tail <= 0.8 * prev)) return false;
    for (int j = 0; j <= K; ++j) o[j] = S[j];
    return true;
}

// o[j] = d^K E_mu_mu / dmu^(K-j) dalpha^j, j = 0..K, K = 0 to 4, with
// E_mu_mu = -1/mu + g, g = 1/r, r = sqrt(mu^2 + alpha^2). With c = mu/r and
// s = alpha/r, d^(K-j,j) g = r^-(K+1) Q(c, s) for j >= 1 (Q a polynomial
// of degree K, using c^2 + s^2 = 1), and for j = 0 d^K g = (-1)^K K! r^-(K+1) P_K(c),
// P_K the Legendre polynomial. The j = 0 component adds d^K(-1/mu) =
// -(-1)^K K! mu^-(K+1), which cancels most of it when alpha << mu; with
// delta = 1 - c = u^2 / (sqrt(1 + u^2) (1 + sqrt(1 + u^2))), u = alpha/mu, the sum is
// (-1)^(K+1) K! mu^-(K+1) [ (1 - c^(K+1)) + c^(K+1) (1 - P_K(1 - delta)) ],
// both terms non-negative, the first by expm1 and the second by a polynomial
// in delta that has no constant term (zero for K = 0).
static inline void pig2_mm_orders(const int K, const double mu, const double al, double* o) {
    const double r = std::hypot(mu, al);
    const double c = mu / r, s = al / r;
    const double u2 = (al / mu) * (al / mu);
    const double lg = std::log1p(u2);                    // 2 log(1 / c)
    double delta = 1.0;                                  // c -> 0 when u2 overflows
    if (R_finite(u2)) {
        const double sc = std::sqrt(1.0 + u2);
        delta = u2 / (sc * (1.0 + sc));
    }
    const double ck = std::exp(-0.5 * (K + 1) * lg);                // c^(K+1)
    const double omc = -std::expm1(-0.5 * (K + 1) * lg);            // 1 - c^(K+1)
    const double rk = std::pow(r, -(K + 1));
    const double c2 = c * c, s2 = s * s;
    if (K == 0) {
        o[0] = -omc / mu;
    } else if (K == 1) {
        o[0] = (omc + ck * delta) / (mu * mu);
        o[1] = -rk * s;
    } else if (K == 2) {
        const double B = delta * (3.0 - 1.5 * delta);
        o[0] = -2.0 * std::pow(mu, -3) * (omc + ck * B);
        o[1] = rk * 3.0 * c * s;
        o[2] = rk * (2.0 * s2 - c2);
    } else if (K == 3) {
        const double B = delta * (6.0 + delta * (-7.5 + delta * 2.5));
        o[0] = 6.0 * std::pow(mu, -4) * (omc + ck * B);
        o[1] = rk * s * (3.0 * s2 - 12.0 * c2);
        o[2] = rk * c * (3.0 * c2 - 12.0 * s2);
        o[3] = rk * s * (9.0 * c2 - 6.0 * s2);
    } else {
        const double B = delta * (10.0 + delta * (-22.5 + delta * (17.5 - 4.375 * delta)));
        o[0] = -24.0 * std::pow(mu, -5) * (omc + ck * B);
        o[1] = rk * c * s * (60.0 * c2 - 45.0 * s2);
        o[2] = rk * (-12.0 * c2 * c2 + 81.0 * c2 * s2 - 12.0 * s2 * s2);
        o[3] = rk * c * s * (-45.0 * c2 + 60.0 * s2);
        o[4] = rk * (9.0 * c2 * c2 - 72.0 * c2 * s2 + 24.0 * s2 * s2);
    }
}

}  // namespace pig2_aa

#endif
