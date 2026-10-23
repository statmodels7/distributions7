#ifndef D7_PT_PIG_H
#define D7_PT_PIG_H

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "pig_asym.h"
#include "pt_sqrt.h"
#include "pt_pig_debye.h"
#include "pt_pig_mass.h"

// The Poisson-inverse Gaussian in both parametrizations: what the kernels of
// pig_hd.cpp and the scalar registry share. The special functions are
// log S_y and its derivatives, from psi_derivs() and pig_logS_w() at one
// response or from the recurrences PigSeq and PigSeqW over the support;
// pig1_slot() and the pig2 component functions below take them as
// arguments and compute one partial each, with the arithmetic of the row
// kernels pig1_row() and pig2_row(), slot for slot.

namespace d7 {

// psi(alpha) = -alpha + log S_y(alpha) and four derivatives in alpha, where
// S_y(alpha) = sum_{k<y} t_k, t_0 = 1, t_{k+1} = t_k w (y + k)(y - k - 1)/(k+1),
// w = 1/(2 alpha), is the Bessel polynomial; S_y = e^alpha sqrt(2 alpha/pi)
// K_{y-1/2}(alpha). Two routes, by the count:
//
//   y < PIG_Y_DEBYE: the sum itself, at most 20 terms, psi_derivs_sum().
//   y >= PIG_Y_DEBYE: the uniform (Debye) expansion of K_nu(alpha), nu = y - 1/2,
//     with every derivative in closed form, psi_derivs_debye(): the cost does
//     not depend on y, which may reach 2^53.
//
// Both write h[0..4] whatever the order, the unread ones as zero, and
// logS = log S_y and A1 = d log S_y / d alpha when asked for.

// The sum, for 1 <= y < 21. With w <= 1 (alpha >= 1/2) the terms are formed
// forward from t_0 = 1; with w > 1 they are formed backward from the last
// term, t_{n-1} = 1 in the units of the anchor, whose log is carried apart,
// so no power of w is formed. The derivatives of log S in alpha are
// combinations with positive coefficients of the cumulants kappa_j of k
// under the weights t_k / S (s = -log alpha):
//   d log S   = -kappa_1 / a
//   d^2       = (kappa_2 + kappa_1) / a^2
//   d^3       = -(kappa_3 + 3 kappa_2 + 2 kappa_1) / a^3
//   d^4       = (kappa_4 + 6 kappa_3 + 11 kappa_2 + 6 kappa_1) / a^4.
// The cumulants of order two to four are taken about the mean, which is
// found first; the rising-factorial moments of the former code cancel at
// order four.
const double PIG_Y_DEBYE = 21.0;

inline void psi_derivs_sum(double y, double alpha, int order, double h[5],
                           double* logS, double* A1out) {
    const int n = (int) y;
    double t[20];
    double anchor = 0.0;
    if (alpha >= 0.5) {
        const double w = 0.5 / alpha;
        t[0] = 1.0;
        for (int k = 0; k < n - 1; ++k)
            t[k + 1] = t[k] * w * (y + k) * (y - k - 1.0) / (k + 1.0);
    } else {
        t[n - 1] = 1.0;
        for (int k = n - 2; k >= 0; --k)
            t[k] = t[k + 1] * (2.0 * alpha * (k + 1.0)) /
                ((y + k) * (y - k - 1.0));
        anchor = R::lgammafn(y + (n - 1)) - R::lgammafn(y - (n - 1)) -
            R::lgammafn((double) n) - (n - 1) * std::log(2.0 * alpha);
    }
    // S - t_0, summed from the smallest terms; forward t_0 = 1, and
    // log S = log1p(S - 1) keeps the digits when S is close to one
    // (a large alpha), where log(S) loses them.
    double Sr = 0.0, s1 = 0.0;
    for (int k = n - 1; k >= 1; --k) {
        Sr += t[k];
        s1 += k * t[k];
    }
    const double S = t[0] + Sr;
    const double k1 = s1 / S;
    double k2 = 0.0, k3 = 0.0, k4 = 0.0;
    if (order >= 2) {
        double c2 = 0.0, c3 = 0.0, c4 = 0.0;
        for (int k = 0; k < n; ++k) {
            const double d = k - k1, d2 = d * d;
            c2 += t[k] * d2;
            if (order < 3) continue;
            c3 += t[k] * d2 * d;
            if (order < 4) continue;
            c4 += t[k] * d2 * d2;
        }
        c2 /= S; c3 /= S; c4 /= S;
        k2 = c2; k3 = c3; k4 = c4 - 3.0 * c2 * c2;
    }
    const double ls = (alpha >= 0.5) ? std::log1p(Sr) : anchor + std::log(S);
    const double ia = 1.0 / alpha;
    if (logS) *logS = ls;
    if (A1out) *A1out = -k1 * ia;
    h[0] = -alpha + ls;
    h[1] = h[2] = h[3] = h[4] = 0.0;
    if (order < 1) return;
    h[1] = -1.0 - k1 * ia;
    if (order < 2) return;
    h[2] = (k2 + k1) * ia * ia;
    if (order < 3) return;
    h[3] = -(k3 + 3.0 * k2 + 2.0 * k1) * ia * ia * ia;
    if (order < 4) return;
    h[4] = (k4 + 6.0 * k3 + 11.0 * k2 + 6.0 * k1) * ia * ia * ia * ia;
}

// The uniform expansion, for y >= 21. With nu = y - 1/2, s = sqrt(nu^2 +
// alpha^2), c = nu / s, d = alpha / s and the Debye polynomial
// P(c) = sum_{k<=K} (-1)^k u_k(c) nu^-k,
//   psi = -(1/4) log(1 + nu^2/alpha^2) - nu eta + log P(c),
//   eta = s/nu + log(alpha / (nu + s)),
// the first two terms being elementary in alpha (the A and B pieces below
// are the derivatives of -nu eta and of the logarithm), and the derivatives
// of log P(c(alpha)) following from the cumulants of P in c and the
// derivatives of c in alpha (Faa di Bruno written out). K is the smallest
// number of terms for which the first omitted one is below 1e-17
// (PIG_UK_NU); it falls from 15 at nu = 20.5 to 1 at nu = 1e9. log S and
// A1 are formed without cancellation: -alpha + nu eta, of size nu, would
// swallow log S, so log S is assembled from -nu^2/(alpha + s), a log1p and
// log1p(P - 1), and A1 from -nu^2/(alpha (s + alpha)).
inline void psi_derivs_debye(double y, double alpha, int order, double h[5],
                             double* logS, double* A1out) {
    const double nu = y - 0.5, x = alpha;
    const double s = std::hypot(nu, x);
    const double c = nu / s, d = x / s;
    const double u = 1.0 / x, v = 1.0 / s;
    const double n2 = nu * nu;
    int K = PIG_KD;
    for (int k = 0; k <= PIG_KD; ++k)
        if (nu > PIG_UK_NU[k]) { K = k; break; }
    // P(c) as a polynomial in c, then P and its derivatives at c by Horner
    double P[3 * PIG_KD + 1];
    for (int e = 0; e <= 3 * K; ++e) P[e] = 0.0;
    const double inu = 1.0 / nu;
    double f = 1.0;
    for (int k = 0; k <= K; ++k) {
        const double* uk = PIG_U + PIG_UOFF[k];
        for (int e = 0; e <= 3 * k; ++e) P[e] += f * uk[e];
        f *= -inu;
    }
    const int nd = (order < 1) ? 1 : order;
    double q[5] = {0.0, 0.0, 0.0, 0.0, 0.0};   // Taylor coefficients in c
    for (int e = 3 * K; e >= 0; --e) {
        for (int i = nd; i >= 1; --i) q[i] = q[i] * c + q[i - 1];
        q[0] = q[0] * c + P[e];
    }
    const double p0 = q[0];
    // the same Horner without the constant term, for log1p(p0 - 1)
    double acc = 0.0;
    for (int e = 3 * K; e >= 1; --e) acc = acc * c + P[e];
    const double q0 = acc * c;
    const double r1 = q[1] / p0;
    const double t1 = -c * d * v;
    const double g1 = r1 * t1;
    // log(1 + nu^2/alpha^2), without overflow of the square
    const double qn = nu * u;
    const double lq = (qn > 1e100) ? 2.0 * std::log(qn) : std::log1p(qn * qn);
    if (logS) {
        const double dd = (nu + n2 / (s + x)) * u;
        *logS = -n2 / (x + s) + nu * std::log1p(dd) - 0.25 * lq +
            std::log1p(q0);
    }
    if (A1out)
        *A1out = -n2 * u / (s + x) + 0.5 * c * c * u + g1;
    const double sn = s / nu;
    const double eta = sn + std::log((x / nu) / (1.0 + sn));
    h[0] = -0.25 * lq - nu * eta + std::log(p0);
    h[1] = h[2] = h[3] = h[4] = 0.0;
    if (order < 1) return;
    const double E1 = -s * u + 0.5 * c * c * u;
    h[1] = E1 + g1;
    if (order < 2) return;
    const double nc = nu * c, d2 = d * d;
    const double r2 = 2.0 * q[2] / p0;
    const double l2 = r2 - r1 * r1;
    const double t2 = c * (3.0 * d2 - 1.0) * v * v;
    const double A1e = nc * u * u;
    const double B2 = -0.5 * c * c * (1.0 + 2.0 * d2) * u * u;
    h[2] = A1e + B2 + (r1 * t2 + l2 * t1 * t1);
    if (order < 3) return;
    const double r3 = 6.0 * q[3] / p0;
    const double l3 = r3 - 3.0 * r1 * r2 + 2.0 * r1 * r1 * r1;
    const double t3 = -3.0 * c * d * (5.0 * d2 - 3.0) * v * v * v;
    const double d4 = d2 * d2;
    const double A2e = -nc * (2.0 + d2) * u * u * u;
    const double B3 = c * c * (1.0 + d2 + 4.0 * d4) * u * u * u;
    h[3] = A2e + B3 + (r1 * t3 + 3.0 * l2 * t1 * t2 + l3 * t1 * t1 * t1);
    if (order < 4) return;
    const double r4 = 24.0 * q[4] / p0;
    const double l4 = r4 - 4.0 * r1 * r3 - 3.0 * r2 * r2 +
        12.0 * r1 * r1 * r2 - 6.0 * r1 * r1 * r1 * r1;
    const double t4 = 3.0 * c * (3.0 - 30.0 * d2 + 35.0 * d4) * v * v * v * v;
    const double A3e = 3.0 * nc * (2.0 + d2 + d4) * u * u * u * u;
    const double B4 = -3.0 * c * c * (1.0 + d2 + 8.0 * d4 * d2) *
        u * u * u * u;
    h[4] = A3e + B4 + (r1 * t4 + l2 * (4.0 * t1 * t3 + 3.0 * t2 * t2) +
        6.0 * l3 * t1 * t1 * t2 + l4 * t1 * t1 * t1 * t1);
}

inline void psi_derivs(double y, double alpha, int order, double h[5],
                       double* logS = nullptr, double* A1out = nullptr) {
    if (y < 1.0) {
        h[0] = -alpha; h[1] = -1.0; h[2] = h[3] = h[4] = 0.0;
        if (logS) *logS = 0.0;
        if (A1out) *A1out = 0.0;
        return;
    }
    if (y < PIG_Y_DEBYE) psi_derivs_sum(y, alpha, order, h, logS, A1out);
    else psi_derivs_debye(y, alpha, order, h, logS, A1out);
}

// The block's column layout, which the wrappers and the guard read:
// order k occupies PIG_LEN[k] columns starting at PIG_OFF[k].
const int PIG_OFF[5] = {0, 1, 3, 6, 10};
const int PIG_LEN[5] = {1, 2, 3, 4, 5};

// The partials of a bivariate function in (mu, sigma), orders 0 to 4, sit in
// fifteen slots: slot PIG_OFF[i + j] + j holds d^{i+j} / dmu^i dsigma^j.
inline int pix(int i, int j) { return PIG_OFF[i + j] + j; }

// h(F(mu, sigma)) for a univariate h with derivatives h[0..4] at F's value
// and F's partials in F[15]: Faa di Bruno written out per component, the
// same expressions the pig2 row composes psi with.
inline void pig_compose(const double h[5], const double* F, int order,
                        double* out) {
    out[0] = h[0];
    if (order < 1) return;
    double m = F[1], u = F[2];
    out[1] = h[1] * m;
    out[2] = h[1] * u;
    if (order < 2) return;
    double mm = F[3], mu = F[4], uu = F[5];
    out[3] = h[2] * m * m + h[1] * mm;
    out[4] = h[2] * m * u + h[1] * mu;
    out[5] = h[2] * u * u + h[1] * uu;
    if (order < 3) return;
    double mmm = F[6], mmu = F[7], muu = F[8], uuu = F[9];
    out[6] = h[3] * m * m * m + 3.0 * h[2] * mm * m + h[1] * mmm;
    out[7] = h[3] * m * m * u + h[2] * (mm * u + 2.0 * mu * m) + h[1] * mmu;
    out[8] = h[3] * m * u * u + h[2] * (uu * m + 2.0 * mu * u) + h[1] * muu;
    out[9] = h[3] * u * u * u + 3.0 * h[2] * uu * u + h[1] * uuu;
    if (order < 4) return;
    out[10] = h[4] * m * m * m * m + 6.0 * h[3] * mm * m * m +
        h[2] * (3.0 * mm * mm + 4.0 * mmm * m) + h[1] * F[10];
    out[11] = h[4] * m * m * m * u +
        h[3] * (3.0 * mu * m * m + 3.0 * mm * m * u) +
        h[2] * (3.0 * mm * mu + 3.0 * mmu * m + mmm * u) + h[1] * F[11];
    out[12] = h[4] * m * m * u * u +
        h[3] * (mm * u * u + uu * m * m + 4.0 * mu * m * u) +
        h[2] * (mm * uu + 2.0 * mu * mu + 2.0 * muu * m + 2.0 * mmu * u) +
        h[1] * F[12];
    out[13] = h[4] * m * u * u * u +
        h[3] * (3.0 * mu * u * u + 3.0 * uu * m * u) +
        h[2] * (3.0 * uu * mu + 3.0 * muu * u + uuu * m) + h[1] * F[13];
    out[14] = h[4] * u * u * u * u + 6.0 * h[3] * uu * u * u +
        h[2] * (3.0 * uu * uu + 4.0 * uuu * u) + h[1] * F[14];
}

// log S_y as a function of w = 1/(2 alpha), with its four derivatives in w,
// L[r] = d^r log S_y / dw^r. S_y(w) = sum_k a_k w^k is a polynomial with
// positive coefficients a_k = (y+k-1)! / ((y-k)! k!), so log S and its
// derivatives are finite as w -> 0 (the Poisson limit). In alpha the same
// quantities carry alpha^-r and the chain onto sigma cancels at a small
// dispersion, which is why the derivatives are taken in w directly. Two
// routes, by the count, as psi_derivs() has for pig2:
//
//   y < PIG_Y_DEBYE: the sum itself, pig_logS_w_sum();
//   y >= PIG_Y_DEBYE: the uniform (Debye) expansion, pig_logS_w_debye(),
//     whose cost does not depend on y, which may reach 2^53.
//
// Only the orders up to `order` are formed; L[0] = log S_y.

// The switch of the sum between its two forms, on X = y (y - 1) w, the mean
// of k under the terms of S_y at small w. Measured against 60 digits for
// y <= 20: the Taylor form is better below it, the cumulants about the mean
// above it, and their errors cross at X between 1 and 1.5.
const double PIG_W_SWITCH = 1.25;

// The sum, for 1 <= y < 21: a_0 = 1, a_{k+1} = a_k (y + k)(y - k - 1)/(k + 1)
// for k < y - 1, held in a fixed array.
//
// X <= PIG_W_SWITCH: the Taylor coefficients c_r = S^(r)/r! at w by Horner's
// scheme with its derivatives, on Q(x) = sum_{k>=1} a_k x^(k-1) and then
// S = 1 + x Q. Every term is positive, no power of w is divided out, and
// log S = log1p(S - 1) keeps the digits at a small w. The derivatives of
// log S are the cumulants of the ratios n_r = r! c_r / c_0,
//   L1 = n1, L2 = n2 - n1^2, L3 = n3 - 3 n1 n2 + 2 n1^3, ...
// which lose about y^(r-1) of the relative accuracy as w -> 0.
//
// X > PIG_W_SWITCH: the cumulants kappa_j of k under the weights t_k / S,
// t_k = a_k w^k, taken about the mean, which is found first:
//   L1 = kappa_1 / w
//   L2 = (kappa_2 - kappa_1) / w^2
//   L3 = (kappa_3 - 3 kappa_2 + 2 kappa_1) / w^3
//   L4 = (kappa_4 - 6 kappa_3 + 11 kappa_2 - 6 kappa_1) / w^4.
// These cancel as w -> 0 and are accurate above the switch. For w > 1 the
// weights are formed from the last term down, t_k = a_k w^-(n-1-k), and
// log S = (n - 1) log w + log(sum t_k), so no power of w is formed.
inline void pig_logS_w_sum(double y, double w, int order, double L[5]) {
    const int n = (int) y;
    double a[20];
    a[0] = 1.0;
    for (int k = 0; k < n - 1; ++k)
        a[k + 1] = a[k] * (y + k) * (y - k - 1.0) / (k + 1.0);
    if (y * (y - 1.0) * w <= PIG_W_SWITCH) {
        double dd[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        for (int k = n - 1; k >= 1; --k) {
            for (int r = order; r >= 1; --r) dd[r] = dd[r] * w + dd[r - 1];
            dd[0] = dd[0] * w + a[k];
        }
        const double Sr = w * dd[0];
        for (int r = order; r >= 1; --r) dd[r] = dd[r] * w + dd[r - 1];
        dd[0] = dd[0] * w + 1.0;
        L[0] = std::log1p(Sr);
        if (order < 1) return;
        const double n1 = dd[1] / dd[0];
        L[1] = n1;
        if (order < 2) return;
        const double n2 = 2.0 * dd[2] / dd[0];
        L[2] = n2 - n1 * n1;
        if (order < 3) return;
        const double n3 = 6.0 * dd[3] / dd[0];
        L[3] = n3 - 3.0 * n1 * n2 + 2.0 * n1 * n1 * n1;
        if (order < 4) return;
        const double n4 = 24.0 * dd[4] / dd[0];
        L[4] = n4 - 4.0 * n1 * n3 - 3.0 * n2 * n2 + 12.0 * n1 * n1 * n2 -
            6.0 * n1 * n1 * n1 * n1;
        return;
    }
    double t[20] = {0.0};
    if (w <= 1.0) {
        double pw = 1.0;
        for (int k = 0; k < n; ++k) { t[k] = a[k] * pw; pw *= w; }
    } else {
        const double v = 1.0 / w;
        double u = 1.0;
        for (int k = n - 1; k >= 0; --k) { t[k] = a[k] * u; u *= v; }
    }
    double Sr = 0.0, s1 = 0.0;
    for (int k = n - 1; k >= 1; --k) {
        Sr += t[k];
        s1 += k * t[k];
    }
    const double S = t[0] + Sr;
    const double k1 = s1 / S;
    L[0] = (w <= 1.0) ? std::log1p(Sr) : (n - 1) * std::log(w) + std::log(S);
    if (order < 1) return;
    const double iw = 1.0 / w;
    double k2 = 0.0, k3 = 0.0, k4 = 0.0;
    if (order >= 2) {
        double c2 = 0.0, c3 = 0.0, c4 = 0.0;
        for (int k = 0; k < n; ++k) {
            const double d = k - k1, d2 = d * d;
            c2 += t[k] * d2;
            if (order < 3) continue;
            c3 += t[k] * d2 * d;
            if (order < 4) continue;
            c4 += t[k] * d2 * d2;
        }
        c2 /= S; c3 /= S; c4 /= S;
        k2 = c2; k3 = c3; k4 = c4 - 3.0 * c2 * c2;
    }
    L[1] = k1 * iw;
    if (order < 2) return;
    L[2] = (k2 - k1) * iw * iw;
    if (order < 3) return;
    L[3] = (k3 - 3.0 * k2 + 2.0 * k1) * iw * iw * iw;
    if (order < 4) return;
    L[4] = (k4 - 6.0 * k3 + 11.0 * k2 - 6.0 * k1) * iw * iw * iw * iw;
}

// The uniform expansion, for y >= 21, the same as psi_derivs_debye() but
// differentiated in w. With nu = y - 1/2, z = 2 nu w = nu / alpha, r =
// sqrt(1 + z^2), c = z / r, d = 1 / r and P(c) the Debye polynomial of
// psi_derivs_debye(),
//   log S = nu f(z) - (1/4) log(1 + z^2) + log P(c),
//   f(z) = asinh z - z / (1 + r),
// a function of z and nu alone, so that d/dw = 2 nu d/dz and
//   L_k = (2 nu)^k [nu f^(k)(z) + h^(k)(z) + (log P(c(z)))^(k)],
// h = -(1/4) log(1 + z^2). The derivatives of f, h and c in z are
// written out in c and d (f' = g = 1 / (1 + r) = d / (1 + d), every term
// free of cancellation at a small z), those of log P follow from the
// cumulants of P in c (Faa di Bruno written out). There is no chain from
// alpha, whose terms cancel as w -> 0.
inline void pig_logS_w_debye(double y, double w, int order, double L[5]) {
    const double nu = y - 0.5;
    const double z = 2.0 * nu * w;
    const double r = std::hypot(1.0, z);
    const double d = 1.0 / r;
    const double c = (z > 1e150) ? 1.0 : z * d;
    int K = PIG_KD;
    for (int k = 0; k <= PIG_KD; ++k)
        if (nu > PIG_UK_NU[k]) { K = k; break; }
    double P[3 * PIG_KD + 1];
    for (int e = 0; e <= 3 * K; ++e) P[e] = 0.0;
    const double inu = 1.0 / nu;
    double f = 1.0;
    for (int k = 0; k <= K; ++k) {
        const double* uk = PIG_U + PIG_UOFF[k];
        for (int e = 0; e <= 3 * k; ++e) P[e] += f * uk[e];
        f *= -inu;
    }
    const int nd = (order < 1) ? 1 : order;
    double q[5] = {0.0, 0.0, 0.0, 0.0, 0.0};   // Taylor coefficients in c
    for (int e = 3 * K; e >= 0; --e) {
        for (int i = nd; i >= 1; --i) q[i] = q[i] * c + q[i - 1];
        q[0] = q[0] * c + P[e];
    }
    const double p0 = q[0];
    // the same Horner without the constant term, for log1p(p0 - 1)
    double acc = 0.0;
    for (int e = 3 * K; e >= 1; --e) acc = acc * c + P[e];
    const double q0 = acc * c;
    const double g = d / (1.0 + d);
    // asinh z and log(1 + z^2) without overflow of the square
    double asz, lq;
    if (z > 1e100) {
        const double lz = std::log(2.0 * nu) + std::log(w);
        asz = 0.6931471805599453 + lz;
        lq = 2.0 * lz;
    } else {
        asz = std::asinh(z);
        lq = std::log1p(z * z);
    }
    L[0] = nu * (asz - c / (1.0 + d)) - 0.25 * lq + std::log1p(q0);
    L[1] = L[2] = L[3] = L[4] = 0.0;
    if (order < 1) return;
    const double two = 2.0 * nu;
    const double d2 = d * d, c2 = c * c;
    const double r1 = q[1] / p0;
    const double t1 = d2 * d;
    const double f1 = g, h1 = -0.5 * c * d;
    L[1] = two * (nu * f1 + h1 + r1 * t1);
    if (order < 2) return;
    const double r2 = 2.0 * q[2] / p0;
    const double l2 = r2 - r1 * r1;
    const double t2 = -3.0 * c * d2 * d2;
    const double f2 = -c * g * g;
    const double h2 = -0.5 * d2 * (d2 - c2);
    L[2] = two * two * (nu * f2 + h2 + (r1 * t2 + l2 * t1 * t1));
    if (order < 3) return;
    const double r3 = 6.0 * q[3] / p0;
    const double l3 = r3 - 3.0 * r1 * r2 + 2.0 * r1 * r1 * r1;
    const double t3 = (12.0 * c2 - 3.0 * d2) * d2 * d2 * d;
    const double f3 = g * g * (2.0 * c2 * g - d2 * d);
    const double h3 = c * d2 * d * (3.0 * d2 - c2);
    L[3] = two * two * two * (nu * f3 + h3 +
        (r1 * t3 + 3.0 * l2 * t1 * t2 + l3 * t1 * t1 * t1));
    if (order < 4) return;
    const double r4 = 24.0 * q[4] / p0;
    const double l4 = r4 - 4.0 * r1 * r3 - 3.0 * r2 * r2 +
        12.0 * r1 * r1 * r2 - 6.0 * r1 * r1 * r1 * r1;
    const double d4 = d2 * d2;
    const double t4 = 15.0 * c * d4 * d2 * (3.0 * d2 - 4.0 * c2);
    const double f4 = c * g * g * (3.0 * d4 + 6.0 * d2 * d * g -
        6.0 * c2 * g * g);
    const double h4 = 3.0 * d4 * (d4 - 6.0 * c2 * d2 + c2 * c2);
    L[4] = two * two * two * two * (nu * f4 + h4 +
        (r1 * t4 + l2 * (4.0 * t1 * t3 + 3.0 * t2 * t2) +
         6.0 * l3 * t1 * t1 * t2 + l4 * t1 * t1 * t1 * t1));
}

inline void pig_logS_w(double y, double w, int order, double L[5]) {
    L[0] = L[1] = L[2] = L[3] = L[4] = 0.0;
    if (y < 1.0) return;
    if (y < PIG_Y_DEBYE) pig_logS_w_sum(y, w, order, L);
    else pig_logS_w_debye(y, w, order, L);
}

// ---------------------------------------------------------------------------
// The expected information, exactly, by one pass over the support.
// ---------------------------------------------------------------------------
//
// E[l_ab] is a sum over y of the mass times the observed component, and
// each observed component reads log S_y and the rising-factorial moments
// m_r of k under the terms of S_y. Computed afresh at each y by
// psi_derivs() those cost y terms apiece, so the sum over the support is
// quadratic in its length. They satisfy a recurrence in y instead. With
// S_y(a) = sqrt(2a/pi) e^a K_{y-1/2}(a), the Bessel recurrence
// K_{v+1} = K_{v-1} + (2v/a) K_v reads
//
//   S_{y+1} = S_{y-1} + ((2y - 1)/a) S_y,          S_0 = S_1 = 1,
//
// and differentiating it r times in a gives, for the moments M_r = (-a)^r
// S^(r) = m_r S, a recurrence of positive terms. The derivatives of log S are
// not formed from those moments: their cumulants cancel, by about y^r, when
// the law of the index k of the terms of S_y is close to a Poisson law (a
// large a), which is the range of the expected information beyond
// PIG_ASYM_CUT at a large mu.
//
// As a function of u = 1/a, S_y is the polynomial that PigSeqW carries in w,
// with w replaced by u/2, since S_{y+1} = S_{y-1} + (2y - 1) u S_y. What is
// carried is therefore the same: the factorial cumulants f_r, r = 1..4, of
// the law of k (the derivatives at x = 1 of g_y(x) = log(S_y(u x) / S_y(u)))
// and their steps a_r = (g_y - g_{y-1})^(r)(1), with
//
//   S_{y+1}(u x) / S_{y+1}(u) = p S_{y-1}(u x) / S_{y-1}(u) +
//                               v x S_y(u x) / S_y(u),
//
// p = S_{y-1}/S_{y+1}, v = (2y - 1) u S_y/S_{y+1}, p + v = 1, so that
// g_{y+1} - g_y = log Z, Z(x) = p exp(-(g_y - g_{y-1})) + v x. The derivatives
// of Z at x = 1 are z_1 = v - p a_1 and z_r = p E_r (r >= 2), E_r being the
// derivatives of exp(-(g_y - g_{y-1})) (Bell polynomials in the a_r); the new
// steps are the derivatives of log Z, and f_r(y + 1) = f_r(y) + a_r(y + 1).
//
// The derivatives of log S in a follow from the f_r. With q_r =
// (-1)^r d^r(log S)/da^r = u^r F_r, F_1 = f_1 and F_{r+1} = r F_r + u d F_r/du,
// where u d f_j/du = j f_j + f_{j+1}, so that
//
//   F_1 = f_1,  F_2 = 2 f_1 + f_2,  F_3 = 6 f_1 + 6 f_2 + f_3,
//   F_4 = 24 f_1 + 36 f_2 + 12 f_3 + f_4
//
// (Lah numbers, all positive). The derivatives of psi = -a + log S are
// -1 - u F_1, u^2 F_2, -u^3 F_3, u^4 F_4. No raw moment enters, and every
// call costs O(1). The steps themselves cancel by about y in the same range.
//
// The ratio S_{y+1}/S_y - 1 is carried directly, the step to log S being a
// log1p: formed as rho + t - 1 it cancels at a large a, where it is of size
// y/a.
struct PigSeq {
    double a, u;           // alpha and 1/alpha
    int y;
    double st[4];          // steps of f_r, r = 1..4, from y - 1 to y
    double f[4];           // f_r at y
    double logS;           // log S_y
    double e;              // S_y / S_{y-1} - 1
    explicit PigSeq(double alpha)
        : a(alpha), u(1.0 / alpha), y(0), logS(0.0), e(0.0) {
        for (int r = 0; r < 4; ++r) st[r] = f[r] = 0.0;
    }
    // advance from y to y + 1
    void next() {
        if (y == 0) { y = 1; e = 0.0; return; }   // S_1 = S_0 = 1
        const double t = (2.0 * y - 1.0) * u;
        const double d = t - e / (1.0 + e);       // S_{y+1}/S_y - 1
        const double R = 1.0 + d;
        const double p = 1.0 / (R * (1.0 + e));   // S_{y-1}/S_{y+1}
        const double v = t / R;                   // (t S_y)/S_{y+1}
        const double a1 = st[0], a2 = st[1], a3 = st[2], a4 = st[3];
        const double a11 = a1 * a1;
        const double z1 = v - p * a1;
        const double z2 = p * (a11 - a2);
        const double z3 = p * (a1 * (3.0 * a2 - a11) - a3);
        const double z4 = p * (a11 * (a11 - 6.0 * a2) + 4.0 * a1 * a3 +
                               3.0 * a2 * a2 - a4);
        const double z11 = z1 * z1;
        const double s2 = z2 - z11;
        const double s3 = z3 - z1 * (3.0 * z2 - 2.0 * z11);
        const double s4 = z4 - 4.0 * z1 * z3 - z2 * (3.0 * z2 - 12.0 * z11) -
            6.0 * z11 * z11;
        st[0] = z1; st[1] = s2; st[2] = s3; st[3] = s4;
        f[0] += z1; f[1] += s2; f[2] += s3; f[3] += s4;
        logS += std::log1p(d);
        e = d;
        ++y;
    }
    // what psi_derivs() hands a row, with the same meaning:
    // h[0..4], then log S, then A1 = d log S / d a
    void pre(double out[7]) const {
        const double u2 = u * u;
        out[0] = -a + logS;
        out[1] = -1.0 - u * f[0];
        out[2] = u2 * (2.0 * f[0] + f[1]);
        out[3] = -u2 * u * (6.0 * f[0] + 6.0 * f[1] + f[2]);
        out[4] = u2 * u2 * (24.0 * f[0] + 36.0 * f[1] + 12.0 * f[2] + f[3]);
        out[5] = logS;
        out[6] = -u * f[0];
    }
};

// The same recurrence in w = 1/(2 alpha), for pig1, whose row reads the
// derivatives of log S in w. S_{y+1} = S_{y-1} + 2(2y - 1) w S_y. The
// derivatives of log S are not formed from the moments N_r = S^(r)(w): the
// cumulants of N_r / S cancel, by about y^r, when the law of the index k of
// the terms of S_y is close to a Poisson law (a small w), which is the range
// of the expected information beyond PIG_ASYM_CUT at a large mu.
//
// What is carried are the factorial cumulants of that law, f_r = w^r L_r
// (r = 1..4), which are the derivatives at x = 1 of g_y(x) = log(S_y(w x) /
// S_y(w)), and their steps a_r = (g_y - g_{y-1})^(r)(1). The recurrence reads
//
//   S_{y+1}(w x) / S_{y+1}(w) = p S_{y-1}(w x) / S_{y-1}(w) +
//                               v x S_y(w x) / S_y(w),
//
// p = S_{y-1}/S_{y+1}, v = 2(2y - 1) w S_y/S_{y+1}, p + v = 1, so that
//
//   g_{y+1} - g_y = log Z,   Z(x) = p exp(-(g_y - g_{y-1})) + v x.
//
// The derivatives of Z at x = 1 are z_1 = v - p a_1, z_r = p E_r (r >= 2),
// with E_r the derivatives of exp(-(g_y - g_{y-1})) (Bell polynomials in the
// a_r), and the new steps are the derivatives of log Z (cumulants of the z_r).
// Then f_r(y + 1) = f_r(y) + a_r(y + 1) and L_r = f_r / w^r. No raw moment
// enters, so the digits of a near-Poisson law are kept; the steps themselves
// cancel by about y in the same range. Every call costs O(1).
struct PigSeqW {
    double w;
    int y;
    double a[4];           // steps of f_r, r = 1..4, from y - 1 to y
    double f[4];           // f_r at y
    double ip[4];          // w^-r
    double logS, e;
    explicit PigSeqW(double ww) : w(ww), y(0), logS(0.0), e(0.0) {
        const double iw = 1.0 / ww;
        for (int r = 0; r < 4; ++r) a[r] = f[r] = 0.0;
        ip[0] = iw; ip[1] = iw * iw; ip[2] = ip[1] * iw; ip[3] = ip[2] * iw;
    }
    void next() {
        if (y == 0) { y = 1; e = 0.0; return; }
        const double two = 2.0 * (2.0 * y - 1.0);
        const double t = two * w;
        const double d = t - e / (1.0 + e);     // S_{y+1}/S_y - 1
        const double R = 1.0 + d;
        const double p = 1.0 / (R * (1.0 + e)), v = t / R;
        const double a1 = a[0], a2 = a[1], a3 = a[2], a4 = a[3];
        const double a11 = a1 * a1;
        const double z1 = v - p * a1;
        const double z2 = p * (a11 - a2);
        const double z3 = p * (a1 * (3.0 * a2 - a11) - a3);
        const double z4 = p * (a11 * (a11 - 6.0 * a2) + 4.0 * a1 * a3 +
                               3.0 * a2 * a2 - a4);
        const double z11 = z1 * z1;
        const double s2 = z2 - z11;
        const double s3 = z3 - z1 * (3.0 * z2 - 2.0 * z11);
        const double s4 = z4 - 4.0 * z1 * z3 - z2 * (3.0 * z2 - 12.0 * z11) -
            6.0 * z11 * z11;
        a[0] = z1; a[1] = s2; a[2] = s3; a[3] = s4;
        f[0] += z1; f[1] += s2; f[2] += s3; f[3] += s4;
        logS += std::log1p(d);
        e = d;
        ++y;
    }
    // log S and its four derivatives in w
    void pre(double out[7]) const {
        out[0] = logS;
        out[1] = f[0] * ip[0];
        out[2] = f[1] * ip[1];
        out[3] = f[2] * ip[2];
        out[4] = f[3] * ip[3];
        out[5] = out[6] = 0.0;
    }
};


// The threshold on tau = (1 + mu) x, x = sigma (pig1) or 1/alpha (pig2),
// below which the series of pig_asym.h replaces the sums. Measured against
// exact sums (stabilita/pig_cmp.R, pig*_dexpected_ref.py): at tau = 0.05 the
// order-20 series is within 2e-15 of every component, where the sums lose up
// to 5e-11; above 0.1 the series converges too slowly.
const double PIG_ASYM_CUT = 0.06;

// slot k of pig_compose(): one partial of h(F), from the partials of F in
// the down-set of k, with the expressions pig_compose() writes
inline double pig_compose_slot(const double h[5], const double* F, int k) {
    switch (k) {
    case 0: return h[0];
    case 1: return h[1] * F[1];
    case 2: return h[1] * F[2];
    case 3: { double m = F[1]; return h[2] * m * m + h[1] * F[3]; }
    case 4: { double m = F[1], u = F[2]; return h[2] * m * u + h[1] * F[4]; }
    case 5: { double u = F[2]; return h[2] * u * u + h[1] * F[5]; }
    case 6: {
        double m = F[1], mm = F[3];
        return h[3] * m * m * m + 3.0 * h[2] * mm * m + h[1] * F[6];
    }
    case 7: {
        double m = F[1], u = F[2], mm = F[3], mu = F[4];
        return h[3] * m * m * u + h[2] * (mm * u + 2.0 * mu * m) + h[1] * F[7];
    }
    case 8: {
        double m = F[1], u = F[2], mu = F[4], uu = F[5];
        return h[3] * m * u * u + h[2] * (uu * m + 2.0 * mu * u) + h[1] * F[8];
    }
    case 9: {
        double u = F[2], uu = F[5];
        return h[3] * u * u * u + 3.0 * h[2] * uu * u + h[1] * F[9];
    }
    case 10: {
        double m = F[1], mm = F[3], mmm = F[6];
        return h[4] * m * m * m * m + 6.0 * h[3] * mm * m * m +
            h[2] * (3.0 * mm * mm + 4.0 * mmm * m) + h[1] * F[10];
    }
    case 11: {
        double m = F[1], u = F[2], mm = F[3], mu = F[4], mmm = F[6],
            mmu = F[7];
        return h[4] * m * m * m * u +
            h[3] * (3.0 * mu * m * m + 3.0 * mm * m * u) +
            h[2] * (3.0 * mm * mu + 3.0 * mmu * m + mmm * u) + h[1] * F[11];
    }
    case 12: {
        double m = F[1], u = F[2], mm = F[3], mu = F[4], uu = F[5],
            mmu = F[7], muu = F[8];
        return h[4] * m * m * u * u +
            h[3] * (mm * u * u + uu * m * m + 4.0 * mu * m * u) +
            h[2] * (mm * uu + 2.0 * mu * mu + 2.0 * muu * m + 2.0 * mmu * u) +
            h[1] * F[12];
    }
    case 13: {
        double m = F[1], u = F[2], mu = F[4], uu = F[5], muu = F[8],
            uuu = F[9];
        return h[4] * m * u * u * u +
            h[3] * (3.0 * mu * u * u + 3.0 * uu * m * u) +
            h[2] * (3.0 * uu * mu + 3.0 * muu * u + uuu * m) + h[1] * F[13];
    }
    default: {
        double u = F[2], uu = F[5], uuu = F[9];
        return h[4] * u * u * u * u + 6.0 * h[3] * uu * u * u +
            h[2] * (3.0 * uu * uu + 4.0 * uuu * u) + h[1] * F[14];
    }
    }
}

inline bool pig_in_support(double yy) {
    return R_finite(yy) && yy >= 0.0 && yy == std::floor(yy);
}

// ---- pig1 (mu, sigma) -------------------------------------------------------

// w = 1/(2 alpha), at which log S is differentiated; pig1_row()'s expression
inline double pig1_w(double m, double sg) {
    double c = 1.0 + 2.0 * sg * m;
    double s = d7::sqrt_cr(c);
    return 0.5 * sg / s;
}

// The partial d^(i+j) l / dmu^i dsigma^j at one response, L being log S_y and
// its derivatives in w (pig_logS_w() or PigSeqW). It runs pig1_row()'s body
// on the slots of the down-set of (i, j) alone.
inline double pig1_slot(int i, int j, double yy, double m, double sg,
                        const double* L) {
    const int k = pix(i, j);
    if (!pig_in_support(yy)) return (k == 0) ? R_NegInf : R_NaN;
    double c = 1.0 + 2.0 * sg * m;
    double s = d7::sqrt_cr(c);
    double w = 0.5 * sg / s;
    double C[15] = {0};
    C[0] = c; C[1] = 2.0 * sg; C[2] = 2.0 * m; C[4] = 2.0;
    double ic = 1.0 / c;
    double hlog[5] = {std::log(c), ic, -ic * ic, 2.0 * ic * ic * ic,
                      -6.0 * ic * ic * ic * ic};
    double hsq[5] = {s, 0.5 / s, -0.25 / (s * c), 0.375 / (s * c * c),
                     -0.9375 / (s * c * c * c)};
    double ip = 1.0 / (1.0 + s);
    double hphi[5] = {ip, -ip * ip, 2.0 * ip * ip * ip,
                      -6.0 * ip * ip * ip * ip, 24.0 * ip * ip * ip * ip * ip};
    double is = 1.0 / s;
    double hpsi[5] = {is, -is * is, 2.0 * is * is * is,
                      -6.0 * is * is * is * is, 24.0 * is * is * is * is * is};
    double S[15] = {0}, PS[15] = {0}, Wp[15] = {0};
    for (int o = 0; o <= i + j; ++o) {
        for (int jj = 0; jj <= o; ++jj) {
            int ii = o - jj;
            if (ii > i || jj > j) continue;
            int kk = pix(ii, jj);
            S[kk] = pig_compose_slot(hsq, C, kk);
            PS[kk] = pig_compose_slot(hpsi, S, kk);
            Wp[kk] = 0.5 * (sg * PS[kk] + (jj > 0 ? jj * PS[pix(ii, jj - 1)] : 0.0));
        }
    }
    Wp[0] = w;
    double PHk = pig_compose_slot(hphi, S, k);
    double PHi = (i > 0) ? pig_compose_slot(hphi, S, pix(i - 1, j)) : 0.0;
    double G = -2.0 * (m * PHk + (i > 0 ? i * PHi : 0.0));
    double LC = pig_compose_slot(hlog, C, k);
    double LS = pig_compose_slot(L, Wp, k);
    double h = -yy / 2.0;
    double val = h * LC + G + LS;
    if (j == 0) {
        double im = 1.0 / m;
        double ym;
        switch (i) {
        case 0: ym = yy * std::log(m); break;
        case 1: ym = yy * im; break;
        case 2: ym = -yy * im * im; break;
        case 3: ym = 2.0 * yy * im * im * im; break;
        default: ym = -6.0 * yy * im * im * im * im;
        }
        val += ym;
    }
    if (k == 0) val -= R::lgammafn(yy + 1.0);
    return val;
}

inline double pig1_value(double yy, double m, double sg, const double* L) {
    return pig1_slot(0, 0, yy, m, sg, L);
}
inline double pig1_score_mu(double yy, double m, double sg, const double* L) {
    return pig1_slot(1, 0, yy, m, sg, L);
}
inline double pig1_score_sigma(double yy, double m, double sg, const double* L) {
    return pig1_slot(0, 1, yy, m, sg, L);
}
inline double pig1_hess_mu_mu(double yy, double m, double sg, const double* L) {
    return pig1_slot(2, 0, yy, m, sg, L);
}
inline double pig1_hess_sigma_sigma(double yy, double m, double sg,
                                    const double* L) {
    return pig1_slot(0, 2, yy, m, sg, L);
}
// the cross partial, which pig1_row() writes at slot pix(1, 1)
inline double pig1_hess_mu_sigma(double yy, double m, double sg,
                                 const double* L) {
    return pig1_slot(1, 1, yy, m, sg, L);
}

// ---- pig2 (mu, alpha) -------------------------------------------------------

// the per-observation pieces pig2_row() forms before its first partial
struct Pig2Base {
    double r, ir, ia, v0, v1, v2, t, b, Q1, Q2, im, r_m, r_a, rm2, ra2, u,
        u_m, T_m, T_a;
};

inline Pig2Base pig2_base(double yy, double m, double aa) {
    Pig2Base B;
    B.r = std::hypot(m, aa);
    B.ir = 1.0 / B.r;
    B.ia = 1.0 / aa;
    B.v0 = B.ia * B.ia;
    B.v1 = -2.0 * B.v0 * B.ia;
    B.v2 = 6.0 * B.v0 * B.v0;
    B.t = m * B.ia;
    double wt = std::hypot(1.0, B.t);
    B.b = aa / (B.t + wt);
    B.Q1 = -yy - B.b;
    B.Q2 = yy + 2.0 * B.b;
    B.im = 1.0 / m;
    B.r_m = m * B.ir;
    B.r_a = aa * B.ir;
    B.rm2 = B.r_m * B.r_m;
    B.ra2 = B.r_a * B.r_a;
    B.u = m + B.r;
    B.u_m = 1.0 + B.r_m;
    double S_m = B.u_m * B.v0;
    double S_a = B.r_a * B.v0 + B.u * B.v1;
    B.T_m = B.b * S_m;
    B.T_a = B.b * S_a;
    return B;
}

// logS the log of S_y; A1 = d log S / d alpha; p2 the second derivative of
// psi(alpha) = -alpha + log S, from psi_derivs() or PigSeq
inline double pig2_value(double yy, double m, double aa, double logS) {
    if (!pig_in_support(yy)) return R_NegInf;
    double r = std::hypot(m, aa);
    double ia = 1.0 / aa;
    double t = m * ia;
    double core = -m * (aa / (m + r)) * (1.0 + m / (aa + r));
    return yy * std::log(m) - yy * std::asinh(t) + core + logS -
        R::lgammafn(yy + 1.0);
}

inline double pig2_score_mu(double yy, double m, const Pig2Base& B) {
    if (!pig_in_support(yy)) return R_NaN;
    double P_m = B.Q1 * B.T_m;
    return yy * B.im + P_m;
}

inline double pig2_score_alpha(double yy, double m, double aa,
                               const Pig2Base& B, double A1) {
    if (!pig_in_support(yy)) return R_NaN;
    double dD_a = -B.r_m * (m / (B.r + aa));
    return yy * m * B.ir * B.ia + A1 + dD_a;
}

inline double pig2_hess_mu_mu(double yy, const Pig2Base& B) {
    if (!pig_in_support(yy)) return R_NaN;
    double r_mm = B.ra2 * B.ir;
    double S_mm = r_mm * B.v0;
    double T_mm = B.b * S_mm;
    double P_mm = B.Q2 * B.T_m * B.T_m + B.Q1 * T_mm;
    return -yy * B.im * B.im + P_mm;
}

inline double pig2_hess_mu_alpha(double yy, const Pig2Base& B) {
    if (!pig_in_support(yy)) return R_NaN;
    double r_ma = -B.r_m * B.r_a * B.ir;
    double S_ma = r_ma * B.v0 + B.u_m * B.v1;
    double T_ma = B.b * S_ma;
    return B.Q2 * B.T_m * B.T_a + B.Q1 * T_ma;
}

inline double pig2_hess_alpha_alpha(double yy, const Pig2Base& B, double p2) {
    if (!pig_in_support(yy)) return R_NaN;
    double r_aa = B.rm2 * B.ir;
    double S_aa = r_aa * B.v0 + 2.0 * B.r_a * B.v1 + B.u * B.v2;
    double T_aa = B.b * S_aa;
    double P_aa = B.Q2 * B.T_a * B.T_a + B.Q1 * T_aa;
    return yy * B.ia * B.ia + p2 + P_aa;
}

// ---- the routers of the scalar registry ------------------------------------

inline void pig1_score_curv(int k, double y, const double* th, double* out) {
    const double m = th[0], sg = th[1];
    double L[5] = {0, 0, 0, 0, 0};
    if (pig_in_support(y)) pig_logS_w(y, pig1_w(m, sg), 2, L);
    if (k == 0) {
        out[0] = pig1_score_mu(y, m, sg, L);
        out[1] = pig1_hess_mu_mu(y, m, sg, L);
    } else {
        out[0] = pig1_score_sigma(y, m, sg, L);
        out[1] = pig1_hess_sigma_sigma(y, m, sg, L);
    }
}

inline void pig2_score_curv(int k, double y, const double* th, double* out) {
    const double m = th[0], aa = th[1];
    double p[5] = {0, 0, 0, 0, 0}, logS = 0.0, A1 = 0.0;
    if (pig_in_support(y)) psi_derivs(y, aa, 2, p, &logS, &A1);
    const Pig2Base B = pig2_base(y, m, aa);
    if (k == 0) {
        out[0] = pig2_score_mu(y, m, B);
        out[1] = pig2_hess_mu_mu(y, B);
    } else {
        out[0] = pig2_score_alpha(y, m, aa, B, A1);
        out[1] = pig2_hess_alpha_alpha(y, B, p[2]);
    }
}

// The (k, k) expected second derivative and its derivative in parameter k,
// by pig_expected_obs()'s pass over the support, the information and its
// derivative in one pass: the information stops at its own rule (pw = 4)
// and the derivative runs on to its rule (pw = 6), each with the sums the
// kernels pig*_expected_cpp and pig*_dexpected1_cpp form. Through the
// second Bartlett identity, E[l_kk] = -E[l_k^2] and
// d_k E[l_kk] = -E[2 l_kk l_k + l_k^3]. Returns in ok[0], ok[1] whether each
// sum was completed.
template <class Seq, class Comp>
inline void pig_diag_obs(double m, double arg, double tail, Comp comp,
                         double* E, double* D, bool* ok) {
    *E = 0.0;
    *D = 0.0;
    ok[0] = ok[1] = false;
    Seq seq(arg);
    double pre[7], cum = 0.0, fprev = 0.0;
    bool e_done = false;
    const long cap = 100000000L;
    for (long yy = 0; yy < cap; ++yy) {
        if (yy > 0) seq.next();
        seq.pre(pre);
        double val, g, h;
        comp((double) yy, pre, &val, &g, &h);
        double f = std::exp(val);
        if (!(f > 0.0)) {
            if ((double) yy > m) {
                ok[0] = ok[1] = true;
                return;
            }
            continue;
        }
        cum += f;
        double gab = g * g;
        if (!e_done) *E -= f * gab;
        *D -= f * (h * g + g * h + gab * g);
        if ((double) yy > m && fprev > 0.0) {
            double r = f / fprev;
            double rt = (r < 1.0) ? r / (1.0 - r) : R_PosInf;
            double base = f * std::max(rt, tail);
            if (!e_done && base * std::pow((yy + 1.0) / (m + 1.0), 4) <=
                1e-17 * cum) {
                e_done = true;
                ok[0] = true;
            }
            if (base * std::pow((yy + 1.0) / (m + 1.0), 6) <= 1e-17 * cum) {
                ok[1] = true;
                return;
            }
        }
        fprev = f;
    }
}

inline void pig1_info_dinfo(int k, double y, const double* th, double* out) {
    (void) y;
    const double m = th[0], x = th[1];
    double E = 0.0, D = 0.0;
    bool ok[2] = {false, false};
    if (R_finite(m) && R_finite(x) && m > 0.0 && x > 0.0) {
        if ((1.0 + m) * x <= PIG_ASYM_CUT) {
            E = (k == 0) ? pig1_asym0_0(m, x) : pig1_asym0_1(m, x);
            D = (k == 0) ? pig1_asym1_0(m, x) : pig1_asym1_3(m, x);
            ok[0] = ok[1] = true;
        } else {
            double w = 0.5 * x / d7::sqrt_cr(1.0 + 2.0 * x * m);
            pig_diag_obs<PigSeqW>(m, w, 2.0 * x * m,
                [&](double yy, const double* pre, double* val, double* g,
                    double* h) {
                    *val = pig1_value(yy, m, x, pre);
                    if (k == 0) {
                        *g = pig1_score_mu(yy, m, x, pre);
                        *h = pig1_hess_mu_mu(yy, m, x, pre);
                    } else {
                        *g = pig1_score_sigma(yy, m, x, pre);
                        *h = pig1_hess_sigma_sigma(yy, m, x, pre);
                    }
                }, &E, &D, ok);
        }
    }
    out[0] = ok[0] ? E : NA_REAL;
    out[1] = ok[1] ? D : NA_REAL;
}

inline void pig2_info_dinfo(int k, double y, const double* th, double* out) {
    (void) y;
    const double m = th[0], x = th[1];
    double E = 0.0, D = 0.0;
    bool ok[2] = {false, false};
    if (R_finite(m) && R_finite(x) && m > 0.0 && x > 0.0) {
        const double xs = 1.0 / x;
        if ((1.0 + m) * xs <= PIG_ASYM_CUT) {
            E = (k == 0) ? pig2_asym0_0(m, xs) : pig2_asym0_1(m, xs);
            D = (k == 0) ? pig2_asym1_0(m, xs) : pig2_asym1_3(m, xs);
            ok[0] = ok[1] = true;
        } else {
            double t = m / x;
            double b = x / (t + std::hypot(1.0, t));
            pig_diag_obs<PigSeq>(m, x, 2.0 * m / b,
                [&](double yy, const double* pre, double* val, double* g,
                    double* h) {
                    *val = pig2_value(yy, m, x, pre[5]);
                    const Pig2Base B = pig2_base(yy, m, x);
                    if (k == 0) {
                        *g = pig2_score_mu(yy, m, B);
                        *h = pig2_hess_mu_mu(yy, B);
                    } else {
                        *g = pig2_score_alpha(yy, m, x, B, pre[6]);
                        *h = pig2_hess_alpha_alpha(yy, B, pre[2]);
                    }
                }, &E, &D, ok);
        }
    }
    out[0] = ok[0] ? E : NA_REAL;
    out[1] = ok[1] ? D : NA_REAL;
}

// ---- the joint registry -----------------------------------------------------
//
// The expected Hessian E[l_ab] (three entries, in the order mu_mu, p2_p2,
// mu_p2) and its derivative in each of the two parameters (six entries,
// D[2 p + c] for pair p and parameter c), by the pass over the support of
// pig_expected_obs() in pig_hd.cpp: E[l_ab] = -E[l_a l_b] and
// d_c E[l_ab] = -E[l_ac l_b + l_a l_bc + l_a l_b l_c], the sums stopping at
// the rules of pw = 4 (E) and pw = 6 (D) of the kernels pig*_expected_cpp
// and pig*_dexpected1_cpp. With WANT_D false the pass ends where the
// information's rule is met, as the kernel of order 0 does. ok[0], ok[1]
// report whether E and D were completed.
template <class Seq, bool WANT_D, class Comp>
inline void pig_joint_obs(double m, double arg, double tail, Comp comp,
                          double* E, double* D, bool* ok) {
    static const int PA[3] = {0, 1, 0};
    static const int PB[3] = {0, 1, 1};
    for (int j = 0; j < 3; ++j) E[j] = 0.0;
    for (int j = 0; j < 6; ++j) D[j] = 0.0;
    ok[0] = ok[1] = false;
    Seq seq(arg);
    double pre[7], cum = 0.0, fprev = 0.0;
    bool e_done = false;
    const long cap = 100000000L;
    for (long yy = 0; yy < cap; ++yy) {
        if (yy > 0) seq.next();
        seq.pre(pre);
        double val, g[2], h[3];
        comp((double) yy, pre, &val, g, h);
        double f = std::exp(val);
        if (!(f > 0.0)) {
            if ((double) yy > m) {
                ok[0] = ok[1] = true;
                return;
            }
            continue;
        }
        cum += f;
        for (int p = 0; p < 3; ++p) {
            int a = PA[p], b = PB[p];
            double ga = g[a], gb = g[b], gab = ga * gb;
            if (!e_done) E[p] -= f * gab;
            if (WANT_D) {
                for (int c = 0; c < 2; ++c) {
                    double gc = g[c];
                    D[2 * p + c] -= f * (h[a + c] * gb + ga * h[b + c] +
                                         gab * gc);
                }
            }
        }
        if ((double) yy > m && fprev > 0.0) {
            double r = f / fprev;
            double rt = (r < 1.0) ? r / (1.0 - r) : R_PosInf;
            double base = f * std::max(rt, tail);
            if (!e_done && base * std::pow((yy + 1.0) / (m + 1.0), 4) <=
                1e-17 * cum) {
                e_done = true;
                ok[0] = true;
                if (!WANT_D) return;
            }
            if (WANT_D &&
                base * std::pow((yy + 1.0) / (m + 1.0), 6) <= 1e-17 * cum) {
                ok[1] = true;
                return;
            }
        }
        fprev = f;
    }
}

// E and D of pig1 at (mu, sigma): the series of pig_asym.h near the Poisson
// limit, the support sum elsewhere, NA where a sum was not completed or the
// parameters are outside the domain (the branches of pig_expected_run())
template <bool WANT_D>
inline void pig1_joint_expected(double m, double x, double* E, double* D) {
    bool ok[2] = {false, false};
    for (int j = 0; j < 3; ++j) E[j] = 0.0;
    for (int j = 0; j < 6; ++j) D[j] = 0.0;
    if (R_finite(m) && R_finite(x) && m > 0.0 && x > 0.0) {
        if ((1.0 + m) * x <= PIG_ASYM_CUT) {
            pig1_asym0(m, x, E);
            if (WANT_D) pig1_asym1(m, x, D);
            ok[0] = ok[1] = true;
        } else {
            double w = 0.5 * x / d7::sqrt_cr(1.0 + 2.0 * x * m);
            pig_joint_obs<PigSeqW, WANT_D>(m, w, 2.0 * x * m,
                [&](double yy, const double* pre, double* val, double* g,
                    double* h) {
                    *val = pig1_value(yy, m, x, pre);
                    g[0] = pig1_score_mu(yy, m, x, pre);
                    g[1] = pig1_score_sigma(yy, m, x, pre);
                    h[0] = pig1_hess_mu_mu(yy, m, x, pre);
                    h[1] = pig1_hess_mu_sigma(yy, m, x, pre);
                    h[2] = pig1_hess_sigma_sigma(yy, m, x, pre);
                }, E, D, ok);
        }
    }
    if (!ok[0]) for (int j = 0; j < 3; ++j) E[j] = NA_REAL;
    if (WANT_D && !ok[1]) for (int j = 0; j < 6; ++j) D[j] = NA_REAL;
}

// E and D of pig2 at (mu, alpha), with the same branches
template <bool WANT_D>
inline void pig2_joint_expected(double m, double x, double* E, double* D) {
    bool ok[2] = {false, false};
    for (int j = 0; j < 3; ++j) E[j] = 0.0;
    for (int j = 0; j < 6; ++j) D[j] = 0.0;
    if (R_finite(m) && R_finite(x) && m > 0.0 && x > 0.0) {
        const double xs = 1.0 / x;
        if ((1.0 + m) * xs <= PIG_ASYM_CUT) {
            pig2_asym0(m, xs, E);
            if (WANT_D) pig2_asym1(m, xs, D);
            ok[0] = ok[1] = true;
        } else {
            double t = m / x;
            double b = x / (t + std::hypot(1.0, t));
            pig_joint_obs<PigSeq, WANT_D>(m, x, 2.0 * m / b,
                [&](double yy, const double* pre, double* val, double* g,
                    double* h) {
                    *val = pig2_value(yy, m, x, pre[5]);
                    const Pig2Base B = pig2_base(yy, m, x);
                    g[0] = pig2_score_mu(yy, m, B);
                    g[1] = pig2_score_alpha(yy, m, x, B, pre[6]);
                    h[0] = pig2_hess_mu_mu(yy, B);
                    h[1] = pig2_hess_mu_alpha(yy, B);
                    h[2] = pig2_hess_alpha_alpha(yy, B, pre[2]);
                }, E, D, ok);
        }
    }
    if (!ok[0]) for (int j = 0; j < 3; ++j) E[j] = NA_REAL;
    if (WANT_D && !ok[1]) for (int j = 0; j < 6; ++j) D[j] = NA_REAL;
}

// The routers of the joint registry, one per order: every component at one
// observation, on the parameter scale (th = mu, then sigma for pig1 and
// alpha for pig2). The Hessian and the expected information fill a full
// 2 x 2 matrix and the derivative of the expected information a 2 x 2 x 2
// array, column-major, out[a + 2 b + 4 c] being the derivative of E[l_ab] in
// parameter c.
inline void pig1_joint_grad(double y, const double* th, double* out) {
    const double m = th[0], sg = th[1];
    double L[5] = {0, 0, 0, 0, 0};
    if (pig_in_support(y)) pig_logS_w(y, pig1_w(m, sg), 1, L);
    out[0] = pig1_score_mu(y, m, sg, L);
    out[1] = pig1_score_sigma(y, m, sg, L);
}

inline void pig1_joint_hess(double y, const double* th, double* out) {
    const double m = th[0], sg = th[1];
    double L[5] = {0, 0, 0, 0, 0};
    if (pig_in_support(y)) pig_logS_w(y, pig1_w(m, sg), 2, L);
    out[0] = pig1_hess_mu_mu(y, m, sg, L);
    out[3] = pig1_hess_sigma_sigma(y, m, sg, L);
    out[1] = out[2] = pig1_hess_mu_sigma(y, m, sg, L);
}

inline void pig1_joint_info(double y, const double* th, double* out) {
    (void) y;
    double E[3], D[6];
    pig1_joint_expected<false>(th[0], th[1], E, D);
    out[0] = E[0];
    out[3] = E[1];
    out[1] = out[2] = E[2];
}

inline void pig1_joint_dinfo(double y, const double* th, double* out) {
    (void) y;
    double E[3], D[6];
    pig1_joint_expected<true>(th[0], th[1], E, D);
    for (int c = 0; c < 2; ++c) {
        out[4 * c] = D[c];
        out[3 + 4 * c] = D[2 + c];
        out[1 + 4 * c] = out[2 + 4 * c] = D[4 + c];
    }
}

inline void pig2_joint_grad(double y, const double* th, double* out) {
    const double m = th[0], aa = th[1];
    double p[5] = {0, 0, 0, 0, 0}, logS = 0.0, A1 = 0.0;
    if (pig_in_support(y)) psi_derivs(y, aa, 1, p, &logS, &A1);
    const Pig2Base B = pig2_base(y, m, aa);
    out[0] = pig2_score_mu(y, m, B);
    out[1] = pig2_score_alpha(y, m, aa, B, A1);
}

inline void pig2_joint_hess(double y, const double* th, double* out) {
    const double m = th[0], aa = th[1];
    double p[5] = {0, 0, 0, 0, 0}, logS = 0.0, A1 = 0.0;
    if (pig_in_support(y)) psi_derivs(y, aa, 2, p, &logS, &A1);
    const Pig2Base B = pig2_base(y, m, aa);
    out[0] = pig2_hess_mu_mu(y, B);
    out[3] = pig2_hess_alpha_alpha(y, B, p[2]);
    out[1] = out[2] = pig2_hess_mu_alpha(y, B);
}

inline void pig2_joint_info(double y, const double* th, double* out) {
    (void) y;
    double E[3], D[6];
    pig2_joint_expected<false>(th[0], th[1], E, D);
    out[0] = E[0];
    out[3] = E[1];
    out[1] = out[2] = E[2];
}

inline void pig2_joint_dinfo(double y, const double* th, double* out) {
    (void) y;
    double E[3], D[6];
    pig2_joint_expected<true>(th[0], th[1], E, D);
    for (int c = 0; c < 2; ++c) {
        out[4 * c] = D[c];
        out[3 + 4 * c] = D[2 + c];
        out[1 + 4 * c] = out[2 + 4 * c] = D[4 + c];
    }
}

} // namespace d7

#endif
